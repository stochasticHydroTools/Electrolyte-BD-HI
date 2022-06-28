/*Raul P. Pelaez 2022.
 */
#include"uammd.cuh"
#include"Integrator/BrownianDynamics.cuh"
#include"misc/TabulatedFunction.cuh"
#include"external/spline.h"
#include <cmath>
#include <memory>

using namespace uammd;

namespace dry_detail{

  auto splineMobility(std::vector<real4> &data){
    //The spline library needs input data in ascending order
    std::sort(data.begin(), data.end(),[](real4 a, real4 b){return a.x<b.x;});
    std::vector<double> Z(data.size());
    auto Mx = Z;
    auto My = Mx;
    auto Mz = My;
    std::transform(data.begin(), data.end(), Z.begin(),  [](real4 a){return a.x;});
    std::transform(data.begin(), data.end(), Mx.begin(), [](real4 a){return a.y;});
    std::transform(data.begin(), data.end(), My.begin(), [](real4 a){return a.z;});
    std::transform(data.begin(), data.end(), Mz.begin(), [](real4 a){return a.w;});
    tk::spline mobilityx, mobilityy, mobilityz;
    mobilityx.set_points(Z, Mx);
    mobilityy.set_points(Z, My);
    mobilityz.set_points(Z, Mz);
    return std::make_tuple(mobilityx, mobilityy, mobilityz);
  }

  //Reads the height vs the three self mobilities from the file
  auto readMobilityFile(std::string fileName){
    std::ifstream in(fileName);
    std::istream_iterator<real4> begin(in), end;
    std::vector<real4> data{begin, end};
    return data;
  }
  //Returns a functor that returns the derivative of "y" at any point via spline interpolation
  template<class Functor>
  auto computeDerivative(Functor y, int N){
    std::vector<double> derivative(N,0);
    auto x = derivative;
    double h = 2.0/(N-1);
    x[0] = -1;
    derivative[0] = (y(-1+h)-y(-1))/h;
    for(int i = 1; i<N-1; i++){
      real z = -1+h*i;
      x[i] = z;
      derivative[i] = (y(z+h)-y(z-h))/(2.0*h);
    }
    x[N-1] = 1;
    derivative[N-1] = (y(1)-y(1-h))/h;
    tk::spline sy;
    sy.set_points(x, derivative);
    return sy;
  }

  auto initializeTable(std::vector<real4> &data, real Lz){
    tk::spline mobilityx, mobilityy, mobilityz;
    std::tie(mobilityx, mobilityy, mobilityz) = dry_detail::splineMobility(data);
    const real hdiff = 1e-3;
    const int ntablePoints = std::max(int(Lz/hdiff), 1<<20);
    auto allM =[&](real z){return make_real4(mobilityx(z),mobilityy(z),mobilityz(z), mobilityz.deriv(1, z));};
    return TabulatedFunction<real4>(ntablePoints, -1, 1, allM);
  }

}

//The () operator of this struct must return the normalized self mobility and its derivative when given a position,
// i.e returning 6*pi*eta*a*{M0(x,y,z), \nabla M0(x,y,z)}
class DryMobilityWithThermalDrift{

  TabulatedFunction<real4> mobilityAndDerivative;
  real Lz;

public:

  DryMobilityWithThermalDrift(std::vector<real4> data, real Lz):Lz(Lz){
    System::log<System::MESSAGE>("[SelfMobility] Initialized, Lz=%g", Lz);
    dry_detail::initializeTable(data, Lz);
  }

  DryMobilityWithThermalDrift(std::string fileName, real Lz):Lz(Lz){
    System::log<System::MESSAGE>("[SelfMobility] Initialized, Lz=%g", Lz);
    auto data = dry_detail::readMobilityFile(fileName);
    this->mobilityAndDerivative = dry_detail::initializeTable(data, Lz);
  }

  //The position received by this function will be folded to the main cell.
  __device__ auto operator()(real3 pos){
    real z = (pos.z - floor(pos.z/Lz + real(0.5))*Lz)/Lz; //Height in the [-0.5,0.5] interval
    const auto MandDiffM = mobilityAndDerivative(real(2.0)*z);
    real3 M = make_real3(MandDiffM);
    real diffMz = MandDiffM.w;
    return thrust::make_pair(M, make_real3(0, 0, real(2.0)*diffMz/Lz));
  }

};

class DryMobility{
  real m0;
public:
  DryMobility(BD::Parameters par, real dryHydrodynamicRadius){
    this->m0 = par.hydrodynamicRadius/dryHydrodynamicRadius;
  }

  __device__ auto operator()(real3 pos){
    return thrust::make_pair(make_real3(m0), make_real3(0, 0, 0));
  }


};

// class WetMobility{
//   using DPStokes = DPStokesSlab_ns::DPStokes;
//   std::shared_ptr<DPStokes> dpstokes;

//   auto createDPStokes(std::shared_ptr<ParticleData> pd){
//     DPStokes::Parameters par;
//     par.nx = 1;
//     par.ny = 1;
//     par.nz = -1;
//     par.dt = 1;
//     par.viscosity = 1;
//     par.Lx = 0;
//     par.Ly = 0;
//     par.H = 0;
//     par.tolerance = 1e-7;
//     par.w = 0;
//     par.w_d = 0;
//     par.hydrodynamicRadius = 0;
//     par.beta = -1;
//     par.beta_d = -1;
//     par.alpha = -1;
//     par.alpha_d = -1;
//     par.mode = DPStokesSlab_ns::WallMode::bottom;

//     return std::make_shared<DPStokes>(pd, par);

//   }


// public:

//   WetMobility(std::shared_ptr<ParticleData> pd){
//     dpstokes = createDPStokes(pd);
//   }

//};




class WetMobility{
  real m0;
public:
  WetMobility(BD::Parameters par, real wetHydrodynamicRadius){
    this->m0 = par.hydrodynamicRadius/wetHydrodynamicRadius;
  }

  __device__ auto operator()(real3 pos){
    return thrust::make_pair(make_real3(m0), make_real3(0, 0, 0));
  }

};


auto createDryMobility(BD::Parameters par, real dryHydrodynamicRadius){
  return std::make_shared<DryMobility>(par, dryHydrodynamicRadius);
}

auto createWetMobility(BD::Parameters par, real wetHydrodynamicRadius){
  return std::make_shared<WetMobility>(par, wetHydrodynamicRadius);
}

class BDWithDryDiffusion: public BD::BaseBrownianIntegrator{
public:
  enum class update_rules{leimkuhler, euler_maruyama};
private:
  std::shared_ptr<DryMobility> dryMobility;
  std::shared_ptr<WetMobility> wetMobility;
  thrust::device_vector<real3> noisePrevious;
  update_rules brownian_rule;
public:
  struct Parameters: public BD::Parameters{
    real wetRadius;
    update_rules brownianUpdateRule;
  };

  BDWithDryDiffusion(shared_ptr<ParticleData> pd,
		     Parameters par):
    BaseBrownianIntegrator(pd, par){
    this->seed = sys->rng().next32();
    this->steps = 0;
    this->brownian_rule = par.brownianUpdateRule;
    sys->log<System::MESSAGE>("[BDWithThermalDrift] Initialized with seed %u", this->seed);
    wetMobility = createWetMobility(par, par.wetRadius);
    real dryRadius = 1/(1/par.wetRadius + 1/par.hydrodynamicRadius);
    dryMobility = createDryMobility(par, dryRadius);
  }

  void forwardTime() override;

private:
  void updatePositions();
};


namespace BDWithThermalDrift_ns{
  using update_rules = BDWithDryDiffusion::update_rules;

  __device__ real3 genNoise(int i, uint stepNum, uint seed){
    Saru rng(i, stepNum, seed);
    return make_real3(rng.gf(0, 1), rng.gf(0, 1).x);
  }

  //This integration scheme allows for a self mobility depending on the position.
  //With the associated non-zero thermal drift term.
  template<update_rules rule, class DryMobility, class WetMobility>
  __global__ void integrateGPU(real4* pos,
			       ParticleGroup::IndexIterator indexIterator,
			       const int* originalIndex,
			       const real4* force,
			       real selfMobility,
			       DryMobility dryMobility,
			       WetMobility wetMobility,
			       real3* noisePrevious,
			       real dt,
			       real temperature,
			       int N,
			       uint stepNum, uint seed){
    uint id = blockIdx.x*blockDim.x+threadIdx.x;
    if(id>=N) return;
    int i = indexIterator[id];
    real3 R = make_real3(pos[i]);
    const real3 F = make_real3(force[i]);
    const auto dry = dryMobility(R);
    const auto wet = wetMobility(R);
    const real3 M = selfMobility*(dry.first+wet.first);
    R += dt*(M*F);
    if(temperature > 0){
      int ori = originalIndex[i];
      const auto Bn = sqrt(real(2.0)*temperature*M*dt);
      const auto dWn = genNoise(ori, stepNum, seed);
      if(rule == update_rules::euler_maruyama){
	R += Bn*dWn + temperature*dt*selfMobility*(dry.second + wet.second);
      }
      else if(rule ==update_rules::leimkuhler){
	R += real(0.5)*(Bn*dWn+noisePrevious[ori]) + temperature*dt*selfMobility*(dry.second+wet.second);
	noisePrevious[ori] = Bn*dWn;
      }
    }
    pos[i].x = R.x;
    pos[i].y = R.y;
    pos[i].z = R.z;
  }

}

void BDWithDryDiffusion::forwardTime(){
  steps++;
  sys->log<System::DEBUG1>("[BD::Leimkuhler] Performing integration step %d", steps);
  updateInteractors();
  computeCurrentForces();
  updatePositions();
}

void BDWithDryDiffusion::updatePositions(){
  int numberParticles = pg->getNumberParticles();
  noisePrevious.resize(numberParticles);
  if(steps==1)
    thrust::fill(noisePrevious.begin(), noisePrevious.end(), real3());
  int BLOCKSIZE = 128;
  uint Nthreads = BLOCKSIZE<numberParticles?BLOCKSIZE:numberParticles;
  uint Nblocks = numberParticles/Nthreads +  ((numberParticles%Nthreads!=0)?1:0);
  auto groupIterator = pg->getIndexIterator(access::location::gpu);
  auto pos = pd->getPos(access::location::gpu, access::mode::readwrite);
  auto force = pd->getForce(access::location::gpu, access::mode::read);
  auto originalIndex = pd->getIdOrderedIndices(access::location::gpu);
  using namespace BDWithThermalDrift_ns;
  auto foo =  integrateGPU<update_rules::euler_maruyama, DryMobility, WetMobility>;
  if(brownian_rule == update_rules::leimkuhler){
    foo =  integrateGPU<update_rules::leimkuhler, DryMobility, WetMobility>;
  }
  foo<<<Nblocks, Nthreads, 0, st>>>(pos.raw(),
				    groupIterator,
				    originalIndex,
				    force.raw(),
				    selfMobility,
				    *dryMobility,
				    *wetMobility,
				    thrust::raw_pointer_cast(noisePrevious.data()),
				    dt,
				    temperature,
				    numberParticles,
				    steps, seed);
}
