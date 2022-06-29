/*Raul P. Pelaez 2022.
 */
#include"uammd.cuh"
#include"Integrator/BrownianDynamics.cuh"
#include"misc/TabulatedFunction.cuh"
#include"external/spline.h"
#include <cmath>
#include <memory>
#include"Integrator/BDHI/FCM/FCM_impl.cuh"
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
  real m0;
public:

  DryMobilityWithThermalDrift(BD::Parameters par, real dryHydrodynamicRadius, std::vector<real4> data, real Lz):Lz(Lz){
    System::log<System::MESSAGE>("[SelfMobility] Initialized, Lz=%g", Lz);
    dry_detail::initializeTable(data, Lz);
    this->m0 = 1.0/(6*M_PI*par.viscosity*dryHydrodynamicRadius);
  }

  DryMobilityWithThermalDrift(BD::Parameters par, real dryHydrodynamicRadius,
			      std::string fileName, real Lz):Lz(Lz){
    System::log<System::MESSAGE>("[SelfMobility] Initialized, Lz=%g", Lz);
    auto data = dry_detail::readMobilityFile(fileName);
    this->mobilityAndDerivative = dry_detail::initializeTable(data, Lz);
    this->m0 = 1.0/(6*M_PI*par.viscosity*dryHydrodynamicRadius);
  }

  //The position received by this function will be folded to the main cell.
  __device__ auto operator()(real3 pos){
    real z = (pos.z - floor(pos.z/Lz + real(0.5))*Lz)/Lz; //Height in the [-0.5,0.5] interval
    const auto MandDiffM = mobilityAndDerivative(real(2.0)*z);
    real3 M = make_real3(MandDiffM);
    real diffMz = MandDiffM.w;
    return thrust::make_pair(m0*M, m0*make_real3(0, 0, real(2.0)*diffMz/Lz));
  }

};

class DryMobility{
  real m0;
public:
  DryMobility(BD::Parameters par, real dryHydrodynamicRadius){
    this->m0 = 1.0/(6*M_PI*par.viscosity*dryHydrodynamicRadius);
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

enum class update_rules { leimkuhler, euler_maruyama };

struct DryParameters: public BD::Parameters{
  real wetRadius;
  update_rules brownianUpdateRule;
  std::string dryMobilityFile;
  real Lxy;
  real H;
};

class WetMobilityFCM{
  real m0;
  real temperature;
  real dt;
  uint seed1, seed2;
  using FCM = BDHI::FCM_impl<BDHI::FCM_ns::Kernels::Gaussian, BDHI::FCM_ns::Kernels::GaussianTorque>;
  std::shared_ptr<FCM> fcm;
  std::shared_ptr<ParticleData> pd;
  BDHI::cached_vector<real3> hydrodynamicDisplacements;
public:
  WetMobilityFCM(BD::Parameters par, std::shared_ptr<ParticleData> pd, real wetHydrodynamicRadius):
    pd(pd){
    this->temperature = par.temperature;
    this->dt = par.dt;
    this->m0 = 1.0/(6*M_PI*par.viscosity*wetHydrodynamicRadius);
    FCM::Parameters fcm_par;
    fcm_par.hydrodynamicRadius = wetHydrodynamicRadius;
    fcm_par.temperature = par.temperature;
    fcm_par.dt = par.dt;
    fcm_par.box = Box(128);
    fcm_par.viscosity = par.viscosity;
    fcm_par.tolerance = 1e-5;
    fcm = std::make_shared<FCM>(fcm_par);
  }

  void update(){
    auto pos = pd->getPos(access::gpu, access::read);
    auto force = pd->getForce(access::gpu, access::read);
    auto disp = fcm->computeHydrodynamicDisplacements(pos.raw(), force.raw(), nullptr,
						      pos.size(), temperature, 1.0, 0);
    hydrodynamicDisplacements = disp.first;
  }

  real3* getDeterministicVelocities(){
    return thrust::raw_pointer_cast(hydrodynamicDisplacements.data());
  }

  real3* getStochasticVelocities(){
    return nullptr; //thrust::raw_pointer_cast(hydrodynamicDisplacements.data());
  }

};

class WetMobilityDPStokes{
  real m0;
  real temperature;
  real dt;
  uint seed1, seed2;
  std::shared_ptr<DPStokes> dpstokes;
  std::shared_ptr<ParticleData> pd;
  BDHI::cached_vector<real3> hydrodynamicDisplacements;
public:
  WetMobilityDPStokes(DryParameters par, std::shared_ptr<ParticleData> pd, real wetHydrodynamicRadius):
    pd(pd){
    this->temperature = par.temperature;
    this->dt = par.dt;
    this->m0 = 1.0/(6*M_PI*par.viscosity*wetHydrodynamicRadius);
    FCM::Parameters fcm_par;
    fcm_par.hydrodynamicRadius = wetHydrodynamicRadius;
    fcm_par.temperature = par.temperature;
    fcm_par.dt = par.dt;
    fcm_par.box = Box(128);
    fcm_par.viscosity = par.viscosity;
    fcm_par.tolerance = 1e-5;
    fcm = std::make_shared<FCM>(fcm_par);
  }

  void update(){
    auto pos = pd->getPos(access::gpu, access::read);
    auto force = pd->getForce(access::gpu, access::read);
    auto disp = fcm->computeHydrodynamicDisplacements(pos.raw(), force.raw(), nullptr,
						      pos.size(), temperature, 1.0, 0);
    hydrodynamicDisplacements = disp.first;
  }

  real3* getDeterministicVelocities(){
    return thrust::raw_pointer_cast(hydrodynamicDisplacements.data());
  }

  real3* getStochasticVelocities(){
    return nullptr; //thrust::raw_pointer_cast(hydrodynamicDisplacements.data());
  }

};



class BDWithDryDiffusion: public BD::BaseBrownianIntegrator{
public:
  using update_rules = update_rules;
  using DryMobility = DryMobilityWithThermalDrift;
  using WetMobility = WetMobilityFCM;
  using Parameters = DryParameters;
private:
  std::shared_ptr<DryMobility> dryMobility;
  std::shared_ptr<WetMobility> wetMobility;
  thrust::device_vector<real3> noisePrevious;
  update_rules brownian_rule;
public:
  BDWithDryDiffusion(shared_ptr<ParticleData> pd,
		     Parameters par):
    BaseBrownianIntegrator(pd, par){
    this->seed = sys->rng().next32();
    this->steps = 0;
    this->brownian_rule = par.brownianUpdateRule;
    sys->log<System::MESSAGE>("[BDWithThermalDrift] Initialized with seed %u", this->seed);
    wetMobility = std::make_shared<WetMobility>(par, pd, par.wetRadius);
    real dryRadius = 1/( 1/par.hydrodynamicRadius - 1/par.wetRadius);
    dryMobility = std::make_shared<DryMobility>(par, dryRadius, par.dryMobilityFile, par.H);
  }

  void forwardTime() override;

private:
  void updatePositions();
};


namespace BDWithThermalDrift_ns{
  using update_rules = BDWithDryDiffusion::update_rules;

  __device__ real3 genNoise(Saru &rng){
    return make_real3(rng.gf(0, 1), rng.gf(0, 1).x);
  }

  //This integration scheme allows for a self mobility depending on the position.
  //With the associated non-zero thermal drift term.
  template<update_rules rule, class DryMobility>
  __global__ void integrateGPU(real4* pos,
			       ParticleGroup::IndexIterator indexIterator,
			       const int* originalIndex,
			       const real4* force,
			       DryMobility dryMobility,
			       real3 *wetMF,
			       real3 *wetBdW,
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
    R += dt*(dry.first*F + wetMF[i]);
    if(temperature > 0){
      int ori = originalIndex[i];
      Saru rng(ori, stepNum, seed);
      const auto wetnoise = wetBdW?(wetBdW[i]*dt):real3();
      const auto Bn = sqrt(real(2.0)*temperature*dry.first*dt);
      const auto dWn = genNoise(rng);
      const auto noise =  wetnoise + Bn*dWn;
      if(rule == update_rules::euler_maruyama){
	R += noise + temperature*dt*dry.second;
      }
      else if(rule ==update_rules::leimkuhler){
	R += real(0.5)*(noise+noisePrevious[ori]) + temperature*dt*dry.second;
	noisePrevious[ori] = noise;
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
  wetMobility->update();
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
  auto foo =  integrateGPU<update_rules::euler_maruyama, DryMobility>;
  if(brownian_rule == update_rules::leimkuhler){
    foo =  integrateGPU<update_rules::leimkuhler, DryMobility>;
  }
  foo<<<Nblocks, Nthreads, 0, st>>>(pos.raw(),
				    groupIterator,
				    originalIndex,
				    force.raw(),
				    *dryMobility,
				    wetMobility->getDeterministicVelocities(),
				    wetMobility->getStochasticVelocities(),
				    thrust::raw_pointer_cast(noisePrevious.data()),
				    dt,
				    temperature,
				    numberParticles,
				    steps, seed);
}
