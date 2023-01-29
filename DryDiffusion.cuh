/*Raul P. Pelaez 2022.  This  file provides an UAMMD Integrator called
  DryWetBD.  The Integrator separates the dynamics into two parts, one
  is dry (i.e  neglecting hydrodynamic interactions) and  the wet part
  has  an  increased  hydrodynamic radius  but  includes  hydrodynamic
  interactions.

  In particular the DryWetBD Integrator encodes a slit channel geometry.

  The Wet part uses the DPStokes  UAMMD module, while the Dry part has
  a selfmobility that depends on the height (thus  including thermal
  drift).

  The Dry part works by either reading a file with mobilities (see the
  main code  for instructions)  or by  precomputing the  self mobility
  DPStokes  with the  dry radius  and then  spline interpolating  when
  necessary.

  The resulting  dynamics neglect near  field hydrodynamics up  to the
  dry  radius,  but  otherwise  represent a  BDHI  simulation  with  a
  hydrodynamic radius given by:

  a_total^{-1} = a_dry^{-1} + a_wet^{-1}

  This Integrator has the following parameters:

  real hydrodynamicRadius = -1.0: Total hydrodynamic radius
  real wetRadius: The hydrodynamic radius of the wet part update_rules
  brownianUpdateRule:  Can be  either update_rules::euler_maruyama  or
  update_rules::leimkuhler

  std::string dryMobilityFile: Optional. The name of a file containing
  the mobility information  for the mobility. If  present the mobility
  will depend on  the height of the particle according  to the data in
  this file.This file must have two  columns with a list of normalized
  heights  (so Z  must  go from  -1 to  1)  and normalized  mobilities
  (i.e. 6*pi*eta*a*M0)  in X, Y  and Z.  The values for  each particle
  will  be  linearly  interpolated  from  the  data  provided  in  the
  file. The order of the values does not matter. Example:

  --- mobility.dat---
  -1.0 1.0 1.0 1.0
  0.0 1.0 1.0 1.0
  1.0 1.0 1.0 1.0
  -------------------
  If the option is not present the mobility will be autocomputed using
  DPStokes.

  real Lxy: Size of the domain in the plane
  real H: Height of the domain
  real temperature = 0: Temperature in kT
  real viscosity = 1
  real dt = 0: Time step


  Notes: The  DPStokes module  requires partitioning  the domain  in a
  grid with a  size that depends on the hydrodynamic  radius, thus the
  radius cannot be enforced exactly in  general (it can be enforced by
  modifying the box dimensions to be a multiple of the cell size for a
  given hydrodynamic  radius).  Alternatively, the grid  and ES kernel
  parameters can be provided for DPStokes, which allow a fine external
  control of this. This functionality is not currently exposed by this
  interface,  but  this   is  an  easy  modification,   which  can  be
  implemented  by  providing  the  parameters  to  DPStokes  (see  the
  constructor of the WetMobilityDPStokes class).
 */
#include"uammd.cuh"
#include"Integrator/BrownianDynamics.cuh"
#include"misc/TabulatedFunction.cuh"
#include"external/spline.h"
#include <cmath>
#include <memory>
#include "misc/LanczosAlgorithm.cuh"
#include "Integrator/BDHI/FCM/FCM_impl.cuh"
#include"Integrator/BDHI/DoublyPeriodic/DPStokesSlab.cuh"
#include "utils/execution_policy.cuh"
#include<fstream>

using namespace uammd;

namespace dry_detail{

  //This function receives  as input a vector  containing an arbitrary
  //number of points (3 minimum) storing data for the self mobility at
  //different heights.  It returns  three splines, allowing  to sample
  //the mobility  at any  point.  The  data must  contain data  for at
  //least the two domain limits.  The format of the data is: data[i] =
  //{z[i],  Mxx(z[i]), Myy(z[i]),  Mzz(z[i])}
  //The mobilities must be normalized by 6*pi*eta*a The height must be
  //normalized by H/2 (meaning that the walls are at +-1).
  auto splineMobility(std::vector<real4> &data){
    if(data.size() <3) throw std::runtime_error("splineMobility requires at least three points");
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

  //Reads the height vs the three self mobilities from the file, returns a vector of real4 values, each one containing a line of the file, which must comply with the following format:
  //The mobilities must be normalized by 6*pi*eta*a. The height must be
  //normalized by H/2 (meaning that the walls are at +-1).
  //Each line of the file must contain: z Mxx(z) Myy(z) Mzz(z).
  //The format of the data is: data[i] = {z[i],  Mxx(z[i]), Myy(z[i]),  Mzz(z[i])}
  auto readMobilityFile(std::string fileName){
    std::ifstream in(fileName);
    std::istream_iterator<real4> begin(in), end;
    std::vector<real4> data{begin, end};
    return data;
  }

  //Returns a spline functor that returns the derivative of "y" at any
  //point  via spline  interpolation The  input is  any callable,  for
  //instance  a spline.  N  is the  number of  sample  points used  to
  //construct the spline for the derivative.
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

  //Constructs a  tabulated function for  GPU use, allowing  to sample
  //Mxx, Myy, Mzz,  div(Mzz) for any height, stored as  a simple real4
  //value.  The  data must contain  data for  at least the  two domain
  //limits.  The  format of the  data is: data[i] =  {z[i], Mxx(z[i]),
  //Myy(z[i]),  Mzz(z[i])}  The  mobilities   must  be  normalized  by
  //6*pi*eta*a The height must be  normalized by H/2 (meaning that the
  //walls are at +-1).
  auto initializeTable(std::vector<real4> &data, real Lz){
    tk::spline mobilityx, mobilityy, mobilityz;
    std::tie(mobilityx, mobilityy, mobilityz) = dry_detail::splineMobility(data);
    const real hdiff = 1e-3;
    const int ntablePoints = std::max(int(Lz/hdiff), 1<<20);
    auto allM =[&](real z){return make_real4(mobilityx(z),mobilityy(z),mobilityz(z), mobilityz.deriv(1, z));};
    return TabulatedFunction<real4>(ntablePoints, -1, 1, allM);
  }

}

//This class provides the dry part of the mobility.
//The  () operator  of this  class returns  the normalized  self
// mobility and  its derivative when  given a position,  i.e returning
// 6*pi*eta*a*{Mxx(z), Myy(z), Mzz(z), \nabla Mzz(z)}
class DryMobilityWithThermalDrift{
  TabulatedFunction<real4> mobilityAndDerivative;
  real Lz;
  real m0;
public:
  //The default constructor
  DryMobilityWithThermalDrift(){}

  //This constructor requires a  viscosity, a dry hydrodynamic radius,
  //the height  of the domain  and a  vector with mobility  data.  The
  //data must  contain data for at  least the two domain  limits.  The
  //format  of the  data is:  data[i] =  {z[i], Mxx(z[i]),  Myy(z[i]),
  //Mzz(z[i])}  The mobilities  must be  normalized by  6*pi*eta*a The
  //height must  be normalized by H/2  (meaning that the walls  are at
  //+-1).
  DryMobilityWithThermalDrift(real viscosity, real dryHydrodynamicRadius, std::vector<real4> data, real Lz):Lz(Lz){
    System::log<System::MESSAGE>("[SelfMobility] Initialized, Lz=%g", Lz);
    dry_detail::initializeTable(data, Lz);
    this->m0 = 1.0/(6*M_PI*viscosity*dryHydrodynamicRadius);
  }

  //This constructor requires a  viscosity, a dry hydrodynamic radius,
  //the  height of  the  domain  and the  name  of  a file  containing
  //mobility data as required by the readMobilityFile function.
  DryMobilityWithThermalDrift(real viscosity, real dryHydrodynamicRadius,
			      std::string fileName, real Lz):Lz(Lz){
    System::log<System::MESSAGE>("[SelfMobility] Initialized, Lz=%g", Lz);
    auto data = dry_detail::readMobilityFile(fileName);
    this->mobilityAndDerivative = dry_detail::initializeTable(data, Lz);
    this->m0 = 1.0/(6*M_PI*viscosity*dryHydrodynamicRadius);
  }

  //Given a  position this  function returns  the self  mobility (Mxx,
  //Myy,  Mzz) and  the  divergence of  Mzz with  respect  to z.   The
  //position received by this function will be folded to the main cell
  //(-0.5, 0.5) in the third direction.
  __device__ auto operator()(real3 pos){
    real z = (pos.z - floor(pos.z/Lz + real(0.5))*Lz)/Lz; //Height in the [-0.5,0.5] interval
    const auto MandDiffM = mobilityAndDerivative(real(2.0)*z);
    real3 M = make_real3(MandDiffM);
    real diffMz = MandDiffM.w;
    return thrust::make_pair(m0*M, m0*make_real3(0, 0, real(2.0)*diffMz/Lz));
  }

};

//An enumerator  to distinguish  between the different  update schemes
//accepted by WetDryBD
enum class update_rules { leimkuhler, euler_maruyama };

//Parameters for the WetDryBD Integrator
struct WetDryParameters: BD::Parameters{
  real wetRadius;
  update_rules brownianUpdateRule;
  std::string dryMobilityFile;
  real Lxy;
  real H;
  // int w;
  // int nxy_stokes;
  // int nz_stokes;
  // real beta;
};

namespace dry_detail{

  //Returns  a  vector  of  size  size  filled  with  random  gaussian
  //distributed  numbers with  mean  0  and std  1.  Using the  second
  //argument as seed.
  BDHI::cached_vector<real3> createGaussianNoise(int size, uint seed){
    BDHI::cached_vector<real3> dW(size);
    auto cit = thrust::make_counting_iterator(0);
    thrust::transform(cit, cit+ size,
		      dW.begin(),
		      [=]__device__(int i){
			Saru rng(seed, i);
			return make_real3(rng.gf(0,1), rng.gf(0,1).x);
		      });
    return dW;
  }
}

//From the box size (Lxy in the  plane and H in height), viscosity and
//a hydrodynamic radius this function returns the parameters needed by
//DPStokes. Parameters for a support of w=6 according to the paper
//Note that this function does not set the temperature nor the time step.
auto getDPStokesParamtersOnlyForce(real Lxy, real H, real viscosity, real hydrodynamicRadius){
  real h = hydrodynamicRadius/1.554;
  int nxy = int(Lxy/h +0.5);
  DPStokesSlab_ns::DPStokes::Parameters par;
  par.nx = nxy;
  par.ny = par.nx;
  par.nz = int(M_PI*H/(2*h));
  par.w = 6;
  par.beta = 1.714*par.w;
  par.alpha = par.w*0.5;
  par.mode = DPStokesSlab_ns::WallMode::slit;
  par.viscosity = viscosity;
  par.Lx = par.Ly = Lxy;
  par.H = H;
  par.tolerance = 1e-4;
  return par;
}

//This  boilerplate  function  only  exists because  the  two  classes
//exposing   the   DPStokes   algorithm   in   UAMMD   (DPStokes   and
//DPStokesIntegrator)  use   slightly  different  structs   for  their
//parameters. This function calls getDPStokesParameters and adapts the
//output.
using DPStokes = DPStokesSlab_ns::DPStokesIntegrator;
auto getDPStokesIntegratorParamtersOnlyForce(real Lxy, real H, real viscosity, real hydrodynamicRadius){
  auto par = getDPStokesParamtersOnlyForce(Lxy, H, viscosity, hydrodynamicRadius);
  DPStokes::Parameters pari;
  pari.nx = par.nx;
  pari.ny = par.ny;
  pari.nz = par.nz;
  pari.w = par.w;
  pari.beta = par.beta;
  pari.alpha = par.alpha;
  pari.mode = DPStokesSlab_ns::WallMode::slit;
  pari.viscosity = par.viscosity;
  pari.Lx = par.Lx;
  pari.Ly= par.Ly;
  pari.H = H;
  return pari;
}

//This class provides the wet part of the mobility. It exposes functions to get the different terms in the BDHI equation (mainly MF, BdW and kT div(M)).
class WetMobilityDPStokes{
  std::shared_ptr<DPStokes> dpstokes;
  std::shared_ptr<ParticleData> pd;
  BDHI::cached_vector<real3> hydrodynamicDisplacements, fluctuations, thermalDrift;
  real temperature = 0;
public:
  //This constructor requires a ParticleData instance and a WetDryParameters
  WetMobilityDPStokes(WetDryParameters par, std::shared_ptr<ParticleData> pd):
    pd(pd){
    auto dpstokes_par = getDPStokesIntegratorParamtersOnlyForce(par.Lxy, par.H,
								par.viscosity, par.wetRadius);
    dpstokes_par.dt = par.dt;
    dpstokes_par.temperature = par.temperature;
    this->temperature = par.temperature;
    dpstokes = std::make_shared<DPStokes>(pd, dpstokes_par);
  }

  //Update the DPStokesIntegrator with the latest particle positions and forces
  void update(){
    hydrodynamicDisplacements = dpstokes->computeDeterministicDisplacements();
    if(temperature){
      fluctuations = dpstokes->computeFluctuations();
      thermalDrift = dpstokes->computeThermalDrift();
    }
  }

  //Provides a pointer to the deterministic displacements: MF
  real3* getDeterministicVelocities(){
    return thrust::raw_pointer_cast(hydrodynamicDisplacements.data());
  }

  //Provides a pointer to the stochastic displacements: BdW
  real3* getStochasticVelocities(){
    return thrust::raw_pointer_cast(fluctuations.data());
  }

  //Provides a pointer to the thermal drift term: kT*dt div_z(M)
  real3* getThermalDrift(){
    return thrust::raw_pointer_cast(thermalDrift.data());
  }
};

//Computes a vector with mobility data  as required by the Dry part of
//the  mobility. Uses  the  DPStokes algorithm  with the  hydrodynamic
//radius  in  the second  argument  and  samples the  selfmobility  at
//several heights.
auto computeMobilityDataForDryDiffusion(WetDryParameters par,
					real hydrodynamicRadius){
  auto dppar = getDPStokesParamtersOnlyForce(par.Lxy, par.H, par.viscosity, hydrodynamicRadius);
  auto dpstokes = std::make_shared<DPStokesSlab_ns::DPStokes>(dppar);
  constexpr int nsamples = 1000;
  std::vector<real4> mobilityData(nsamples);
  for(int i = 0; i<nsamples; i++){
    real z = -par.H*0.5 + par.H*(i/real(nsamples-1));
    mobilityData[i].x = z;
    auto pos = thrust::make_constant_iterator<real4>({0,0,z,0});
    //Mxx
    auto disp = dpstokes->Mdot(pos, thrust::make_constant_iterator<real4>({1,0,0,0}), 1, 0);
    //Myy
    mobilityData[i].y = 6*M_PI*par.viscosity*hydrodynamicRadius*real3(disp[0]).x;
    disp = dpstokes->Mdot(pos, thrust::make_constant_iterator<real4>({0,1,0,0}), 1, 0);
    //Mzz
    mobilityData[i].z = 6*M_PI*par.viscosity*hydrodynamicRadius*real3(disp[0]).y;
    disp = dpstokes->Mdot(pos, thrust::make_constant_iterator<real4>({0,0,1,0}), 1, 0);
    mobilityData[i].w = 6*M_PI*par.viscosity*hydrodynamicRadius*real3(disp[0]).z;
  }
  return mobilityData;
}


//This class  provides an Integrator for  the BDHI equation in  a slit
// channel by separating each term into two contributions: A wet and a
// dry part.  Both the dry and  wet parts are provided  by the classes
// above    (DryMobilityWithThermalDrift    and    WetMobilityDPStokes
// respectively)
class DryWetBD: public BD::BaseBrownianIntegrator{
public:
  using update_rules = update_rules;
  using DryMobility = DryMobilityWithThermalDrift;
  using WetMobility = WetMobilityDPStokes;
  using Parameters = WetDryParameters;
private:
  std::shared_ptr<DryMobility> dryMobility;
  std::shared_ptr<WetMobility> wetMobility;
  thrust::device_vector<real3> noisePrevious; //Used to store the noise of the previous step in Leimkuhler
  update_rules brownian_rule;
  bool isFullWet = false;
  bool isFullDry = false;

public:
  DryWetBD(shared_ptr<ParticleData> pd,
	   Parameters par):
    BaseBrownianIntegrator(pd, par){
    this->seed = sys->rng().next32();
    this->steps = 0;
    this->brownian_rule = par.brownianUpdateRule;
    sys->log<System::MESSAGE>("[BDWithThermalDrift] Initialized with seed %u", this->seed);
    real dryRadius = 0;
    if(par.wetRadius <= 0){
      dryRadius = par.hydrodynamicRadius;
      this->isFullDry = true;
      sys->log<System::MESSAGE>("[BDWithThermalDrift] Enabling full dry mode");
    }
    else if(par.wetRadius <= par.hydrodynamicRadius){
      par.wetRadius = par.hydrodynamicRadius;
      dryRadius = 0;
    }
    else{
      dryRadius = 1/( 1/par.hydrodynamicRadius - 1/par.wetRadius);
    }
    if(dryRadius > 0){
      sys->log<System::MESSAGE>("[BDWithThermalDrift] Dry radius is: %g", dryRadius);
      sys->log<System::MESSAGE>("[BDWithThermalDrift] Wet radius is: %g", par.wetRadius);
      if(par.dryMobilityFile.empty()){
	sys->log<System::MESSAGE>("[BDWithThermalDrift] Computing self mobility data using DPStokes ");
	auto mobilityData = computeMobilityDataForDryDiffusion(par, dryRadius);
	dryMobility = std::make_shared<DryMobility>(par.viscosity, dryRadius, mobilityData, par.H);
      }
      else{
	sys->log<System::MESSAGE>("[BDWithThermalDrift] Reading self mobility from" +
				  par.dryMobilityFile);
	dryMobility = std::make_shared<DryMobility>(par.viscosity, dryRadius, par.dryMobilityFile, par.H);
      }
    }
    else{
      sys->log<System::MESSAGE>("[BDWithThermalDrift] Enabling full wet mode");
      this->isFullWet = true;
      dryMobility = std::make_shared<DryMobility>(); //A dummy instance
    }
    if(not isFullDry)
      wetMobility = std::make_shared<WetMobility>(par, pd);
  }

  void forwardTime() override;

private:
  void updatePositions();
};


namespace BDWithThermalDrift_ns{
  using update_rules = DryWetBD::update_rules;
  //Draws three Gaussian random numbers from the rng
  __device__ real3 genNoise(Saru &rng){
    return make_real3(rng.gf(0, 1), rng.gf(0, 1).x);
  }

  //Adds to each particle the displacement coming from the dry and wet parts of the BDHI eq.
  //It encodes two update schemes.
  //Euler Maruyama
  // dX^n = dt*( (MF + BdW + kTdiv_z(M))^n_wet + (MF + BdW + kTdiv_z(M))^n_dry)
  //Leimkuhler
  // dX^n = dt*( (MF + kTdiv_z(M))^n_wet + (MF + kTdiv_z(M))^n_dry
  //+ 0.5((BdW)^n_dry+(BdW)^{n-1}_dry)dt
  //+ 0.5((BdW)^{n-1}_wet+(BdW)^{n-1}_wet)dt
  //)
  template<update_rules rule, class DryMobility>
  __global__ void integrateGPU(real4* pos,
			       ParticleGroup::IndexIterator indexIterator,
			       const int* originalIndex,
			       const real4* force,
			       DryMobility dryMobility,
			       real3 *wetMF,
			       real3 *wetBdW,
			       real3 *wetThermalDrift,
			       real3* noisePrevious,
			       bool isFullWet,
			       real dt,
			       real temperature,
			       int N,
			       uint stepNum, uint seed){
    uint id = blockIdx.x*blockDim.x+threadIdx.x;
    if(id>=N) return;
    int i = indexIterator[id];
    real3 R = make_real3(pos[i]);
    const real3 F = make_real3(force[i]);
    const auto dry = isFullWet?thrust::make_pair(real3(), real3()):dryMobility(R);
    auto wetMF_ = wetMF?wetMF[i]:real3();
    R += dt*(dry.first*F + wetMF_);
    if(temperature > 0){
      int ori = originalIndex[i];
      Saru rng(ori, stepNum, seed);
      const auto wetnoise = wetBdW?(wetBdW[i]*dt):real3();
      const auto Bn = sqrt(real(2.0)*temperature*dry.first*dt);
      const auto dWn = genNoise(rng);
      const auto noise =  wetnoise + Bn*dWn;
      if(rule == update_rules::euler_maruyama){
	R += noise;
      }
      else if(rule ==update_rules::leimkuhler){
	R += real(0.5)*(noise+noisePrevious[ori]);
	noisePrevious[ori] = noise;
      }
      //Add thermal drift
      auto wetTD = wetThermalDrift?wetThermalDrift[i]:real3();
      R += wetTD + temperature*dt*dry.second;
    }
    pos[i].x = R.x;
    pos[i].y = R.y;
    pos[i].z = R.z;
  }

}

//Takes the particles to the next time step
void DryWetBD::forwardTime(){
  steps++;
  sys->log<System::DEBUG1>("[BD::Leimkuhler] Performing integration step %d", steps);
  updateInteractors();
  computeCurrentForces();
  if(not isFullDry){
    wetMobility->update();
  }
  updatePositions();
}

void DryWetBD::updatePositions(){
  int numberParticles = pg->getNumberParticles();
  noisePrevious.resize(numberParticles);
  if(steps==1)
    thrust::fill(uammd::cached_device_execution_policy.on(st),
		 noisePrevious.begin(), noisePrevious.end(), real3());
  auto groupIterator = pg->getIndexIterator(access::location::gpu);
  auto pos = pd->getPos(access::location::gpu, access::mode::readwrite);
  auto force = pd->getForce(access::location::gpu, access::mode::read);
  auto originalIndex = pd->getIdOrderedIndices(access::location::gpu);
  using namespace BDWithThermalDrift_ns;
  auto foo =  integrateGPU<update_rules::euler_maruyama, DryMobility>;
  if(brownian_rule == update_rules::leimkuhler){
    foo =  integrateGPU<update_rules::leimkuhler, DryMobility>;
  }
  int BLOCKSIZE = 128;
  uint Nthreads = BLOCKSIZE<numberParticles?BLOCKSIZE:numberParticles;
  uint Nblocks = numberParticles/Nthreads +  ((numberParticles%Nthreads!=0)?1:0);
  real3* wetMF   = wetMobility?wetMobility->getDeterministicVelocities():nullptr;
  real3* wetBdW  = wetMobility?wetMobility->getStochasticVelocities():nullptr;
  real3* wetDivM = wetMobility?wetMobility->getThermalDrift():nullptr;
  foo<<<Nblocks, Nthreads, 0, st>>>(pos.raw(),
				    groupIterator,
				    originalIndex,
				    force.raw(),
				    *dryMobility,
				    wetMF, wetBdW, wetDivM,
				    thrust::raw_pointer_cast(noisePrevious.data()),
				    this->isFullWet,
				    dt,
				    temperature,
				    numberParticles,
				    steps, seed);
}
