/*Raul P.  Pelaez 2020-2023. Dry and  wet diffusion in a slit channel
  with electrostatic interactions.

  About this code
  ------------------

  Besides including many functionalities  from UAMMD, this source also
  includes two headers that are part of this particular project. These
  are:
  1.  RepulsivePotential:  The  definition   of  a  LJ-like  potential
  compatible with UAMMD
  2.  DryDiffusion:  An UAMMD  Integrator  that  mixes UAMMD's  Doubly
  Periodic Stokes  Integrator with  a bare-diffusion  Integrator. This
  allows   to  separate   the  mobility   in  two   parts,  one   with
  hydro. interactions  for the far  field (expensive) and one  for the
  self mobility (fast).

  This source is organized as a series of functions and utility structures.

  Many  of  this  functions   are  called  createSomething.  In  these
  instances, Something refers to some UAMMD-related structure, such as
  the Integrator or the class that takes care of electrostatics.

  The physical domain being simulated is  walled, and as such is prone
  to a particular error where a  particle teleports outside a wall due
  to numerical error.  This source includes a functionality consisting
  of saving the  simulation state every few steps and  rewinding it in
  the  case of  something  like this  happening.  These functions  are
  checkOverlap, saveConfiguration, restoreLastSavedConfiguration,...

  The initialization of  the UAMMD structures and  the main simulation
  loop is carried  out in main().  The length of  this function is due
  to  the implementation  of the  rewind functionality,  which in  its
  current state suffers from a lot of code repetition.

  Basicaly, initialization comes  by calling createIntegrator followed
  by  a series  of calls  to  addInteractor (which  adds a  particular
  interaction to the Integrator) depending  on the options.  Then, the
  simulation is taken  to the next step by calling  forwardTime on the
  Integrator.   Every now  and then  the  state of  the simulation  is
  probed via sim.pd (an instance  of ParticleData, the UAMMD structure
  that holds particle states) and written to disk.

  This code  is GPU enabled  via CUDA, which sometimes  requires doing
  things in  a kind of  convoluted manner. For instance,  the function
  checkOverlap  requires  defining  a structure  called  CheckOverlap,
  whose operator() member  has a fancy __device__  decorator.  Then, a
  thrust  algorithm is  used  to find  possible  overlaps.  You  might
  wonder  why a  simple loop  cannot be  used here.  We could  do that
  following the  rules of CUDA,  but the  boiler plate would  be worse
  and, you'll  have to believe  me, the  resulting code would  be even
  more obtuse  (and inefficient, for  what is worth).

  In CUDA, and in general when writting parallel code in C++, thinking
  in   terms  of   these  "algorithms"   (such  as   copy,  transform,
  find_if,...) is  many times "the  way", in terms  of expressiveness,
  efficiency  and   versatility.   For  instance,  we   can  make  the
  checkOverlap   function   run  in   the   CPU   just  by   replacing
  thrust::cuda::par by thrust::cuda::host. Or the other way around. So
  once you learn to think in  data-driven algorithms, you can get away
  with  not learning  to code  in  a GPU,  since switching  to a  CUDA
  implementation consists of changing a simple word.

*/
#include"uammd.cuh"
#include"RepulsivePotential.cuh"
#include"Interactor/PairForces.cuh"
#include"Interactor/SpectralEwaldPoisson.cuh"
#include"Interactor/ExternalForces.cuh"
#include"Integrator/BrownianDynamics.cuh"
#include"DryDiffusion.cuh"
#include"utils/InputFile.h"
#include"Interactor/DoublyPeriodic/DPPoissonSlab.cuh"
#include<fstream>
#include<limits>
#include<random>
using namespace uammd;
using std::make_shared;
using std::endl;

class RepulsiveWall{
  RepulsivePotentialFunctor::PairParameters params;
  real H;
  real imageDistanceMultiplier = 2.0; //Controls the distance of the image ---->   0  |  0, if this parameter is 2, particles interact with images, if 1, image particles are located at the wall height.
public:
  RepulsiveWall(real H, RepulsivePotentialFunctor::PairParameters ip, real imageDistanceMultiplier):
    H(H),params(ip),imageDistanceMultiplier(imageDistanceMultiplier){}

  __device__ ForceEnergyVirial sum(Interactor::Computables comp, real4 pos /*, real mass */){
    real distanceToImage = abs(abs(pos.z) - H * real(0.5))*imageDistanceMultiplier;
    real fz = RepulsivePotentialFunctor::force(distanceToImage * distanceToImage, params) * distanceToImage;
    ForceEnergyVirial fev;
    fev.force = make_real3(0, 0, fz*(pos.z<0?real(-1.0):real(1.0)));
    return fev;
  }

  auto getArrays(ParticleData *pd){
    auto pos = pd->getPos(access::gpu, access::read);
    return std::make_tuple(pos.raw());
  }
};

class ExternalField{
  real3 externalField;
public:
  ExternalField(real3 externalField):externalField(externalField){}

  __device__ ForceEnergyVirial sum(Interactor::Computables comp, real charges){
    real3 externalForce = charges*externalField;
    ForceEnergyVirial result;
    result.force = externalForce;
    return result;
  }

  auto getArrays(ParticleData *pd){
    auto charges = pd->getCharge(access::gpu, access::read); // a number
    return charges.begin();
  }
};


struct Parameters{
  int numberParticles;
  real Lxy, H;
  int Nxy = -1;
  int support = 10;
  real numberStandardDeviations = 4;
  real upsampling = 1.2;
  real tolerance = 1e-4;
  real temperature;
  real permitivity, permitivityBottom, permitivityTop;

  real bottomWallSurfaceValue = 0;

  int numberSteps, printSteps, relaxSteps;
  real dt, viscosity, hydrodynamicRadius, wetHydrodynamicRadius;

  real gw;
  real U0, sigma, r_m, p, cutOff;
  real wall_U0, wall_sigma, wall_r_m, wall_p, wall_cutOff;
  real imageDistanceMultiplier;

  std::string outfile, readFile, forcefile, fieldfile;
  std::string mobilityFile;


  std::string brownianUpdateRule = "EulerMaruyama";
  bool idealParticles=false;
  bool noElectrostatics=false;
  int w = 6;
  real beta = 10.13641758;
  int nxy_stokes;
  int nz_stokes;
  real hxy_stokes;
  
  real3 externalField;
  int fold;
};

struct UAMMD{
  std::shared_ptr<ParticleData> pd;
  std::shared_ptr<thrust::device_vector<real4>> savedPositions;
  Parameters par;
};

Parameters readParameters(std::string fileName);

void initializeParticles(UAMMD sim){
  auto pos = sim.pd->getPos(access::location::cpu, access::mode::write);
  auto charge = sim.pd->getCharge(access::location::cpu, access::mode::write);
  if(sim.par.readFile.empty()){
    std::generate(pos.begin(), pos.end(),
		  [&](){
		    real Lxy = sim.par.Lxy;
		    real H = sim.par.H;
		    real3 p;
		    real pdf;
		    do{
		      p = make_real3(sim.pd->getSystem()->rng().uniform3(-0.5, 0.5))*make_real3(Lxy, Lxy, H-2*sim.par.gw);
		      pdf = 1.0;
		    }while(sim.pd->getSystem()->rng().uniform(0, 1) > pdf);
		    return make_real4(p, 0);
		  });
    fori(0, sim.par.numberParticles){
      charge[i] = ((i%2)-0.5)*2;
    }
  }
  else{
    std::ifstream in(sim.par.readFile);
    fori(0, sim.par.numberParticles){
      in>>pos[i].x>>pos[i].y>>pos[i].z>>charge[i];
      pos[i].w = 0;
    }
  }
  thrust::copy(pos.begin(), pos.end(), sim.savedPositions->begin());
}

UAMMD initialize(int argc, char *argv[]){
  UAMMD sim;
  auto sys = std::make_shared<System>(argc, argv);
  std::random_device r;
  auto now = static_cast<long long unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  sys->rng().setSeed(now);
  std::string datamain = argc>1?argv[1]:"data.main";
  sim.par = readParameters(datamain);
  sim.pd = std::make_shared<ParticleData>(sim.par.numberParticles, sys);
  sim.savedPositions = std::make_shared<thrust::device_vector<real4>>();
  sim.savedPositions->resize(sim.par.numberParticles);
  initializeParticles(sim);
  return sim;
}

auto string2BrownianRule(std::string str) {
  if(str == "EulerMaruyama") return DryWetBD::update_rules::euler_maruyama;
  else if(str =="Leimkuhler")   return DryWetBD::update_rules::leimkuhler;
  else
    throw std::runtime_error("Invalid brownian rule");
}

auto createIntegrator(UAMMD sim){
  using BD = DryWetBD;
  BD::Parameters par;
  par.temperature = sim.par.temperature;
  par.viscosity = sim.par.viscosity;
  par.hydrodynamicRadius = sim.par.hydrodynamicRadius;
  par.dt = sim.par.dt;
  par.wetRadius = sim.par.wetHydrodynamicRadius;
  par.brownianUpdateRule = string2BrownianRule(sim.par.brownianUpdateRule);
  par.dryMobilityFile = sim.par.mobilityFile;
  par.H = sim.par.H;
  par.Lxy = sim.par.Lxy;
  //par.hxy_stokes = sim.par.hxy_stokes
  // par.w = sim.par.w;
  // par.nxy_stokes = sim.par.nxy_stokes;
  // par.nz_stokes = sim.par.nz_stokes;
  return std::make_shared<BD>(sim.pd, par);
}

auto createDoublyPeriodicElectrostaticInteractor(UAMMD sim){
  DPPoissonSlab::Parameters par;
  par.Lxy = make_real2(sim.par.Lxy);
  par.H = sim.par.H;
  DPPoissonSlab::Permitivity perm;
  perm.inside = sim.par.permitivity;
  perm.top = sim.par.permitivityTop;
  if(sim.par.permitivityTop<0){ //Metallic boundary
    perm.top = std::numeric_limits<real>::infinity();
  }
  perm.bottom = sim.par.permitivityBottom;
  if(sim.par.permitivityBottom<0){ //Metallic boundary
    perm.bottom = std::numeric_limits<real>::infinity();
  }
  par.permitivity = perm;
  par.gw = sim.par.gw;
  par.tolerance = sim.par.tolerance;
  if(sim.par.upsampling > 0){
    par.upsampling=sim.par.upsampling;
  }
  if(sim.par.numberStandardDeviations > 0){
    par.numberStandardDeviations=sim.par.numberStandardDeviations;
  }
  if(sim.par.support > 0){
    par.support=sim.par.support;
  }
  if(sim.par.Nxy > 0){
    par.Nxy = sim.par.Nxy;
  }
  par.support = sim.par.support;
  par.numberStandardDeviations = sim.par.numberStandardDeviations;
  auto dppoisson = std::make_shared<DPPoissonSlab>(sim.pd, par);
  if(sim.par.bottomWallSurfaceValue){
    System::log<System::MESSAGE>("[DPPoisson] Setting the bottom wall zero mode Fourier value to %g", sim.par.bottomWallSurfaceValue);
    dppoisson->setSurfaceValuesZeroModeFourier({0, 0, sim.par.bottomWallSurfaceValue, 0});
  }
  return dppoisson;
}

auto createWallRepulsionInteractor(UAMMD sim){
  RepulsivePotentialFunctor::PairParameters potpar;
  potpar.cutOff2 = sim.par.wall_cutOff*sim.par.wall_cutOff;
  potpar.sigma = sim.par.wall_sigma;
  potpar.U0 = sim.par.wall_U0;
  potpar.r_m = sim.par.wall_r_m;
  potpar.p = sim.par.wall_p;
  return make_shared<ExternalForces<RepulsiveWall>>(sim.pd, make_shared<RepulsiveWall>(sim.par.H, potpar, sim.par.imageDistanceMultiplier));
}

auto createExternalFieldInteractor(UAMMD sim){
  return make_shared<ExternalForces<ExternalField>>(sim.pd, make_shared<ExternalField>(sim.par.externalField));
}

auto createPotential(UAMMD sim){
  auto pot = std::make_shared<RepulsivePotential>();
  RepulsivePotential::InputPairParameters ppar;
  ppar.cutOff = sim.par.cutOff;
  ppar.U0 = sim.par.U0;
  ppar.sigma = sim.par.sigma;
  ppar.r_m = sim.par.r_m;
  ppar.p = sim.par.p;
  System::log<System::MESSAGE>("Repulsive rcut: %g", ppar.cutOff);
  pot->setPotParameters(0, 0, ppar);
 return pot;
}

template<class UsePotential> auto createShortRangeInteractor(UAMMD sim){
  auto pot = createPotential(sim);
  using SR = PairForces<UsePotential>;
  typename SR::Parameters params;
  real Lxy = sim.par.Lxy;
  real H = sim.par.H;
  params.box = Box(make_real3(Lxy, Lxy, H));
  params.box.setPeriodicity(1,1,0);
  auto pairForces = std::make_shared<SR>(sim.pd, params, pot);
  return pairForces;
}

void writeSimulation(UAMMD sim, std::vector<real4> fieldAtParticles){
  auto pos = sim.pd->getPos(access::location::cpu, access::mode::read);
  auto charge = sim.pd->getCharge(access::location::cpu, access::mode::read);
  auto force = sim.pd->getForce(access::location::cpu, access::mode::read);
  static std::ofstream out(sim.par.outfile);
  static std::ofstream outf(sim.par.forcefile);
  static std::ofstream outfield(sim.par.fieldfile);
  Box box(make_real3(sim.par.Lxy, sim.par.Lxy, sim.par.H));
  box.setPeriodicity(1,1,0);
  real3 L = box.boxSize;
  out<<"#Lx="<<L.x*0.5<<";Ly="<<L.y*0.5<<";Lz="<<L.z*0.5<<";"<<std::endl;
  if(outf.good())outf<<"#"<<std::endl;
  if(outfield.good())outfield<<"#"<<std::endl;
  fori(0, sim.par.numberParticles){
    real3 p;
    if(sim.par.fold == 1){
      p = box.apply_pbc(make_real3(pos[i]));
    } else {
      p = make_real3(pos[i]);
    }
    // real3 p = make_real3(pos[i]);
    real q = charge[i];
    out<<std::setprecision(2*sizeof(real))<<p<<" "<<q<<"\n";
    if(outf.good()){
      outf<<std::setprecision(2*sizeof(real))<<force[i]<<"\n";
    }
    if(outfield.good() and fieldAtParticles.size()>0){
      outfield<<std::setprecision(2*sizeof(real))<<fieldAtParticles[i]<<"\n";
    }
  }
  out<<std::flush;
}

struct CheckOverlap {
  real H;
  CheckOverlap(real H):H(H){

  }

  __device__ bool operator()(real4 p){
    return abs(p.z) >= (real(0.5)*H);
  }

};

bool checkWallOverlap(UAMMD sim){
  auto pos = sim.pd->getPos(access::location::gpu, access::mode::read);
  //int overlappingCharges = thrust::count_if(thrust::cuda::par, pos.begin(), pos.end(), CheckOverlap(sim.par.H));
  //return overlappingCharges > 0;
  auto overlappingPos = thrust::find_if(thrust::cuda::par, pos.begin(), pos.end(), CheckOverlap(sim.par.H));
  return overlappingPos != pos.end();
}

void restoreLastSavedConfiguration(UAMMD sim) {
  auto pos = sim.pd->getPos(access::location::gpu, access::mode::write);
  thrust::copy(thrust::cuda::par, sim.savedPositions->begin(), sim.savedPositions->end(), pos.begin());
}

void saveConfiguration(UAMMD sim) {
  auto pos = sim.pd->getPos(access::location::gpu, access::mode::read);
  thrust::copy(thrust::cuda::par, pos.begin(), pos.end(), sim.savedPositions->begin());
}

int main(int argc, char *argv[]){
  auto sim = initialize(argc, argv);
  auto bd = createIntegrator(sim);
  std::shared_ptr<DPPoissonSlab> dpslab;
  if(not sim.par.idealParticles){
    if(not sim.par.noElectrostatics){
      dpslab = createDoublyPeriodicElectrostaticInteractor(sim);
      bd->addInteractor(dpslab);
    }
    if(sim.par.U0 > 0){
      bd->addInteractor(createShortRangeInteractor<RepulsivePotential>(sim));
    }
  }
  bd->addInteractor(createWallRepulsionInteractor(sim));
  if(sim.par.externalField.x != 0 or sim.par.externalField.y != 0 or sim.par.externalField.z != 0){
    bd->addInteractor(createExternalFieldInteractor(sim));
  }
  int numberRetries=0;
  int numberRetriesThisStep=0;
  int lastStepSaved=0;
  constexpr int saveRate = 100;
  constexpr int maximumRetries = 1e6;
  constexpr int maximumRetriesPerStep=1e4;
  forj(0, sim.par.relaxSteps){
    bd->forwardTime();
    if(checkWallOverlap(sim)){
      numberRetries++;
      if(numberRetries>maximumRetries){
	throw std::runtime_error("Too many steps with wall overlapping charges detected, aborting run");
      }
      numberRetriesThisStep++;
      if(numberRetriesThisStep>maximumRetriesPerStep){
	throw std::runtime_error("Cannot recover from configuration with wall overlapping charges, aborting run");
      }
      j=lastStepSaved;
      restoreLastSavedConfiguration(sim);
      continue;
    }
    if(j%saveRate==0){
      numberRetriesThisStep = 0;
      lastStepSaved=j;
      saveConfiguration(sim);
    }

  }
  Timer tim;
  tim.tic();
  lastStepSaved=0;
  forj(0, sim.par.numberSteps){
    bd->forwardTime();
    if(checkWallOverlap(sim)){
      numberRetries++;
      if(numberRetries>maximumRetries){
	throw std::runtime_error("Too many steps with wall overlapping charges detected, aborting run");
      }
      numberRetriesThisStep++;
      if(numberRetriesThisStep>maximumRetriesPerStep){
	throw std::runtime_error("Cannot recover from configuration with wall overlapping charges, aborting run");
      }
      j=lastStepSaved;
      restoreLastSavedConfiguration(sim);
      continue;
    }
    if(j%saveRate==0){
      numberRetriesThisStep=0;
      lastStepSaved=j;
      saveConfiguration(sim);
    }
    if(sim.par.printSteps > 0 and j%sim.par.printSteps==0){
      std::vector<real4> fieldAtParticles;
      if(not sim.par.fieldfile.empty() and dpslab){
	// System::log<System::ERROR>("This functionality is not available");
	auto d_field = dpslab->computeFieldAtParticles();
	fieldAtParticles.resize(d_field.size());
	thrust::copy(d_field.begin(), d_field.end(), fieldAtParticles.begin());
      }
      writeSimulation(sim, fieldAtParticles);
      numberRetriesThisStep=0;
      lastStepSaved=j;
      saveConfiguration(sim);
    }
  }
  System::log<System::MESSAGE>("Number of rejected configurations: %d (%g%% of total)", numberRetries, (double)numberRetries/(sim.par.numberSteps + sim.par.relaxSteps)*100.0);
  auto totalTime = tim.toc();
  System::log<System::MESSAGE>("mean FPS: %.2f", sim.par.numberSteps/totalTime);
  return 0;
}

Parameters readParameters(std::string datamain){
  InputFile in(datamain);
  Parameters par;
  in.getOption("fold", InputFile::Required)>>par.fold;
  in.getOption("Lxy", InputFile::Required)>>par.Lxy;
  in.getOption("H", InputFile::Required)>>par.H;
  in.getOption("numberSteps", InputFile::Required)>>par.numberSteps;
  in.getOption("printSteps", InputFile::Required)>>par.printSteps;
  in.getOption("relaxSteps", InputFile::Required)>>par.relaxSteps;
  in.getOption("dt", InputFile::Required)>>par.dt;
  in.getOption("numberParticles", InputFile::Required)>>par.numberParticles;
  in.getOption("temperature", InputFile::Required)>>par.temperature;
  in.getOption("viscosity", InputFile::Required)>>par.viscosity;
  in.getOption("hydrodynamicRadius", InputFile::Required)>>par.hydrodynamicRadius;
  in.getOption("outfile", InputFile::Required)>>par.outfile;
  in.getOption("useMobilityFromFile", InputFile::Optional)>>par.mobilityFile;
  in.getOption("forcefile", InputFile::Optional)>>par.forcefile;
  in.getOption("fieldfile", InputFile::Optional)>>par.fieldfile;
  in.getOption("U0", InputFile::Required)>>par.U0;
  in.getOption("r_m", InputFile::Required)>>par.r_m;
  in.getOption("p", InputFile::Required)>>par.p;
  in.getOption("sigma", InputFile::Required)>>par.sigma;
  in.getOption("readFile", InputFile::Optional)>>par.readFile;
  in.getOption("wetHydrodynamicRadius", InputFile::Required)>>par.wetHydrodynamicRadius;
  in.getOption("gw", InputFile::Required)>>par.gw;
  in.getOption("tolerance", InputFile::Optional)>>par.tolerance;
  in.getOption("permitivity", InputFile::Required)>>par.permitivity;
  in.getOption("permitivityTop", InputFile::Required)>>par.permitivityTop;
  in.getOption("permitivityBottom", InputFile::Required)>>par.permitivityBottom;
  in.getOption("externalField", InputFile::Required)>>par.externalField;

  in.getOption("Nxy", InputFile::Required)>>par.Nxy;

  in.getOption("wall_U0", InputFile::Required)>>par.wall_U0;
  in.getOption("wall_r_m", InputFile::Required)>>par.wall_r_m;
  in.getOption("wall_p", InputFile::Required)>>par.wall_p;
  in.getOption("wall_sigma", InputFile::Required)>>par.wall_sigma;
  in.getOption("imageDistanceMultiplier", InputFile::Required)>>par.imageDistanceMultiplier;
  par.wall_cutOff = par.wall_sigma*pow(2,1.0/par.wall_p);
  par.cutOff = par.sigma*pow(2,1.0/par.p);
  in.getOption("BrownianUpdateRule", InputFile::Optional)>>par.brownianUpdateRule;
  if(in.getOption("idealParticles", InputFile::Optional))
    par.idealParticles = true;
  if(in.getOption("noElectrostatics", InputFile::Optional))
    par.noElectrostatics = true;

  in.getOption("bottomWallSurfaceValue", InputFile::Optional)>>par.bottomWallSurfaceValue;

  // in.getOption("hxy_stokes", InputFile::Required)>>par.hxy_stokes;
  // in.getOption("w", InputFile::Required)>>par.w;
  // in.getOption("beta", InputFile::Required)>>par.beta;
  // in.getOption("nxy_stokes", InputFile::Required)>>par.nxy_stokes;
  // in.getOption("nz_stokes", InputFile::Required)>>par.nz_stokes;

  return par;
}
