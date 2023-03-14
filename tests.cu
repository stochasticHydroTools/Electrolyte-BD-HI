/* Raul P. Pelaez and Aref Hashemi 2023. Tests for the dry wet slit channel Integrator

 */
#include <gtest/gtest.h>
#include "gmock/gmock.h"
#include <random>
#include <fstream>
#include "DryDiffusion.cuh"
#include"Interactor/DoublyPeriodic/DPPoissonSlab.cuh"
using namespace uammd;

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

  real3 externalField;
  int fold;
};


struct UAMMD{
  std::shared_ptr<ParticleData> pd;
  std::shared_ptr<thrust::device_vector<real4>> savedPositions;
  Parameters par;
};

//Adds a force to the first particle, the force defaults to 1,0,0 if not provided
struct miniInteractor: public Interactor{
  real3 f;
public:
  miniInteractor(std::shared_ptr<ParticleData> pd, real3 f = {1,0,0}):
    Interactor(pd), f(f){
  }

  void sum(Computables comp, cudaStream_t st =0) override{
    pd->getForce(access::cpu, access::write)[0] = make_real4(f);
  }
};

// ############## Tests by Aref ############## //

// // Add a force to a particle and recompute it using the Intercator
// TEST(Playing, ReadWriteParticleForce){
//   auto pd = std::make_shared<ParticleData>(2);
//   auto external = std::make_shared<miniInteractor>(pd);
//   external->sum({.force=true});
//   real3 F0 = make_real3(pd->getForce(access::location::cpu, access::mode::read)[0]);
//   std::cout << "force on particle #1 = " << F0 << std::endl;
//   EXPECT_THAT(F0.x, ::testing::DoubleNear(1, 1e-5));
// }

// // Reading and writing particle positions
// TEST(Playing, ReadWriteParticlePosition){
//   auto pd = std::make_shared<ParticleData>(2);
//   pd->getPos(access::location::cpu, access::mode::write)[0] = {1,0,0};
//   pd->getPos(access::location::cpu, access::mode::write)[1] = {1.5,0,0};
//   real x0 = pd->getPos(access::location::cpu, access::mode::read)[0].x;
//   real x1 = pd->getPos(access::location::cpu, access::mode::read)[1].x;
//   std::cout << "particle #1 position = " << x0 << std::endl;
//   std::cout << "particle #2 position = " << x1 << std::endl;
//   EXPECT_THAT(x0, ::testing::DoubleNear(1, 1e-5));
//   EXPECT_THAT(x1, ::testing::DoubleNear(1.5, 1e-5));
// }

// // Compute electrostatic fields for 2 particles
// TEST(FULLDRY, ComputeElectrostaticField){
//   using DPP = DPPoissonSlab;
//   DPP::Parameters par;
//   // needs a real2 for Lxy (see DPPoissonSlab.cuh)
//   par.Lxy = {32, 32};
//   par.H = 20;
//   par.gw = 0.25;
//   DPP::Permitivity perm;
//   perm.inside = 1.0;
//   perm.top = 1.0;
//   perm.bottom = 0.05;
//   par.permitivity = perm;
//   par.Nxy = 72;
  
//   auto pd = std::make_shared<ParticleData>(2);
//   pd->getPos(access::location::cpu, access::mode::write)[0] = {6.001,4.02,0.05};
//   pd->getPos(access::location::cpu, access::mode::write)[1] = {5.11,5.24,-2.1};
//   pd->getCharge(access::location::cpu, access::mode::write)[0] = 1;
//   pd->getCharge(access::location::cpu, access::mode::write)[1] = -1;

//   real c0 = pd->getCharge(access::location::cpu, access::mode::read)[0];
//   std::cout << "charge of particle #1 = " << c0 << std::endl;
//   real c1 = pd->getCharge(access::location::cpu, access::mode::read)[1];
//   std::cout << "charge of particle #2 = " << c1 << std::endl;
  
//   auto poisson = std::make_shared<DPPoissonSlab>(pd, par);
//   thrust::device_vector<real4> field = poisson->computeFieldAtParticles();
//   std::vector<real4> fieldAtParticles;
//   fieldAtParticles.resize(field.size());
//   thrust::copy(field.begin(), field.end(), fieldAtParticles.begin());
//   real E0x = fieldAtParticles[0].x;
//   std::cout << "x field at particle #1 (electrostatics) = " << E0x << std::endl;
//   real magicalValue = -0.003902211008208;//given by the MATLB code
//   EXPECT_THAT(E0x, ::testing::DoubleNear(magicalValue, 1e-5));
// }

// // Test if a simple integration is working properly
// TEST(FULLDRY, SimpleIntegration){
//   using BD = DryWetBD;
//   BD::Parameters par;
//   par.temperature = 0;
//   // \mu = 1/(6\pi\eta a) = 1
//   par.viscosity = 1.0/(6*M_PI);
//   par.hydrodynamicRadius = 1.0;
//   par.dt = 1.0;
//   par.wetRadius = -1;
//   par.brownianUpdateRule = DryWetBD::update_rules::euler_maruyama;
//   par.H = 32;
//   par.Lxy = 64;
//   auto pd = std::make_shared<ParticleData>(1);
//   pd->getPos(access::cpu, access::write)[0] = {0,0,0};
//   auto bd = std::make_shared<BD>(pd, par);
//   bd->addInteractor(std::make_shared<miniInteractor>(pd));
//   bd->forwardTime();
//   // F = [1, 0, 0]; \mu = 1; U = 1; \Delta t = 1. \Rightarrow \Delta x = 1
//   real dx = pd->getPos(access::cpu, access::write)[0].x;
//   EXPECT_THAT(dx, ::testing::DoubleNear(1, 1e-5));
// }

// Test if an integration works for a pair of particles interacting electrostatically
TEST(FULLDRY, Integration){
  UAMMD sim;
  
  sim.par.numberParticles = 2;
  sim.pd = std::make_shared<ParticleData>(sim.par.numberParticles);
  
  sim.pd->getPos(access::location::cpu, access::mode::write)[0] = {6.001,4.02,0.05};
  sim.pd->getPos(access::location::cpu, access::mode::write)[1] = {5.11,5.24,-2.1};
  sim.pd->getCharge(access::location::cpu, access::mode::write)[0] = 1;
  sim.pd->getCharge(access::location::cpu, access::mode::write)[1] = -1;

  sim.par.Lxy = 32;
  sim.par.H = 20;
  sim.par.gw = 0.25;
  sim.par.permitivity = 1;
  sim.par.permitivityTop = 1;
  sim.par.permitivityBottom = 0.05;
  sim.par.Nxy = 72;
  sim.par.temperature = 0;
  sim.par.viscosity = 1.0/(6*M_PI);
  sim.par.hydrodynamicRadius = 1.0;
  sim.par.wetHydrodynamicRadius = -1;
  sim.par.dt = 0.01;
  sim.par.brownianUpdateRule = "EulerMaruyama";

  using DPP = DPPoissonSlab;
  DPP::Parameters par;
  par.Lxy = make_real2(sim.par.Lxy);
  par.H = sim.par.H;
  par.gw = sim.par.gw;
  DPP::Permitivity perm;
  perm.inside = sim.par.permitivity;
  perm.top = sim.par.permitivityTop;
  perm.bottom = sim.par.permitivityBottom;
  par.permitivity = perm;
  par.Nxy = sim.par.Nxy;
  auto poisson = std::make_shared<DPPoissonSlab>(sim.pd, par);
  
  using BD = DryWetBD;
  BD::Parameters parBD;
  parBD.temperature = sim.par.temperature;
  parBD.viscosity = sim.par.viscosity;
  parBD.hydrodynamicRadius = sim.par.hydrodynamicRadius;
  parBD.wetRadius = sim.par.wetHydrodynamicRadius;
  parBD.dt = sim.par.dt;
  parBD.brownianUpdateRule = DryWetBD::update_rules::euler_maruyama;
  parBD.Lxy = sim.par.Lxy;
  parBD.H = sim.par.H;
  auto bd = std::make_shared<BD>(sim.pd, parBD);
  bd->addInteractor(poisson);
  bd->forwardTime();

  real dx = sim.pd->getPos(access::cpu, access::write)[0].x-6.001;
  std::cout << "displacement is " << dx << std::endl;
  real magicalValue = -0.0000390221100821;//given by the MATLB code
  EXPECT_THAT(dx, ::testing::DoubleNear(magicalValue, 1e-5));
}


// // ############## Tests by Raul ############## //
// TEST(DryWetMobility, CanBeCreated){
//   using BD = DryWetBD;
//   BD::Parameters par;
//   par.temperature = 1.0;
//   par.viscosity = 1.0;
//   par.hydrodynamicRadius = 1.0;
//   par.dt = 1.0;
//   par.wetRadius = 0.9;
//   par.brownianUpdateRule = DryWetBD::update_rules::euler_maruyama;
//   //par.dryMobilityFile = sim.par.mobilityFile;
//   par.H = 16;
//   par.Lxy = 32;
//   auto pd = std::make_shared<ParticleData>(1);
//   auto bd = std::make_shared<BD>(pd, par);
// }


// TEST(FullDryMobility, SelfMobilityIsCorrect){
//   using BD = DryWetBD;
//   writeDefaultMobilityFile();
//   BD::Parameters par;
//   par.temperature = 0;
//   par.viscosity = 1.0/(6*M_PI);
//   par.hydrodynamicRadius = 1.0;
//   par.dt = 1.0;
//   par.wetRadius = -1;
//   par.brownianUpdateRule = DryWetBD::update_rules::euler_maruyama;
//   par.dryMobilityFile = "mob.dat";
//   par.H = 32;
//   par.Lxy = 64;
//   auto pd = std::make_shared<ParticleData>(1);
//   pd->getPos(access::cpu, access::write)[0] = real4();
//   auto bd = std::make_shared<BD>(pd, par);
//   bd->addInteractor(std::make_shared<miniInteractor>(pd));
//   bd->forwardTime();
//   real M0 = pd->getPos(access::cpu, access::write)[0].x;
//   ASSERT_THAT(M0, ::testing::DoubleNear(1, 1e-5));
// }

// TEST(FullWetMobility, SelfMobilityIsCorrectAtMiddlePlaneForLargeDomain){
//   using BD = DryWetBD;
//   writeDefaultMobilityFile();
//   BD::Parameters par;
//   par.temperature = 0;
//   par.viscosity = 1.0/(6*M_PI);
//   par.hydrodynamicRadius = 1.0;
//   par.dt = 1.0;
//   par.wetRadius = par.hydrodynamicRadius; //0<wetRadius<=hydrodynamicRadius means full wet
//   par.brownianUpdateRule = DryWetBD::update_rules::euler_maruyama;
//   par.dryMobilityFile = "mob.dat";
//   par.H = 128;
//   par.Lxy = 64;
//   auto pd = std::make_shared<ParticleData>(1);
//   pd->getPos(access::cpu, access::write)[0] = real4();
//   auto bd = std::make_shared<BD>(pd, par);
//   bd->addInteractor(std::make_shared<miniInteractor>(pd));
//   bd->forwardTime();
//   real M0 = pd->getPos(access::cpu, access::write)[0].x;
//   ASSERT_THAT(M0, ::testing::DoubleNear(1, 1e-1));
// }


// //Asserts the correctness of the self mobility for a certain wet radius
// //All other parameters are hardcoded, see the function.
// //The total hydrodynamic radius is 1 (meaning that wetRadius >=1 is full wet)
// void computeSelfMobilityWithWetRadius(real wetRadius){
//   using BD = DryWetBD;
//   writeDefaultMobilityFile();
//   BD::Parameters par;
//   par.temperature = 0;
//   par.viscosity = 1.0/(6*M_PI);
//   par.hydrodynamicRadius = 1.0;
//   par.dt = 1.0;
//   par.wetRadius = wetRadius;
//   par.brownianUpdateRule = DryWetBD::update_rules::euler_maruyama;
//   par.dryMobilityFile = "mob.dat";
//   par.H = 64;
//   par.Lxy = 64;
//   auto pd = std::make_shared<ParticleData>(1);
//   pd->getPos(access::cpu, access::write)[0] = real4();
//   auto bd = std::make_shared<BD>(pd, par);
//   bd->addInteractor(std::make_shared<miniInteractor>(pd));
//   bd->forwardTime();
//   real M0 = pd->getPos(access::cpu, access::write)[0].x;
//   ASSERT_THAT(M0, ::testing::DoubleNear(1, 1e-1))<<"Failed with wet radius "<<wetRadius;
// }

// TEST(DryWetMobility, SelfMobilityIsCorrectForAnyWetRadius){
//   real minWetRadius = 2;
//   real maxWetRadius = 16;
//   int Ntest = 4;
//   fori(0, Ntest){
//     real wetRadius = minWetRadius + i*(maxWetRadius - minWetRadius)/(Ntest-1);
//     computeSelfMobilityWithWetRadius(wetRadius);
//   }
// }
