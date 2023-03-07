/* Raul P. Pelaez 2023. Tests for the dry wet slit channel Integrator

 */
#include <gtest/gtest.h>
#include "gmock/gmock.h"
#include <random>
#include <fstream>
#include "DryDiffusion.cuh"
#include"Interactor/DoublyPeriodic/DPPoissonSlab.cuh"
using namespace uammd;

//Writes a mobility file with constant mobility accross the domain
void writeDefaultMobilityFile(){
  std::ofstream out("mob.dat");
  out<<"-1.0 1.0 1.0 1.0"<<std::endl;
  out<<"0.0 1.0 1.0 1.0"<<std::endl;
  out<<"1.0 1.0 1.0 1.0"<<std::endl;
}

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

// Compute electrostatic fields for 2 particles
TEST(Playing, ComputeElectrostaticForce){
  using DPP = DPPoissonSlab;
  DPP::Parameters par;
  // needs a real2 for Lxy (see DPPoissonSlab.cuh)
  par.Lxy = {32, 32};
  par.H = 20;
  par.gw = 0.25;
  DPP::Permitivity perm;
  perm.inside = 1.0;
  perm.top = 1.0;
  perm.bottom = 0.05;
  par.permitivity = perm;
  par.Nxy = 72;
  
  auto pd = std::make_shared<ParticleData>(2);
  pd->getPos(access::location::cpu, access::mode::write)[0] = {6.001,4.02,0.05};
  pd->getPos(access::location::cpu, access::mode::write)[1] = {5.11,5.24,-2.1};
  pd->getCharge(access::location::cpu, access::mode::write)[0] = 1;
  pd->getCharge(access::location::cpu, access::mode::write)[1] = -1;

  real c0 = pd->getCharge(access::location::cpu, access::mode::read)[0];
  std::cout << "charge of particle #1 = " << c0 << std::endl;
  real c1 = pd->getCharge(access::location::cpu, access::mode::read)[1];
  std::cout << "charge of particle #2 = " << c1 << std::endl;
  
  auto poisson = std::make_shared<DPPoissonSlab>(pd, par);
  thrust::device_vector<real4> field = poisson->computeFieldAtParticles();
  std::vector<real4> fieldAtParticles;
  fieldAtParticles.resize(field.size());
  thrust::copy(field.begin(), field.end(), fieldAtParticles.begin());
  real E0x = fieldAtParticles[0].x;
  std::cout << "x field at particle #1 (electrostatics) = " << E0x << std::endl;
  real magicalValue = -0.003902211008208;//given by the MATLB code
  EXPECT_THAT(E0x, ::testing::DoubleNear(magicalValue, 1e-5));
}

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
//   // Runs forever if mob.dat is not given
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
