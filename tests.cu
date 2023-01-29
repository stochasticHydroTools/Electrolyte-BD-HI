#include <gtest/gtest.h>
#include "gmock/gmock.h"
#include <random>
#include <fstream>
#include "DryDiffusion.cuh"
using namespace uammd;

//Writes a mobility file with constant mobility accross the domain
void writeDefaultMobilityFile(){
  std::ofstream out("mob.dat");
  out<<"-1.0 1.0 1.0 1.0"<<std::endl;
  out<<"0.0 1.0 1.0 1.0"<<std::endl;
  out<<"1.0 1.0 1.0 1.0"<<std::endl;
}

TEST(DryWetMobility, CanBeCreated){
  using BD = DryWetBD;
  BD::Parameters par;
  par.temperature = 1.0;
  par.viscosity = 1.0;
  par.hydrodynamicRadius = 1.0;
  par.dt = 1.0;
  par.wetRadius = 0.9;
  par.brownianUpdateRule = DryWetBD::update_rules::euler_maruyama;
  //par.dryMobilityFile = sim.par.mobilityFile;
  par.H = 16;
  par.Lxy = 32;
  auto pd = std::make_shared<ParticleData>(1);
  auto bd = std::make_shared<BD>(pd, par);
}

//Adds a force to the first particle, the force defaults to 1,0,0 of not provided
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


TEST(FullDryMobility, SelfMobilityIsCorrect){
  using BD = DryWetBD;
  writeDefaultMobilityFile();
  BD::Parameters par;
  par.temperature = 0;
  par.viscosity = 1.0/(6*M_PI);
  par.hydrodynamicRadius = 1.0;
  par.dt = 1.0;
  par.wetRadius = -1;
  par.brownianUpdateRule = DryWetBD::update_rules::euler_maruyama;
  par.dryMobilityFile = "mob.dat";
  par.H = 32;
  par.Lxy = 64;
  auto pd = std::make_shared<ParticleData>(1);
  pd->getPos(access::cpu, access::write)[0] = real4();
  auto bd = std::make_shared<BD>(pd, par);
  bd->addInteractor(std::make_shared<miniInteractor>(pd));
  bd->forwardTime();
  real M0 = pd->getPos(access::cpu, access::write)[0].x;
  ASSERT_THAT(M0, ::testing::DoubleNear(1, 1e-5));
}



TEST(FullWetMobility, SelfMobilityIsCorrectAtMiddlePlaneForLargeDomain){
  using BD = DryWetBD;
  writeDefaultMobilityFile();
  BD::Parameters par;
  par.temperature = 0;
  par.viscosity = 1.0/(6*M_PI);
  par.hydrodynamicRadius = 1.0;
  par.dt = 1.0;
  par.wetRadius = par.hydrodynamicRadius;
  par.brownianUpdateRule = DryWetBD::update_rules::euler_maruyama;
  par.dryMobilityFile = "mob.dat";
  par.H = 128;
  par.Lxy = 64;
  auto pd = std::make_shared<ParticleData>(1);
  pd->getPos(access::cpu, access::write)[0] = real4();
  auto bd = std::make_shared<BD>(pd, par);
  bd->addInteractor(std::make_shared<miniInteractor>(pd));
  bd->forwardTime();
  real M0 = pd->getPos(access::cpu, access::write)[0].x;
  ASSERT_THAT(M0, ::testing::DoubleNear(1, 1e-2));
}
