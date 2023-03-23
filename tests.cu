/* Raul P. Pelaez and Aref Hashemi 2023. Tests for the dry wet slit channel Integrator

 */
#include <gtest/gtest.h>
#include "gmock/gmock.h"
#include <random>
#include <fstream>
#include "DryDiffusion.cuh"
#include"Interactor/DoublyPeriodic/DPPoissonSlab.cuh"
using namespace uammd;

using scalar = double;

// A normalized measure for accuracy
scalar tolerance(scalar correctVal, scalar mult){
  scalar measure;
  if (abs(correctVal) > 1e-15){
    measure = mult*abs(correctVal);
  } else {
    measure = 1e-15;
  }

  return measure;
}

//Writes a mobility file with constant mobility accross the domain
void writeDefaultMobilityFile(){
  std::ofstream out("uniformMob.dat");
  out<<"-1.0 1.0 1.0 1.0"<<std::endl;
  out<<"0.0 1.0 1.0 1.0"<<std::endl;
  out<<"1.0 1.0 1.0 1.0"<<std::endl;
}

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

// Add a force to a particle and recompute it using the Intercator
TEST(Playing, ReadWriteParticleForce){
  auto pd = std::make_shared<ParticleData>(2);
  auto external = std::make_shared<miniInteractor>(pd);
  external->sum({.force=true});
  real3 F0 = make_real3(pd->getForce(access::location::cpu, access::mode::read)[0]);
  std::cout << "force on particle #1 = " << F0 << std::endl;
  EXPECT_THAT(F0.x, ::testing::DoubleNear(1, 1e-5));
}

// Reading and writing particle positions
TEST(Playing, ReadWriteParticlePosition){
  auto pd = std::make_shared<ParticleData>(2);
  pd->getPos(access::location::cpu, access::mode::write)[0] = {1,0,0};
  pd->getPos(access::location::cpu, access::mode::write)[1] = {1.5,0,0};
  real x0 = pd->getPos(access::location::cpu, access::mode::read)[0].x;
  real x1 = pd->getPos(access::location::cpu, access::mode::read)[1].x;
  std::cout << "particle #1 position = " << x0 << std::endl;
  std::cout << "particle #2 position = " << x1 << std::endl;
  EXPECT_THAT(x0, ::testing::DoubleNear(1, 1e-5));
  EXPECT_THAT(x1, ::testing::DoubleNear(1.5, 1e-5));
}

// full dry mode: Compute electrostatic fields for 2 particles
TEST(FULLDRY, ComputeElectrostaticField){
  using DPP = DPPoissonSlab;
  DPP::Parameters par;
  // needs a real2 for Lxy (see DPPoissonSlab.cuh)
  par.Lxy = {76.8, 76.8};
  par.H = 19.2;
  par.gw = 0.25;
  DPP::Permitivity perm;
  perm.inside = 1.0;
  perm.top = 1.0;
  perm.bottom = 0.05;
  par.permitivity = perm;
  par.Nxy = 72;
  
  auto pd = std::make_shared<ParticleData>(2);
  pd->getPos(access::location::cpu, access::mode::write)[0] = {7,6,-4.8};
  pd->getPos(access::location::cpu, access::mode::write)[1] = {11,6,-4.8};
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
  real E0y = fieldAtParticles[0].y;
  real E0z = fieldAtParticles[0].z;
  std::cout << "x field at particle #1 (electrostatics) = " << E0x << std::endl;
  std::cout << "y field at particle #1 (electrostatics) = " << E0y << std::endl;
  std::cout << "z field at particle #1 (electrostatics) = " << E0z << std::endl;
  real magicalValuex = 0.005223629199127;//given by the MATLB code
  real magicalValuey = 0                ;
  real magicalValuez = 0.000166699457402;
  EXPECT_THAT(E0x, ::testing::DoubleNear(magicalValuex, tolerance(magicalValuex,1e-5)));
  EXPECT_THAT(E0y, ::testing::DoubleNear(magicalValuey, tolerance(magicalValuey,1e-5)));
  EXPECT_THAT(E0z, ::testing::DoubleNear(magicalValuez, tolerance(magicalValuez,1e-5)));
}

// full dry mode: Test if a simple integration is working properly
TEST(FULLDRY, SimpleIntegration){
  using BD = DryWetBD;
  BD::Parameters par;
  writeDefaultMobilityFile();
  par.dryMobilityFile = "uniformMob.dat";
  par.temperature = 0;
  // \mu = 1/(6\pi\eta a) = 1
  par.viscosity = 1.0/(6*M_PI);
  par.hydrodynamicRadius = 1.0;
  par.dt = 1.0;
  par.wetRadius = -1;
  par.brownianUpdateRule = DryWetBD::update_rules::euler_maruyama;
  par.H = 19.2;
  par.Lxy = 76.8;
  auto pd = std::make_shared<ParticleData>(1);
  pd->getPos(access::cpu, access::write)[0] = {0,0,0};
  auto bd = std::make_shared<BD>(pd, par);
  bd->addInteractor(std::make_shared<miniInteractor>(pd));
  bd->forwardTime();
  // F = [1, 0, 0]; \mu = 1; U = 1; \Delta t = 1. \Rightarrow \Delta x = 1
  real dx = pd->getPos(access::cpu, access::write)[0].x;
  EXPECT_THAT(dx, ::testing::DoubleNear(1, 1e-5));
}

// full dry mode: Test if an integration works for a pair of particles interacting electrostatically
// self mobility given by Einstein equation
TEST(FULLDRY, IntegrationFlatMobility){
  UAMMD sim;
  
  sim.par.numberParticles = 2;
  sim.pd = std::make_shared<ParticleData>(sim.par.numberParticles);
  
  sim.pd->getPos(access::location::cpu, access::mode::write)[0] = {7,6,-4.8};
  sim.pd->getPos(access::location::cpu, access::mode::write)[1] = {11,6,-4.8};
  sim.pd->getCharge(access::location::cpu, access::mode::write)[0] = 1;
  sim.pd->getCharge(access::location::cpu, access::mode::write)[1] = -1;

  // {
  //   auto pos = sim.pd->getPos(access::cpu, access::write);
  //   pos[0] = {1,2,3};
  //   pos[1] = {2,3,4};
  // }
  
  sim.par.Lxy = 76.8;
  sim.par.H = 19.2;
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
  thrust::device_vector<real4> field = poisson->computeFieldAtParticles();
  std::vector<real4> fieldAtParticles;
  fieldAtParticles.resize(field.size());
  thrust::copy(field.begin(), field.end(), fieldAtParticles.begin());
  scalar F0x = fieldAtParticles[0].x;
  scalar F0y = fieldAtParticles[0].y;
  scalar F0z = fieldAtParticles[0].z;
  
  using BD = DryWetBD;
  BD::Parameters parBD;
  writeDefaultMobilityFile();
  parBD.dryMobilityFile = "uniformMob.dat";
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

  scalar dx = sim.pd->getPos(access::cpu, access::write)[0].x-7;
  scalar dy = sim.pd->getPos(access::cpu, access::write)[0].y-6;
  scalar dz = sim.pd->getPos(access::cpu, access::write)[0].z+4.8;
  std::cout << "x displacement is " << dx << std::endl;
  std::cout << "y displacement is " << dy << std::endl;
  std::cout << "z displacement is " << dz << std::endl;
  scalar magicalValuedx = F0x*1*sim.par.dt;//0.0000522362919913 given by the MATLB code  
  scalar magicalValuedy = F0y*1*sim.par.dt;//0
  scalar magicalValuedz = F0z*1*sim.par.dt;//0.0000016669945740;
  EXPECT_THAT(dx, ::testing::DoubleNear(magicalValuedx, tolerance(magicalValuedx,1e-5)));
  EXPECT_THAT(dy, ::testing::DoubleNear(magicalValuedy, tolerance(magicalValuedy,1e-5)));
  EXPECT_THAT(dz, ::testing::DoubleNear(magicalValuedz, tolerance(magicalValuedz,1e-5)));
}


// test is the self mobility is computed correctly
TEST(FULLDRY, selfMobility){
  auto pd = std::make_shared<ParticleData>(1);//dummy instance
  using BD = DryWetBD;
  BD::Parameters parBD;
  parBD.viscosity = 1.0/(6*M_PI);
  parBD.hydrodynamicRadius = 1;
  parBD.Lxy = 76.8;
  parBD.H = 19.2;
  auto bd = std::make_shared<BD>(pd, parBD);

  real z = 4;
  real4 mobilityData = computeSelfMobility(parBD, z);
  scalar muxx = mobilityData.y;
  scalar muyy = mobilityData.z;
  scalar muzz = mobilityData.w;
  std::cout << muxx << std::endl;
  std::cout << muyy << std::endl;
  std::cout << muzz << std::endl;
  scalar magicalValuemuxx = 0.862495177070900;// computed by the DPStokes python code at z = 4R_h above the bottom wall
  scalar magicalValuemuyy = 0.862458866786545;
  scalar magicalValuemuzz = 0.725839526229897;
  std::cout << "expected xx mobility is " << magicalValuemuxx << std::endl;
  std::cout << "expected yy mobility is " << magicalValuemuyy << std::endl;
  std::cout << "expected zz mobility is " << magicalValuemuzz << std::endl;
  EXPECT_THAT(muxx, ::testing::DoubleNear(magicalValuemuxx, tolerance(magicalValuemuxx,1e-2)));
  EXPECT_THAT(muyy, ::testing::DoubleNear(magicalValuemuyy, tolerance(magicalValuemuyy,1e-2)));
  EXPECT_THAT(muzz, ::testing::DoubleNear(magicalValuemuzz, tolerance(magicalValuemuzz,1e-2)));
}


// full dry mode: Test if an integration works for a pair of particles interacting electrostatically
// self mobility precomputed by the Stokes solver
TEST(FULLDRY, Integration){
  UAMMD sim;
  
  sim.par.numberParticles = 2;
  sim.pd = std::make_shared<ParticleData>(sim.par.numberParticles);
  
  sim.par.Lxy = 76.8;
  sim.par.H = 19.2;
  scalar yp0  = 0.5*sim.par.Lxy;
  scalar yp1  = yp0;
  scalar zp0 = -0.5*sim.par.H+4;
  scalar zp1 = -0.5*sim.par.H+7;
  scalar xp0 = 0.5*sim.par.Lxy+2;
  scalar xp1 = 0.5*sim.par.Lxy-2;
  
  sim.pd->getPos(access::location::cpu, access::mode::write)[0] = {xp0,yp0,zp0};
  sim.pd->getPos(access::location::cpu, access::mode::write)[1] = {xp1,yp1,zp1};
  sim.pd->getCharge(access::location::cpu, access::mode::write)[0] = 1;
  sim.pd->getCharge(access::location::cpu, access::mode::write)[1] = -1;

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
  thrust::device_vector<real4> field = poisson->computeFieldAtParticles();
  std::vector<real4> fieldAtParticles;
  fieldAtParticles.resize(field.size());
  thrust::copy(field.begin(), field.end(), fieldAtParticles.begin());
  scalar F0x = fieldAtParticles[0].x;
  scalar F0y = fieldAtParticles[0].y;
  scalar F0z = fieldAtParticles[0].z;

  
  using BD = DryWetBD;
  BD::Parameters parBD;
  // writeDefaultMobilityFile();
  // parBD.dryMobilityFile = "uniformMob.dat";
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
  
  scalar dx = sim.pd->getPos(access::cpu, access::write)[0].x-xp0;
  scalar dy = sim.pd->getPos(access::cpu, access::write)[0].y-yp0;
  scalar dz = sim.pd->getPos(access::cpu, access::write)[0].z-zp0;
  std::cout << "x displacement is " << dx << std::endl;
  std::cout << "y displacement is " << dy << std::endl;
  std::cout << "z displacement is " << dz << std::endl;
  scalar mu_xx = 0.862495177070900;//Computed by DPStokes solver
  scalar mu_yy = 0.862458866786545;
  scalar mu_zz = 0.725839526229897;
  scalar magicalValuedx = mu_xx*F0x*sim.par.dt;
  scalar magicalValuedy = mu_yy*F0y*sim.par.dt;
  scalar magicalValuedz = mu_zz*F0z*sim.par.dt;
  std::cout << "expected x displacement is " << magicalValuedx << std::endl;
  std::cout << "expected y displacement is " << magicalValuedy << std::endl;
  std::cout << "expected z displacement is " << magicalValuedz << std::endl;
  EXPECT_THAT(dx, ::testing::DoubleNear(magicalValuedx, tolerance(magicalValuedx,1e-2)));
  EXPECT_THAT(dy, ::testing::DoubleNear(magicalValuedy, tolerance(magicalValuedy,1e-2)));
  EXPECT_THAT(dz, ::testing::DoubleNear(magicalValuedz, tolerance(magicalValuedz,1e-2)));
}


// full wet mode: Test if an integration works for a pair of particles interacting electrostatically
TEST(FULLWET, Integration){
  UAMMD sim;
  
  sim.par.numberParticles = 2;
  sim.pd = std::make_shared<ParticleData>(sim.par.numberParticles);
  sim.par.Lxy = 76.8;
  sim.par.H = 19.2;
  scalar yp0  = 0.5*sim.par.Lxy;
  scalar yp1  = yp0;
  scalar zp0 = -0.5*sim.par.H+4;
  scalar zp1 = -0.5*sim.par.H+7;
  scalar xp0 = 0.5*sim.par.Lxy+2;
  scalar xp1 = 0.5*sim.par.Lxy-2;
  
  sim.pd->getPos(access::location::cpu, access::mode::write)[0] = {xp0,yp0,zp0};
  sim.pd->getPos(access::location::cpu, access::mode::write)[1] = {xp1,yp1,zp1};
  sim.pd->getCharge(access::location::cpu, access::mode::write)[0] = 1;
  sim.pd->getCharge(access::location::cpu, access::mode::write)[1] = -1;

  real c0 = sim.pd->getCharge(access::location::cpu, access::mode::read)[0];
  real c1 = sim.pd->getCharge(access::location::cpu, access::mode::read)[1];
   
  sim.par.gw = 0.25;
  sim.par.permitivity = 1;
  sim.par.permitivityTop = 1;
  sim.par.permitivityBottom = 0.05;
  sim.par.Nxy = 135;
  sim.par.temperature = 0;
  sim.par.viscosity = 1.0/(6*M_PI);
  sim.par.hydrodynamicRadius = 1.0;
  sim.par.wetHydrodynamicRadius = 1.0;
  sim.par.dt = 1;
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
  thrust::device_vector<real4> field = poisson->computeFieldAtParticles();
  std::vector<real4> fieldAtParticles;
  fieldAtParticles.resize(field.size());
  thrust::copy(field.begin(), field.end(), fieldAtParticles.begin());
  scalar F0x = c0*fieldAtParticles[0].x;
  scalar F0y = c0*fieldAtParticles[0].y;
  scalar F0z = c0*fieldAtParticles[0].z;
  scalar F1x = c1*fieldAtParticles[1].x;
  scalar F1y = c1*fieldAtParticles[1].y;
  scalar F1z = c1*fieldAtParticles[1].z;
  std::cout << "F0x = " << F0x << std::endl;
  std::cout << "F0y = " << F0y << std::endl;
  std::cout << "F0z = " << F0z << std::endl;
  std::cout << "F1x = " << F1x << std::endl;
  std::cout << "F1y = " << F1y << std::endl;
  std::cout << "F1z = " << F1z << std::endl;
  
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

  scalar dx = sim.pd->getPos(access::cpu, access::write)[1].x-xp1;
  scalar dy = sim.pd->getPos(access::cpu, access::write)[1].y-yp1;
  scalar dz = sim.pd->getPos(access::cpu, access::write)[1].z-zp1;
  std::cout << "x displacement is " << dx << std::endl;
  std::cout << "y displacement is " << dy << std::endl;
  std::cout << "z displacement is " << dz << std::endl;
  
  //given by the DPStokes solver code (pair_mobility.py)
  scalar nu_xx = 0.156615987796052 ;
  scalar nu_xy = 0                 ;
  scalar nu_xz = -0.036746481747739;
  scalar nu_yx = 0                 ;
  scalar nu_yy = 0.066428530938183 ;
  scalar nu_yz = 0                 ;
  scalar nu_zx = -0.085289117508094;
  scalar nu_zy = 0                 ;
  scalar nu_zz = 0.041363902842424 ;
  scalar mu_xx = 0.924166257588772 ;
  scalar mu_yy = 0.924130139892737 ;
  scalar mu_zz = 0.840468139888785 ;
  scalar magicalValuedx = (nu_xx*F0x+nu_xy*F0y+nu_xz*F0z+mu_xx*F1x)*sim.par.dt;
  scalar magicalValuedy = (nu_yx*F0x+nu_yy*F0y+nu_yz*F0z+mu_yy*F1y)*sim.par.dt;
  scalar magicalValuedz = (nu_zx*F0x+nu_zy*F0y+nu_zz*F0z+mu_zz*F1z)*sim.par.dt;
  std::cout << "expected x displacement is " << magicalValuedx << std::endl;
  std::cout << "expected y displacement is " << magicalValuedy << std::endl;
  std::cout << "expected z displacement is " << magicalValuedz << std::endl;
  
  EXPECT_THAT(dx, ::testing::DoubleNear(magicalValuedx, tolerance(magicalValuedx,1e-2)));
  EXPECT_THAT(dy, ::testing::DoubleNear(magicalValuedy, tolerance(magicalValuedy,1e-2)));
  EXPECT_THAT(dz, ::testing::DoubleNear(magicalValuedz, tolerance(magicalValuedz,1e-2)));
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
//   par.dryMobilityFile = "uniformMob.dat";
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
//   par.dryMobilityFile = "uniformMob.dat";
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
//   par.dryMobilityFile = "uniformMob.dat";
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
