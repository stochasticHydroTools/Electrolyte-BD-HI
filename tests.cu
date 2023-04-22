/* Raul P. Pelaez and Aref Hashemi 2023. Tests for the dry wet slit channel Integrator

 */
#include <gtest/gtest.h>
#include "gmock/gmock.h"
#include <random>
#include <fstream>
#include "DryDiffusion.cuh"
#include "Interactor/DoublyPeriodic/DPPoissonSlab.cuh"
#include "uammd/src/utils/complex.cuh"
#include <thrust/device_vector.h>
#include <random>
#include "/usr/include/eigen3/Eigen/Dense"
using namespace uammd;

using scalar = double;

// A normalized measure for accuracy
bool tolerance(scalar val, scalar expectedval, scalar numdigits = 1.9){
  // std::cout << val << " " << expectedval << std::endl;
  if (val == expectedval){
    return true;
  } else {
    if (abs(expectedval) > 1e-15){
      scalar numCorrectDigits = log10(abs(expectedval/(val-expectedval)));
      if (numCorrectDigits > numdigits){
	return true;
      } else {
	return false;
      }
    } else {
      if (-log10(abs(val-expectedval))>14){
	return true;
      } else {
	return false;
      }
    }
  }

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
  real dt, viscosity, hydrodynamicRadius, wetHydrodynamicRadius, wetFraction;

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

// Printing average velocity velocity (in the xy plane) as a function of z
TEST(FLUIDVELOCITY,CanPrint){
  DPStokesSlab_ns::DPStokes::Parameters par;
  par.viscosity = 1.0/(6*M_PI);
  par.hydrodynamicRadius = 1;
  par.Lx = 76.8;
  par.H = 19.2;
  real hxy_stokes = 0.64;

  auto dppar = getDPStokesParamtersOnlyForce(par.Lx, par.H, par.viscosity, par.hydrodynamicRadius, hxy_stokes);
  auto dpstokes = std::make_shared<DPStokesSlab_ns::DPStokes>(dppar);

  real z = 0;
  auto pos = thrust::make_constant_iterator<real4>({0,0,z,0});
  auto force = thrust::make_constant_iterator<real4>({1,0,0,0});
  auto MdotOut = dpstokes->Mdot(pos, force, 1);
  real muxx = 6*M_PI*par.viscosity*par.hydrodynamicRadius*real3(MdotOut[0]).x;
  std::cout << "mobility = " << muxx << std::endl;

  std::vector<double> averageVelocity = dpstokes->computeAverageVelocity(pos, force, 1, 0);// 0 denotes x direction
  for (int i=0;i<averageVelocity.size();i++){
    std::cout << averageVelocity[i] << std::endl;
  }
  // a dummy check!
  double a = 1;
  EXPECT_THAT(a, ::testing::DoubleNear(1, 1e-5));
}

// validation of the average velocity versus continuum
// Many particles are randomly placed in the simulation box with an ad-hoc force (in the x direction) on each of them that correspond to a force density; uammd solutions will be compared to the analytical solution to the 1D Stokes equation. Data will be written in 'fluidVelocity.dat' as z, uammd vel, analytical vel.
TEST(FLUIDVELOCITY,Validation){
  DPStokesSlab_ns::DPStokes::Parameters par;
  par.viscosity = 1.0/(6*M_PI);
  par.hydrodynamicRadius = 1;
  real hxy_stokes = 0.64*par.hydrodynamicRadius;
  par.Lx = 40;
  par.H = 10;

  auto dppar = getDPStokesParamtersOnlyForce(par.Lx, par.H, par.viscosity, par.hydrodynamicRadius, hxy_stokes);
  auto dpstokes = std::make_shared<DPStokesSlab_ns::DPStokes>(dppar);

  real forceDensity = 1;//f_x in the Stokes equarion: \mu\frac{d^2u}{dz^2}=-f_x
  real packFactor = 10;
  int numberParticles = (int)(pow(par.Lx,2)*par.H/pow(hxy_stokes,3)*packFactor);
  std::cout << "number of farticles: " << numberParticles << std::endl;
  real forceParticle = forceDensity*forceDensity*pow(par.Lx,2)*par.H/((real)numberParticles);
  real loz = -0.5*par.H, rangez = par.H;
  real loxy = 0, rangexy = par.Lx;
  real x, y, z;
  srand((unsigned) time(NULL));
  auto pd = std::make_shared<ParticleData>(numberParticles);
  std::ofstream randPosOut("randomPos.dat");
  for (int i=0;i<numberParticles;i++){
    z = loz+(real)rand()/(real)RAND_MAX*rangez;
    x = loxy+(real)rand()/(real)RAND_MAX*rangexy;
    y = loxy+(real)rand()/(real)RAND_MAX*rangexy;
    randPosOut << x << " " << y << " " << z << std::endl;
    pd->getPos(access::gpu, access::write)[i] = {x,y,z,0};
    pd->getForce(access::gpu, access::write)[i] = {forceParticle,0,0,0};
  }
  randPosOut.close();

  auto pos = pd->getPos(access::gpu, access::read);
  auto force = pd->getForce(access::gpu, access::read);
  
  std::vector<double> averageVelocity = dpstokes->computeAverageVelocity(pos, force, numberParticles, 0);
  std::ofstream out("fluidVelocity.dat");
  // Note Chebyshev points are (b+a)/2+(b-a)/2*cos(j*pi/(nz-1)) for j = 0,1,...,nz-1
  real a = -0.5*par.H, b = 0.5*par.H;
  int nz = averageVelocity.size();
  std::vector<double> averageVelocityContinuum(nz);// solution from solving the continuum Stokes sequation with a force density: u = fx*H^2/(8*mu)*(1-4(z/H)^2)
  real zpos;
  real err = 0;
  for (int i=0;i<nz;i++){
    zpos = (b+a)/2+(b-a)/2*cos(i*M_PI/(nz-1));
    averageVelocityContinuum[i] = forceDensity*pow(par.H,2)/(8*par.viscosity)*(1-4*pow(zpos/par.H,2));
    out << zpos/b << " " << averageVelocity[i] << " " << averageVelocityContinuum[i] << std::endl;
    err += averageVelocityContinuum[i]!=0?pow(abs((averageVelocityContinuum[i]-averageVelocity[i])/averageVelocityContinuum[i]),2):abs(averageVelocityContinuum[i]-averageVelocity[i]);
  }
  out.close();
  // One can also plot the data to observe how the two solutions compare.
  std::cout << "mean squared error = " << err/nz << std::endl; 
  real tol = 0.01;
  EXPECT_THAT(err/nz, ::testing::Le(tol));
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

// Full Dry Mode: Compute electrostatic fields for 2 particles
// The computed electric fields at particles will be compared to that obtained from the MATLAB code.
TEST(FULLDRY, ComputeElectrostaticField){
  using DPP = DPPoissonSlab;
  DPP::Parameters par;
  // needs a real2 for Lxy (see DPPoissonSlab.cuh)
  par.Lxy = {76.8, 76.8};
  par.H = 19.2;
  par.gw = 0.25;
  DPP::Permitivity perm;
  perm.inside = 1.0;
  perm.top = 0.05;
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
  real expectedE0x = 0.005240820721856;//given by the MATLB code
  real expectedE0y = 0                ;
  real expectedE0z = 0.000164511031114;
  EXPECT_THAT(tolerance(E0x, expectedE0x, 4), ::testing::IsTrue);
  EXPECT_THAT(tolerance(E0y, expectedE0y, 4), ::testing::IsTrue);
  EXPECT_THAT(tolerance(E0z, expectedE0z, 4), ::testing::IsTrue);
}

// Full Dry Mode: Test if a simple integration is working properly
// An external force [1 0 0] is applied on a particle and we compute its displacement when mobility is 1.
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
  par.wetFraction = 0;
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

// Full Dry Mode: Test if an integration works for a pair of particles interacting electrostatically
// Here we assume that the self mobility is given by the Einstein equation.
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
  sim.par.permitivityTop = 0.05;
  sim.par.permitivityBottom = 0.05;
  sim.par.Nxy = 72;
  sim.par.temperature = 0;
  sim.par.viscosity = 1.0/(6*M_PI);
  sim.par.hydrodynamicRadius = 1.0;
  sim.par.wetFraction = 0;
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
  parBD.wetFraction = sim.par.wetFraction;
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
  scalar expecteddx = F0x*1*sim.par.dt;//0.0000524082072186 given by the MATLB code  
  scalar expecteddy = F0y*1*sim.par.dt;//0
  scalar expecteddz = F0z*1*sim.par.dt;//0.0000016451103111;
  EXPECT_THAT(tolerance(dx, expecteddx, 4), ::testing::IsTrue);
  EXPECT_THAT(tolerance(dy, expecteddy, 4), ::testing::IsTrue);
  EXPECT_THAT(tolerance(dz, expecteddz, 4), ::testing::IsTrue);
}


// Test if the self mobility is computed correctly
TEST(FULLDRY, selfMobility){
  auto pd = std::make_shared<ParticleData>(1);//dummy instance
  using BD = DryWetBD;
  BD::Parameters parBD;
  parBD.viscosity = 1.0/(6*M_PI);
  parBD.hydrodynamicRadius = 1;
  parBD.Lxy = 76.8;
  parBD.H = 19.2;
  parBD.hxy_stokes = 0.64;
  parBD.wetFraction = 0;
  auto bd = std::make_shared<BD>(pd, parBD);

  real z = 4;
  real4 mobilityData = computeSelfMobility(parBD, z);
  scalar muxx = mobilityData.y;
  scalar muyy = mobilityData.z;
  scalar muzz = mobilityData.w;
  std::cout << muxx << std::endl;
  std::cout << muyy << std::endl;
  std::cout << muzz << std::endl;
  scalar expectedmuxx = 0.849327755959755;// computed by the DPStokes python code at z = 4R_h above the bottom wall
  scalar expectedmuyy = 0.849291424072988;
  scalar expectedmuzz = 0.724165665045591;
  std::cout << "expected xx mobility is " << expectedmuxx << std::endl;
  std::cout << "expected yy mobility is " << expectedmuyy << std::endl;
  std::cout << "expected zz mobility is " << expectedmuzz << std::endl;
  EXPECT_THAT(tolerance(muxx, expectedmuxx), ::testing::IsTrue);
  EXPECT_THAT(tolerance(muyy, expectedmuyy), ::testing::IsTrue);
  EXPECT_THAT(tolerance(muzz, expectedmuzz), ::testing::IsTrue);
}


// Full Dry Mode: Test if an integration works for a pair of particles interacting electrostatically
// Here self mobility is precomputed by the Stokes solver.
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
  sim.par.permitivityTop = 0.05;
  sim.par.permitivityBottom = 0.05;
  sim.par.Nxy = 72;
  sim.par.temperature = 0;
  sim.par.viscosity = 1.0/(6*M_PI);
  sim.par.hydrodynamicRadius = 1.0;
  sim.par.wetFraction = 0;
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
  parBD.wetFraction = sim.par.wetFraction;
  parBD.dt = sim.par.dt;
  parBD.brownianUpdateRule = DryWetBD::update_rules::euler_maruyama;
  parBD.Lxy = sim.par.Lxy;
  parBD.H = sim.par.H;
  parBD.hxy_stokes = 0.64;
  auto bd = std::make_shared<BD>(sim.pd, parBD);
  bd->addInteractor(poisson);
  bd->forwardTime();
  
  scalar dx = sim.pd->getPos(access::cpu, access::write)[0].x-xp0;
  scalar dy = sim.pd->getPos(access::cpu, access::write)[0].y-yp0;
  scalar dz = sim.pd->getPos(access::cpu, access::write)[0].z-zp0;
  std::ofstream out("FULLDRYdisplacements.dat");// Will be needed for the next test (periodicity)
  out << dx << " " << dy << " " << dz << std::endl;
  std::cout << "x displacement is " << dx << std::endl;
  std::cout << "y displacement is " << dy << std::endl;
  std::cout << "z displacement is " << dz << std::endl;
  scalar mu_xx = 0.849327755959755;// computed by the DPStokes python code at z = 4R_h above the bottom wall
  scalar mu_yy = 0.849291424072988;
  scalar mu_zz = 0.724165665045591;
  scalar expecteddx = mu_xx*F0x*sim.par.dt;
  scalar expecteddy = mu_yy*F0y*sim.par.dt;
  scalar expecteddz = mu_zz*F0z*sim.par.dt;
  std::cout << "expected x displacement is " << expecteddx << std::endl;
  std::cout << "expected y displacement is " << expecteddy << std::endl;
  std::cout << "expected z displacement is " << expecteddz << std::endl;
  EXPECT_THAT(tolerance(dx, expecteddx), ::testing::IsTrue);
  EXPECT_THAT(tolerance(dy, expecteddy), ::testing::IsTrue);
  EXPECT_THAT(tolerance(dz, expecteddz), ::testing::IsTrue);
}


// Full Dry Mode: periodicity
// This test follows the previous test.
TEST(FULLDRY, Periodicity){
  UAMMD sim;
  
  sim.par.numberParticles = 2;
  sim.pd = std::make_shared<ParticleData>(sim.par.numberParticles);
  
  sim.par.Lxy = 76.8;
  sim.par.H = 19.2;
  scalar yp0  = 0.5*sim.par.Lxy;
  scalar yp1  = yp0;
  scalar zp0 = -0.5*sim.par.H+4;
  scalar zp1 = -0.5*sim.par.H+7;
  scalar xp0 = 1*sim.par.Lxy+2;// Moved to the next block
  scalar xp1 = 1*sim.par.Lxy-2;
  
  sim.pd->getPos(access::location::cpu, access::mode::write)[0] = {xp0,yp0,zp0};
  sim.pd->getPos(access::location::cpu, access::mode::write)[1] = {xp1,yp1,zp1};
  sim.pd->getCharge(access::location::cpu, access::mode::write)[0] = 1;
  sim.pd->getCharge(access::location::cpu, access::mode::write)[1] = -1;

  sim.par.gw = 0.25;
  sim.par.permitivity = 1;
  sim.par.permitivityTop = 0.05;
  sim.par.permitivityBottom = 0.05;
  sim.par.Nxy = 72;
  sim.par.temperature = 0;
  sim.par.viscosity = 1.0/(6*M_PI);
  sim.par.hydrodynamicRadius = 1.0;
  sim.par.wetFraction = 0;
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
  parBD.wetFraction = sim.par.wetFraction;
  parBD.dt = sim.par.dt;
  parBD.brownianUpdateRule = DryWetBD::update_rules::euler_maruyama;
  parBD.Lxy = sim.par.Lxy;
  parBD.H = sim.par.H;
  parBD.hxy_stokes = 0.64;
  auto bd = std::make_shared<BD>(sim.pd, parBD);
  bd->addInteractor(poisson);
  bd->forwardTime();
  
  scalar dx = sim.pd->getPos(access::cpu, access::write)[0].x-xp0;
  scalar dy = sim.pd->getPos(access::cpu, access::write)[0].y-yp0;
  scalar dz = sim.pd->getPos(access::cpu, access::write)[0].z-zp0;
  std::ifstream in("FULLDRYdisplacements.dat");// Should be available from the previous test
  double s;
  std::vector<double> expectedvals;
  while(in >> s){
    expectedvals.push_back(s);
  }
  
  std::cout << "x displacement is " << dx << std::endl;
  std::cout << "y displacement is " << dy << std::endl;
  std::cout << "z displacement is " << dz << std::endl;
  
  EXPECT_THAT(tolerance(dx, expectedvals[0]), ::testing::IsTrue);
  EXPECT_THAT(tolerance(dy, expectedvals[1]), ::testing::IsTrue);
  EXPECT_THAT(tolerance(dz, expectedvals[2]), ::testing::IsTrue);
}

// Full Wet Mode: Test if an integration works for a pair of particles
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
  sim.par.permitivityTop = 0.05;
  sim.par.permitivityBottom = 0.05;
  sim.par.Nxy = 135;
  sim.par.temperature = 0;
  sim.par.viscosity = 1.0/(6*M_PI);
  sim.par.hydrodynamicRadius = 1.0;
  sim.par.wetFraction = 1;
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
  parBD.wetFraction = sim.par.wetFraction;
  parBD.dt = sim.par.dt;
  parBD.brownianUpdateRule = DryWetBD::update_rules::euler_maruyama;
  parBD.Lxy = sim.par.Lxy;
  parBD.H = sim.par.H;
  parBD.hxy_stokes = 0.64;
  auto bd = std::make_shared<BD>(sim.pd, parBD);
  bd->addInteractor(poisson);
  bd->forwardTime();

  scalar dx = sim.pd->getPos(access::cpu, access::write)[1].x-xp1;
  scalar dy = sim.pd->getPos(access::cpu, access::write)[1].y-yp1;
  scalar dz = sim.pd->getPos(access::cpu, access::write)[1].z-zp1;
  std::ofstream out("FULLWETdisplacements.dat");// Will be needed for the next test (periodicity)
  out << dx << " " << dy << " " << dz << std::endl;
  std::cout << "x displacement is " << dx << std::endl;
  std::cout << "y displacement is " << dy << std::endl;
  std::cout << "z displacement is " << dz << std::endl;
  
  //given by the DPStokes solver code (pair_mobility.py)
  scalar nu_xx = 0.138291767457352;
  scalar nu_xy = 0;
  scalar nu_xz = -0.038706105653531;
  scalar nu_yx = 0;
  scalar nu_yy = 0.047304644896516;
  scalar nu_yz = 0;
  scalar nu_zx = -0.081720265183282;
  scalar nu_zy = 0;
  scalar nu_zz = 0.037554341201021;
  scalar mu_xx = 0.893993039433745;
  scalar mu_yy = 0.893956879655325;
  scalar mu_zz = 0.829029667444467;
  scalar expecteddx = (nu_xx*F0x+nu_xy*F0y+nu_xz*F0z+mu_xx*F1x)*sim.par.dt;
  scalar expecteddy = (nu_yx*F0x+nu_yy*F0y+nu_yz*F0z+mu_yy*F1y)*sim.par.dt;
  scalar expecteddz = (nu_zx*F0x+nu_zy*F0y+nu_zz*F0z+mu_zz*F1z)*sim.par.dt;
  std::cout << "expected x displacement is " << expecteddx << std::endl;
  std::cout << "expected y displacement is " << expecteddy << std::endl;
  std::cout << "expected z displacement is " << expecteddz << std::endl;
  EXPECT_THAT(tolerance(dx, expecteddx), ::testing::IsTrue);
  EXPECT_THAT(tolerance(dy, expecteddy), ::testing::IsTrue);
  EXPECT_THAT(tolerance(dz, expecteddz), ::testing::IsTrue);
}

// Full Wet Mode: periodicity
// This test follows the previous test.
TEST(FULLWET, Periodicity){
  UAMMD sim;
  
  sim.par.numberParticles = 2;
  sim.pd = std::make_shared<ParticleData>(sim.par.numberParticles);
  sim.par.Lxy = 76.8;
  sim.par.H = 19.2;
  scalar yp0  = 0.5*sim.par.Lxy;
  scalar yp1  = yp0;
  scalar zp0 = -0.5*sim.par.H+4;
  scalar zp1 = -0.5*sim.par.H+7;
  scalar xp0 = 1*sim.par.Lxy+2;// Moved to the next block
  scalar xp1 = 1*sim.par.Lxy-2;
  
  sim.pd->getPos(access::location::cpu, access::mode::write)[0] = {xp0,yp0,zp0};
  sim.pd->getPos(access::location::cpu, access::mode::write)[1] = {xp1,yp1,zp1};
  sim.pd->getCharge(access::location::cpu, access::mode::write)[0] = 1;
  sim.pd->getCharge(access::location::cpu, access::mode::write)[1] = -1;

  real c0 = sim.pd->getCharge(access::location::cpu, access::mode::read)[0];
  real c1 = sim.pd->getCharge(access::location::cpu, access::mode::read)[1];
   
  sim.par.gw = 0.25;
  sim.par.permitivity = 1;
  sim.par.permitivityTop = 0.05;
  sim.par.permitivityBottom = 0.05;
  sim.par.Nxy = 135;
  sim.par.temperature = 0;
  sim.par.viscosity = 1.0/(6*M_PI);
  sim.par.hydrodynamicRadius = 1.0;
  sim.par.wetFraction = 1;
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
  parBD.wetFraction = sim.par.wetFraction;
  parBD.dt = sim.par.dt;
  parBD.brownianUpdateRule = DryWetBD::update_rules::euler_maruyama;
  parBD.Lxy = sim.par.Lxy;
  parBD.H = sim.par.H;
  parBD.hxy_stokes = 0.64;
  auto bd = std::make_shared<BD>(sim.pd, parBD);
  bd->addInteractor(poisson);
  bd->forwardTime();

  scalar dx = sim.pd->getPos(access::cpu, access::write)[1].x-xp1;
  scalar dy = sim.pd->getPos(access::cpu, access::write)[1].y-yp1;
  scalar dz = sim.pd->getPos(access::cpu, access::write)[1].z-zp1;
  std::ifstream in("FULLWETdisplacements.dat");// Should be available from the previous test
  double s;
  std::vector<double> expectedvals;
  while(in >> s){
    expectedvals.push_back(s);
  }
  
  std::cout << "x displacement is " << dx << std::endl;
  std::cout << "y displacement is " << dy << std::endl;
  std::cout << "z displacement is " << dz << std::endl;
  
  EXPECT_THAT(tolerance(dx, expectedvals[0]), ::testing::IsTrue);
  EXPECT_THAT(tolerance(dy, expectedvals[1]), ::testing::IsTrue);
  EXPECT_THAT(tolerance(dz, expectedvals[2]), ::testing::IsTrue);
}


// 50% Wet, 50% Dry: Test if an integration works for a pair of particles
TEST(DRYWET, Integration){
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
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
  sim.par.permitivityTop = 0.05;
  sim.par.permitivityBottom = 0.05;
  sim.par.Nxy = 135;
  sim.par.temperature = 0;
  sim.par.viscosity = 1.0/(6*M_PI);
  sim.par.hydrodynamicRadius = 1.0;
  sim.par.wetFraction = 0.5;
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

  VectorXd F(6);
  F << c0*fieldAtParticles[0].x, c0*fieldAtParticles[0].y, c0*fieldAtParticles[0].z,
       c1*fieldAtParticles[1].x, c1*fieldAtParticles[1].y, c1*fieldAtParticles[1].z;
  
  using BD = DryWetBD;
  BD::Parameters parBD;
  parBD.temperature = sim.par.temperature;
  parBD.viscosity = sim.par.viscosity;
  parBD.hydrodynamicRadius = sim.par.hydrodynamicRadius;
  parBD.wetFraction = sim.par.wetFraction;
  parBD.dt = sim.par.dt;
  parBD.brownianUpdateRule = DryWetBD::update_rules::euler_maruyama;
  parBD.Lxy = sim.par.Lxy;
  parBD.H = sim.par.H;
  parBD.hxy_stokes = 1.28;
  auto bd = std::make_shared<BD>(sim.pd, parBD);
  bd->addInteractor(poisson);
  bd->forwardTime();

  scalar dx = sim.pd->getPos(access::cpu, access::write)[1].x-xp1;
  scalar dy = sim.pd->getPos(access::cpu, access::write)[1].y-yp1;
  scalar dz = sim.pd->getPos(access::cpu, access::write)[1].z-zp1;
  std::cout << "x displacement is " << dx << std::endl;
  std::cout << "y displacement is " << dy << std::endl;
  std::cout << "z displacement is " << dz << std::endl;

  // mobility matrix computed by the DPStokes solver (python interface)
  MatrixXd mWet(6,6);
  mWet.row(0) << 0.708884362108862, 0, 0, 0.257442239887106, 0, -0.126448728742101;
  mWet.row(1) << 0, 0.708894111245968, 0, 0, 0.120242666875890, 0;
  mWet.row(2) << 0, 0, 0.489019735544522, -0.051068001579525, 0, 0.081803173383461;
  mWet.row(3) << 0.257442991316777, 0, -0.051065753076216, 0.791052032682486, 0, 0;
  mWet.row(4) << 0, 0.120242407250183, 0, 0, 0.791061727610104, 0;
  mWet.row(5) << -0.126448405742272, 0, 0.081803264694953, 0, 0, 0.666364258639729;
  std::cout << "M_wet = " << std::endl;
  std::cout << mWet << std::endl;

  MatrixXd mDry(6,6);
  mDry.row(0) << 0.708884362108862, 0, 0, 0, 0, 0;
  mDry.row(1) << 0, 0.708894111245968, 0, 0, 0, 0;
  mDry.row(2) << 0, 0, 0.489019735544522, 0, 0, 0;
  mDry.row(3) << 0, 0, 0, 0.791052032682486, 0, 0;
  mDry.row(4) << 0, 0, 0, 0, 0.791061727610104, 0;
  mDry.row(5) << 0, 0, 0, 0, 0, 0.666364258639729;
  std::cout << "M_dry = " << std::endl;
  std::cout << mDry << std::endl;

  VectorXd displacement(6);
  displacement = 0.5*(mWet+mDry)*F*sim.par.dt;

  scalar expecteddx = displacement(3);
  scalar expecteddy = displacement(4);
  scalar expecteddz = displacement(5);
  std::cout << "expected x displacement is " << expecteddx << std::endl;
  std::cout << "expected y displacement is " << expecteddy << std::endl;
  std::cout << "expected z displacement is " << expecteddz << std::endl;
  EXPECT_THAT(tolerance(dx, expecteddx), ::testing::IsTrue);
  EXPECT_THAT(tolerance(dy, expecteddy), ::testing::IsTrue);
  EXPECT_THAT(tolerance(dz, expecteddz), ::testing::IsTrue);
}

// ############## Tests by Raul ############## //
TEST(DryWetMobility, CanBeCreated){
  using BD = DryWetBD;
  BD::Parameters par;
  par.temperature = 1.0;
  par.viscosity = 1.0;
  par.hydrodynamicRadius = 1.0;
  par.dt = 1.0;
  par.wetFraction = 0.7;
  par.brownianUpdateRule = DryWetBD::update_rules::euler_maruyama;
  par.H = 16;
  par.Lxy = 32;
  par.hxy_stokes = -1.0;
  auto pd = std::make_shared<ParticleData>(1);
  auto bd = std::make_shared<BD>(pd, par);
}


TEST(FullDryMobility, SelfMobilityIsCorrect){
  using BD = DryWetBD;
  writeDefaultMobilityFile();
  BD::Parameters par;
  par.temperature = 0;
  par.viscosity = 1.0/(6*M_PI);
  par.hydrodynamicRadius = 1.0;
  par.dt = 1.0;
  par.wetFraction = 0;
  par.brownianUpdateRule = DryWetBD::update_rules::euler_maruyama;
  par.dryMobilityFile = "uniformMob.dat";
  par.H = 32;
  par.Lxy = 64;
  par.hxy_stokes = -1.0;
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
  par.wetFraction = 1;
  par.brownianUpdateRule = DryWetBD::update_rules::euler_maruyama;
  par.dryMobilityFile = "uniformMob.dat";
  par.H = 128;
  par.Lxy = 64;
  par.hxy_stokes = -1.0;
  auto pd = std::make_shared<ParticleData>(1);
  pd->getPos(access::cpu, access::write)[0] = real4();
  auto bd = std::make_shared<BD>(pd, par);
  bd->addInteractor(std::make_shared<miniInteractor>(pd));
  bd->forwardTime();
  real M0 = pd->getPos(access::cpu, access::write)[0].x;
  ASSERT_THAT(M0, ::testing::DoubleNear(1, 1e-1));
}


//Asserts the correctness of the self mobility for a certain wet radius
//All other parameters are hardcoded, see the function.
//The total hydrodynamic radius is 1
void computeSelfMobilityWithWetRadius(real wetFraction){
  using BD = DryWetBD;
  writeDefaultMobilityFile();
  BD::Parameters par;
  par.temperature = 0;
  par.viscosity = 1.0/(6*M_PI);
  par.hydrodynamicRadius = 1.0;
  par.dt = 1.0;
  par.wetFraction = wetFraction;
  par.brownianUpdateRule = DryWetBD::update_rules::euler_maruyama;
  par.dryMobilityFile = "uniformMob.dat";
  par.H = 64;
  par.Lxy = 64;
  par.hxy_stokes = -1.0;
  auto pd = std::make_shared<ParticleData>(1);
  pd->getPos(access::cpu, access::write)[0] = real4();
  auto bd = std::make_shared<BD>(pd, par);
  bd->addInteractor(std::make_shared<miniInteractor>(pd));
  bd->forwardTime();
  real M0 = pd->getPos(access::cpu, access::write)[0].x;
  ASSERT_THAT(M0, ::testing::DoubleNear(1, 1e-1))<<"Failed with wet fraction "<<wetFraction;
}

TEST(DryWetMobility, SelfMobilityIsCorrectForAnyWetRadius){
  real minwetFraction = 0;
  real maxwetFraction = 1;
  int Ntest = 4;
  fori(0, Ntest){
    real wetFraction = minwetFraction + i*(maxwetFraction - minwetFraction)/(Ntest-1);
    computeSelfMobilityWithWetRadius(wetFraction);
  }
}

