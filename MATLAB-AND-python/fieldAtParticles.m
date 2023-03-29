% this MATLAB script, together with run.m calculates the electric field for a pair of particles
% it needs stochasticHydroTools/DoublyPeriodic/Poisson stuff to run 
clc;clear
addpath('./DPFunctions')
addpath('./Poisson')

% Params
H  = 19.2000 ;
L  = 0.5*76.8;
gw = 0.25;
dt = 0.01;
viscosity = 1/(6*pi);

Hew = 6.4;
nxy = 72;
Nz = 85;

hexy = 2*L/nxy;
nC = 2;
pts = [7 6 4.8;11 6 4.8];
charges = [1;-1];
epsilon = 1;
epsilonTop = 0.05;
epsilonBottom = 0.05;

pargiven = 1;

% epsratiot: dielectric jump ratio at the top (epsilon_t=1/epsratiot)
% epsratiob: dielectric jump ratio at the bottom (epsilon_b=1/epsratiob)
% Use 'M' for metallic
if epsilonTop < 0
    epsratiot = 'M';
else
    epsratiot = 1/epsilonTop;
end
if epsilonBottom < 0
    epsratiob = 'M';
else
    epsratiob = 1/epsilonBottom;
end
 
run
mobility = 1;
displacements_MATLAB = ECharges.*charges*mobility*dt;
disp('displacements:')
fprintf('%20.16f %20.16f %20.16f %d\n', displacements_MATLAB(1,1:3),1)
fprintf('%20.16f %20.16f %20.16f %d\n', displacements_MATLAB(2,1:3),-1)
