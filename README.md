


# Introduction

This repo encodes the following simulation:  

A group of Gaussian sources of charge with width gw evolve via Brownian Dynamics inside a doubly periodic domain of dimensions Lxy, Lxy, H (periodic in XY and open in Z).  
Optionally, a flag can be passed in the input to enable triply periodic electrostatics.  

Sources interact via an electrostatic potential and repell each other via a repulsive LJ-like potential.  

There are two (potentially) charged walls on $z=-H/2$ and $z=H/2$, this walls have the same charge equal to half of the sum of all charges in the system but with opposite sign (so the overall system is always electroneutral). This is a quirk of the algorithm which is present unless the walls are metallic. The potential (when walls are metallic) or total charge (otherwise) on the bottom wall can be set to avoid this.  

Charges are repelled by the wall via the same repulsive potential they have between each other. See option **imageDistanceMultiplier** below.  
The three domains demarcated by the walls (above, between and below) may have different, arbitrary permittivities.  

The repulsive potential can be found in RepulsivePotential.cuh (should be available along this file) and has the following form:  

$$U_{LJ}(r)=4U_0 ( (\sigma/r)^{2p} - (\sigma/r)^p ) + U_0$$  

The repulsive force is then defined by parts using  

$$F_{LJ}=-\partial U_{LJ}/\partial r$$  

$$F(r\le r_m) = F_{LJ}(r=r_m)$$  

$$F(r_m\lt r\le 2^{1/p}\sigma) = F_{LJ}(r)$$  

$$F(r\gt 2^{1/p}\sigma) = 0$$  

Hydrodynamic interactions in the slit channel can be incorporated via a wet hydrodynamic radius (using DPStokes).  
The special wet radius equal to zero means neglecting hydrodynamics (pure dry diffusion).  
On the contrary, making the wet radius equal to the total radius results in a pure wet diffusion simulation.  


# Compilation

This code depends on:  

1.  CUDA
2.  MKL or Lapack+blas

If you conda available, you can install everything you need to compile with:  

    conda install -c nvidia cuda mkl-devel cmake make 

Two methods of compilation are available: Make and CMake  
Tests can only be compiled using CMake  


## Make

To use Make, you will have to tweak the Makefile according to your system and the simply run:  

    $ make

This will compile slab.cu into slab.  


## CMake

CMake will try to autodetect the necessary paths in your system. You can compile both slab.cu and the tests with:  

    $ mkdir -p build
    $ cd build
    $ cmake ..
    $ make

It is possible that CMake complains about some missing dependency (maybe MKL, or Lapack/blas), in this case you have to install said dependency and then run cmake again.  
Note that the CMake configuration will compile in debug mode by default. This is much slower to both compile and run.  
You can disable debugging by running cmake as:  

    $cmake -DCMAKE_BUILD_TYPE=Release ..


# Tests

Several unit tests are available in tests.cu. They use [Gtest](https://github.com/google/googletest). CMake will autodownload it for you. Compile with CMake and then run ./tests  
Gtest provides a lot of useful functionality for dealing with tests (for instance, running only a subset of them). Try ./tests &#x2013;help  


# Electroosmosis

The directory EOSSimulations includes bash scripts to run simulations of electroosmosis in nanochannels. The simulation has two components:  


## Equilibrium

Run by `bash runEquilibrium.bash charged numberSimulations wetFraction`. Here charged can be 0 (uncharged surface) or 1 (charged surface), numberSimuations (int) is the total number of independent simulations, and 0<=wetFraction<=1 is the wet fraction. The results will be stored in DIR="Equilibrium\_\*Surface-longrun"  

## Electroosmosis

This should be run after equilibrium. The corresponding equilibrium distribution results directory (i.e., "Equilibrium\_\*Surface-longrun") should be available in EOSSimulations/. Run by `bash runEOS.bash charged numberSimulations wetFraction`.

## Postprocessing

Run by `bash postprocessing.bash charged numberSimulations`. The equilibrium and all of the electroosmosis (i.e., "EOS_*") results directories should be available in EOSSimulations/.

The subdirectory Analysis includes MATLAB scripts to plot the data. Sample plots are provided.  


# USAGE:

This code expects to find a file called data.main in the folder where it is executed.  
data.main contains a series of parameters that allow to customize the simulation and must have the following  format:  

    ----data.main starts
    #Lines starting with # are ignored
    option [argument] #Anything after the argument is ignored
    flag #An option that does not need arguments
    ----end of data.main


## The following options are available:

**numberParticles**: The number of charges in the simulation  

**gw**: The Gaussian width of the charges  

**H**: The thickness of the domain for periodic and Lz for triply periodic.  

**Lxy**: The dimensions of the box in XY  

**permitivity**: Permittivity inside the slab  

**permitivityBottom** Below z=-H/2. If the value is negative it means metallic boundary (infinite permittivity).  

**permitivityTop** Above z=H/2. If the value is negative it means metallic boundary (infinite permittivity).  

**bottomWallSurfaceValue** The zero mode value of the Fourier transform of the bottom wall surface value (potential when the boundary is metallic and surface charge otherwise).  

**noElectrostatics** Optional. If present electrostatics are not included, every other option related to electrostatics is ignored.  

**temperature**: Temperature for the Brownian Dynamics integrator, the diffusion coefficient will be D=T/(6\*pi\*viscosity\*hydrodynamicRadius). This temperature is therefore given in units of energy.  

**viscosity**: For wet diffusion in BD  

**hydrodynamicRadius**: Total hydrodynamic radius  

**wetFraction**: between 0 and 1; wetRadius = hydrodynamicRaius/wetFrcation; dryRadius = hydrodynamicRadius/(1-wetFraction)  

**dt**: Time step for the BD integrator  

**U0, sigma, r\_m, p**: Parameters for the repulsive interaction. If U0=0 the steric repulsion is turned off.  

**wall\_U0, wall\_sigma, wall\_r\_m, wall\_p** Parameters for the ion-wall repulsive interaction.  

**imageDistanceMultiplier** Multiplies the distance of the particles to the wall by this amount. For instance, if 2, particles interact with their images, if 1, particles are repelled to the wall (as if the image was at the wall's height)  

**numberSteps**: The simulation will run for this many steps  

**printSteps**: If greater than 0, the positions and forces will be printed every printSteps steps  

**relaxSteps**: The simulation will run without printing for this many steps.  

**outfile**: Positions and charge will be written to this file, each snapshot is separated by a #, each line will contain X Y Z Charge. Can be /dev/stdout to print to screen.  

**forcefile**: Optional, if present forces acting on particles will be written to this file.  

**fieldfile**: Optional, if present electric field acting on particles will be written to this file.  

**velocityfile** Average fluid velocity along the X direction will be written to this file.  

**readFile**: Optional, if present charge positions will be read from this file with the format X Y Z Charge. numberParticles lines will be read. Can be /dev/stdin to read from pipe.  

**Nxy**: The number of cells in XY for the DPPoisson algorithm.  

**fold**: 0 or 1; whether to fold the periodic box (1) or not (0).  

**hxy\_stokes**: grid size in the planar xy direction. It will be automatically computed if set to a negative value.  

**externalField** real3: applied external electric field in the x, y, z directions.  

**useMobilityFromFile**: Optional, if this option is present, the mobility will depend on the height of the particle according to the data in this file.This file must have two columns with a list of normalized heights (so Z must go from -1 to 1) and normalized mobilities (i.e. 6\*pi\*eta\*a) in X, Y and Z. The values for each particle will be linearly interpolated from the data provided in the file. The order of the values does not matter.  
Example:  

    --- mobility.dat---
    -1.0 1.0 1.0 1.0
     0.0 1.0 1.0 1.0
     1.0 1.0 1.0 1.0
    -------------------

If the option is not present the mobility will be autocomputed using DPStokes.  

**BrownianUpdateRule**: Optional. Can either be EulerMaruyama (default) or Leimkuhler.  

**idealParticles**: Optional. If this flag is present particles will not interact between them in any way.  

