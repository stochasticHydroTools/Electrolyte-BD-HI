

CXX=g++
NVCC=nvcc
#You might have to change this if you want to use MKL instead of lapacke. Only one is needed.
LAPACKE_FLAGS=-llapacke -I/usr/include/lapacke -lcblas
VERBOSITY=5
#MKL_FLAGS=-DUSE_MKL -DMKL_ILP64 -m64 -I${MKLROOT}/include -L${MKLROOT}/lib/intel64  -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl


UAMMD_ROOT=uammd/
CUDA_ROOT:=$(shell dirname `which nvcc`)/..
INCLUDEFLAGS= -I$(CUDA_ROOT)/include -I $(UAMMD_ROOT)/src -I $(UAMMD_ROOT)/src/third_party
BASIC_LINE= $(NVCC) -ccbin=$(CXX) -O3 -std=c++14 -x cu $(INCLUDEFLAGS) --expt-relaxed-constexpr $(MKL_FLAGS) $(LAPACKE_FLAGS) -DMAXLOGLEVEL=$(VERBOSITY) --expt-extended-lambda #-DDOUBLE_PRECISION


all: slab

slab: slab.cu RepulsivePotential.cuh DryDiffusion.cuh
	$(BASIC_LINE) $<  -o $@ -lcufft -lcublas

clean:
	rm -f poisson
