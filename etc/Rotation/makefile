#CPU
CXX      = g++
CC       = gcc

#GPU
CUDA_DIR      = /usr/local/cuda3.2/cuda/
NVCC          = $(CUDA_DIR)/bin/nvcc
CUDA_LIBS     = -L$(CUDA_DIR)/lib64 -lcublas -lcudart	
CUDA_INCLUDE  = -I$(CUDA_DIR)/include

#Flags
CXXFLAGS  = -fopenmp -fPIC -fno-exceptions -w -g -O0		 \
	    -DGPU_ACTIVE -DVERBOSE -DPROFILING 
INCLUDE   = $(CUDA_INCLUDE)
#LDFLAGS   = -L/usr/lib64/ /usr/lib64/liblapack.so.3.0.3 -lblas $(CUDA_LIBS)
LDFLAGS   = -L/usr/lib64/ -L/opt/intel/mkl/10.2.5.035/lib/em64t -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread $(CUDA_LIBS)
CUDAFLAGS = -arch=sm_20

#rules
all: SyntheticRotation.exe

SyntheticRotation.exe: SyntheticRotation.o Logger.o CPUKernels.o CudaKernels.o 
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

%.o: %.cc 
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

%.o: %.cu
	$(NVCC) $(CUDAFLAGS) $(CUDA_INCLUDE) -c -o $@ $<

clean: 
	-$(RM) *.o *.exe *.out 
