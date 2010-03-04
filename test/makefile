CXX  = icpc
CC   = icpc

MKLPATH  = /opt/intel/Compiler/11.1/059/mkl/lib/em64t
CUDA_DIR = /usr/local/cuda
NVCC = $(CUDA_DIR)/bin/nvcc

LDFLAGS  = -L$(CUDA_DIR)/lib64 -lcublas -lcudart -Wl,-rpath           \
$(CUDA_DIR)/lib64 -L$(MKLPATH) $(MKLPATH)/libmkl_solver_lp64.a        \
-Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core      \
-Wl,--end-group -openmp -lpthread 

CPPFLAGS = -g -O0 -I$(CUDA_DIR)/include                               \
-I/usr/local/apps/intel/mkl/10.0.3.020/include -I../include/          \
-I../src/ -fPIC -DNDEBUG #-DCUDA_SHT_PROFILING 

CUDAFLAGS =  -I../src/ -I../include/ -I$(CUDA_DIR)/include

all: 

DeviceTest: DeviceGPU.o 
	$(CXX) $(CPPFLAGS) $(LDFLAGS) $@.cc $^ -o $@.out

DeviceGPU.o: ../include/DeviceGPU.cu
	$(NVCC) $(CUDAFLAGS) -c ../include/DeviceGPU.cu

ScalarsTest:

VectorsTest:

SHTransTest: cuda_sht.o 

driver: cuda_sht.o blas_sht.o 

cuda_sht.o: ../src/cuda_sht.C
	$(CXX) $(CPPFLAGS) -c ../src/cuda_sht.C

blas_sht.o: ../src/blas_sht.C
	$(CXX) $(CPPFLAGS) -c ../src/blas_sht.C

SurfaceTest: cuda_sht.o

.PHONY: clean

clean: 	
	@-rm -f *.o *.out  

# CurvFlow: cuda_sht.o 

# CUTest: CUvec.o 
# 	$(CXX) $(CPPFLAGS) $(LDFLAGS) $(INCLUDE) $@.cc $^ -o $@.out

CUvec.o: ../etc/CUvec.cu
	$(NVCC) -c ../etc/CUvec.cu