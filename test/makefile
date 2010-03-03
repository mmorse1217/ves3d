
CXX      = icpc
CC       = icpc
#CUDA     = /usr/local/cuda/bin/nvcc

CUDA_DIR = /usr/local/cuda

# INCLUDE  = -I../include/ -I../src/ -I$(CUDA_DIR)/include
# CXXFLAGS = -g -O0 -fPIC  #-DNDEBUG 
# LDFLAGS  = -L$(CUDA_DIR)/lib64 -lcuda -lcudart -Wl,-rpath $(CUDA_DIR)/lib64  

MKLPATH = /opt/intel/Compiler/11.1/059/mkl/lib/em64t

LDFLAGS = -L$(CUDA_DIR)/lib64 -lcublas -lcudart -Wl,-rpath            \
$(CUDA_DIR)/lib64 -L$(MKLPATH) $(MKLPATH)/libmkl_solver_lp64.a        \
-Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core      \
-Wl,--end-group -openmp -lpthread 

CXXFLAGS = -g -O0 -I$(CUDA_DIR)/include                               \
-I/usr/local/apps/intel/mkl/10.0.3.020/include -I../include/          \
-I../src/ -fPIC  #-DNDEBUG #-DCUDA_SHT_PROFILING 


# SHTFLAGS = -g -O0 -I$(CUDA_DIR)/include					\
# -I/usr/local/apps/intel/mkl/10.0.3.020/include -I../include/ -I../src/	\
# -fPIC -L$(CUDA_DIR)/lib64 -lcublas -lcudart -Wl,-rpath			\
# $(CUDA_DIR)/lib64 -L$(MKLPATH) $(MKLPATH)/libmkl_solver_lp64.a		\
# -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core	\
# -Wl,--end-group -openmp -lpthread



all: 

DeviceTest: #DeviceGPU.o 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCLUDE) $@.cc $^ -o $@.out

DeviceGPU.o: ../include/DeviceGPU.cu
	$(CUDA) $(INCLUDE) -c ../include/DeviceGPU.cu

ScalarsTest:
	$(CXX) $(CXXFLAGS) $(INCLUDE) $@.cc -o $@.out

VectorsTest:
	$(CXX) $(CXXFLAGS) $(INCLUDE) $@.cc -o $@.out

SHTransTest: cuda_sht.o blas_sht.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) $@.cc $^ -o $@.out

driver: cuda_sht.o blas_sht.o driver.o 

cuda_sht.o: ../src/cuda_sht.C
	$(CXX) $(SHTFLAGS) -c ../src/cuda_sht.C

blas_sht.o: ../src/blas_sht.C
	$(CXX) $(SHTFLAGS) -c ../src/blas_sht.C

.PHONY: clean

clean: 	
	@-rm -f *.o *.out  

# MKLPATH = /opt/intel/Compiler/11.1/059/mkl/lib/em64t

# LDFLAGS = -L$(CUDA_DIR)/lib64 -lcublas -lcudart -Wl,-rpath            \
# $(CUDA_DIR)/lib64 -L$(MKLPATH) $(MKLPATH)/libmkl_solver_lp64.a        \
# -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core      \
# -Wl,--end-group -openmp -lpthread 

# CPPFLAGS = -g -O0 -I$(CUDA_DIR)/include                               \
# -I/usr/local/apps/intel/mkl/10.0.3.020/include -I../include/          \
# -I../src/ -fPIC  #-DNDEBUG #-DCUDA_SHT_PROFILING 

# SurfaceTest: cuda_sht.o

# CurvFlow: cuda_sht.o 

# CUTest: CUvec.o 
# 	$(CXX) $(CPPFLAGS) $(LDFLAGS) $(INCLUDE) $@.cc $^ -o $@.out

# CUvec.o: ../src/CUvec.cu
# 	$(CUDA) -c ../src/CUvec.cu

# MEXSUFFIX = mexa64	
# MEX       = mex
# MEXFLAGS  = -cxx -I../include/ -I../src/
# CurvFlowMex: CurvFlowMex.cc cuda_sht.o CurvFlow.o 
# 	$(MEX) $(MEXFLAGS) $(LDFLAGS) -output CurvFlowMex CurvFlowMex.cc CurvFlow.o cuda_sht.o

