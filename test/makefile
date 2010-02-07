MKLPATH = /opt/intel/Compiler/11.1/059/mkl/lib/em64t
CXX = icpc
CC = icpc
CUDA_DIR = /usr/local/cuda
LDFLAGS = -L$(CUDA_DIR)/lib64 -lcublas -lcudart -Wl,-rpath              \
$(CUDA_DIR)/lib64 -L$(MKLPATH) $(MKLPATH)/libmkl_solver_lp64.a          \
-Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core        \
-Wl,--end-group -openmp -lpthread

CPPFLAGS = -g -O0 -I$(CUDA_DIR)/include                                 \
           -I/usr/local/apps/intel/mkl/10.0.3.020/include               \
           -I../include/                                                \
           -I../src/                                                    \
           -DCUDA_SHT_PROFILING

all: CurvFlow

SHScalarsTest:

SHVectorsTest:

SphHarmTest: cuda_sht.o blas_sht.o SphHarmTest.o

SurfaceTest: cuda_sht.o SurfaceTest.o

CurvFlow: cuda_sht.o CurvFlow.o

driver: cuda_sht.o blas_sht.o driver.o 

cuda_sht.o: ../src/cuda_sht.C
	$(CXX) $(CPPFLAGS) -c ../src/cuda_sht.C

blas_sht.o: ../src/blas_sht.C
	$(CXX) $(CPPFLAGS) -c ../src/blas_sht.C

.PHONY: clean

clean: 	
	@-rm -f *.o *.out CurvFlow

