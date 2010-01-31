MKLPATH = /opt/intel/Compiler/11.1/059/mkl/lib/em64t
CXX = icpc
CC = icpc
CUDA_DIR = /usr/local/cuda
LDFLAGS = -L$(CUDA_DIR)/lib64 -lcublas -lcudart -Wl,-rpath              \
$(CUDA_DIR)/lib64 -L$(MKLPATH) $(MKLPATH)/libmkl_solver_lp64.a          \
-Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core        \
-Wl,--end-group -openmp -lpthread

CPPFLAGS = -I$(CUDA_DIR)/include                                        \
           -I/usr/local/apps/intel/mkl/10.0.3.020/include               \
           -I../include/                                                \
           -I../src/                                                    \
           -DCUDA_SHT_PROFILING

CXXFLAGS = -g -O0


#all: driver

SHScalars:
        $(CXX) $(CFLAGS) $(CPPFLAGS) -o $@.out $@Test.cc 

SHVectors:
        $(CXX) $(CFLAGS) $(CPPFLAGS) -o $@.out $@Test.cc 

SphHarm:
        $(CXX) $(CFLAGS) $(CPPFLAGS) -o $@.out $@Test.cc 

driver: cuda_sht.o driver.o blas_sht.o

clean: 
        rm *.o *.out


