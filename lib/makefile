include ${VES3D_DIR}/makefile.in.files/makefile.in
LIB_SRC= ${VES3D_DIR}/src/CPUKernels.cc ${VES3D_DIR}/src/Logger.cc

ARFILES=libves3d.a

ifeq ($(USE_GPU),yes)
	LIB_SRC_GPU= ${VES3D_DIR}/src/CudaKernels.cu ${VES3D_DIR}/src/transpose_kernel.cu
	LIB_OBJ_GPU=$(LIB_SRC_GPU:.cu=.o)	
	ARFILES+= libves3dgpu.a 
endif

DEP = $(LIB_SRC:.cpp=.d)
LIB_OBJ=$(LIB_SRC:.cc=.o)	

all:$(ARFILES)

libves3d.a: $(LIB_OBJ) $(MakeFiles)
	$(AR) $(ARFLAGS) ${VES3D_DIR}/lib/libves3d.a $(LIB_OBJ)
	$(RANLIB) ${VES3D_DIR}/lib/libves3d.a

libves3dgpu.a: $(LIB_OBJ_GPU) $(MakeFiles)
	$(AR) $(ARFLAGS) ${VES3D_DIR}/lib/libves3dgpu.a $(LIB_OBJ_GPU)
	$(RANLIB) ${VES3D_DIR}/lib/libves3dgpu.a