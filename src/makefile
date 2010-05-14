include ${VES3D_DIR}/makefile.in.files/makefile.in
LIB_SRC= BlasSht.cc CPUKernels.cc Logger.cc

ARFILES=libves3.a

ifeq ($(USE_GPU),yes)
	LIB_SRC+= CudaSht.cc
	LIB_SRC_GPU= CudaKernels.cu transpose_kernel.cu
	LIB_OBJ_GPU=$(LIB_SRC_GPU:.cu=.o)	
	ARFILES+= libves3dgpu.a 
endif

DEP = $(LIB_SRC:.cpp=.d)
LIB_OBJ=$(LIB_SRC:.cc=.o)	


all:$(ARFILES)

libves3.a: $(LIB_OBJ)
	$(AR) $(ARFLAGS) $(VES3D_DIR)/lib/libves3d.a $(LIB_OBJ)
	$(RANLIB) $(VES3D_DIR)/lib/libves3d.a
	/bin/rm -rf *.o

libves3dgpu.a: $(LIB_OBJ_GPU)
	$(AR) $(ARFLAGS) $(VES3D_DIR)/lib/libves3dgpu.a $(LIB_OBJ_GPU)
	$(RANLIB) $(VES3D_DIR)/lib/libves3dgpu.a
	/bin/rm -rf *.o