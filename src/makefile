all:libves3.a

include ${VES3D_DIR}/makefile.in
LIB_SRC= BlasSht.cc 

ifeq ($(USE_GPU),yes)
	LIB_OBJ=$(LIB_SRC:.cc=.o) CudaSht.o  DeviceGPUSrc.o
else
	LIB_OBJ=$(LIB_SRC:.cc=.o)
endif

DEP = $(LIB_SRC:.cpp=.d)

libves3.a: $(LIB_OBJ)
	$(AR) $(ARFLAGS) $(VES3D_DIR)/lib/libves3d.a $(LIB_OBJ)
	$(RANLIB) $(VES3D_DIR)/lib/libves3d.a
	/bin/rm -rf *.o
