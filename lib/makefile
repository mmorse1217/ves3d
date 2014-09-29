include $(VES3D_DIR)/makefile.in.files/makefile.in

LIB_SRC = $(VES3D_SRC)/CPUKernels.cc \
	  $(VES3D_SRC)/Logger.cc	 \
	  $(VES3D_SRC)/Enums.cc      \
	  $(VES3D_SRC)/Error.cc      \
	  $(VES3D_SRC)/DataIO.cc	 \
	  $(VES3D_SRC)/anyoption.cc

ARFILES = libves3d.a

ifeq ($(USE_GPU),yes)
  LIB_SRC_GPU = $(VES3D_DIR)/src/CudaKernels.cu
  LIB_SRC    +=	$(VES3D_DIR)/src/CudaApiGlobals.cc
  LIB_OBJ_GPU = $(LIB_SRC_GPU:.cu=.o)
  ARFILES    += libves3dgpu.a
endif

LIB_OBJ = $(LIB_SRC:.cc=.o)

all: $(ARFILES)

libves3d.a: $(LIB_OBJ) $(MAKEIN)
	$(AR) $(ARFLAGS) $(VES3D_DIR)/lib/libves3d.a $(LIB_OBJ)
	$(RANLIB) $(VES3D_DIR)/lib/libves3d.a

libves3dgpu.a: $(LIB_OBJ_GPU) $(MAKEIN)
	$(AR) $(ARFLAGS) $(VES3D_DIR)/lib/libves3dgpu.a $(LIB_OBJ_GPU)
	$(RANLIB) $(VES3D_DIR)/lib/libves3dgpu.a

clean:
	-$(RM) *.a
