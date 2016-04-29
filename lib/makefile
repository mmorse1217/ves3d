########################################################
### Do not add any customization to this file.       ###
### put customizations in the platform or cxx files. ###
########################################################

ifndef VES3D_DIR
$(error "$${VES3D_DIR} environment variable is not set.")
endif

# include rules and flag for compiler/host
VES3D_MKDIR ?= ${VES3D_DIR}/makefile.in.files
include ${VES3D_MKDIR}/makefile.in

LIB_SRC = ${VES3D_SRCDIR}/CPUKernels.cc 	\
	  ${VES3D_SRCDIR}/Logger.cc	 	\
	  ${VES3D_SRCDIR}/Enums.cc      	\
	  ${VES3D_SRCDIR}/Error.cc      	\
	  ${VES3D_SRCDIR}/DataIO.cc 		\
	  ${VES3D_SRCDIR}/anyoption.cc		\
	  ${VES3D_SRCDIR}/legendre_rule.cc

LIB_SRC_GPU = ${VES3D_SRCDIR}/CudaKernels.cu
ifeq (${VES3D_USE_GPU},yes)
  LIB_SRC += ${VES3D_SRCDIR}/CudaApiGlobals.cc
endif

LIB_OBJ     = $(LIB_SRC:.cc=.o)
LIB_OBJ_GPU = $(LIB_SRC_GPU:.cu=.o)

all: ${ARFILES}

libves3d.a: ${LIB_OBJ} ${MAKE_DEP}
	${AR} ${ARFLAGS} ${VES3D_LIBDIR}/libves3d.a ${LIB_OBJ}
	${RANLIB} ${VES3D_LIBDIR}/libves3d.a

libves3dgpu.a: ${LIB_OBJ_GPU} ${MAKE_DEP}
	${AR} ${ARFLAGS} ${VES3D_LIBDIR}/libves3dgpu.a ${LIB_OBJ_GPU}
	${RANLIB} ${VES3D_LIBDIR}/libves3dgpu.a

clean:

distclean:
	-$(RM) *.a
