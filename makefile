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

# TODO:add to makefile.col
#initialize include and lib paths
SURFACE_DIR = ${BLENDSURF_DIR} ${FACEMAP_DIR}
SURFACE_INC = ${SURFACE_DIR} ${FACEMAP_DIR}/src
SURFACE_LIB = ${SURFACE_DIR} ${FACEMAP_DIR}/lib

#p4est
P4EST_INC = ${P4EST_DIR}/include 
SC_INC = ${P4EST_DIR}/sc/include
P4EST_LIB = ${P4EST_DIR}/lib 
SC_LIB = ${P4EST_DIR}/lib

#FFTW
FFTW_INC = ${FFTW_DIR}/include
FFTW_LIB = ${FFTW_DIR}/lib

#Petsc
PETSC_BUILD_INC = ${PETSC_DIR}/include ${PETSC_DIR}/${PETSC_ARCH}/include
PETSC_LIB = ${PETSC_DIR}/${PETSC_ARCH}/lib

#PvFMM
#PVFMM_INC = ${PVFMM_DIR}/include utils/pvfmm-utils/
#PVFMM_LIB = ${PVFMM_DIR}/lib/pvfmm utils/pvfmm-utils/
PVFMM_INC = /usr/local/include/pvfmm ${HEDGEHOG_DIR}/utils/pvfmm-utils/
PVFMM_LIB = /usr/local/lib/pvfmm ${HEDGEHOG_DIR}/utils/pvfmm-utils/

MPI_INC = ${MPI_HOME}/include
MPI_LIB = ${MPI_HOME}/lib

VTK_INC = ${VTK_DIR}/${VTK_INC_DIR}
VTK_LIB = ${VTK_DIR}/lib


VTK_LIBS_TO_LINK = -lvtkIOXML-8.1 -lvtkChartsCore-8.1 -lvtkCommonDataModel-8.1 -lvtkCommonCore-8.1 -lvtkFiltersCore-8.1 -lvtkCommonExecutionModel-8.1
#aggregate
DEPS_INC = ${SURFACE_INC} ${SC_INC} ${P4EST_INC} ${FFTW_INC} ${PETSC_BUILD_INC} ${PVFMM_INC} ${MPI_INC} ${VTK_INC}
DEPS_LIB = ${SURFACE_LIB} ${SC_LIB} ${P4EST_LIB} ${FFTW_LIB} ${PETSC_LIB} ${PVFMM_LIB} ${MPI_LIB} ${FORTRAN_LIB} ${VTK_LIB}
DEPS_LIB_FLAGS = -lhedgehog -lpatchwork -lblend -lsc -lp4est -lfftw3 -lpvfmmwrap -lpvfmm -lpetsc -lexpat -lblas -llapack -lmpi -lm -fopenmp ${VTK_LIBS_TO_LINK}
HEDGEHOG_INC = -I/usr/include -I${HEDGEHOG_DIR}/src/ ${DEPS_INC:%=-I%}
HEDGEHOG_LIB = -L${HEDGEHOG_DIR}/lib ${DEPS_LIB:%=-L%} ${DEPS_LIB_FLAGS}
#

LDFLAGS += -L/opt/intel/compilers_and_libraries_2018.2.164/mac/compiler/lib
LDFLAGS += -L/opt/intel/compilers_and_libraries_2018.2.164/mac/mkl/lib
LDFLAGS += -L/Users/libinlu/Documents/Projects/ves3d/contact3d/lib
LDFLAGS += -L/usr/local/Cellar/fftw/3.3.5/lib
LDFLAGS += -L/usr/local/Cellar/cgal/4.9/lib
LDFLAGS += -L/usr/local/lib

LDLIBS += -lCGAL -lgmp -llapack -liagm ${HEDGEHOG_LIB} # -liomp5

EIGEN_INC := /opt/local/include/eigen3
CONTACT_INC := /Users/libinlu/Documents/Projects/ves3d/contact3d/inc
CXXFLAGS += -I$(EIGEN_INC) -I$(CONTACT_INC) ${HEDGEHOG_INC} -I${MKL_INC}
#

# targets of install
VES3D_BINS = ves3d

print-%:
	@echo '$*=$($*)'

all: install

lib:
	${MAKE} -C  ${VES3D_LIBDIR}

install: $(addprefix ${VES3D_BINDIR}/,${VES3D_BINS})

test:
	${MAKE} -C  ${VES3D_TSTDIR}

check:
	${MAKE} -C ${VES3D_TSTDIR} check

doc: ${MAKE_DEP}
	${DOX} ${VES3D_DOCDIR}/Doxyfile

tag: ${MAKE_DEP}
	${TAGS} ${VES3D_SRCDIR}/*.cc  ${VES3D_INCDIR}/*.h ${VES3D_TSTDIR}/*.cc ${VES3D_TSTDIR}/*.h

all-tag: ${MAKE_DEP}
	${TAGS} ${VES3D_SRCDIR}/*  ${VES3D_INCDIR}/* ${VES3D_TSTDIR}/*

${VES3D_BINDIR}/%: ${VES3D_SRCDIR}/%.o ${MAKE_DEP} lib
	-@${MKDIR} $(dir $@)
	${CXX} ${CXXFLAGS} ${CPPFLAGS} ${LDFLAGS} $< ${LDLIBS} -o $@
	${CXX} -MM -MT ${VES3D_SRCDIR}/$*.o ${CXXFLAGS} ${CPPFLAGS} -c -o ${VES3D_SRCDIR}/$*.d ${VES3D_SRCDIR}/$*.cc

clean:
	${MAKE} -C ${VES3D_TSTDIR} clean
	${MAKE} -C ${VES3D_LIBDIR} clean
	-${RM} *.o ${VES3D_SRCDIR}/*.o ${VES3D_SRCDIR}/*.d
	-${RM} -r ${VES3D_DOCDIR}/latex ${VES3D_DOCDIR}/html TAGS
	-${RM} *.${VES3D_CHKEXT} ${VES3D_EXPRDIR}/*.${VES3D_CHKEXT} *.vtp *.pvtp

distclean:
	${MAKE} clean
	${MAKE} -C ${VES3D_LIBDIR} distclean
	${MAKE} -C ${VES3D_TSTDIR} distclean
	-${RM} $(addprefix ${VES3D_BINDIR}/,${VES3D_BINS})

# include dependency file for templates (generated when making bin files)
-include $(addprefix ${VES3D_SRCDIR}/, $(addsuffix .d, ${VES3D_BINS}))

help:
	@echo "    Useful compilation flags are:"
	@echo "        VES3D_DEBUG      [yes|no] affects only compiler flags (not printing verbosity)"
	@echo "        VES3D_PROFILE    [yes|no] turns profiling macros on"
	@echo "        VES3D_VERBOSE    [1|2|3]  1=most verbose, 2=normal, 3=quiet"
	@echo "        VES3D_TESTING    [yes|no] equal to DEBUG and PROFILE"
	@echo "        VES3D_VERSION    default to 'hg id -n' (local to repo) can be set as argument if hg is not working on the node"

.PHONY: all lib install test check doc tags all-tags clean distclean
.SECONDARY:
