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

LDFLAGS += -L/opt/intel/compilers_and_libraries_2018.1.126/mac/compiler/lib
LDFLAGS += -L/opt/intel/compilers_and_libraries_2018.1.126/mac/mkl/lib
LDFLAGS += -L/Users/libinlu/Documents/Projects/ves3d/contact3d/lib
LDFLAGS += -L/usr/local/Cellar/fftw/3.3.5/lib
LDFLAGS += -L/usr/local/Cellar/cgal/4.9/lib
LDFLAGS += -L/usr/local/lib

LDLIBS += -lCGAL -lgmp -llapack -liagm ${HEDGEHOG_LIB} # -liomp5

EIGEN_INC := /opt/local/include/eigen3
CONTACT_INC := /Users/libinlu/Documents/Projects/ves3d/contact3d/inc
CXXFLAGS += -I$(EIGEN_INC) -I$(CONTACT_INC) ${HEDGEHOG_INC} -std=c++11 -fabi-version=6 -fpermissive
#


TEST = RigidTest.exe 
       #ExtenTestData.exe \
       #ShearTestData3.exe \

RUNEXE   = 0
RUNLOG = test_results.txt
LASTRUN  = test_results_last.txt
EXECLOG  = tests.log
EXECNAME =
EXECS   := $(shell ls *.exe)
SHELL=/bin/bash -o pipefail

all: ${TEST}

#build the dependencies for next time. We can have a rule for
#dependency files (%.d: %.cc) but that builds the dependency everytime
%.exe: %.o $(prefix ${VES3D_LIBDIR}, ${ARFILES}) ${MAKE_DEP}
	${CXX} ${CXXFLAGS} ${CPPFLAGS} ${LDFLAGS} $< ${LDLIBS} -o $@
	${CXX} -MM ${CXXFLAGS} ${CPPFLAGS} -c -o $*.d $*.cc
	@if [ "${RUNEXE}" -ne "0" ]; then ./$@; fi

check:  ${TEST}
	${MAKE} run

run-failed: EXECS=$(shell cat ${LASTRUN}|tr '\n' ' ')
run-failed:
	@echo -n '' > ${LASTRUN}
	@${MAKE} run EXECS="${EXECS}"


run: ${EXECS}
	@echo -- TEST TIME: `date` -- list of executed test --------------- >> ${EXECLOG}
	tail -1 ${EXECLOG} >> ${RUNLOG}
	@for exec in ${EXECS}; do				\
	     echo Running $$exec;				\
	     ${MAKE} EXECNAME=$$exec run-exec;			\
	     sleep 1;						\
	done;
	@tail -$(words ${EXECS} "+1") ${RUNLOG}
	@echo -- Failed test to be run with 'make run-failed' -----------------------------------
	@cat ${LASTRUN}

run-exec:
	./${EXECNAME} | tee -a ${EXECLOG}; (test $$? -eq 0 && (echo "passed: ${EXECNAME}" >> ${RUNLOG})) || \
	(echo "failed: ${EXECNAME}" >> ${RUNLOG} && echo "${EXECNAME}" >> ${LASTRUN})


clean:
	-${RM} *.o *.out *.d *.${VES3D_CHKEXT} ${RUNLOG}  ${EXECLOG} *.vtp *.pvtp *.exe

distclean:
	${MAKE} clean
	-${RM} *.exe ${LASTRUN}

.PHONY: clean check run
.SECONDARY:
.SHELLFLAGS = -o pipefail

-include $(TEST:.exe=.d)
