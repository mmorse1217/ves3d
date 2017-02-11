ifndef VES3D_DIR
$(error "$${VES3D_DIR} environment variable is not set.")
endif

# include rules and flag for compiler/host
VES3D_MKDIR ?= ${VES3D_DIR}/makefile.in.files
include ${VES3D_MKDIR}/makefile.in

TEST =		ArrayTest.exe			\
	BiCGStabTest.exe		\
	BlasToyTest.exe			\
	DataIOTest.exe			\
	DeviceTest.exe			\
	EnumsTest.exe			\
	ErrorTest.exe			\
	EvolveSurfaceTest.exe		\
	LoggerTest.exe			\
	MovePoleTest.exe		\
	ParametersTest.exe		\
	ParsingTest.exe			\
	SHTransTest.exe			\
	ScalarsTest.exe			\
	SimulationTest.exe		\
	StokesDoubleLayerTest.exe	\
	StokesTest.exe			\
	StreamableTest.exe 		\
        SurfaceTest.exe			\
        Tr1Test.exe			\
        VectorsTest.exe			\
#	MemoryManagerTest.exe		\

ifeq (${VES3D_USE_PVFMM},yes)
  TEST += PVFMMInterfaceTest.exe	\
	  NearSingularTest.exe
endif

ifeq (${VES3D_USE_PETSC},yes)
  TEST += ParallelLinSolverPetscTest.exe \
          InterfacialVelocityTest.exe
endif

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
	-${RM} *.o *.out *.d *.${VES3D_CHKEXT} ${RUNLOG}  ${EXECLOG} *.vtp *.pvtp

distclean:
	${MAKE} clean
	-${RM} *.exe ${LASTRUN}

.PHONY: clean check run
.SECONDARY:
.SHELLFLAGS = -o pipefail

-include $(TEST:.exe=.d)
