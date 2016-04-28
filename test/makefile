ifndef VES3D_DIR
$(error "$${VES3D_DIR} environment variable is not set.")
endif

# include rules and flag for compiler/host
VES3D_MKDIR ?= ${VES3D_DIR}/makefile.in.files
include ${VES3D_MKDIR}/makefile.in

TEST = 	ArrayTest.exe			\
	BiCGStabTest.exe		\
	BlasToyTest.exe			\
	DataIOTest.exe			\
	DeviceTest.exe			\
	EnumsTest.exe			\
	ErrorTest.exe			\
	EvolveSurfaceTest.exe		\
	EvolveSurfaceTestMultThread.exe \
	LoggerTest.exe			\
	MovePoleTest.exe		\
	ParametersTest.exe		\
	ParsingTest.exe			\
	RepartitionTest.exe		\
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

ifeq (${VES3D_DEBUG},yes)
  TEST +=BenchmarkInterfacialVel.exe
endif

RUNEXE   = 0
FAILFILE = test_results.txt
EXECLOG  = tests.log
EXECNAME =
SHELL = /bin/bash

all: ${TEST}

#build the dependencies for next time. We can have a rule for
#dependency files (%.d: %.cc) but that builds the dependency everytime
%.exe: %.o $(prefix ${VES3D_LIBDIR}, ${ARFILES}) ${MAKE_DEP}
	${CXX} ${CXXFLAGS} ${CPPFLAGS} ${LDFLAGS} $< ${LDLIBS} -o $@
	${CXX} -MM ${CXXFLAGS} ${CPPFLAGS} -c -o $*.d $*.cc
	@if [ "${RUNEXE}" -ne "0" ]; then ./$@; fi

clean:
	-${RM} *.o *.out *.d *.${VES3D_CHKEXT} ${FAILFILE}  ${EXECLOG}

distclean:
	${MAKE} clean
	-${RM} *.exe

check:  ${TEST}
	${MAKE} run

run:
	@echo -- TEST TIME: `date` --------------- >> ${EXECLOG}
	tail -1 ${EXECLOG} >> ${FAILFILE}
	@for exec in *.exe; do					\
	    echo Running $$exec;				\
	    ${MAKE} EXECNAME=$$exec run-exec;			\
	    sleep 1;						\
	done;
	@tail -20 ${FAILFILE}

run-exec:
	./${EXECNAME} | tee -a ${EXECLOG}; test $$? -eq 0 || echo "${EXECNAME} failed" >> ${FAILFILE}

.PHONY: clean check run
.SECONDARY:
.SHELLFLAGS = -o pipefail

-include $(TEST:.exe=.d)
