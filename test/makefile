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
	ParsingTest.exe			\
	RepartitionTest.exe		\
	SHTransTest.exe			\
	ScalarsTest.exe			\
	StokesTest.exe			\
	StokesDoubleLayerTest.exe	\
        SurfaceTest.exe			\
        Tr1Test.exe			\
        VectorsTest.exe			\

ifeq (${VES3D_USE_PVFMM},yes)
  TEST +=PVFMMInterfaceTest.exe \
				 NearSingularTest.exe
endif

ifeq (${VES3D_TESTING},yes)
  TEST +=BenchmarkInterfacialVel.exe
endif

all: ${TEST}

#build the dependencies for next time. We can have a rule for
#dependency files (%.d: %.cc) but that builds the dependency everytime
%.exe: %.o $(prefix ${VES3D_LIBDIR}, ${ARFILES}) ${MAKE_DEP}
	${CXX} ${CXXFLAGS} ${CPPFLAGS} ${LDFLAGS} $< ${LDLIBS} -o $@
	${CXX} -MM ${CXXFLAGS} ${CPPFLAGS} -c -o $*.d $*.cc

clean:
	-${RM} *.o *.exe *.out *.d

check:  ${TEST}
	${MAKE} run

run:
	@for exec in *.exe; do echo "Running " $$exec; ./$$exec;sleep 1;done

.PHONY: clean check run
.SECONDARY:

-include $(TEST:.exe=.d)
