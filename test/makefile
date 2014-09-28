EXEC = 	ArrayTest.exe			\
	BenchmarkInterfacialVel.exe	\
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
	SpharmTest.exe			\
	StokesTest.exe			\
        SurfaceTest.exe			\
        Tr1Test.exe			\
        VectorsTest.exe			\


ARFILES = $(VES3D_DIR)/lib/libves3d.a

ifeq ($(USE_GPU),yes)
  ARFILES += $(VES3D_DIR)/lib/libves3dgpu.a
endif

ifeq ($(USE_PVFMM),yes)
  EXEC +=PVFMMInterfaceTest.exe
endif

MAKEIN = $(VES3D_DIR)/makefile.in.files/makefile.in
include $(MAKEIN)

all: $(EXEC)

#build the dependencies for next time. We can have a rule for
#dependency files (%.d: %.cc) but that builds the dependency everytime
%.exe: %.o $(ARFILES) $(MAKEIN)
	$(CXX) $(CPPFLAGS) $< -o $@ $(LDFLAGS)
	$(CXX) -MM $(CPPFLAGS) $*.cc -o $*.d

.PHONY: clean check run
.SECONDARY:

clean:
	-$(RM) *.o *.exe *.out *.d

check:  $(EXEC)
	$(MAKE) run

run:
	@for exec in *.exe; do echo "\n\nRunning " $$exec; ./$$exec; sleep 1;done

-include $(EXEC:.exe=.d)
