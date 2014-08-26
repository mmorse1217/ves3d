EXEC = 	BenchmarkInterfacialVel.exe	\
	BiCGStabTest.exe		\
	BlasToyTest.exe			\
	DeviceTest.exe			\
	EnumsTest.exe			\
	EvalVelocityTest.exe		\
	EvolveSurfaceTest.exe		\
	EvolveSurfaceTestMultThread.exe \
	LoggerTest.exe			\
	MovePoleTest.exe		\
	ParsingTest.exe			\
	RepartitionTest.exe		\
	SHTransTest.exe			\
	ScalarsTest.exe			\
	StokesTest.exe			\
        SurfaceTest.exe			\


ARFILES = $(VES3D_DIR)/lib/libves3d.a

ifeq ($(USE_GPU),yes)
  ARFILES += $(VES3D_DIR)/lib/libves3dgpu.a
endif

MAKEIN = $(VES3D_DIR)/makefile.in.files/makefile.in
include $(MAKEIN)

all: $(EXEC)

%.exe: %.o $(ARFILES) %.d $(MAKEIN)
	$(CXX) $(CPPFLAGS) $< -o $@ $(LDFLAGS)

%.d: %.cc
	$(CXX) -MM $(CPPFLAGS) $*.cc -o $@

.PHONY: clean check

clean:
	-$(RM) *.o *.exe *.out *.d

check:  $(EXEC)
	@for exec in *.exe; do echo "\n\nRunning " $$exec; ./$$exec; done

-include $(EXEC:.exe=.d)
