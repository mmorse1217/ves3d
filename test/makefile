EXEC =  LoggerTest.exe        \
	DeviceTest.exe        \
	BlasToyTest.exe       \
	ScalarsTest.exe       \
	BiCGStabTest.exe      \
	SurfaceTest.exe       \
	EvolveSurfaceTest.exe \
	EvolveSurfaceTestMultThread.exe \
#BenchmarkInterfacialVel.exe \
#MovePoleTest.exe #StokesTest.exe 

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
	@for exec in *.exe; do ./$$exec; done

-include $(EXEC:.exe=.d)
