EXEC =  LoggerTest.exe  \
	DeviceTest.exe  \
	BlasToyTest.exe \
	ScalarsTest.exe #BiCGStabTest.exe #EvolveSurfaceTest.exe SurfaceTest.exe SHTransTest.exe #StokesTest.exe TimeStepperTest.exe VectorsTest.exe        

include $(VES3D_DIR)/makefile.in
-include $(EXEC:.exe=.d)

ARFILES = libves3d.a

ifeq ($(USE_GPU),yes)
  ARFILES += libves3dgpu.a 
endif

all: $(EXEC)

%.exe: %.o $(ARFILES) %.d $(MAKEDEP)
	$(CXX) $(CPPFLAGS) $< -o $@ $(LDFLAGS)

%.d: %.cc
	$(CXX) -MM $(CPPFLAGS) $*.cc -o $@ $(LDFLAGS)

.PHONY: clean check

clean: 
	-$(RM) -rf *.o *.exe *.out *.d

check:  $(EXEC)
#call the executables