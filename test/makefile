EXEC =  LoggerTest.exe     \
	DeviceTest.exe     \
	BlasToyTest.exe    \
	ScalarsTest.exe    \
	BiCGStabTest.exe   \
	SurfaceTest.exe #EvolveSurfaceTest.exe SHTransTest.exe #StokesTest.exe TimeStepperTest.exe VectorsTest.exe        

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
	@for exec in *.exe; do ./$$exec; done


include $(VES3D_DIR)/makefile.in
-include $(EXEC:.exe=.d)
