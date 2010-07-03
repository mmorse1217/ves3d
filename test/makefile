EXEC =  LoggerTest.exe     \
	DeviceTest.exe     \
	BlasToyTest.exe    \
	ScalarsTest.exe    \
	BiCGStabTest.exe   \
	SurfaceTest.exe #EvolveSurfaceTest.exe #StokesTest.exe VectorsTest.exe        

ARFILES = $(VES3D_DIR)/lib/libves3d.a

ifeq ($(USE_GPU),yes)
  ARFILES += $(VES3D_DIR)/lib/libves3dgpu.a 
endif

include ${VES3D_DIR}/makefile.in.files/makefile.in

all: $(EXEC)

%.exe: %.o $(ARFILES) %.d $(MAKEIN)
	$(CXX) $(CPPFLAGS) $< -o $@ $(LDFLAGS)

%.d: %.cc
	$(CXX) -MM $(CPPFLAGS) $*.cc -o $@ $(LDFLAGS)

.PHONY: clean check

clean: 
	-$(RM) -rf *.o *.exe *.out *.d

check:  $(EXEC)
	@for exec in *.exe; do ./$$exec; done

-include $(EXEC:.exe=.d)
