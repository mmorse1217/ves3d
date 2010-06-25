#include ${VES3D_DIR}/makefile.in.files/makefile.in
-include ${VES3D_DIR}/lib/makefile

EXEC=BiCGStabTest.exe BlasToyTest.exe DeviceTest.exe EvolveSurfaceTest.exe \
    LoggerTest.exe ScalarsTest.exe #SHTransTest.exe StokesTest.exe         \
    SurfaceTest.exe TimeStepperTest.exe VectorsTest.exe

-include $(EXEC:.exe=.d)

all: $(EXEC)

%.exe: %.o libves3d.a %.d $(MakeFiles)
	$(CXX) $(CPPFLAGS) $< -o $@ $(LDFLAGS)

%.d: 
	$(CXX) -MM $(CPPFLAGS) $*.cc -o $@ $(LDFLAGS)

