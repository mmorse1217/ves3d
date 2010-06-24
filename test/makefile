#include ${VES3D_DIR}/makefile.in.files/makefile.in
include ${VES3D_DIR}/lib/makefile

all:BiCGStabTest.exe BlasToyTest.exe DeviceTest.exe EvolveSurfaceTest.exe \
    LoggerTest.exe ScalarsTest.exe #SHTransTest.exe StokesTest.exe         \
    SurfaceTest.exe TimeStepperTest.exe VectorsTest.exe

%.exe: %.cc libves3d.a $(MakeFiles)
	$(CXX) $(CPPFLAGS) $< -o $@ $(LDFLAGS)

