include ${VES3D_DIR}/makefile.in.files/makefile.in

all: executable

executable: BiCGStabTest.exe BlasToyTest.exe DeviceTest.exe \
	EvolveSurfaceTest.exe LoggerTest.exe ScalarsTest.exe\
	SHTransTest.exe StokesTest.exe SurfaceTest.exe      \
	TimeStepperTest.exe VectorsTest.exe 

%.exe: %.cc
	$(CXX) $(CPPFLAGS) $*.cc -o $*.exe $(LDFLAGS)

.PHONY: clean

clean: 	
	@-rm -f *.o *.exe *.out


