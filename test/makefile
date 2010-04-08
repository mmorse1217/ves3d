all:ScalarsTest ExplicitStepper VectorsTest SurfaceTest #DeviceTest 

include ${VES3D_DIR}/makefile.in.files/makefile.in

ScalarsTest: 
	$(CXX)  $(CPPFLAGS) $@.cc $^ -o $@.exe  $(LDFLAGS) 

VectorsTest: 
	$(CXX)  $(CPPFLAGS) $@.cc $^ -o $@.exe  $(LDFLAGS) 

SurfaceTest: 
	$(CXX)  $(CPPFLAGS) $@.cc $^ -o $@.exe  $(LDFLAGS) 

SurfaceTestGPU: 
	$(CXX)  $(CPPFLAGS) $@.cc $^ -o $@.exe  $(LDFLAGS) 

DeviceTest: 
	$(CXX)  $(CPPFLAGS) $@.cc $^ -o $@.exe  $(LDFLAGS) 

ExplicitStepper:
	$(CXX) $(CPPFLAGS)  $@.cc $^ -o $@.exe  $(LDFLAGS) 

ExplicitStepperPara:
	$(CXX) $(CPPFLAGS)  $@.cc $^ -o $@.exe  $(LDFLAGS) 

ExplicitStepperGPU:
	$(CXX) $(CPPFLAGS)  $@.cc $^ -o $@.exe  $(LDFLAGS) 

ExplicitStepperAll:
	$(CXX) $(CPPFLAGS)  $@.cc $^ -o $@.exe  $(LDFLAGS) 

SHTransTest:
	$(CXX) $(CPPFLAGS)  $@.cc $^ -o $@.exe  $(LDFLAGS) 

StokesTest:
	$(CXX) $(CPPFLAGS)  $@.cc $^ -o $@.exe  $(LDFLAGS) 

BlasToyTest:
	$(CXX) $(CPPFLAGS)  $@.cc $^ -o $@.exe  $(LDFLAGS) 

.PHONY: clean

clean: 	
	@-rm -f *.o *.exe


#all:ScalarsTest  VectorsTest  DeviceTest AllTests SHTransTest ExplicitStepper
