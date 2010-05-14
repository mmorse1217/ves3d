all:

include ${VES3D_DIR}/makefile.in.files/makefile.in

DeviceTest: 
	$(CXX)  $(CPPFLAGS) $@.cc $^ -o $@.exe  $(LDFLAGS) 

LoggerTest:
	$(CXX)  $(CPPFLAGS) $@.cc $^ -o $@.exe  $(LDFLAGS) 

BlasToyTest:
	$(CXX) $(CPPFLAGS)  $@.cc $^ -o $@.exe  $(LDFLAGS) 

ScalarsTest: 
	$(CXX)  $(CPPFLAGS) $@.cc $^ -o $@.exe  $(LDFLAGS) 

VectorsTest: 
	$(CXX)  $(CPPFLAGS) $@.cc $^ -o $@.exe  $(LDFLAGS) 

SHTransTest:
	$(CXX) $(CPPFLAGS)  $@.cc $^ -o $@.exe  $(LDFLAGS) 

.PHONY: clean

clean: 	
	@-rm -f *.o *.exe

###########################################################
ArrayTest: 
	$(CXX)  $(CPPFLAGS) $@.cc $^ -o $@.exe  $(LDFLAGS) 


SurfaceTest: 
	$(CXX)  $(CPPFLAGS) $@.cc $^ -o $@.exe  $(LDFLAGS) 

SurfaceTestGPU: 
	$(CXX)  $(CPPFLAGS) $@.cc $^ -o $@.exe  $(LDFLAGS) 

ExplicitStepper:
	$(CXX) $(CPPFLAGS)  $@.cc $^ -o $@.exe  $(LDFLAGS) 

ExplicitStepperPara:
	$(CXX) $(CPPFLAGS)  $@.cc $^ -o $@.exe  $(LDFLAGS) 

ExplicitStepperGPU:
	$(CXX) $(CPPFLAGS)  $@.cc $^ -o $@.exe  $(LDFLAGS) 

ExplicitStepperAll:
	$(CXX) $(CPPFLAGS)  $@.cc $^ -o $@.exe  $(LDFLAGS) 

fileOfVes:
	$(CXX) $(CPPFLAGS)  $@.cc $^ -o $@.exe  $(LDFLAGS) 

StokesTest:
	$(CXX) $(CPPFLAGS)  $@.cc $^ -o $@.exe  $(LDFLAGS) 
