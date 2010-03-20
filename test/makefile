all:ScalarsTest DeviceTest ExplicitStepper VectorsTest

include ${VES3D_DIR}/makefile.in

ScalarsTest: 
	$(CXX)  $(CPPFLAGS) $@.cc $^ -o $@.exe  $(LDFLAGS) 

VectorsTest: 
	$(CXX)  $(CPPFLAGS) $@.cc $^ -o $@.exe  $(LDFLAGS) 

DeviceTest: 
	$(CXX)  $(CPPFLAGS) $@.cc $^ -o $@.exe  $(LDFLAGS) 

ExplicitStepper:
	$(CXX) $(CPPFLAGS)  $@.cc $^ -o $@.exe  $(LDFLAGS) 

.PHONY: clean

clean: 	
	@-rm -f *.o *.out  


#all:ScalarsTest  VectorsTest  DeviceTest AllTests SHTransTest ExplicitStepper
