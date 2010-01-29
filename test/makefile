CFLAGS	= 
CPPFLAGS= -I../include/ -I../src/
CXX 	= g++
CFLAGS  = -Wall 


SHScalars:
	$(CXX) $(CFLAGS) $(CPPFLAGS) -o $@.out $@Test.cc 

SHVectors:
	$(CXX) $(CFLAGS) $(CPPFLAGS) -o $@.out $@Test.cc 

clean: 
	rm *.o *.out
