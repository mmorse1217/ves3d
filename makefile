CXX       = g++
CFLAGS    = -Wall -fPIC -fopenmp -O3 -funroll-loops  
INCLUDE   = -I ./include/
DOX       = doxygen

./examples/ex1.o: ./examples/ex1.cpp 
	$(CXX) $(CFLAGS) $(INCLUDE) ./examples/ex1.cpp 

doc: 
	doxygen ./docs/Doxyfile

clean: 
	rm -rf ./src/*.o ./src/*.h ./docs/latex ./docs/html 
