CXX       = g++
CFLAGS    = -Wall -fPIC -fopenmp -O3 -funroll-loops -I${PETSC_DIR}/src/dm/mesh/sieve
INCLUDE   = -I ./include/ -I${PETSC_DIR}/include/
DOX       = doxygen

doc: 
	${DOX} ./docs/Doxyfile

clean: 
	rm -rf ./src/*.o ./src/*.h ./docs/latex ./docs/html 

