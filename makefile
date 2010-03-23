DOX = doxygen

.PHONY: docs clean

docs: 
	$(DOX) ./docs/Doxyfile

clean: 
	-rm -rf *.o ./src/*.o  ./docs/latex ./docs/html 

include makefile.in.files/makefile.in
include makefile.in.files/makefile.in.common

