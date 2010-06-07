DOX = doxygen
TAGS = etags

all: docs tags

.PHONY: docs clean tags

docs: 
	$(DOX) ./docs/Doxyfile
tags:
	$(TAGS) src/* include/*

clean: 
	-rm -rf *.o ./src/*.o  ./docs/latex ./docs/html  lib/*.a

include makefile.in.files/makefile.in
include makefile.in.files/makefile.in.common

