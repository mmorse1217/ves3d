include makefile.in.files/makefile.in

DOX = doxygen
TAGS = etags

all: docs tags lib test

.PHONY: docs clean tags lib test

docs: $(MakeFiles)
	$(DOX) ./docs/Doxyfile

tags: $(MakeFiles)
	$(TAGS) src/* include/*

lib:
	$(MAKE) -C lib/ all

test:
	$(MAKE) -C test/ all

clean: 
	-@rm -rf *.o ./src/*.o  ./docs/latex ./docs/html  \
	lib/*.a test/*.exe test/*.out TAGS




