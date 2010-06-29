include ${VES3D_DIR}/makefile.in

all: docs tags lib test

.PHONY: docs tags lib test clean

docs: $(MAKEDEP)
	-$(DOX) ./docs/Doxyfile

tags: $(MAKEDEP)
	-$(TAGS) src/* include/*

lib:
	$(MAKE) -C  ${VES3D_DIR}/lib/ all

test:
	$(MAKE) -C  ${VES3D_DIR}/test/ all

clean: 
	-$(RM) -rf *.o ./src/*.o  ./docs/latex ./docs/html  \
	lib/*.a test/*.exe test/*.out TAGS


