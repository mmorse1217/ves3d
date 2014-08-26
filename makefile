include $(VES3D_DIR)/makefile.in.files/makefile.in

all: lib

.PHONY: docs tags all-tags lib test experiment clean

docs: $(MAKEDEP)
	-$(DOX) ./docs/Doxyfile

tags: $(MAKEDEP)
	-$(TAGS) src/*.cc  include/*.h test/*.cc test/*.h

all-tags: $(MAKEDEP)
	-$(TAGS) src/* include/* tests/* experiments/* etc/*

lib:
	$(MAKE) -C  ${VES3D_DIR}/lib/

test:
	$(MAKE) -C  ${VES3D_DIR}/test/

experiment:
	$(MAKE) -C  ${VES3D_DIR}/experiment/

clean:
	$(MAKE) -C ${VES3D_DIR}/lib clean
	$(MAKE) -C ${VES3D_DIR}/test clean
	$(MAKE) -C ${VES3D_DIR}/experiment clean
	-$(RM) -rf *.o ./src/*.o  ./docs/latex ./docs/html TAGS
