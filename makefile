include $(VES3D_DIR)/makefile.in.files/makefile.in

all: lib

.PHONY: docs TAGS lib test clean

docs: $(MAKEDEP)
	-$(DOX) ./docs/Doxyfile

TAGS: $(MAKEDEP)
	-$(TAGS) src/* include/* 

tags: $(MAKEDEP)
	-$(TAGS) src/*.cc  include/*.h test/*.cc test/*.h

lib:
	$(MAKE) -C  ${VES3D_DIR}/lib/ all

test:
	$(MAKE) -C  ${VES3D_DIR}/test/ all

clean: 
	-$(RM) -rf *.o ./src/*.o  ./docs/latex ./docs/html TAGS \
	./test/*.exe ./test/*.o ./examples/*.exe ./examples/*.o


