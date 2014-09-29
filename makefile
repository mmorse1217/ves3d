MAKEIN = $(VES3D_DIR)/makefile.in.files/makefile.in
include $(MAKEIN)

BINS = ves3d_seq_direct ves3d_pvfmm

all: lib

.PHONY: docs tags all-tags lib install test experiment clean check

docs: $(MAKEDEP)
	-$(DOX) ./docs/Doxyfile

tags: $(MAKEDEP)
	-$(TAGS) $(VES3D_SRC)/*.cc  $(VES3D_INC)/*.h test/*.cc test/*.h

all-tags: $(MAKEDEP)
	-$(TAGS) $(VES3D_SRC)/* $(VES3D_INC)/* tests/* experiments/* etc/*

lib:
	$(MAKE) -C  ${VES3D_DIR}/lib/

install: $(addprefix $(VES3D_BIN)/,$(BINS))

$(VES3D_BIN)/%: $(VES3D_SRC)/%.o
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CPPFLAGS) $< -o $@ $(LDFLAGS)
	$(CXX) -MM $(CPPFLAGS) -MT $(VES3D_SRC)/$*.o -o $(VES3D_SRC)/$*.d  $(VES3D_SRC)/$*.cc

test:
	$(MAKE) -C  ${VES3D_DIR}/test/

check:
	$(MAKE) -C ${VES3D_DIR}/test check

experiment:
	$(MAKE) -C  ${VES3D_DIR}/experiment/

clean:
	$(MAKE) -C ${VES3D_DIR}/lib clean
	$(MAKE) -C ${VES3D_DIR}/test clean
	-$(RM) *.o $(VES3D_SRC)/*.o ./docs/latex ./docs/html TAGS
	-$(RM) $(VES3D_SRC)/*.d $(addprefix $(VES3D_BIN)/,$(BINS))

# include dependency file for templates
-include $(addprefix $(VES3D_SRC)/, $(addsuffix .d, $(BINS)))
