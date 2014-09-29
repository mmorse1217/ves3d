MAKEIN = $(VES3D_DIR)/makefile.in.files/makefile.in
include $(MAKEIN)

BINS = ves3d_seq_direct

# include dependency file for templates
-include $(addprefix $(VES3D_SRC)/, $(addsuffix .d, $(BINS)))

all: lib

.PHONY: docs tags all-tags lib test experiment clean bin
.SECONDARY:

docs: $(MAKEDEP)
	-$(DOX) ./docs/Doxyfile

tags: $(MAKEDEP)
	-$(TAGS) $(VES3D_SRC)/*.cc  $(VES3D_INC)/*.h test/*.cc test/*.h

all-tags: $(MAKEDEP)
	-$(TAGS) $(VES3D_SRC)/* $(VES3D_INC)/* tests/* experiments/* etc/*

lib:
	$(MAKE) -C  ${VES3D_DIR}/lib/

bin: $(addprefix $(VES3D_BIN)/,$(BINS))

$(VES3D_BIN)/%: $(VES3D_SRC)/%.o
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CPPFLAGS) $< -o $@ $(LDFLAGS)
	$(CXX) -MM $(CPPFLAGS) -MT $(VES3D_SRC)/$*.o -o $(VES3D_SRC)/$*.d  $(VES3D_SRC)/$*.cc

test:
	$(MAKE) -C  ${VES3D_DIR}/test/

experiment:
	$(MAKE) -C  ${VES3D_DIR}/experiment/

clean:
	$(MAKE) -C ${VES3D_DIR}/lib clean
	$(MAKE) -C ${VES3D_DIR}/test clean
	$(MAKE) -C ${VES3D_DIR}/experiment clean
	-$(RM) *.o $(VES3D_SRC)/*.o ./docs/latex ./docs/html TAGS
	-$(RM) $(VES3D_SRC)/*.d $(addprefix $(VES3D_BIN)/,$(BINS))
