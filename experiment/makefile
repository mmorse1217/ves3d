ARFILES = $(VES3D_DIR)/lib/libves3d.a

ifeq ($(USE_GPU),yes)
  ARFILES += $(VES3D_DIR)/lib/libves3dgpu.a
endif

MAKEIN = $(VES3D_DIR)/makefile.in.files/makefile.in
include $(MAKEIN)

#build the dependencies for next time. We can have a rule for
#dependency files (%.d: %.cc) but that builds the dependency everytime
%.exe: %.o $(ARFILES)
	$(CXX) $(CPPFLAGS) $< -o $@ $(LDFLAGS)
	$(CXX) -MM $(CPPFLAGS) $*.cc -o $*.d

.PHONY: clean
.SECONDARY:

clean:
	-$(RM) *.o *.exe *.out
