########################################################
### Do not add any customization to this file.       ###
### put customizations in the platform or cxx files. ###
########################################################

ifndef VES3D_DIR
$(error "$${VES3D_DIR} environment variable is not set.")
endif

# include rules and flag for compiler/host
VES3D_MKDIR ?= ${VES3D_DIR}/makefile.in.files
include ${VES3D_MKDIR}/makefile.in

# targets of install
VES3D_BINS = ves3d_seq_direct
ifeq (${VES3D_USE_PVFMM},yes)
  VES3D_BINS += ves3d_pvfmm
endif

all: lib

lib:
	${MAKE} -C  ${VES3D_LIBDIR}

install: $(addprefix ${VES3D_BINDIR}/,${VES3D_BINS})

test:
	${MAKE} -C  ${VES3D_TSTDIR}

check:
	${MAKE} -C ${VES3D_TSTDIR} check

docs: ${MAKE_DEP}
	${DOX} ${VES3D_DOCDIR}/Doxyfile

tags: ${MAKE_DEP}
	${TAGS} ${VES3D_SRCDIR}/*.cc  ${VES3D_INCDIR}/*.h ${VES3D_TSTDIR}/*.cc ${VES3D_TSTDIR}/*.h

all-tags: ${MAKE_DEP}
	${TAGS} ${VES3D_SRCDIR}/*  ${VES3D_INCDIR}/* ${VES3D_TSTDIR}/*

${VES3D_BINDIR}/%: ${VES3D_SRCDIR}/%.o ${MAKE_DEP}
	-@${MKDIR} $(dir $@)
	${CXX} ${CXXFLAGS} ${CPPFLAGS} ${LDFLAGS} $< ${LDLIBS} -o $@
	${CXX} -MM -MT ${VES3D_SRCDIR}/$*.o ${CXXFLAGS} ${CPPFLAGS} -c -o ${VES3D_SRCDIR}/$*.d ${VES3D_SRCDIR}/$*.cc

clean:
	${MAKE} -C ${VES3D_LIBDIR} clean
	${MAKE} -C ${VES3D_TSTDIR} clean
	-${RM} *.o ${VES3D_SRCDIR}/*.o ${VES3D_SRCDIR}/*.d ${VES3D_DOCDIR}/latex ${VES3D_DOCDIR}/html TAGS

distclean:
	${MAKE} clean
	-${RM} $(addprefix ${VES3D_BINDIR}/,${VES3D_BINS})

# include dependency file for templates (generated in when making bin files)
-include $(addprefix ${VES3D_SRCDIR}/, $(addsuffix .d, ${VES3D_BINS}))

.PHONY: lib install config test check experiment docs tags all-tags clean distclean
.SECONDARY: