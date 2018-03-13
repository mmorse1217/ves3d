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

# TODO:add to makefile.col
LDFLAGS += -L/opt/intel/compilers_and_libraries_2018.1.126/mac/compiler/lib
LDFLAGS += -L/opt/intel/compilers_and_libraries_2018.1.126/mac/mkl/lib
LDFLAGS += -L/Users/libinlu/Documents/Projects/ves3d/contact3d/lib
LDFLAGS += -L/usr/local/Cellar/fftw/3.3.5/lib
LDFLAGS += -L/usr/local/Cellar/cgal/4.9/lib
LDFLAGS += -L/usr/local/lib

LDLIBS += -lCGAL -lgmp -llapack -liomp5 -liagm

EIGEN_INC := /opt/local/include/eigen3
CONTACT_INC := /Users/libinlu/Documents/Projects/ves3d/contact3d/inc
CXXFLAGS += -I$(EIGEN_INC) -I$(CONTACT_INC) -std=c++11 -fabi-version=6 -fpermissive
#

# targets of install
VES3D_BINS = ves3d

all: install

lib:
	${MAKE} -C  ${VES3D_LIBDIR}

install: $(addprefix ${VES3D_BINDIR}/,${VES3D_BINS})

test:
	${MAKE} -C  ${VES3D_TSTDIR}

check:
	${MAKE} -C ${VES3D_TSTDIR} check

doc: ${MAKE_DEP}
	${DOX} ${VES3D_DOCDIR}/Doxyfile

tag: ${MAKE_DEP}
	${TAGS} ${VES3D_SRCDIR}/*.cc  ${VES3D_INCDIR}/*.h ${VES3D_TSTDIR}/*.cc ${VES3D_TSTDIR}/*.h

all-tag: ${MAKE_DEP}
	${TAGS} ${VES3D_SRCDIR}/*  ${VES3D_INCDIR}/* ${VES3D_TSTDIR}/*

${VES3D_BINDIR}/%: ${VES3D_SRCDIR}/%.o ${MAKE_DEP} lib
	-@${MKDIR} $(dir $@)
	${CXX} ${CXXFLAGS} ${CPPFLAGS} ${LDFLAGS} $< ${LDLIBS} -o $@
	${CXX} -MM -MT ${VES3D_SRCDIR}/$*.o ${CXXFLAGS} ${CPPFLAGS} -c -o ${VES3D_SRCDIR}/$*.d ${VES3D_SRCDIR}/$*.cc

clean:
	${MAKE} -C ${VES3D_TSTDIR} clean
	${MAKE} -C ${VES3D_LIBDIR} clean
	-${RM} *.o ${VES3D_SRCDIR}/*.o ${VES3D_SRCDIR}/*.d
	-${RM} -r ${VES3D_DOCDIR}/latex ${VES3D_DOCDIR}/html TAGS
	-${RM} *.${VES3D_CHKEXT} ${VES3D_EXPRDIR}/*.${VES3D_CHKEXT} *.vtp *.pvtp

distclean:
	${MAKE} clean
	${MAKE} -C ${VES3D_LIBDIR} distclean
	${MAKE} -C ${VES3D_TSTDIR} distclean
	-${RM} $(addprefix ${VES3D_BINDIR}/,${VES3D_BINS})

# include dependency file for templates (generated when making bin files)
-include $(addprefix ${VES3D_SRCDIR}/, $(addsuffix .d, ${VES3D_BINS}))

help:
	@echo "    Useful compilation flags are:"
	@echo "        VES3D_DEBUG      [yes|no] affects only compiler flags (not printing verbosity)"
	@echo "        VES3D_PROFILE    [yes|no] turns profiling macros on"
	@echo "        VES3D_VERBOSE    [1|2|3]  1=most verbose, 2=normal, 3=quiet"
	@echo "        VES3D_TESTING    [yes|no] equal to DEBUG and PROFILE"
	@echo "        VES3D_VERSION    default to 'hg id -n' (local to repo) can be set as argument if hg is not working on the node"

.PHONY: all lib install test check doc tags all-tags clean distclean
.SECONDARY:
