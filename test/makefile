ifndef VES3D_DIR
$(error "$${VES3D_DIR} environment variable is not set.")
endif

# include rules and flag for compiler/host
VES3D_MKDIR ?= ${VES3D_DIR}/makefile.in.files
include ${VES3D_MKDIR}/makefile.in

LDFLAGS += -L/opt/intel/compilers_and_libraries_2017.1.126/mac/compiler/lib
LDFLAGS += -L/Users/libinlu/Documents/Projects/ves3d/contact3d/lib
LDFLAGS += -L/usr/local/Cellar/fftw/3.3.5/lib
LDFLAGS += -L/usr/local/Cellar/cgal/4.9/lib
LDFLAGS += -L/usr/local/lib
LDLIBS += -lCGAL -lgmp -llapack -liomp5 -liagm

EIGEN_INC := /opt/local/include/eigen3
CONTACT_INC := /Users/libinlu/Documents/Projects/ves3d/contact3d/inc
CXXFLAGS += -I$(EIGEN_INC) -I$(CONTACT_INC) -std=c++11 -DQUIET


TEST = ExtenTestData.exe \
       #ShearTestData3.exe \

RUNEXE   = 0
RUNLOG = test_results.txt
LASTRUN  = test_results_last.txt
EXECLOG  = tests.log
EXECNAME =
EXECS   := $(shell ls *.exe)
SHELL=/bin/bash -o pipefail

all: ${TEST}

#build the dependencies for next time. We can have a rule for
#dependency files (%.d: %.cc) but that builds the dependency everytime
%.exe: %.o $(prefix ${VES3D_LIBDIR}, ${ARFILES}) ${MAKE_DEP}
	${CXX} ${CXXFLAGS} ${CPPFLAGS} ${LDFLAGS} $< ${LDLIBS} -o $@
	${CXX} -MM ${CXXFLAGS} ${CPPFLAGS} -c -o $*.d $*.cc
	@if [ "${RUNEXE}" -ne "0" ]; then ./$@; fi

check:  ${TEST}
	${MAKE} run

run-failed: EXECS=$(shell cat ${LASTRUN}|tr '\n' ' ')
run-failed:
	@echo -n '' > ${LASTRUN}
	@${MAKE} run EXECS="${EXECS}"


run: ${EXECS}
	@echo -- TEST TIME: `date` -- list of executed test --------------- >> ${EXECLOG}
	tail -1 ${EXECLOG} >> ${RUNLOG}
	@for exec in ${EXECS}; do				\
	     echo Running $$exec;				\
	     ${MAKE} EXECNAME=$$exec run-exec;			\
	     sleep 1;						\
	done;
	@tail -$(words ${EXECS} "+1") ${RUNLOG}
	@echo -- Failed test to be run with 'make run-failed' -----------------------------------
	@cat ${LASTRUN}

run-exec:
	./${EXECNAME} | tee -a ${EXECLOG}; (test $$? -eq 0 && (echo "passed: ${EXECNAME}" >> ${RUNLOG})) || \
	(echo "failed: ${EXECNAME}" >> ${RUNLOG} && echo "${EXECNAME}" >> ${LASTRUN})


clean:
	-${RM} *.o *.out *.d *.${VES3D_CHKEXT} ${RUNLOG}  ${EXECLOG} *.vtp *.pvtp

distclean:
	${MAKE} clean
	-${RM} *.exe ${LASTRUN}

.PHONY: clean check run
.SECONDARY:
.SHELLFLAGS = -o pipefail

-include $(TEST:.exe=.d)
