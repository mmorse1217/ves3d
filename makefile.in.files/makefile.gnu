###################################################################
## Should only define/adjust CXXFLAGS, CPPFLAGS, LDFLAGS, LDLIBS ##
###################################################################

MAKE_DEP += ${VES3D_MKDIR}/makefile.gnu

CXXFLAGS += -fopenmp -fPIC -fno-exceptions -w

ifeq ($(VES3D_DEBUG),yes)
  CXXFLAGS += -O0 -g -pedantic #-gstabs+
else
  CXXFLAGS += -O3 -finline-functions -funroll-loops -w	\
	      -funsafe-loop-optimizations -msse3 #-malign-double
endif

## for ref, these are the implicit rule gnu make uses for compile/link
# f: f.cc
# 	[CXX] [CXXFLAGS] [CPPFLAGS] [LDFLAGS]  f.cc [LOAD] [LDLIBS] -o f
#
# f.o : n.cc
#       [CXX] [CXXFLAGS] [CPPFLAGS]  -c -o f.o f.cc
#
# f : f.o
# 	[CXX] [LDFLAGS]  f.o [LOAD] [LDLIBS] -o f
