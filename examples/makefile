CFLAGS	        = 
CPPFLAGS        = -I/home/abtin/school/projects/Ves3D/include/
CXX 		= mpicxx

include ${PETSC_DIR}/conf/base

ex1: ex1.o   chkopts
	-${CXX} -o ex1 ex1.o  ${PETSC_SYS_LIB}
	${RM} -f ex1.o

ex2: ex2.o   chkopts
	-${CLINKER} -o ex2 ex2.o  ${PETSC_SYS_LIB}
	${RM} -f ex2.o

include ${PETSC_DIR}/conf/test