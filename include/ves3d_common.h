/**
 * @file
 * @author Rahimian, Abtin <arahimian@acm.org> 
 * @revision $Revision$
 * @tags $Tags$
 * @date $Date$
 *
 * @brief
 */

/*
 * Copyright (c) 2014, Abtin Rahimian 
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef _VES3D_COMMON_H_
#define _VES3D_COMMON_H_

#include <stdlib.h> //getenv
#include <string>

#define VES3D_INITIALIZE(argc,argv,file,help)
#define VES3D_FINALIZE()

#ifdef HAS_MPI
#include "mpi.h"
#define VES3D_COMM_WORLD  MPI_COMM_WORLD
#define VES3D_COMM_SELF   MPI_COMM_SELF
#define VES3D_COMM_STDOUT MPI_COMM_WORLD

#undef VES3D_INITIALIZE
#undef VES3D_FINALIZE
#define VES3D_INITIALIZE(argc,argv,file,help) MPI_Init(argc,argv);
#define VES3D_FINALIZE() MPI_Finalize();
#endif //HAS_MPI

#ifdef HAS_PETSC
#include "petscksp.h"
#undef VES3D_COMM_WORLD
#undef VES3D_COMM_SELF
#undef VES3D_COMM_STDOUT
#define VES3D_COMM_WORLD  PETSC_COMM_WORLD
#define VES3D_COMM_SELF   PETSC_COMM_SELF
#define VES3D_COMM_STDOUT PETSC_COMM_WORLD

#define CHK_PETSC(expr) expr;
#define PRINTF(comm, ... ) PetscPrintf(comm, __VA_ARGS__);
#define PRINTF_ERR(comm, ... ) PetscFPrintf(comm, PETSC_STDERR, __VA_ARGS__);
#define ABORT(ierr, msg) SETERRABORT(PETSC_COMM_WORLD, ierr, msg)

#undef VES3D_INITIALIZE
#undef VES3D_FINALIZE
#define VES3D_INITIALIZE(argc,argv,file,help) PetscInitialize(argc, argv, file, help);
#define VES3D_FINALIZE() PetscFinalize();
#endif //HAS_PETSC

#define STR_EXPAND(tok) #tok
#define STR(tok) STR_EXPAND(tok)

const std::string VES3D_PATH(getenv("VES3D_DIR"));
const std::string VERSION(STR(VES3D_VERSION));

#ifdef DOUBLE_PREC
typedef double real_t;
#endif //DOUBLE_PREC

#ifdef SINGLE_PREC
typedef float real_t;
#endif //SINGLE_PREC

#endif //_VES3D_COMMON_H_
