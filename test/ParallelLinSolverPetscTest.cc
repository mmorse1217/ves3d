/**
 * @file
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @revision $Revision$
 * @tags $Tags$
 * @date $Date$
 *
 * @brief unit test
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

#include "ParallelLinSolverPetscTest.h"
#include "ParallelLinSolver_Petsc.h"

Error_t matvec(const ParallelLinOp<real_t> *M, const real_t *x, real_t *y)
{
    size_t lrsz, lcsz, grsz, gcsz, mn, mx;
    CHK(M->GetSizes(lrsz, lcsz, grsz, gcsz));
    int rank;
    MPI_Comm_rank(*M->MPIComm(), &rank);
    COUTDEBUG("my matvec local row="<<lrsz<<", local column size="<<lcsz);

    mn = (lrsz < lcsz ) ? lrsz : lcsz;
    mx = (lrsz > lcsz ) ? lrsz : lcsz;
    for (int i=0; i<mn; ++i)
    	y[i] = (rank*mx+i+1)*x[i];

    return ErrorEvent::Success;
}

Error_t precond(const ParallelLinSolver<real_t> *KSP, const real_t *x, real_t *y)
{
    typedef ParallelLinSolver<real_t> PSol;
    typename PSol::matvec_type *M;
    CHK(KSP->Operator(&M));

    size_t lrsz, lcsz, grsz, gcsz, mn, mx;
    CHK(M->GetSizes(lrsz, lcsz, grsz, gcsz));
    int rank;
    MPI_Comm_rank(*M->MPIComm(), &rank);
    COUTDEBUG("my precond local row="<<lrsz<<", local column size="<<lcsz);

    mn = (lrsz < lcsz ) ? lrsz : lcsz;
    mx = (lrsz > lcsz ) ? lrsz : lcsz;
    for (int i=0; i<mn; ++i)
    	y[i] = x[i]/(rank*mx+i+1);

    return ErrorEvent::Success;
}

int main(int argc, char **argv){
    PetscInitialize(&argc, &argv, NULL, NULL);

    {   // The vector
	ParallelVecPetsc<real_t> x(VES3D_COMM_WORLD);
	TestParallelVec(&x, 4);
    }

    {   // The matrix
        ParallelLinOpPetsc<real_t> A(VES3D_COMM_WORLD);
        TestParallelOp(&A, matvec, 6);
    }

    {   // The solver
        ParallelLinSolverPetsc<real_t> ksp(VES3D_COMM_WORLD);
        TestParallelSolver(&ksp, matvec, precond, 10);
    }

    PetscFinalize();
}
