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

#include "ParallelLinSolverInterface.h"
#include "Logger.h"
#include "mpi.h"
#include "TestTools.h"

template<typename T>
void TestParallelVec(ParallelVec<T> *pv, size_t sz=10){
    INFO("Testing an instance of ParallelVec");

    typedef ParallelVec<T> PVec;
    typedef typename PVec::size_type  size_type;
    typedef typename PVec::value_type value_type;

    INFO("Testing vec sizes");
    CHK(pv->SetSizes(sz));
    CHK(pv->SetName("Original"));
    CHK(pv->Configure());

    int wsz, rank;
    MPI_Comm_size(*pv->MPIComm(), &wsz);
    MPI_Comm_rank(*pv->MPIComm(), &rank);

    size_type lsz, gsz;
    CHK(pv->GetSizes(lsz, gsz));
    INFO("local size="<<lsz<<", global size="<<gsz<<", world="<<wsz);
    testtools::AssertTrue(lsz==sz, "correct local size", "incorrect local size");
    testtools::AssertTrue(gsz==wsz*sz, "correct global size", "incorrect global size");

    INFO("Testing get array");
    lsz = 0;
    value_type *arr(NULL), *ref(NULL);
    CHK(pv->GetArray(arr, lsz));
    testtools::AssertTrue(arr!=NULL, "non-null array", "null array returned");
    testtools::AssertTrue(lsz==sz, "correct local size", "wrong local size");

    CHK(pv->RestoreArray(arr));
    testtools::AssertTrue(arr==NULL, "null restored array", "non-null restored array");

    INFO("Testing set value");
    value_type *xv = ref = new value_type[lsz];
    for (int i = 0; i<lsz; ++i) xv[i] = rank*lsz+i;
    pv->SetValuesLocal(xv);
    testtools::AssertTrue(xv==ref, "pointer location unmodified", "pointer location modified");

    CHK(pv->GetArray(arr, lsz));
    bool pass(true);
    for (int i=0; i<lsz; ++i)
	pass &= arr[i]==xv[i];

    testtools::AssertTrue(pass, "values set correctly", "incorrect values");
    delete[] ref;
    CHK(pv->RestoreArray(arr));

    CHK(pv->GetArray(arr, lsz));
    pass=true;
    for (int i=0; i<lsz; ++i)
	pass &= arr[i]==rank*lsz+i;
    testtools::AssertTrue(pass, "values set correctly", "incorrect values");
    CHK(pv->RestoreArray(arr));

    INFO("Testing Replication");
    PVec *pvp(NULL);
    pv->ReplicateTo(&pvp);
    pvp->SetName("Copy");

    size_type rlsz, rgsz;
    pvp->GetSizes(rlsz, rgsz);

    testtools::AssertTrue(rlsz==lsz, "correct local size", "incorrect local size");
    testtools::AssertTrue(rgsz==gsz, "correct golbal size", "incorrect global size");
    testtools::AssertTrue(pv->MPIComm()==pvp->MPIComm(), "same mpi communicator", "incorrect mpi communicator");

    INFO("Testing Assembly");
    CHK(pv->AssemblyBegin());
    CHK(pv->AssemblyEnd());

    INFO("Testing utlity functions");
    value_type nrm;
    CHK(pv->Norm(nrm, PVec::NORM_INFINITY));
    INFO("inf norm="<<nrm);
    testtools::AssertTrue(nrm==wsz*lsz-1,"correct inf norm", "incorrect inf norm");

    WHENVERBOSE(pv->View());
    WHENVERBOSE(pvp->View());

    MPI_Barrier(*pv->MPIComm());
    COUTMASTER(emph<<"Prallel vector test passed\n"
	<<"----------------------------------------------------------------------------"<<emph<<std::endl);
}

template<typename T>
Error_t test_matvec_wrapper(const ParallelLinOp<T> *M, const T *x, T *y)
{
    INFO("The test matvec wrapper");
    testtools::AssertTrue(true, "matvec called properly", "");
    typename ParallelLinOp<T>::apply_type ctx(NULL);
    CHK(M->Context((const void**) &ctx));
    testtools::AssertTrue(*ctx!=NULL, "well defined context", "bad context");
    ctx(M,x,y);

    return ErrorEvent::Success;
}

template<typename T>
void TestParallelOp(ParallelLinOp<T> *A,
    typename ParallelLinOp<T>::apply_type matvec,
    size_t sz=10)
{
    INFO("Testing an instance of ParallelLinOp");

    typedef ParallelLinOp<T> POp;
    typedef typename ParallelLinOp<T>::vec_type PVec;
    typedef typename PVec::size_type  size_type;
    typedef typename PVec::value_type value_type;

    // set up and sizes
    INFO("Testing setup");
    size_t rsz(sz), csz(2*sz);
    CHK(A->SetSizes(rsz,csz));
    CHK(A->SetName("A"));
    CHK(A->SetContext((void*) matvec));
    CHK(A->SetApply(test_matvec_wrapper));
    CHK(A->Configure());

    size_type lrsz, lcsz, grsz, gcsz;
    int wsz, rank;

    CHK(A->GetSizes(lrsz, lcsz, grsz, gcsz));
    INFO("local sizes="<<lrsz<<"/"<<lcsz
	<<", global sizes="<<grsz<<"/"<<gcsz);
    MPI_Comm_size(*A->MPIComm(), &wsz);
    MPI_Comm_rank(*A->MPIComm(), &rank);

    testtools::AssertTrue(lrsz==rsz    , "correct local row size"           , "incorrect local row size");
    testtools::AssertTrue(lcsz==csz    , "correct local column size"        , "incorrect local column size");
    testtools::AssertTrue(grsz==wsz*rsz, "correct local global row size"    , "incorrect global row size");
    testtools::AssertTrue(gcsz==wsz*csz, "correct local global column size" , "incorrect global column size");

    const void *ctx;
    CHK(A->Context(&ctx));
    testtools::AssertTrue(ctx==(void*) matvec, "correct context", "incorrect context");

    // Setting up the vec
    PVec *x(NULL), *b(NULL);
    CHK(A->VecFactory(&x));
    CHK(x->SetSizes(csz));
    CHK(x->SetName("x"));
    CHK(x->Configure());
    value_type *xv = new value_type[csz];
    for (size_type i = 0; i<csz; ++i) xv[i] = rank*csz+i;
    x->SetValuesLocal(xv);

    CHK(A->VecFactory(&b));
    CHK(b->SetName("rhs"));
    CHK(b->SetSizes(rsz));
    CHK(b->Configure());

    INFO("Testing Apply");
    CHK(A->Apply(x, b));
    WHENVERBOSE(x->View());
    WHENVERBOSE(b->View());
    MPI_Barrier(*A->MPIComm()); //for nicer print log

    size_type n;
    value_type *arr, *ref;
    CHK(b->GetArray(arr, n));
    testtools::AssertTrue(arr!=NULL && n==lrsz, "correct pointer and size", "bad pointer or size");

    ref = new value_type[n];
    matvec(A,xv,ref);
    bool pass(true);
    for (int i=0; i<n; ++i)
	pass &= arr[i]==ref[i];
    testtools::AssertTrue(pass,"correct matvec", "bad matvec");
    CHK(x->RestoreArray(arr));

    delete[] xv;
    delete[] ref;
    MPI_Barrier(*A->MPIComm()); //for nicer print log

    COUTMASTER(emph<<"Prallel linear operator test passed\n"
	<<"----------------------------------------------------------------------------"<<emph<<std::endl);
}

template<typename T>
void TestParallelSolver(
    ParallelLinSolver<T> *KSP,
    typename ParallelLinSolver<T>::apply_type matvec,
    typename ParallelLinSolver<T>::precond_type precond,
    size_t sz=10,
    int repeat=10)
{
    INFO("Testing an instance of ParallelLinSolver");

    typedef ParallelLinSolver<T> PSol;
    typedef typename PSol::matvec_type POp;
    typedef typename PSol::vec_type PVec;
    typedef typename PVec::size_type size_type;
    typedef typename PVec::value_type value_type;

    int rank;
    MPI_Comm_rank(*KSP->MPIComm(), &rank);

    // Setting up the operator
    POp *A(NULL);
    CHK(KSP->LinOpFactory(&A));
    CHK(A->SetSizes(sz,sz));
    CHK(A->SetName("A"));
    CHK(A->SetApply(matvec));
    CHK(A->Configure());

    // setting up the rhs
    PVec *x(NULL);
    CHK(KSP->VecFactory(&x));
    CHK(x->SetSizes(sz));
    CHK(x->SetName("reference"));

    CHK(x->Configure());
    value_type *xv = new value_type[sz];
    for (size_type i = 0; i<sz; ++i) xv[i] = 1.0;
    x->SetValuesLocal(xv);

    PVec *b(NULL), *u(NULL);
    CHK(x->ReplicateTo(&b));
    CHK(b->SetName("rhs"));
    CHK(A->Apply(x, b));

    CHK(x->ReplicateTo(&u));
    CHK(b->SetName("solution"));

    // setting up the solver
    value_type rtol(1e-12);
    CHK(KSP->SetTolerances(rtol));
    CHK(KSP->SetOperator(A));
    CHK(KSP->Configure());

    // slove
    // preconditioned
    CHK(KSP->UpdatePrecond(precond));
    for (int i=0;i<repeat; ++i)
    {
	CHK(KSP->Solve(b, u));
	size_t niter;
	CHK(KSP->IterationNumber(niter));
	testtools::AssertTrue(niter==1, "one iter w/ exact preconditioner", "bad preconditioner");
    }

    // unpreconditioned
    CHK(KSP->UpdatePrecond(NULL));
    value_type res, err;
    for (int i=0;i<repeat; ++i)
    {
	CHK(KSP->Solve(b, u));
	size_t niter;
	CHK(KSP->IterationNumber(niter));
	INFO("Try"<<i<<": KSP num_iterations="<<niter);
	WHENVERBOSE(KSP->ViewReport());

	CHK(KSP->OperatorResidual(b, u, res));
	INFO("KSP residual="<<res);
    }
    testtools::AssertTrue(res<1e3*rtol, "KSP converged", "KSP failed");

    // error
    CHK(u->axpy(-1.0, x));
    CHK(u->Norm(err));
    INFO("norm="<<err);
    testtools::AssertTrue(err<1e3*rtol, "acceptable error", "large failed");

    MPI_Barrier(*KSP->MPIComm()); //for nicer print log
    COUTMASTER(emph<<"Prallel linear solver test passed\n"
	    <<"----------------------------------------------------------------------------"<<emph<<std::endl);

}
