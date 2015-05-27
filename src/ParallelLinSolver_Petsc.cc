/**
 * @file
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @revision $Revision$
 * @tags $Tags$
 * @date $Date$
 *
 * @brief Wrapper implementation for PETSC
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

////////////////////////////////////////////////////////////////////////////////
// Parallel Vector Implementation //////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<typename T>
ParallelVecPetsc<T>::ParallelVecPetsc(MPI_Comm &comm) :
    comm_(&comm)
{
    COUTDEBUG("Creating a parallel vector");
    ierr = VecCreate(*comm_, &pv_);
    CHK_PETSC(ierr);
}

template<typename T>
ParallelVecPetsc<T>::ParallelVecPetsc(Vec &pv) :
    pv_(pv)
{
    COUTDEBUG("Creating a parallel vector using existing petsc vec");
    comm_ = new MPI_Comm;
    ierr = PetscObjectGetComm((PetscObject) pv, comm_);
    CHK_PETSC(ierr);
}

template<typename T>
ParallelVecPetsc<T>::~ParallelVecPetsc(){
    const char *name;
    ierr = PetscObjectGetName((PetscObject) pv_, &name); CHK_PETSC(ierr);
    COUTDEBUG("Destroying a parallel vector (name: "<<name<<")");
    ierr = VecDestroy(&pv_);
    CHK_PETSC(ierr);
}

template<typename T>
Error_t ParallelVecPetsc<T>::SetSizes(size_type lsz, size_type gsz)
{
    PetscInt n(lsz), N(PETSC_DECIDE);
    if(gsz>0 || n==0) N=gsz;

    COUTDEBUG("Setting parallel vector sizes to "<<n<<" and "<<N);
    ierr = VecSetSizes(pv_, n, N);
    CHK_PETSC(ierr);
    return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelVecPetsc<T>::GetSizes(size_type &lsz, size_type &gsz) const
{
    PetscInt n, N;
    ierr = VecGetLocalSize(pv_, &n); CHK_PETSC(ierr);
    ierr = VecGetSize(pv_, &N); CHK_PETSC(ierr);
    lsz  = n;
    gsz  = N;
    return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelVecPetsc<T>::Configure()
{
    ierr = VecSetFromOptions(pv_); CHK_PETSC(ierr);
    return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelVecPetsc<T>::SetName(const char *name)
{
    ierr = PetscObjectSetName((PetscObject) pv_, name);
    CHK_PETSC(ierr);

   return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelVecPetsc<T>::GetArray(iterator &i, size_type &lsz)
{
    ierr = VecGetArray(pv_, static_cast<PetscScalar**>(&i));
    CHK_PETSC(ierr);

    PetscInt n;
    ierr = VecGetLocalSize(pv_, &n); CHK_PETSC(ierr);
    lsz = n;
    return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelVecPetsc<T>::GetArray(const_iterator &i, size_type &lsz) const
{
    PetscScalar *arr;
    ierr = VecGetArray(pv_, &arr);
    CHK_PETSC(ierr);
    i = arr;

    PetscInt n;
    ierr = VecGetLocalSize(pv_, &n); CHK_PETSC(ierr);
    lsz = n;
    return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelVecPetsc<T>::RestoreArray(iterator &i)
{
    ierr = VecRestoreArray(pv_, static_cast<PetscScalar**>(&i));
    CHK_PETSC(ierr);
    ASSERT(i==NULL, "iterator isn't reset");

    return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelVecPetsc<T>::RestoreArray(const_iterator &i) const
{
    ierr = VecRestoreArray(pv_, const_cast<PetscScalar**>(&i));
    CHK_PETSC(ierr);
    ASSERT(i==NULL, "iterator isn't reset");

    return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelVecPetsc<T>::Norm(value_type &nrm, const enum base_type::NormType &type) const
{
    ::NormType ptype(static_cast< ::NormType>(type)); //petsc norm types
    ierr = VecNorm(pv_, ptype, static_cast<PetscReal*>(&nrm));
    CHK_PETSC(ierr);
    return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelVecPetsc<T>::View() const
{
    ierr = VecView(pv_, PETSC_VIEWER_STDOUT_WORLD);
    CHK_PETSC(ierr);
    return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelVecPetsc<T>::axpy(value_type a, const base_type *x)
{
    const ParallelVecPetsc* pvp = static_cast<const ParallelVecPetsc*>(x);
    PetscScalar alpha(a);
    ierr = VecAXPY(this->PetscVec(), alpha, pvp->PetscVec());
    CHK_PETSC(ierr);
    return ErrorEvent::Success;
}


template<typename T>
Error_t ParallelVecPetsc<T>::SetValuesLocal(const_iterator const arr)
{
    PetscInt n;
    ierr = VecGetLocalSize(pv_, &n); CHK_PETSC(ierr);

    PetscScalar *viter;
    ierr = VecGetArray(pv_, &viter); CHK_PETSC(ierr);
    for (PetscInt i(0); i<n; ++i) viter[i]=arr[i];
    ierr = VecRestoreArray(pv_, &viter); CHK_PETSC(ierr);

    return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelVecPetsc<T>::AssemblyBegin()
{
    ierr = VecAssemblyBegin(pv_);
    CHK_PETSC(ierr);
    return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelVecPetsc<T>::AssemblyEnd()
{
    ierr = VecAssemblyEnd(pv_);
    CHK_PETSC(ierr);
    return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelVecPetsc<T>::VecFactory(base_type **newvec) const
{
    COUTDEBUG("Making a vec of the same type");
    ASSERT(*newvec==NULL, "Can only replicate to empty pointer");
    ParallelVecPetsc<T> *nv = new ParallelVecPetsc<T>(*comm_);
    *newvec = nv;

    return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelVecPetsc<T>::ReplicateTo(base_type **newvec) const
{
    COUTDEBUG("Replicating the vec and returning the pointer");
    ASSERT(*newvec==NULL, "Can only replicate to empty pointer");
    ParallelVecPetsc<T> *nv = new ParallelVecPetsc<T>(*comm_);
    ierr = VecDuplicate(pv_,&(nv->pv_));CHK_PETSC(ierr);
    *newvec = nv;

    return ErrorEvent::Success;
}

template<typename T>
const MPI_Comm* ParallelVecPetsc<T>::MPIComm() const
{
    return comm_;
}

template<typename T>
Vec& ParallelVecPetsc<T>::PetscVec()
{
    return pv_;
}

template<typename T>
const Vec& ParallelVecPetsc<T>::PetscVec() const
{
    return pv_;
}

////////////////////////////////////////////////////////////////////////////////
// Parallel Operator Implementation ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<typename T>
ParallelLinOpPetsc<T>::ParallelLinOpPetsc(MPI_Comm &comm) :
    comm_(&comm)
{
    COUTDEBUG("Creating a parallel linear operator");
    ierr = MatCreate(*comm_, &pm_);
    CHK_PETSC(ierr);
}

template<typename T>
ParallelLinOpPetsc<T>::~ParallelLinOpPetsc(){
    const char *name;
    ierr = PetscObjectGetName((PetscObject) pm_, &name); CHK_PETSC(ierr);
    COUTDEBUG("Destroying a parallel linear operator (name: "<<name<<")");
    ierr = MatDestroy(&pm_);
    CHK_PETSC(ierr);
}

template<typename T>
Error_t ParallelLinOpPetsc<T>::SetSizes(size_type lrsz, size_type lcsz,
    size_type grsz, size_type gcsz){

    PetscInt m(lrsz), n(lcsz), M(PETSC_DECIDE), N(PETSC_DECIDE);
    if(grsz>0 || m==0) M=grsz;
    if(gcsz>0 || n==0) N=grsz;

    COUTDEBUG("Setting operator sizes to local "
	<<m<<"/"<<n<<" and global "
	<<M<<"/"<<N);
    ierr = MatSetSizes(pm_, m, n, M, N);
    CHK_PETSC(ierr);
    return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelLinOpPetsc<T>::SetName(const char *name)
{
    ierr = PetscObjectSetName((PetscObject) pm_, name);
    CHK_PETSC(ierr);
    return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelLinOpPetsc<T>::SetContext(const void *ctx)
{
    COUTDEBUG("Setting the context to "<<ctx);
    ctx_ = ctx;
    return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelLinOpPetsc<T>::SetApply(apply_type app)
{
    COUTDEBUG("Setting apply function to "<<(void*) app);
    apply_ = app;
    return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelLinOpPetsc<T>::Configure()
{
    COUTDEBUG("Configuring the object");
    ierr = MatSetType(pm_, MATSHELL); CHK_PETSC(ierr);
    ierr = MatShellSetContext(pm_, static_cast<void*>(this)); CHK_PETSC(ierr);
    ierr = MatSetUp(pm_); CHK_PETSC(ierr);

    typedef void(*void_fptr)(void);
    ierr = MatShellSetOperation(pm_, MATOP_MULT, (void_fptr) PetscMatvecWrapper<T>); CHK_PETSC(ierr);
    return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelLinOpPetsc<T>::GetSizes( size_type &lrsz, size_type &lcsz,
    size_type &grsz, size_type &gcsz) const
{
    PetscInt m, n, M, N;
    ierr = MatGetLocalSize(pm_, &m, &n); CHK_PETSC(ierr);
    ierr = MatGetSize(pm_, &M, &N); CHK_PETSC(ierr);
    lrsz = m; lcsz = n;
    grsz = M; gcsz = N;
    return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelLinOpPetsc<T>::Context(const void **ctx) const
{
    *ctx = ctx_;
    return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelLinOpPetsc<T>::Apply(const_iterator x, iterator y) const
{
    return apply_(this, x, y);
}

template<typename T>
Error_t ParallelLinOpPetsc<T>::Apply(const vec_type *x, vec_type *y) const
{
    size_t xsz, ysz;
    size_t lrsz, lcsz, grsz, gcsz;
    const_iterator xi;
    iterator yi;

    CHK(x->GetArray(xi, xsz));
    CHK(y->GetArray(yi, ysz));
    ASSERT(comm_ == x->MPIComm(), "vector x should have the same communicator");
    ASSERT(comm_ == y->MPIComm(), "vector y should have the same communicator");

    ASSERT(!GetSizes(lrsz, lcsz, grsz, gcsz), "");
    ASSERT(lcsz == xsz, "incompatible vector");
    ASSERT(lrsz == ysz, "incompatible vector");

    Error_t retval(ErrorEvent::Success);
    retval = apply_(this, xi, yi);
    x->RestoreArray(xi);
    y->RestoreArray(yi);

    return retval;
}

template<typename T>
Error_t ParallelLinOpPetsc<T>::VecFactory(vec_type **newvec) const
{
    COUTDEBUG("Making a new compatible vec");
    ASSERT(*newvec==NULL, "Can only replicate to empty pointer");
    ParallelVecPetsc<T> *o = new ParallelVecPetsc<T>(*comm_);
    *newvec = o;

    return ErrorEvent::Success;
}

template<typename T>
Error_t  ParallelLinOpPetsc<T>::LinOpFactory(base_type **newop) const
{
    COUTDEBUG("Making a new linear operator");
    ASSERT(*newop==NULL, "Can only replicate to empty pointer");
    ParallelLinOpPetsc<T> *o = new ParallelLinOpPetsc<T>(*comm_);
    *newop = o;

    return ErrorEvent::Success;
}

template<typename T>
const MPI_Comm* ParallelLinOpPetsc<T>::MPIComm() const
{
    return comm_;
}

template<typename T>
Mat& ParallelLinOpPetsc<T>::PetscMat()
{
    return pm_;
}

template<typename T>
const Mat& ParallelLinOpPetsc<T>::PetscMat() const
{
    return pm_;
}

template<typename T>
PetscErrorCode PetscMatvecWrapper(Mat A, Vec x, Vec y){
    PetscErrorCode ierr;
    ParallelLinOpPetsc<T> *ctx;
    typedef typename ParallelLinOpPetsc<T>::iterator iter;

    ierr = MatShellGetContext(A, (void**) &ctx); CHK_PETSC(ierr);

    iter xi, yi;
    ierr = VecGetArray(x,&xi); CHK_PETSC(ierr);
    ierr = VecGetArray(y,&yi); CHK_PETSC(ierr);
    ctx->Apply(xi, yi);

    ierr = VecRestoreArray(x, &xi);
    CHK_PETSC(ierr);
    ASSERT(xi==NULL, "pointer is nulled");

    ierr = VecRestoreArray(y, &yi);
    CHK_PETSC(ierr);
    ASSERT(yi==NULL, "pointer is nulled");

    return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Parallel Solver Implementation //////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<typename T>
ParallelLinSolverPetsc<T>::ParallelLinSolverPetsc(MPI_Comm &comm) :
    comm_(&comm),
    precond_(NULL)
{
    COUTDEBUG("Creating a parallel linear solver");
    ierr = KSPCreate(*comm_, &ps_); CHK_PETSC(ierr);
}

template<typename T>
ParallelLinSolverPetsc<T>::~ParallelLinSolverPetsc()
{
    const char *name;
    ierr = PetscObjectGetName((PetscObject) ps_, &name); CHK_PETSC(ierr);
    COUTDEBUG("Destroying a parallel linear solver (name: "<<name<<")");
    ierr = KSPDestroy(&ps_); CHK_PETSC(ierr);
}

template<typename T>
Error_t  ParallelLinSolverPetsc<T>::SetTolerances(value_type rtol, value_type abstol,
	value_type dtol, int maxits)
{
    COUTDEBUG("Setting linear solver tolerances to rtol="<<rtol
	<<", abstol="<<abstol
	<<", dtol="<<dtol
	<<", maxits="<<maxits
	      );

    ierr = KSPSetTolerances(ps_,
	(rtol	< 0) ? PETSC_DEFAULT : rtol,
	(abstol < 0) ? PETSC_DEFAULT : abstol,
	(dtol   < 0) ? PETSC_DEFAULT : dtol,
	(maxits < 0) ? PETSC_DEFAULT : maxits
			    );

    CHK_PETSC(ierr);
    return ErrorEvent::Success;
}

template<typename T>
Error_t  ParallelLinSolverPetsc<T>::SetOperator(matvec_type *mv)
{
    COUTDEBUG("Setting up the operator for the linear solver");
    ASSERT(comm_==mv->MPIComm(), "Operator should have the same MPI communicator");

    mv_ = static_cast<petsc_matvec_type*>(mv);
#if PETSC_VERSION<35
    ierr = KSPSetOperators(ps_, mv_->PetscMat(), mv_->PetscMat(), SAME_PRECONDITIONER); CHK_PETSC(ierr);
#else
    ierr = KSPSetOperators(ps_, mv_->PetscMat(), mv_->PetscMat()); CHK_PETSC(ierr);
#endif

    return ErrorEvent::Success;
}

template<typename T>
Error_t  ParallelLinSolverPetsc<T>::Operator(matvec_type **mv) const
{
    *mv  = mv_;
    return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelLinSolverPetsc<T>::SetPrecondContext(const void *ctx)
{
    precond_ctx_ = ctx;
    return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelLinSolverPetsc<T>::PrecondContext(const void **ctx) const
{
    *ctx = precond_ctx_;
    return ErrorEvent::Success;
}

template<typename T>
Error_t  ParallelLinSolverPetsc<T>::UpdatePrecond(precond_type precond)
{
    if (precond_ == NULL){
	PC pc;
	ierr = KSPGetPC(ps_, &pc); CHK_PETSC(ierr);
	ierr = PCSetType(pc, PCSHELL); CHK_PETSC(ierr);
	ierr = PCShellSetContext(pc, static_cast<void*>(this)); CHK_PETSC(ierr);
	ierr = PCShellSetApply(pc, PetscPrecondWrapper<T>); CHK_PETSC(ierr);
    } else {
	if (precond == NULL){
	    PC pc;
	    ierr = KSPGetPC(ps_, &pc); CHK_PETSC(ierr);
	    ierr = PCSetType(pc, PCNONE); CHK_PETSC(ierr);
	    ierr = PCShellSetContext(pc, static_cast<void*>(NULL)); CHK_PETSC(ierr);
	    ierr = PCShellSetApply(pc, NULL); CHK_PETSC(ierr);
	}
    }

    precond_ = precond;
    return ErrorEvent::Success;
}

template<typename T>
Error_t  ParallelLinSolverPetsc<T>::Configure()
{
    COUTDEBUG("Configuring the linear solver");
    ierr = KSPSetFromOptions(ps_); CHK_PETSC(ierr);
    ierr = KSPMonitorSet(ps_, PetscKSPMonitor<T>,NULL, NULL);  CHK_PETSC(ierr);
    return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelLinSolverPetsc<T>::InitialGuessNonzero(bool flg) const
{
    ierr = KSPSetInitialGuessNonzero(ps_, (PetscBool) flg);
    return ErrorEvent::Success;
}

template<typename T>
Error_t ParallelLinSolverPetsc<T>::VecFactory(vec_type **newvec) const
{
    COUTDEBUG("Making a new compatible vec");
    ASSERT(*newvec==NULL, "Can only replicate to empty pointer");
    ParallelVecPetsc<T> *o = new ParallelVecPetsc<T>(*comm_);
    *newvec = o;

    return ErrorEvent::Success;
}

template<typename T>
Error_t  ParallelLinSolverPetsc<T>::LinOpFactory(matvec_type **newop) const
{
    COUTDEBUG("Making a new compatible linear operator");
    ASSERT(*newop==NULL, "Can only replicate to empty pointer");
    ParallelLinOpPetsc<T> *o = new ParallelLinOpPetsc<T>(*comm_);
    *newop = o;

    return ErrorEvent::Success;
}

template<typename T>
Error_t  ParallelLinSolverPetsc<T>::LinSolverFactory(base_type **newsol) const
{
    COUTDEBUG("Making a new linear solver");
    ASSERT(*newsol==NULL, "Can only replicate to empty pointer");
    ParallelLinSolverPetsc<T> *o = new ParallelLinSolverPetsc<T>(*comm_);
    *newsol = o;

    return ErrorEvent::Success;
}

template<typename T>
Error_t  ParallelLinSolverPetsc<T>::Solve(const vec_type *rhs, vec_type *x) const
{
    COUTDEBUG("Solving the linear system");
    const petsc_vec_type* rp = static_cast<const petsc_vec_type*>(rhs);
    petsc_vec_type* xp = static_cast<petsc_vec_type*>(x);
    ierr = KSPSolve(ps_, rp->PetscVec(), xp->PetscVec()); CHK_PETSC(ierr);
    return ErrorEvent::Success;
}

template<typename T>
Error_t  ParallelLinSolverPetsc<T>::IterationNumber(size_t &niter) const
{
    PetscInt ni;
    ierr = KSPGetIterationNumber(ps_, &ni); CHK_PETSC(ierr);
    niter = ni;
    return ErrorEvent::Success;
}

template<typename T>
Error_t  ParallelLinSolverPetsc<T>::ViewReport() const
{
    KSPConvergedReason reason;

    ierr = KSPGetConvergedReason(ps_, &reason); CHK_PETSC(ierr);

    if (reason>0){
	WHENCHATTY(PRINTF(*comm_, "KSP converged with reason=%d\n", reason));
    } else {
	PRINTF_ERR(*comm_, "KSP diverged with reason=%d\n", reason);
        return ErrorEvent::DivergenceError;
	ABORT(ierr, "KSP diverged");
    }
    return ErrorEvent::Success;
}

template<typename T>
Error_t  ParallelLinSolverPetsc<T>::OperatorResidual(const vec_type *rhs,
    const vec_type *x, value_type &res) const
{
    const petsc_vec_type* rp = static_cast<const petsc_vec_type*>(rhs);

    petsc_vec_type *val(NULL);
    rp->ReplicateTo( (vec_type**) &val);
    ASSERT(!val->SetName("Temporary matvec holder"), "");
    mv_->Apply(x, val);

    ierr = VecAXPY(val->PetscVec(), -1.0, rp->PetscVec()); CHK_PETSC(ierr);
    ierr = val->Norm(res); CHK_PETSC(ierr);
    delete val;

    return ErrorEvent::Success;
}

template<typename T>
const MPI_Comm* ParallelLinSolverPetsc<T>::MPIComm() const
{
    return comm_;
}

template<typename T>
KSP& ParallelLinSolverPetsc<T>::PetscKSP()
{
    return ps_;
}

template<typename T>
const KSP& ParallelLinSolverPetsc<T>::PetscKSP() const
{
    return ps_;
}

template<typename T>
PetscErrorCode PetscPrecondWrapper(PC P, Vec x, Vec y){
    PetscErrorCode ierr;
    ParallelLinSolverPetsc<T> *ctx;
    typedef typename ParallelLinSolverPetsc<T>::iterator iter;

    ierr = PCShellGetContext(P, (void**) &ctx); CHKERRQ(ierr);
    ASSERT(ctx!=NULL, "context can't be null");
    ASSERT(ctx->precond_!=NULL, "context precond can't be null");

    iter xi, yi;
    ierr = VecGetArray(x,&xi); CHK_PETSC(ierr);
    ierr = VecGetArray(y,&yi); CHK_PETSC(ierr);

    ctx->precond_(ctx, xi, yi);

    ierr = VecRestoreArray(x, &xi);
    CHK_PETSC(ierr);
    ASSERT(xi==NULL, "pointer is nulled");

    ierr = VecRestoreArray(y, &yi);
    CHK_PETSC(ierr);
    ASSERT(yi==NULL, "pointer is nulled");

    return 0;
}


template<typename T>
PetscErrorCode PetscKSPMonitor(KSP K,PetscInt n, PetscReal rnorm, void *dummy){
    WHENCHATTY(PetscPrintf(PETSC_COMM_WORLD,
	    "KSP Residual norm at iteration %D: %14.12e \n",
	    n,
	    rnorm));
  return 0;
}
