/**
 * @file
 * @author Rahimian, Abtin <arahimian@acm.org> 
 * @revision $Revision$
 * @tags $Tags$
 * @date $Date$
 *
 * @brief Parallel linear solvers wrapper for PETSC
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

#ifndef _PARALLELLINSOLVERINTERFACE_PETSC_H_
#define _PARALLELLINSOLVERINTERFACE_PETSC_H_

#include "ParallelLinSolverInterface.h"
#include "ves3d_common.h"
#include "petscksp.h"
#include "Logger.h"

template<typename T>
class ParallelVecPetsc : public ParallelVec<T>
{
  public:
    typedef std::size_t	size_type;
    typedef T value_type;
    typedef T* iterator;
    typedef const T* const_iterator;
    typedef ParallelVec<T> base_type;

    // setup
    explicit ParallelVecPetsc(MPI_Comm &comm);
    explicit ParallelVecPetsc(Vec &pv);

    Error_t SetSizes(size_type lsz, size_type gsz = 0);
    Error_t Configure();
    Error_t SetName(const char *name);

    // management
    Error_t GetSizes(size_type &lsz, size_type &gsz) const;

    Error_t GetArray(iterator &i, size_t &lsz);
    Error_t GetArray(const_iterator &i, size_t &lsz) const;

    Error_t RestoreArray(iterator &i);
    Error_t RestoreArray(const_iterator &i) const;

    Error_t SetValuesLocal(const_iterator const arr);
    Error_t AssemblyBegin();
    Error_t AssemblyEnd();

    // factories
    Error_t VecFactory(base_type **newvec) const;
    Error_t ReplicateTo(base_type **newvec) const;

    // utility
    Error_t Norm(value_type &nrm, const enum base_type::NormType &type=base_type::NORM_2) const;
    Error_t View() const;
    Error_t axpy(value_type a, const base_type *x);

    const MPI_Comm* MPIComm() const;

    // destruction
    ~ParallelVecPetsc();

    // non-virtual
    Vec& PetscVec();
    const Vec& PetscVec() const;

  private:
    mutable PetscErrorCode   ierr;
    MPI_Comm                 *comm_;
    Vec 	             pv_;
};

template<typename T>
class ParallelLinOpPetsc : public ParallelLinOp<T>
{
  public:
    typedef std::size_t	size_type;
    typedef T value_type;
    typedef T* iterator;
    typedef const T* const_iterator;
    typedef ParallelLinOp<T> base_type;
    typedef typename base_type::vec_type vec_type;
    typedef typename base_type::apply_type apply_type;
    typedef ParallelVecPetsc<T>  petsc_vec_type;

    // setup
    explicit ParallelLinOpPetsc(MPI_Comm &comm);
    Error_t SetSizes(size_type lrsz, size_type lcsz, size_type grsz=0, size_type gcsz=0);
    Error_t SetName(const char *name);
    Error_t SetContext(const void *ctx);
    Error_t SetApply(apply_type app);
    Error_t Configure();

    // management
    Error_t GetSizes(size_type &lrsz, size_type &lcsz, size_type &grsz, size_type &gcsz) const;
    Error_t Context(const void **ctx) const;
    Error_t Apply(const_iterator x, iterator y) const;
    Error_t Apply(const vec_type *x, vec_type *y) const;

    // factories
    Error_t VecFactory(vec_type **newvec) const;
    Error_t LinOpFactory(base_type **newop) const;

    // utility
    const MPI_Comm* MPIComm() const;

    // destruction
    ~ParallelLinOpPetsc();

    // non-virtual
    Mat& PetscMat();
    const Mat& PetscMat() const;

  private:
    mutable PetscErrorCode	 ierr;
    MPI_Comm			*comm_;
    Mat				 pm_;
    const void			*ctx_;
    apply_type			 apply_;
};

template<typename T>
PetscErrorCode PetscMatvecWrapper(Mat A, Vec x, Vec y);

template<typename T>
PetscErrorCode PetscPrecondWrapper(PC A, Vec x, Vec y);

template<typename T>
PetscErrorCode PetscKSPMonitor(KSP K,PetscInt n, PetscReal rnorm, void *dummy);

template<typename T>
class ParallelLinSolverPetsc : public ParallelLinSolver<T>
{
  public:
    typedef std::size_t size_type;
    typedef T value_type;
    typedef ParallelLinSolver<T> base_type;
    typedef typename base_type::matvec_type matvec_type;
    typedef typename base_type::vec_type vec_type;
    typedef typename base_type::apply_type apply_type;
    typedef typename base_type::precond_type precond_type;
    typedef ParallelVecPetsc<T>    petsc_vec_type;
    typedef ParallelLinOpPetsc<T>  petsc_matvec_type;

    // setup
    explicit ParallelLinSolverPetsc(MPI_Comm &comm);
    Error_t SetTolerances( value_type rtol, value_type abstol, value_type dtol, int maxits);

    Error_t SetOperator(matvec_type *mv);
    Error_t Operator(matvec_type **mv) const;

    Error_t SetPrecondContext(const void *ctx);
    Error_t PrecondContext(const void **ctx) const;
    Error_t UpdatePrecond(precond_type precond);

    Error_t Configure();
    Error_t InitialGuessNonzero(bool flg) const;

    // factories
    Error_t VecFactory(vec_type **newvec) const;
    Error_t LinOpFactory(matvec_type **newop) const;
    Error_t LinSolverFactory(base_type **newsol) const;

    // application
    Error_t Solve(const vec_type *rhs, vec_type *x) const;

    // utility
    Error_t IterationNumber(size_t &niter) const;
    Error_t ViewReport() const;
    Error_t OperatorResidual(const vec_type *rhs, const vec_type *x, value_type &res) const;
    const MPI_Comm* MPIComm() const;

    // destruction
    ~ParallelLinSolverPetsc();

    // non-virtual
    KSP& PetscKSP();
    const KSP& PetscKSP() const;

  private:
    mutable PetscErrorCode	ierr;
    MPI_Comm               *comm_;
    KSP                     ps_;
    petsc_matvec_type      *mv_;
    const void             *precond_ctx_;
    precond_type            precond_;

    friend PetscErrorCode PetscPrecondWrapper<T>(PC A, Vec x, Vec y);
};

#include "ParallelLinSolver_Petsc.cc"

#endif /* _PARALLELLINSOLVERINTERFACE_PETSC_H_ */
