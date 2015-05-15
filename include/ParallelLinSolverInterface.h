/**
 * @file
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @revision $Revision$
 * @tags $Tags$
 * @date $Date$
 *
 * @brief Base class for parallel linear solvers
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

#ifndef _PARALLELLINSOLVERINTERFACE_H_
#define _PARALLELLINSOLVERINTERFACE_H_

#include <cstddef> //size_t
#include "Error.h"

#ifdef HAS_MPI
#include "mpi.h"
#define  comm_t MPI_Comm
#else
#define  comm_t int
#endif /* HAS_MPI */

template<typename T>
class ParallelVec
{
  public:
    typedef std::size_t size_type;
    typedef T value_type;
    typedef T* iterator;
    typedef const T* const_iterator;
    enum NormType {NORM_1=0, NORM_2=1, NORM_FROBENIUS=2, NORM_INFINITY=3, NORM_1_AND_2=4};

    // setup
    explicit ParallelVec();
    virtual Error_t SetSizes(size_type lsz, size_type gsz=0) = 0;
    virtual Error_t Configure() = 0;
    virtual Error_t SetName(const char *name) = 0;

    // management
    virtual Error_t GetSizes(size_type &lsz, size_type &gsz) const = 0;

    virtual Error_t GetArray(iterator       &i, size_t &lsz) = 0;
    virtual Error_t GetArray(const_iterator &i, size_t &lsz) const = 0;

    virtual Error_t RestoreArray(iterator       &i) = 0;
    virtual Error_t RestoreArray(const_iterator &i) const = 0;

    virtual Error_t SetValuesLocal(const_iterator const arr) = 0;
    virtual Error_t AssemblyBegin() = 0;
    virtual Error_t AssemblyEnd() = 0;

    // factories
    virtual Error_t VecFactory(ParallelVec **newvec) const = 0;
    virtual Error_t ReplicateTo(ParallelVec **other) const = 0;

    // utility
    virtual Error_t Norm(value_type &nrm, const enum NormType &type = NORM_2) const = 0;
    virtual Error_t View() const = 0;
    virtual Error_t axpy(value_type a, const ParallelVec *x) = 0;
    virtual const comm_t* MPIComm() const = 0;

    // destruction
    virtual ~ParallelVec() = 0;

  private:
    //forbid copy and assignment
    ParallelVec(ParallelVec const& rhs);
    ParallelVec& operator=(const ParallelVec& rhs);
};

template<typename T>
ParallelVec<T>::ParallelVec(){}

template<typename T>
ParallelVec<T>::~ParallelVec(){}

template<typename T>
class ParallelLinOp
{
  public:
    typedef std::size_t size_type;
    typedef T value_type;
    typedef T* iterator;
    typedef const T* const_iterator;
    typedef ParallelVec<T> vec_type;
    typedef Error_t (*apply_type)(const ParallelLinOp*, const_iterator, iterator);

    // setup
    explicit ParallelLinOp();
    virtual Error_t SetSizes(size_type lrsz, size_type lcsz, size_type grsz=0, size_type gcsz=0) = 0;
    virtual Error_t SetName(const char *name) = 0;
    virtual Error_t SetContext(const void *ctx) = 0;
    virtual Error_t SetApply(apply_type app) = 0;
    virtual Error_t Configure() = 0;

    // management
    virtual Error_t GetSizes(size_type &lrsz, size_type &lcsz, size_type &grsz, size_type &gcsz) const = 0;
    virtual Error_t Context(const void **ctx) const = 0;
    virtual Error_t Apply(const_iterator x, iterator y) const = 0;
    virtual Error_t Apply(const vec_type *x, vec_type *y) const = 0;

    // factories
    virtual Error_t VecFactory(vec_type **newvec) const = 0;
    virtual Error_t LinOpFactory(ParallelLinOp **newop) const = 0;

    // utility
    virtual const comm_t* MPIComm() const = 0;

    // destruction
    virtual ~ParallelLinOp() = 0;

  private:
    //forbid copy and assignment
    ParallelLinOp(ParallelLinOp const& rhs);
    ParallelLinOp& operator=(const ParallelLinOp& rhs);
};

template<typename T>
ParallelLinOp<T>::ParallelLinOp(){}

template<typename T>
ParallelLinOp<T>::~ParallelLinOp(){}

template<typename T>
class ParallelLinSolver{
  public:
    typedef std::size_t size_type;
    typedef T value_type;
    typedef T* iterator;
    typedef const T* const_iterator;
    typedef ParallelLinOp<T> matvec_type;
    typedef typename ParallelLinOp<T>::vec_type vec_type;
    typedef typename ParallelLinOp<T>::apply_type apply_type;
    static const value_type PLS_DEFAULT = -2.0;
    typedef Error_t (*precond_type)(const ParallelLinSolver*, const_iterator, iterator);

    // setup
    ParallelLinSolver();
    virtual Error_t SetTolerances( value_type rtol = PLS_DEFAULT, value_type abstol = PLS_DEFAULT,
	value_type dtol = PLS_DEFAULT, int maxits = PLS_DEFAULT) = 0;

    virtual Error_t SetOperator(matvec_type *mv) = 0;
    virtual Error_t Operator(matvec_type **mv) const = 0;

    virtual Error_t SetPrecondContext(const void *ctx) = 0;
    virtual Error_t PrecondContext(const void **ctx) const = 0;
    virtual Error_t UpdatePrecond(precond_type precond) = 0;

    virtual Error_t Configure() = 0;
    virtual Error_t InitialGuessNonzero(bool flg) const = 0;

    // factories
    virtual Error_t VecFactory(vec_type **newvec) const = 0;
    virtual Error_t LinOpFactory(matvec_type **newop) const = 0;
    virtual Error_t LinSolverFactory(ParallelLinSolver **newsol) const = 0;

    // application
    virtual Error_t Solve(const vec_type *rhs, vec_type *x) const = 0;

    // utility
    virtual Error_t IterationNumber(size_t &niter) const = 0;
    virtual Error_t ViewReport() const = 0;
    virtual Error_t OperatorResidual(const vec_type *rhs, const vec_type *x, value_type &res) const = 0;
    virtual const comm_t* MPIComm() const = 0;

    // destruction
    virtual ~ParallelLinSolver() = 0;

  private:
    //forbid copy and assignment
    ParallelLinSolver(ParallelLinSolver const& rhs);
    ParallelLinSolver& operator=(const ParallelLinSolver& rhs);
};

template<typename T>
ParallelLinSolver<T>::ParallelLinSolver(){}

template<typename T>
ParallelLinSolver<T>::~ParallelLinSolver(){}

#endif /* _PARALLELLINSOLVERINTERFACE_H_ */
