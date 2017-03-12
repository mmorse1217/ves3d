/**
 * @file GMRESLinSolver.h
 * @author Lu, Libin <ll1488@nyu.edu>
 * @date Fri Feb 24 11:37:30 2017
 */

#ifndef _GMRESLINSOLVER_H_
#define _GMRESLINSOLVER_H_

#include <stdio.h>
#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"
#include "mkl_service.h"
#include "Error.h"

#define GMRESParSize 128

template<typename T>
class GMRESLinSolver
{
    public:
        // if we use passing function as the callback as following, 
        // typedef int (*apply_type)(const double*, double*)
        // we should save the context of the object of the function
        // const void *ctx_ = (static_cast<const void*>)(this)
        // where this the object containing the function int func(const double*, double*)
        // the other way is use template MatVec, pass the object(InterfacialVelocity)
        // and define operator()(const double*, double*) in the object
        
        typedef Error_t (*JacobiImpApply)(const GMRESLinSolver*, const double*, double*);

        Error_t SetContext(const void *ctx);
        Error_t Context(const void **ctx) const;

        int operator()(JacobiImpApply MV, JacobiImpApply PC, T *computed_solution, T *rhs, 
                T reltol, T abstol, size_t N, int maxIters, int restartIters) const;

    private:
        const void *ctx_;
};

#include "GMRESLinSolver.cc"

#endif // _GMRESLINSOLVER_H_
