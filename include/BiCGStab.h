/**
 * @file   BiCGStab.h
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @date   Tue May 18 15:23:36 2010
 * 
 * @brief BiCGSTAB solves the non-symmetric linear system Ax = b using
 * the Preconditioned BiConjugate Gradient Stabilized method, the
 * template is based on the algorithm given in "Templates for the
 * Solution of Linear Systems: Building Blocks for Iterative Methods,
 * 2nd Edition, SIAM".
 */

#ifndef _BICGSTAB_H_
#define _BICGSTAB_H_

#include<iostream>

/// The return type of the BiCGStab function
enum BiCGSReturn {BiCGSSuccess       = 100,
                  MaxIterReached     = 101,
                  BreakDownRhoZero   = 102,
                  BreakDownOmegaZero = 103};
     
template<typename Container>
void Identity(Container &in, Container &out);

template<typename Container>
typename Container::value_type Norm(const Container &x);

template<typename Container>
enum BiCGSReturn BiCGStab(void (*MatVec)(Container &, Container &),
    Container &x, const Container &b, int &max_iter, 
    typename Container::value_type &tol,
    void (*Precond)(Container &, Container &) = &Identity);

#include "BiCGStab.cc"

#endif //_BICGSTAB_H_
