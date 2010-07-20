/**
 * @file   BiCGStab.h
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @date   Tue May 18 15:23:36 2010
 */

#ifndef _BICGSTAB_H_
#define _BICGSTAB_H_

#include<iostream>
#include "Logger.h"

template<typename Container>
class Identity
{
    void operator()(Container &in, Container &out) const;
};

/// The return type of the BiCGStab function
enum BiCGSReturn {BiCGSSuccess       = 100,
                  MaxIterReached     = 101,
                  BreakDownRhoZero   = 102,
                  BreakDownOmegaZero = 103,
                  RelresIsNan        = 104};
/** 
 * BiCGSTAB solves the non-symmetric linear system Ax = b using the
 * Preconditioned BiConjugate Gradient Stabilized method; the template
 * is based on the algorithm given in "Templates for the Solution of
 * Linear Systems: Building Blocks for Iterative Methods, 2nd Edition,
 * SIAM".
 * 
 * <b>*No*</b> compatibility check is performed and it is assumed that
 * the container has a <code>replicate(const Container& )</code>
 * method, and (friend) functions <code>axpy()</code> and
 * <code>AlgebraicDot()</code>.
 */
template<typename Container, typename MatVec, 
         typename Precond = Identity<Container> >
class BiCGStab
{
  private:    
    typename Container::value_type Norm(const Container &x) const;

    mutable Container p, s, t, v, r, rtilde, shat, phat;
  
  public:
    
    /** 
     * Preconditioned version (uses 8n memory).
     *
     * @param MatVec      The Matvec                                    
     * @param x           The answer and also the initial guess
     * @param b           The right hand side                                          
     * @param max_iter    The maximum number of iterations, in return it 
     *                    holds the number of iterations taken
     * @param tol         Desired tolerance, in return it holds the achieved 
     *                    relative residual           
     * @param Precond     The preconditioner (default is the identity)
     * 
     * @return enum type of type BiCGSReturn.
     */
    BiCGSReturn operator()(const MatVec &A, Container &x, 
        const Container &b, int &max_iter, 
        typename Container::value_type &tol, 
        const Precond &P) const;   
    
    /**
     * Not preconditioned version (We could as well pass Identity as
     * the preconditioner, but this version uses 6n memory, compared
     * to 8n).
     *
     * @param MatVec      The Matvec                                    
     * @param x           The answer and also the initial guess
     * @param b           The right hand side                                          
     * @param max_iter    The maximum number of iterations, in return it 
     *                    holds the number of iterations taken
     * @param tol         Desired tolerance, in return it holds the achieved 
     *                    relative residual           
     * 
     * @return enum type of type BiCGSReturn.
     */
    BiCGSReturn operator()(const MatVec &A, Container &x, 
        const Container &b, int &max_iter, 
        typename Container::value_type &tol) const;   
};

#include "BiCGStab.cc"

#endif //_BICGSTAB_H_
