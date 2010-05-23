template<typename Container>
void Identity(Container &in, Container &out)
{
    axpy((typename Container::value_type) 0.0, in, in, out);
}

template<typename Container>
typename Container::value_type Norm(const Container &x)
{
    return(sqrt(Dot(x,x)));
}

std::ostream& operator<<(std::ostream& output, const enum BiCGSReturn &ret)
{
    switch (ret)
    {
        case BiCGSSuccess:
            output<<"BiCGSSucess";
            break;
            
        case MaxIterReached:
            output<<"MaxIterReached";
            break;
            
        case BreakDownRhoZero:
            output<<"BreakDownRhoZero";
            break;

        case BreakDownOmegaZero:
            output<<"BreakDownOmegaZero";
            break;
    }
    
    return output;
}    

/** 
 * BiCGSTAB solves the non-symmetric linear system Ax = b using the
 * Preconditioned BiConjugate Gradient Stabilized method, the template
 * is based on the algorithm given in "Templates for the Solution of
 * Linear Systems: Building Blocks for Iterative Methods, 2nd Edition,
 * SIAM".
 * 
 * *No* compatibility check is performed and it is assumed that the
 * container has a copy ctor Container(const Container& ), and friend
 * functions axpy() and Dot().
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
template<typename Container>
enum BiCGSReturn BiCGStab(void (*MatVec)(Container &in, Container &out),
    Container &x, const Container &b, int &max_iter, 
    typename Container::value_type &tol,
    void (*Precond)(Container &in, Container &out))
{
    typename Container::value_type resid, beta, rho_1, rho_2(1), alpha(1), omega(1);
    Container p(x), s(x), t(x), v(x), r(x), rtilde(x), shat(x), phat(x);

    typename Container::value_type normb = Norm(b);
    MatVec(x, r);
    axpy((typename Container::value_type) -1.0, r, b, r);
    axpy((typename Container::value_type)  0.0, r, r, rtilde);

    normb = (normb == 0.0) ? 1.0 : normb;
  
    if ((resid = Norm(r) / normb) <= tol) {
        tol = resid;
        max_iter = 0;
        return BiCGSSuccess;
    }
    for (int i = 1; i <= max_iter; i++) {
        rho_1 = Dot(rtilde, r);
        if (rho_1 == 0)
        {
            tol = Norm(r) / normb;
            max_iter = i;
            return BreakDownRhoZero;
        }

        if (i == 1)
            axpy((typename Container::value_type) 0.0,r, r, p);
        else {
            beta = (rho_1/rho_2) * (alpha/omega);
            
            axpy(-omega, v, p, p);
            axpy(beta, p, r, p);
        }
        
        Precond(p, phat);

        MatVec(phat, v);
        alpha = rho_1 / Dot(rtilde, v);
        axpy(-alpha, v, r, s);

        if ((resid = Norm(s)/normb) < tol) {
            axpy(alpha, phat, x , x);
            tol = resid;
            max_iter = i;
            return BiCGSSuccess;
        }

        Precond(s, shat);

        MatVec(shat, t);
        omega = Dot(t,s) / Dot(t,t);
        axpy(alpha, phat, x, x);
        axpy(omega, shat, x, x);
        
        axpy(-omega, t, s, r);

        rho_2 = rho_1;
        if ((resid = Norm(r) / normb) < tol) {
            tol = resid;
            max_iter = i;
            return BiCGSSuccess;
        }
        if (omega == 0) {
            tol = Norm(r) / normb;
            max_iter = i;
            return BreakDownOmegaZero;
        }
    }

    tol = resid;
    return MaxIterReached;
}
