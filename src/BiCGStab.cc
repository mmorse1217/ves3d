template<typename Container>
void Identity<Container>::operator()(Container &in, Container &out) const
{
    axpy(1.0, in, out);
}

template<typename Container, typename MatVec, typename Precond>
typename Container::value_type BiCGStab<Container, MatVec, Precond>::Norm(const Container &x) const
{
    return(sqrt(AlgebraicDot(x,x)));
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
        
        case RelresIsNan:
            output<<"RelresIsNan";
            break;

        default:
            CERR("\n The BiCGSReturn type is not recognized. Please update the"
                <<"\n overloaded insertion operator",endl, NULL);
    }
    return output;
}


template<typename Container, typename MatVec, typename Precond>
enum BiCGSReturn BiCGStab<Container, MatVec, Precond>::operator()(const MatVec &A, 
    Container &x, const Container &b, int &max_iter, 
    typename Container::value_type &tol, const Precond &P) const
{
    typename Container::value_type resid, beta, rho_1;
    typename Container::value_type rho_2(1), alpha(1), omega(1);

    p.replicate(x);
    s.replicate(x);
    t.replicate(x);
    v.replicate(x);
    r.replicate(x);
    rtilde.replicate(x);
    shat.replicate(x);
    phat.replicate(x);
    
    typename Container::value_type normb = Norm(b);
    normb = (normb == 0.0) ? 1.0 : normb;
    
    A(x, r);
    axpy((typename Container::value_type) -1.0, r, b, r);
    axpy((typename Container::value_type)  0.0, r, r, rtilde);

    if ((resid = Norm(r) / normb) <= tol) {
            tol = resid;
            max_iter = 0;
            return BiCGSSuccess;
        }
        
        for (int i = 1; i <= max_iter; i++) {
        
            if ( resid != resid )
                return RelresIsNan;

        rho_1 = AlgebraicDot(rtilde, r);
        if (rho_1 == 0)
        {
            tol = Norm(r) / normb;
            max_iter = i;
            return BreakDownRhoZero;
        }

        if (i == 1)
            axpy((typename Container::value_type) 0.0, r, r, p);
        else {
            beta = (rho_1/rho_2) * (alpha/omega);
            
            axpy(-omega, v, p, p);
            axpy(beta, p, r, p);
        }
        
        P(p, phat);

        A(phat, v);
        alpha = rho_1 / AlgebraicDot(rtilde, v);
        axpy(-alpha, v, r, s);

        if ((resid = Norm(s)/normb) < tol) {
            axpy(alpha, phat, x , x);
            tol = resid;
            max_iter = i;
            COUT("\n  BiCGStab:   Iteration = "<<i
                <<"\n                 Relres = "<<scientific<<setprecision(4)<<resid<<endl);
            return BiCGSSuccess;
        }

        P(s, shat);

        A(shat, t);
        omega = AlgebraicDot(t,s) / AlgebraicDot(t,t);
        axpy(alpha, phat, x, x);
        axpy(omega, shat, x, x);
        
        axpy(-omega, t, s, r);

        rho_2 = rho_1;
        if ((resid = Norm(r) / normb) < tol) {
            tol = resid;
            max_iter = i;
            COUT("\n  BiCGStab:   Iteration = "<<i
                <<"\n                 Relres = "<<scientific<<setprecision(4)<<resid<<endl);
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

template<typename Container, typename MatVec, typename Precond>
enum BiCGSReturn BiCGStab<Container, MatVec, Precond>::operator()(const MatVec &A, 
    Container &x, const Container &b, int &num_restart, int &iter_per_restart,
    typename Container::value_type &tol) const   
{
    const int restart_bound(num_restart);
    const int iter_bound(iter_per_restart);
    iter_per_restart = 0;
    const typename Container::value_type tol_in(tol);

    p.replicate(x);
    s.replicate(x);
    t.replicate(x);
    v.replicate(x);
    r.replicate(x);
    rtilde.replicate(x);
        
    typename Container::value_type normb = Norm(b);
    normb = (normb == 0.0) ? 1.0 : normb;

    for (num_restart=0; num_restart<=restart_bound; ++num_restart)
    {
        typename Container::value_type beta, rho_1;
        typename Container::value_type rho_2(1), alpha(1), omega(1);
        A(x, r); 

        axpy((typename Container::value_type) -1.0, r, b, r);
        axpy((typename Container::value_type)  0.0, r, r, rtilde);
        
        if ((tol = Norm(r) / normb) <= tol_in) 
            return BiCGSSuccess;
        
        for (int ii = 0; ii < iter_bound; ++ii) 
        {
            ++iter_per_restart;

            if ( tol != tol )
                return RelresIsNan;

            rho_1 = AlgebraicDot(rtilde, r);
            
            if (rho_1 == 0)
            {
                tol = Norm(r) / normb;
                return BreakDownRhoZero;
            }

            if (ii == 0)
                axpy((typename Container::value_type) 0.0, r, r, p);
            else {
                beta = (rho_1/rho_2) * (alpha/omega);
            
                axpy(-omega, v, p, p);
                axpy(beta, p, r, p);
            }
        
            A(p, v);
            alpha = rho_1 / AlgebraicDot(rtilde, v);
            axpy(-alpha, v, r, s);

            if ((tol = Norm(s)/normb) < tol_in) {
                axpy(alpha, p, x , x);
                COUT("\n  BiCGStab:   Iteration = "<<iter_per_restart<<" ("<<num_restart<<")"
                    <<"\n                 Relres = "<<scientific<<setprecision(4)<<tol<<endl);
                return BiCGSSuccess;
            }
        
            A(s, t);
            omega = AlgebraicDot(t,s) / AlgebraicDot(t,t);
            axpy(alpha, p, x, x);
            axpy(omega, s, x, x);
        
            axpy(-omega, t, s, r);

            rho_2 = rho_1;
            
            COUT("\n  BiCGStab:   Iteration = "<<iter_per_restart<<" ("<<num_restart<<")"
                    <<"\n                 Relres = "<<scientific<<setprecision(4)<<tol<<endl);
            
            if ((tol = Norm(r) / normb) < tol_in)
                return BiCGSSuccess;
            
            if (omega == 0) 
                return BreakDownOmegaZero;
        }
    }
    return MaxIterReached;
}
