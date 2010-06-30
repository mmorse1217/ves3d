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

    Container p, s, t, v, r, rtilde, shat, phat;
    p.replicate(x);
    s.replicate(x);
    t.replicate(x);
    v.replicate(x);
    r.replicate(x);
    rtilde.replicate(x);
    shat.replicate(x);
    phat.replicate(x);

    typename Container::value_type normb = Norm(b);

    A(x, r);
    axpy((typename Container::value_type) -1.0, r, b, r);
    axpy((typename Container::value_type)  0.0, r, r, rtilde);

    normb = (normb == 0.0) ? 1.0 : normb;
  
    if ((resid = Norm(r) / normb) <= tol) {
        tol = resid;
        max_iter = 0;
        return BiCGSSuccess;
    }
    for (int i = 1; i <= max_iter; i++) {

        COUT("\n BiCGStab: iteration = "<<i
            <<"\n           relres    = "<<scientific<<setprecision(4)<<resid<<endl);

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
    Container &x, const Container &b, int &max_iter, typename Container::value_type &tol) const   
{
    typename Container::value_type resid, beta, rho_1;
    typename Container::value_type rho_2(1), alpha(1), omega(1);
    
    Container p, s, t, v, r, rtilde, shat, phat;
    p.replicate(x);
    s.replicate(x);
    t.replicate(x);
    v.replicate(x);
    r.replicate(x);
    rtilde.replicate(x);
        
    typename Container::value_type normb = Norm(b);
    A(x, r);
    axpy((typename Container::value_type) -1.0, r, b, r);
    axpy((typename Container::value_type)  0.0, r, r, rtilde);

    normb = (normb == 0.0) ? 1.0 : normb;
   
    if ((resid = Norm(r) / normb) <= tol) {
        tol = resid;
        max_iter = 0;
        return BiCGSSuccess;
    }
    for (int i = 1; i <= max_iter; i++) {

        COUT("\n BiCGStab: iteration = "<<i
            <<"\n           relres    = "<<scientific<<setprecision(4)<<resid<<endl);
        
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
        
        A(p, v);
        alpha = rho_1 / AlgebraicDot(rtilde, v);
        axpy(-alpha, v, r, s);

        if ((resid = Norm(s)/normb) < tol) {
            axpy(alpha, p, x , x);
            tol = resid;
            max_iter = i;
            return BiCGSSuccess;
        }

        A(s, t);
        omega = AlgebraicDot(t,s) / AlgebraicDot(t,t);
        axpy(alpha, p, x, x);
        axpy(omega, s, x, x);
        
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
