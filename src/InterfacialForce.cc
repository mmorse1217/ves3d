template<typename SurfContainer>
void InterfacialForce<SurfContainer>::bendingForce(const SurfContainer &S, 
    Vec &Fb) const
{ 
    Sca t1, t2;
    
    Fb.replicate(S.getPosition());
    t1.replicate(S.getPosition());
    t2.replicate(S.getPosition());

    S.grad(S.getMeanCurv(), Fb);
    S.div(Fb, t1);
    axpy(static_cast<typename SurfContainer::value_type>(-1), t1, t1);
    
    xy(S.getMeanCurv(),S.getMeanCurv(), t2);
    axpy(static_cast<typename SurfContainer::value_type>(-1), 
        t2, S.getGaussianCurv(), t2);
    
    xy(t2, S.getMeanCurv(), t2);
    axpy(static_cast<typename SurfContainer::value_type>(2), t2, t1, t1);
    xv(t1, S.getNormal(), Fb);

    axpy(Parameters<typename SurfContainer::value_type>::
        getInstance().bending_modulus, Fb, Fb);
}

template<typename SurfContainer>
void InterfacialForce<SurfContainer>::linearBendingForce(const SurfContainer &S, 
    const Vec &x_new, Vec &Fb) const
{ 
    Sca h_lin, tmp;
    
    Fb.replicate(S.getPosition());
    h_lin.replicate(S.getPosition());
    tmp.replicate(S.getPosition());

    xy(S.getMeanCurv(), S.getMeanCurv(), tmp);
    axpy(static_cast<typename SurfContainer::value_type>(-1), tmp, 
        S.getGaussianCurv(), tmp);
    
    S.linearizedMeanCurv(x_new, h_lin);
    xy(h_lin, tmp, tmp);

    S.grad(h_lin, Fb);
    S.div(Fb, h_lin);
    axpy(static_cast<typename SurfContainer::value_type>(-1), h_lin, h_lin);
    
    axpy(static_cast<typename SurfContainer::value_type>(2), 
        tmp, h_lin, h_lin);

    xv(h_lin, S.getNormal(), Fb);
    
    axpy(Parameters<typename SurfContainer::value_type>::
        getInstance().bending_modulus, Fb, Fb);
}

template<typename SurfContainer>
void InterfacialForce<SurfContainer>::tensileForce(const SurfContainer &S, 
    const Sca &tension, Vec &Fs) const
{
    Vec temp;
    
    temp.replicate(S.getPosition());
    Fs.replicate(S.getPosition());

    xv(S.getMeanCurv(), S.getNormal(), Fs);
    xv(tension, Fs, Fs);
    S.grad(tension, temp);
    axpy(static_cast<typename SurfContainer::value_type>(2), Fs, temp, Fs);
}

