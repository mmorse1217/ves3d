template<typename SurfContainer>
InterfacialForce<SurfContainer>::InterfacialForce(const Parameters<value_type> &params) :
    bending_modulus_(params.bending_modulus)
{}

template<typename SurfContainer>
InterfacialForce<SurfContainer>::~InterfacialForce()
{}

template<typename SurfContainer>
void InterfacialForce<SurfContainer>::bendingForce(const SurfContainer &S, 
    Vec_t &Fb) const
{ 
    Fb.replicate(S.getPosition());
    s1.replicate(S.getPosition());
    s2.replicate(S.getPosition());

    S.grad(S.getMeanCurv(), Fb);
    S.div(Fb, s1);
    axpy(static_cast<typename SurfContainer::value_type>(-1), s1, s1);
    
    xy(S.getMeanCurv(),S.getMeanCurv(), s2);
    axpy(static_cast<typename SurfContainer::value_type>(-1), 
        s2, S.getGaussianCurv(), s2);
    
    xy(s2, S.getMeanCurv(), s2);
    axpy(static_cast<typename SurfContainer::value_type>(2), s2, s1, s1);
    xv(s1, S.getNormal(), Fb);

    axpy(bending_modulus_, Fb, Fb);
}

template<typename SurfContainer>
void InterfacialForce<SurfContainer>::linearBendingForce(const SurfContainer &S, 
    const Vec_t &x_new, Vec_t &Fb) const
{     
    Fb.replicate(S.getPosition());
    s1.replicate(S.getPosition());
    s2.replicate(S.getPosition());

    xy(S.getMeanCurv(), S.getMeanCurv(), s2);
    axpy(static_cast<typename SurfContainer::value_type>(-1), 
        S.getGaussianCurv(), s2, s2);
    
    S.linearizedMeanCurv(x_new, s1);
    xy(s1, s2, s2);
    
    S.grad(s1, Fb);
    S.div(Fb, s1);    
    axpy(static_cast<typename SurfContainer::value_type>(2), 
        s2, s1, s1);

    xv(s1, S.getNormal(), Fb);
    axpy(-bending_modulus_, Fb, Fb);
}

template<typename SurfContainer>
void InterfacialForce<SurfContainer>::tensileForce(const SurfContainer &S, 
    const Sca_t &tension, Vec_t &Fs) const
{    
    v1.replicate(S.getPosition());
    Fs.replicate(S.getPosition());

    xv(S.getMeanCurv(), S.getNormal(), Fs);
    xv(tension, Fs, Fs);
    S.grad(tension, v1);
    axpy(static_cast<typename SurfContainer::value_type>(2), Fs, v1, Fs);
}

