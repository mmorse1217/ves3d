template<typename SurfContainer>
void InterfacialForce<SurfContainer>::BendingForce(const SurfContainer &S, 
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
}

template<typename SurfContainer>
void InterfacialForce<SurfContainer>::TensileForce(const SurfContainer &S, 
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

