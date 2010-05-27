template<typename SurfContainer>
void InterfacialForce<SurfContainer>::BendingForce(const SurfContainer &S, 
    Vec &Fb) const
{
    ///@todo resize Fb to match the size of input.
    ///@todo This instantiation is not correct; it assumes to much about the container.
    Sca t1(Fb.getNumSubs(),Fb.getShOrder());
    Sca t2(Fb.getNumSubs(),Fb.getShOrder());
 
    S.Grad(S.h_, Fb);
    S.Div(Fb, t1);
    axpy( -1.0, t1, t1);
    
    xy(S.h_,S.h_,t2);
    axpy(-1.0, t2, S.k_, t2);
    xy(t2, S.h_, t2);
    axpy(2.0, t2, t1, t1);
    xv(t1, S.normal_, Fb);
}

template<typename SurfContainer>
void InterfacialForce<SurfContainer>::TensileForce(const SurfContainer &S, 
    const Sca &tension, Vec &Fs) const
{
    ///@todo resize Fs to match with input
    ///@todo This instantiation is not correct; it assumes to much about the container.
    Vec temp(Fs.getNumSubs(), Fs.getShOrder());
    xv(S.h_, S.normal_, Fs);
    xv(tension, Fs, Fs);
    S.Grad(tension, temp);
    axpy(2.0, Fs, temp, Fs);
}

