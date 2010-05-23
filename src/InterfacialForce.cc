template<typename ScalarContainer, typename VectorContainer,
         template<typename SC, typename VC> class SurfContainer>
void InterfacialForce<ScalarContainer, VectorContainer, 
                      SurfContainer>::BendingForce(const 
                          SurfContainer<ScalarContainer, VectorContainer> &S, 
                          VectorContainer &Fb)
{
    ///@todo resize Fb to match the size of input.
    ///@todo This instantiation is not correct; it assumes to much about the container.
    ScalarContainer t1(Fb.GetShOrder(), Fb.GetNumFuns());
    ScalarContainer t2(Fb.GetShOrder(), Fb.GetNumFuns());
 
    S.Grad(S.h_, Fb);
    S.Div(Fb, t1);
    axpy( -1.0, t1, t1);
    
    xy(S.h_,S.h_,t2);
    axpy(-1.0, t2, S.k_, t2);
    xy(t2, S.h_, t2);
    axpy(2.0, t2, t1, t1);
    xv(t1, S.normal_, Fb);
}

template<typename ScalarContainer, typename VectorContainer,
         template<typename SC, typename VC> class SurfContainer>
void InterfacialForce<ScalarContainer, VectorContainer, 
                      SurfContainer>::TensileForce(const 
                          SurfContainer<ScalarContainer, VectorContainer> &S, 
                          const ScalarContainer &tension, VectorContainer &Fs)
{
    ///@todo resize Fs to match with input
    ///@todo This instantiation is not correct; it assumes to much about the container.
    VectorContainer temp(Fs.GetShOrder(), Fs.GetNumFuns());
    xv(S.h_, S.normal_, Fs);
    xv(tension, Fs, Fs);
    S.Grad(tension, temp);
    axpy(2.0, Fs, temp, Fs);
}

