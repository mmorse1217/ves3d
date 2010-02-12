/**
 * @file   Surface.cc
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @date   Tue Feb  2 13:37:21 2010
 * 
 * @brief  The implementation of the surface class.
 */

template <typename ScalarType> 
Surface<ScalarType>::Surface() :
    p_(0), number_of_surfs_(0){}

template <typename ScalarType> 
Surface<ScalarType>::Surface(int p_in, int number_of_surfs_in) :
    p_(p_in), number_of_surfs_(number_of_surfs_in),
    vector_diff_(p_,3*number_of_surfs_), scalar_diff_(p_,number_of_surfs_)
{
    InitializeAll();
}

template <typename ScalarType> 
Surface<ScalarType>::Surface(int p_in, int number_of_surfs_in, 
    const SHVectors<ScalarType> &x_in) :
    p_(p_in), number_of_surfs_(number_of_surfs_in),
    vector_diff_(p_,3*number_of_surfs_), scalar_diff_(p_,number_of_surfs_)
{
    InitializeAll();
    SetX(x_in);
}

template <typename ScalarType> 
void Surface<ScalarType>::SetX(const SHVectors<ScalarType>& x_in)
{
    x_.SetData(x_in.data_);
    UpdateProps();
}

template <typename ScalarType> 
void Surface<ScalarType>::UpdateProps()
{
    // Calculating the derivative
    vector_diff_.AllDerivatives(x_,xu,xv,xuu,xuv,xvv);

    // First fundamental coefficients
    DotProduct(xu,xu,E);
    DotProduct(xu,xv,F);
    DotProduct(xv,xv,G);
    CrossProduct(xu,xv,normal_);
    
    // Area element
    DotProduct(normal_,normal_,W2);
    w_.SetData(W2.data_);
    w_.Sqrt();

    // Normalizing the normal vector
    xDy(normal_,w_,normal_);
    
    // Second fundamental coefficients 
    DotProduct(xuu,normal_,L);
    DotProduct(xuv,normal_,M);
    DotProduct(xvv,normal_,N);
    
    // Mean curvature
    xTy(E,N,h_);
    xTy(F,M,temp);
    AxPy((ScalarType)-2.0,temp,h_,h_);
    xTy(G,L, temp);
    AxPy((ScalarType)1.0,temp,h_,h_);
    AxPy((ScalarType)0.5,h_,(ScalarType)0.0,h_);
    xDy(h_,W2,h_);

    // Gaussian curvature
    xTy(L,N,k_);
    xTy(M,M,temp);
    AxPy((ScalarType)-1.0,temp,h_,h_);
    xDy(k_,W2,k_);

    //Div and Grad coefficients
    AxPy(F,xv,(ScalarType) 0.0,cu_);
    AxPy((ScalarType) -1.0,cu_, (ScalarType) 0.0, cu_);
    AxPy(G,xu, cu_, cu_);
    xDy(cu_,W2,cu_);
    
    AxPy(F,xu,(ScalarType) 0.0,cv_);
    AxPy((ScalarType) -1.0,cv_, (ScalarType) 0.0, cv_);
    AxPy(E,xv, cv_, cv_);
    xDy(cv_,W2,cv_);
}

template <typename ScalarType> 
void Surface<ScalarType>::InitializeAll()
{
    x_.Resize(p_,number_of_surfs_);
    normal_.Resize(p_,number_of_surfs_); 
    h_.Resize(p_,number_of_surfs_); 
    w_.Resize(p_,number_of_surfs_);
    k_.Resize(p_,number_of_surfs_); 
    cu_.Resize(p_,number_of_surfs_); 
    cv_.Resize(p_,number_of_surfs_);
    xu.Resize(p_,number_of_surfs_);
    xv.Resize(p_,number_of_surfs_);
    xuv.Resize(p_,number_of_surfs_); 
    xvv.Resize(p_,number_of_surfs_); 
    xuu.Resize(p_,number_of_surfs_);
    E.Resize(p_,number_of_surfs_);
    F.Resize(p_,number_of_surfs_);
    G.Resize(p_,number_of_surfs_);
    L.Resize(p_,number_of_surfs_); 
    M.Resize(p_,number_of_surfs_); 
    N.Resize(p_,number_of_surfs_); 
    W2.Resize(p_,number_of_surfs_);
    temp.Resize(p_,number_of_surfs_);
}

template <typename ScalarType> 
void Surface<ScalarType>::SurfGrad(
    const SHScalars<ScalarType>& f_in, 
    SHVectors<ScalarType>& grad_f_out)
{
    SHScalars<ScalarType> fu(p_,number_of_surfs_), fv(p_,number_of_surfs_);
    scalar_diff_.FirstDerivatives(f_in,fu,fv);
    
    AxPy(fu,cu_, (ScalarType) 0.0, grad_f_out);
    AxPy(fv,cv_, grad_f_out, grad_f_out);
}

template <typename ScalarType> 
void Surface<ScalarType>::SurfDiv(
    const SHVectors<ScalarType>& f_in, 
    SHScalars<ScalarType>& div_f_out) 
{
    SHVectors<ScalarType> fu(p_,number_of_surfs_), fv(p_,number_of_surfs_);
    vector_diff_.FirstDerivatives(f_in,fu,fv);
    
    DotProduct(fu,cu_, temp);
    DotProduct(fv,cv_, div_f_out);
    AxPy((ScalarType)1.0,div_f_out,temp, div_f_out);
}
