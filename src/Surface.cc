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
    p_(p_in), 
    number_of_surfs_(number_of_surfs_in), 
    x_(p_,number_of_surfs_),
    normal_(p_,number_of_surfs_), 
    h_(p_,number_of_surfs_), 
    w_(p_,number_of_surfs_),
    k_(p_,number_of_surfs_), 
    cu_(p_,number_of_surfs_), 
    cv_(p_,number_of_surfs_){}

template <typename ScalarType> 
Surface<ScalarType>::Surface(int p_in, int number_of_surfs_in, 
    const SHVectors<ScalarType> &x_in) :
    p_(p_in), 
    number_of_surfs_(number_of_surfs_in), 
    x_(p_,number_of_surfs_),
    normal_(p_,number_of_surfs_), 
    h_(p_,number_of_surfs_), 
    w_(p_,number_of_surfs_),
    k_(p_,number_of_surfs_),
    cu_(p_,number_of_surfs_),
    cv_(p_,number_of_surfs_)
{
    SetX(x_in);
}

template <typename ScalarType> 
void Surface<ScalarType>::SetX(const SHVectors<ScalarType>& x_in)
{
    x_.SetData(x_in.data_);
    
    /// Setting up space for the derivatives (TO BE REMOVED)
    SHVectors<ScalarType> xu (p_,number_of_surfs_);
    SHVectors<ScalarType> xv (p_,number_of_surfs_);
    SHVectors<ScalarType> xuv(p_,number_of_surfs_);
    SHVectors<ScalarType> xvv(p_,number_of_surfs_);
    SHVectors<ScalarType> xuu(p_,number_of_surfs_);

    // Calculating the derivative
    SphHarm<ScalarType> diffOpt(p_,3*number_of_surfs_);
    diffOpt.Derivatives(x_in,xu,xv,xuu,xuv,xvv);

    // First fundamental coefficients
    SHScalars<ScalarType> E(p_,number_of_surfs_); DotProduct(xu,xu,E);
    SHScalars<ScalarType> F(p_,number_of_surfs_); DotProduct(xu,xv,F);
    SHScalars<ScalarType> G(p_,number_of_surfs_); DotProduct(xv,xv,G);
    CrossProduct(xu,xv,normal_);
    
    // Area element
    SHScalars<ScalarType> W2(p_,number_of_surfs_); DotProduct(normal_,normal_,W2);

    int W_len = W2.GetDataLength();
    ScalarType *W_data = new ScalarType [W_len];
    for(int idx=0;idx<W_len; ++idx)
        W_data[idx] = sqrt(W2.data_[idx]);
    w_.SetData(W_data);
    delete W_data;
    W_data = NULL;

    // Normalizing the normal vector
    xDy(normal_,w_,normal_);
    
    // Second fundamental coefficients 
    SHScalars<ScalarType> L(p_,number_of_surfs_); DotProduct(xuu,normal_,L);
    SHScalars<ScalarType> M(p_,number_of_surfs_); DotProduct(xuv,normal_,M);
    SHScalars<ScalarType> N(p_,number_of_surfs_); DotProduct(xvv,normal_,N);
    
    // Mean curvature
    SHScalars<ScalarType> temp(p_,number_of_surfs_);
    xTy(E,N,h_);
    xTy(F,M,temp);
    AxPy(-(ScalarType) 2.0,temp,h_,h_);
    xTy(G,L, temp);
    AxPy( (ScalarType) 1.0,temp,h_,h_);
    AxPy( (ScalarType) 0.5,h_, (ScalarType) 0.0,h_);
    xDy(h_,W2,h_);
    
    // Gaussian curvature
    xTy(L,N,k_);
    xTy(M,M,temp);
    AxPy( (ScalarType) -1.0, temp,h_,h_);
    xDy(k_,W2,k_);
    
    // Cu = (G.*Xu - F.*Xv)./W2;
    // Cv = (E.*Xv - F.*Xu)./W2;
}
