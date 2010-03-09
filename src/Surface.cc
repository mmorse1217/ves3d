/**
 * @file   Surface.cc
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @date   Tue Feb  2 13:37:21 2010
 * 
 * @brief  The implementation of the surface class.
 */

template <typename T> 
Surface<T>::Surface(Device<T> &device_in) :
    device_(device_in), 
    p_(0), 
    n_surfs_(0), 
    x_(device_),
    normal_(device_), 
    h_(device_), 
    w_(device_),
    k_(device_), 
    cu_(device_), 
    cv_(device_),
    xu(device_),
    xv(device_),
    xuv(device_), 
    xvv(device_), 
    xuu(device_),
    bending_force(device_),
    E(device_),
    F(device_),
    G(device_),
    L(device_), 
    M(device_), 
    N(device_), 
    W2(device_),
    temp(device_)
{}

template <typename T> 
Surface<T>::Surface(Device<T> &device_in, int p_in, int n_surfs_in) :
    device_(device_in), 
    p_(p_in), 
    n_surfs_(n_surfs_in),
    x_(device_,p_,n_surfs_),
    normal_(device_, p_,n_surfs_), 
    h_(device_,p_,n_surfs_), 
    w_(device_,p_,n_surfs_),
    k_(device_,p_,n_surfs_), 
    cu_(device_,p_,n_surfs_), 
    cv_(device_,p_,n_surfs_),
    xu(device_,p_,n_surfs_),
    xv(device_,p_,n_surfs_),
    xuv(device_,p_,n_surfs_), 
    xvv(device_,p_,n_surfs_), 
    xuu(device_,p_,n_surfs_),
    bending_force(device_,p_,n_surfs_),
    E(device_,p_,n_surfs_),
    F(device_,p_,n_surfs_),
    G(device_,p_,n_surfs_),
    L(device_,p_,n_surfs_), 
    M(device_,p_,n_surfs_), 
    N(device_,p_,n_surfs_), 
    W2(device_,p_,n_surfs_),
    temp(device_,p_,n_surfs_)
{}

template <typename T> 
Surface<T>::Surface(Device<T> &device_in, int p_in, int n_surfs_in,  
    const Vectors<T> &x_in) :
    device_(device_in), 
    p_(p_in), 
    n_surfs_(n_surfs_in),
    x_(device_,p_,n_surfs_),
    normal_(device_, p_,n_surfs_), 
    h_(device_,p_,n_surfs_), 
    w_(device_,p_,n_surfs_),
    k_(device_,p_,n_surfs_), 
    cu_(device_,p_,n_surfs_), 
    cv_(device_,p_,n_surfs_),
    xu(device_,p_,n_surfs_),
    xv(device_,p_,n_surfs_),
    xuv(device_,p_,n_surfs_), 
    xvv(device_,p_,n_surfs_), 
    xuu(device_,p_,n_surfs_),
    bending_force(device_,p_,n_surfs_),
    E(device_,p_,n_surfs_),
    F(device_,p_,n_surfs_),
    G(device_,p_,n_surfs_),
    L(device_,p_,n_surfs_), 
    M(device_,p_,n_surfs_), 
    N(device_,p_,n_surfs_), 
    W2(device_,p_,n_surfs_),
    temp(device_,p_,n_surfs_)
{
    SetX(x_in);
}

template <typename T> 
void Surface<T>::SetX(const Vectors<T> &x_in)
{
    x_.SetData(x_in.data_);
    UpdateProps();
}

template <typename T> 
void Surface<T>::UpdateProps()
{
    // Calculating the derivative
    T *shc, *work_arr;
    shc = (T*) malloc(6 * p_ * (p_ + 1) * n_surfs_ * sizeof(T));
    work_arr = (T*) malloc(12 * p_ * (p_ + 1) * n_surfs_ * sizeof(T));
    
    device_.AllDerivatives(x_.data_, work_arr, 3*n_surfs_, shc, xu.data_, 
        xv.data_, xuu.data_, xuv.data_, xvv.data_);
    
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
    uyInv(normal_,w_,normal_);
    
    // Second fundamental coefficients 
    DotProduct(xuu,normal_,L);
    DotProduct(xuv,normal_,M);
    DotProduct(xvv,normal_,N);
    
    // Mean curvature
    xy(E,N,h_);
    xy(F,M,temp);
    axpy((T)-2.0,temp,h_,h_);
    xy(G,L, temp);

    axpy((T)1.0,temp,h_,h_);
    axpb((T)0.5,h_,(T)0.0,h_);
    xyInv(h_,W2,h_);

    // Gaussian curvature
    xy(L,N,k_);
    xy(M,M,temp);
    axpy((T)-1.0,temp,k_,k_);
    xyInv(k_,W2,k_);
    
    //Div and Grad coefficients
    xvpb(F,xv,(T) 0.0,cu_);
    axpb((T) -1.0,cu_, (T) 0.0, cu_);
    xvpw(G,xu, cu_, cu_);
    uyInv(cu_,W2,cu_);
    
    xvpb(F,xu,(T) 0.0,cv_);
    axpb((T) -1.0,cv_, (T) 0.0, cv_);
    xvpw(E,xv, cv_, cv_);
    uyInv(cv_,W2,cv_);

     //Bending force
//     SurfGrad(h_, bending_force);
//     SurfDiv(bending_force, E);
//     axpb((T) -1.0,E, (T) 0.0, E);
//     xy(h_,h_,F);
//     axpy((T) -1.0, F, k_, F);
//     xy(F, h_, F);
//     axpy((T) 2.0, F, E, E);
//     xvpb(E, normal_, (T) 0.0, bending_force);
}

template <typename T>
void Surface<T>::StokesMatVec(const Vectors<T> &density_in, Vectors<T> &velocity_out)
{
    int len = density_in.GetDataLength();
    for(int ii=0;ii<len; ii++)
        velocity_out.data_[ii] = density_in.data_[ii];
}

template <typename T> 
void Surface<T>::SurfGrad(
    const Scalars<T> &f_in, Vectors<T> &grad_f_out)
{
    ///@todo The allocation should be moved outside
    Scalars<T> fu(device_, p_, n_surfs_), fv(device_, p_, n_surfs_);
    T *shc, *work_arr;
    shc = (T*) malloc(6 * p_ * (p_ + 1) * n_surfs_ * sizeof(T));
    work_arr = (T*) malloc(12 * p_ * (p_ + 1) * n_surfs_ * sizeof(T));
    
    device_.FirstDerivatives(f_in.data_, work_arr, n_surfs_, shc, fu.data_, fv.data_);
    
    xvpb(fu,cu_, (T) 0.0, grad_f_out);
    xvpw(fv,cv_, grad_f_out, grad_f_out);
}

template <typename T> 
void Surface<T>::SurfDiv(
    const Vectors<T> &f_in, Scalars<T> &div_f_out) 
{
    ///@todo The allocation should be moved outside
    Vectors<T> fu(device_, p_, n_surfs_), fv(device_, p_, n_surfs_);
    T *shc, *work_arr;
    shc = (T*) malloc(6 * p_ * (p_ + 1) * n_surfs_ * sizeof(T));
    work_arr = (T*) malloc(12 * p_ * (p_ + 1) * n_surfs_ * sizeof(T));
    
    device_.FirstDerivatives(f_in.data_, work_arr, 3*n_surfs_, shc, fu.data_, fv.data_);
    
    DotProduct(fu,cu_, temp);
    DotProduct(fv,cv_, div_f_out);
    axpy((T)1.0,div_f_out,temp, div_f_out);
}
