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
    bending_force_(device_),
    tensile_force_(device_),
    kappa_(.01),
    rep_ts_(.1),
    max_vel_(.01),
    iter_max_(10),
    E(device_),
    F(device_),
    G(device_),
    L(device_), 
    M(device_), 
    N(device_),
    Fu(device_),
    Fv(device_),
    shc(0),
    work_arr(0)
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
    bending_force_(device_,p_,n_surfs_),
    tensile_force_(device_,p_,n_surfs_),
    kappa_(.01),
    rep_ts_(.1),
    max_vel_(.01),
    iter_max_(10),
    E(device_,p_,n_surfs_),
    F(device_,p_,n_surfs_),
    G(device_,p_,n_surfs_),
    L(device_,p_,n_surfs_), 
    M(device_,p_,n_surfs_), 
    N(device_,p_,n_surfs_),
    Fu(device_,p_,n_surfs_),
    Fv(device_,p_,n_surfs_)
{
    shc      = device_.Malloc(6  * p_ *(p_ + 1) * n_surfs_);
    work_arr = device_.Malloc(12 * p_ *(p_ + 1) * n_surfs_);
}

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
    bending_force_(device_,p_,n_surfs_),
    tensile_force_(device_,p_,n_surfs_),
    kappa_(.01),
    rep_ts_(.1),
    max_vel_(.01),
    iter_max_(10),
    E(device_,p_,n_surfs_),
    F(device_,p_,n_surfs_),
    G(device_,p_,n_surfs_),
    L(device_,p_,n_surfs_), 
    M(device_,p_,n_surfs_), 
    N(device_,p_,n_surfs_),
    Fu(device_,p_,n_surfs_),
    Fv(device_,p_,n_surfs_)
{
    shc      = device_.Malloc(6  * p_ *(p_ + 1) * n_surfs_);
    work_arr = device_.Malloc(12 * p_ *(p_ + 1) * n_surfs_);
    SetX(x_in);
}

template<typename T>
Surface<T>::~Surface()
{
    device_.Free(shc);
    device_.Free(work_arr);
}

template <typename T> 
void Surface<T>::SetX(const Vectors<T> &x_in)
{
    x_.SetData(x_in.data_);
    UpdateProps();
}

template <typename T> 
void Surface<T>::UpdateFirstForms()
{
    // Spherical harmonic coefficient
    device_.ShAna(x_.data_, work_arr, 3*n_surfs_, shc);

    // First derivatives
    device_.ShSynDu(shc, work_arr, 3*n_surfs_, Fu.data_);//Fu is holding xu
    device_.ShSynDv(shc, work_arr, 3*n_surfs_, Fv.data_);//Fv is holding xv

    // First fundamental coefficients
    DotProduct(Fu,Fu, E);
    DotProduct(Fu,Fv, F);
    DotProduct(Fv,Fv, G);
    CrossProduct(Fu,Fv,normal_);

    // Area element
    DotProduct(normal_, normal_, w_);//w = W^2 for now
    
    // Dividing EFG by W^2
    xyInv(E, w_, E);
    xyInv(F, w_, F);
    xyInv(G, w_, G);
    w_.Sqrt();
    uyInv(normal_,w_,normal_);
}

template <typename T> 
void Surface<T>::UpdateAll()
{
    UpdateFirstForms();///@todo This may lead to a bug. Since there is no guarantee that E,F,G hold their values.

    //Div and Grad coefficients
    xvpb(F, Fv, (T) 0.0, cu_);//F*xv
    axpb((T) -1.0,cu_, (T) 0.0, cu_);
    xvpw(G, Fu, cu_, cu_);
     
    xvpb(F, Fu, (T) 0.0, cv_);//F*xu
    axpb((T) -1.0,cv_, (T) 0.0, cv_);
    xvpw(E, Fv, cv_, cv_);

    // Second derivatives
    device_.ShSynDuu(shc, work_arr, 3*n_surfs_, Fu.data_);//Fu is holding xuu
    DotProduct(Fu, normal_, L);
    
    device_.ShSynDuv(shc, work_arr, 3*n_surfs_, Fu.data_);//Fu is holding xuv
    DotProduct(Fu, normal_, M);

    device_.ShSynDvv(shc, work_arr, 3*n_surfs_, Fu.data_);//Fu is holding xvv
    DotProduct(Fu, normal_, N);
    
    // Gaussian curvature
    xy(L,N,k_);
    xy(M,M,h_);
    axpy((T)-1.0,h_,k_,k_);
    xyInv(k_,w_,k_);
    xyInv(k_,w_,k_);

    // Mean curvature
    xy(E,N,h_);
    axpb((T) .5,h_, (T) 0.0,h_);
    
    xy(F,M,N);
    axpy((T)-1.0, N, h_, h_);
    
    xy(G,L,N);
    axpb((T) .5,N, (T) 0.0, N);
    
    axpy((T)1.0 ,h_,N, h_);
    
    //Bending force
    SurfGrad(h_, Fu);
    SurfDiv(Fu, E);
    axpb((T) -1.0,E, (T) 0.0, E);
    xy(h_,h_,F);
    axpy((T) -1.0, F, k_, F);
    xy(F, h_, F);
    axpy((T) 2.0, F, E, E);
    xvpb(E, normal_, (T) 0.0, bending_force_);
    axpb(kappa_,Fu, (T) 0.0, bending_force_);
    
    //Tensile force
    xvpb(h_,normal_,(T) 0.0, tensile_force_);
}

template <typename T>
void Surface<T>::StokesMatVec(const Vectors<T> &density_in, Vectors<T> &velocity_out)
{
    device_.Memcpy(velocity_out.data_, density_in.data_, density_in.GetDataLength(), MemcpyDeviceToDevice);
}

template <typename T> 
void Surface<T>::SurfGrad(
    const Scalars<T> &f_in, Vectors<T> &grad_f_out)
{
    device_.FirstDerivatives(f_in.data_, work_arr, n_surfs_, shc, E.data_, G.data_);
    
    xvpb(E,cu_, (T) 0.0, grad_f_out);
    xvpw(G,cv_, grad_f_out, grad_f_out);
}

///@todo this can be done in place so that we'll only need one extra work space.
template <typename T> 
void Surface<T>::SurfDiv(
    const Vectors<T> &f_in, Scalars<T> &div_f_out) 
{
    device_.FirstDerivatives(f_in.data_, work_arr, 3*n_surfs_, shc, Fu.data_, Fv.data_);
    
    DotProduct(Fu,cu_, E);
    DotProduct(Fv,cv_, div_f_out);
    axpy((T)1.0,div_f_out,E, div_f_out);
}

template <typename T> 
void Surface<T>::Reparam()
{
    bool do_iterate = true;
    int iter = 0;
    T vel;

    //upsample
    while(do_iterate && (iter++ < iter_max_))
    {
        UpdateFirstForms();///@todo shc may have changed
        //device_.ShAna(x_.data_, work_arr, 3*n_surfs_, shc);
        device_.Filter(shc, work_arr, 3*n_surfs_, shc);
        device_.ShSyn( shc, work_arr, 3*n_surfs_, Fu.data_); //Fu holds the filtered version
        
        axpy((T) -1.0, x_, Fu, Fu);
        DotProduct(normal_,Fu, E);
        axpb((T) -1.0,E, (T) 0.0,E);
        xvpw(E, normal_, Fu, Fu);//The correction velocity
        axpy(rep_ts_, Fu, x_, x_);

        DotProduct(Fu,Fu,E);//@todo seems extra because we work with singles now
        vel = E.Max();
        do_iterate = ( vel > (max_vel_ * max_vel_) );

#ifndef NDEBUG
        cout<<"Iteration"<<iter<<", reparam max vel: "<<vel<<endl;
#endif
    }
    //down-sample
}
