/**
 * @file   Surface.cc
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @date   Tue Feb  2 13:37:21 2010
 * 
 * @brief  The implementation of the surface class.
 */
template <typename ScalarContainer, typename VectorContainer>  
Surface<ScalarContainer, VectorContainer>::Surface(const Vec& x_in) :
    sht_(x_in.getShOrder()), ///@todo make sht_ autonomous
    position_has_changed_outside_(true),
    first_forms_are_stale_(true),
    second_forms_are_stale_(true)
{
    setPosition(x_in);
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::setPosition(const Vec& x_in)
{
    position_has_changed_outside_ = true;
    first_forms_are_stale_ = true;
    second_forms_are_stale_ = true;

    x_.replicate(x_in);
    axpy(static_cast<value_type>(1), x_in, x_);
}

template <typename ScalarContainer, typename VectorContainer>  
VectorContainer& Surface<ScalarContainer, 
                         VectorContainer>::getPositionModifiable()
{
    position_has_changed_outside_ = true;
    first_forms_are_stale_ = true;
    second_forms_are_stale_ = true;

    return(x_);
}

template <typename ScalarContainer, typename VectorContainer>  
const VectorContainer& Surface<ScalarContainer, 
                               VectorContainer>::getPosition() const
{
    return(x_);
}

template <typename ScalarContainer, typename VectorContainer>  
const VectorContainer& Surface<ScalarContainer, 
                               VectorContainer>::getNormal() const
{
    if(first_forms_are_stale_)
        updateFirstForms();
    return(normal_);
}

template <typename ScalarContainer, typename VectorContainer>  
const ScalarContainer& Surface<ScalarContainer, 
                               VectorContainer>::getAreaElement() const
{
    if(first_forms_are_stale_)
        updateFirstForms();
    return(w_);
}

template <typename ScalarContainer, typename VectorContainer>  
const ScalarContainer& Surface<ScalarContainer, 
                               VectorContainer>::getMeanCurv() const
{
    if(second_forms_are_stale_)
        updateAll();
    return(h_);
}

template <typename ScalarContainer, typename VectorContainer>  
const ScalarContainer& Surface<ScalarContainer, 
                               VectorContainer>::getGaussianCurv() const
{
    if(second_forms_are_stale_)
        updateAll();
    return(k_);
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::updateFirstForms() const
{
    if(position_has_changed_outside_)
        checkContainers();
        
    // Spherical harmonic coefficient
    sht_.FirstDerivatives(x_, work_arr, shc, Xu, Xv);
    
    // First fundamental coefficients
    GeometricDot(Xu, Xu, E);
    GeometricDot(Xu, Xv, F);
    GeometricDot(Xv, Xv, G);
    GeometricCross(Xu, Xv, normal_);
    
    // Area element
    GeometricDot(normal_, normal_, w_);//w_ = W^2 for now
    
    // Dividing EFG by W^2
    xyInv(E, w_, E);
    xyInv(F, w_, F);
    xyInv(G, w_, G);
    Sqrt(w_, w_);
    uyInv(normal_, w_, normal_);

    //Div and Grad coefficients
    xv(F, Xv, cu_);
    axpy(static_cast<value_type>(-1), cu_, cu_);
    xvpw(G, Xu, cu_, cu_);
    
    xv(F, Xu, cv_);
    axpy(static_cast<value_type>(-1), cv_, cv_);
    xvpw(E, Xv, cv_, cv_);

    first_forms_are_stale_ = false;
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::updateAll() const
{
    if(first_forms_are_stale_)
        updateFirstForms();
    
    // Second derivatives
    sht_.forward(x_, work_arr, shc);

    sht_.backward_d2u(shc, work_arr, Xu);
    GeometricDot(Xu, normal_, L);
    
    sht_.backward_duv(shc, work_arr, Xu);
    GeometricDot(Xu, normal_, M);

    sht_.backward_d2v(shc, work_arr, Xu);
    GeometricDot(Xu, normal_, N);
    
    // Gaussian curvature
    xy(L,N,k_);
    xy(M,M,h_);
    axpy(static_cast<value_type>(-1), h_, k_, k_);
    xyInv(k_,w_,k_);
    xyInv(k_,w_,k_);
    
    sht_.Filter(k_,work_arr, shc, k_);
    
    // Mean curvature
    xy(E,N,h_);
    axpy(static_cast<value_type>(.5), h_, h_);
    
    xy(F,M,N);
    axpy(static_cast<value_type>(-1), N, h_, h_);
    
    xy(G,L,N);
    axpy(static_cast<value_type>(.5), N, h_, h_); 

    sht_.Filter(h_, work_arr, shc, h_);

    second_forms_are_stale_ = false;
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::grad(const ScalarContainer 
    &f_in, VectorContainer &grad_f_out) const
{
    if(first_forms_are_stale_)
        updateFirstForms();

    sht_.FirstDerivatives(f_in, work_arr, shc, M, N);

    xv(M, cu_, grad_f_out);
    xvpw(N, cv_, grad_f_out, grad_f_out);
    sht_.Filter(grad_f_out, work_arr, shc, grad_f_out);
}

template <typename ScalarContainer, typename VectorContainer> 
void Surface<ScalarContainer, VectorContainer>::div(const VectorContainer 
    &f_in, ScalarContainer &div_f_out) const
{
    if(first_forms_are_stale_)
        updateFirstForms();

    sht_.FirstDerivatives(f_in, work_arr, shc, Xu, Xv);

    GeometricDot(Xu,cu_, N);
    GeometricDot(Xv,cv_, div_f_out);
    axpy(static_cast<value_type>(1),div_f_out, N, div_f_out);

    sht_.Filter(div_f_out,work_arr, shc, div_f_out);
}

template< typename ScalarContainer, typename VectorContainer >
void Surface<ScalarContainer, VectorContainer>::area(ScalarContainer 
    &area_out) const
{
    if(first_forms_are_stale_)
        updateFirstForms();
    
    integrator_(w_, area_out);
}

template< typename ScalarContainer, typename VectorContainer >
void Surface<ScalarContainer, VectorContainer>::volume(ScalarContainer 
    &vol_out) const
{
    if(first_forms_are_stale_)
        updateFirstForms();
    
    GeometricDot(x_,normal_,N);
    axpy(static_cast<value_type>(1)/3, N, N);

    integrator_(N, w_, vol_out);
}

template< typename ScalarContainer, typename VectorContainer >
void Surface<ScalarContainer, VectorContainer>::getCenters(Vec &centers) const
{
    if(first_forms_are_stale_)
        updateFirstForms();

    GeometricDot(x_, x_, N);
    axpy(static_cast<value_type>(.5), N, N);
    xv(N, normal_, Xu);
    
    integrator_(Xu, w_, centers);
    
    M.replicate(centers);
    volume(M);
    uyInv(centers, M, centers);
    M.replicate(x_);
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::checkContainers() const
{
    normal_.replicate(x_);
    w_.replicate(x_);    
    h_.replicate(x_);    
    k_.replicate(x_);    
    cu_.replicate(x_);    
    cv_.replicate(x_);    
    
    ///@todo change these
    shc.replicate(x_);
    work_arr.resize(2 * x_.getNumSubs(), x_.getShOrder(), x_.getGridDim());
    E.replicate(x_);    
    F.replicate(x_);    
    G.replicate(x_);    
    L.replicate(x_);    
    M.replicate(x_);    
    N.replicate(x_);    
    Xu.replicate(x_);    
    Xv.replicate(x_);

    position_has_changed_outside_ = false;
}
