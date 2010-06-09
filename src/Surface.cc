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
    second_forms_are_stale_(true),
    checked_out_work_sca_(0)
{
    setPosition(x_in);
}

template <typename ScalarContainer, typename VectorContainer>  
Surface<ScalarContainer, VectorContainer>::~Surface()
{
    assert(checked_out_work_sca_);

    cout<<"Total size of Q ("<<checked_out_work_sca_<<") "<<scalar_work_q_.size()<<endl;

    while ( !scalar_work_q_.empty() )
    {
        delete scalar_work_q_.front();
        scalar_work_q_.pop();
    }
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

    sht_.forward(x_, work_arr, shc);
    
    sht_.backward_duv(shc, work_arr, Xu);
    GeometricDot(Xu, normal_, h_);

    xy(h_, h_, k_);
    axpy(static_cast<value_type>(-1), k_, k_);

    xy(F, h_, h_);
    axpy(static_cast<value_type>(-1), h_, h_);

    Sca* L(produceSca(x_));
    sht_.backward_d2u(shc, work_arr, Xu);
    GeometricDot(Xu, normal_, *L);
    
    Sca* N(produceSca(x_));
    xy(G, *L, *N);
    axpy(static_cast<value_type>(.5), *N, h_, h_); 

    
    sht_.backward_d2v(shc, work_arr, Xu);
    GeometricDot(Xu, normal_, *N);

    xy(*L, *N, *L);
    axpy(static_cast<value_type>(1), *L, k_, k_); 

    xyInv(k_,w_,k_);
    xyInv(k_,w_,k_);

    xy(E, *N, *N);     
    axpy(static_cast<value_type>(.5), *N, h_, h_); 

    recycleSca(L);
    recycleSca(N);

    sht_.Filter(k_,work_arr, shc, k_);
    sht_.Filter(h_, work_arr, shc, h_);

    second_forms_are_stale_ = false;
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::grad(const ScalarContainer 
    &f_in, VectorContainer &grad_f_out) const
{
    if(first_forms_are_stale_)
        updateFirstForms();

    Sca* scw1(produceSca(x_));
    Sca* scw2(produceSca(x_));

    sht_.FirstDerivatives(f_in, work_arr, shc, *scw1, *scw2);

    xv(*scw1, cu_, grad_f_out);
    recycleSca(scw1);

    xvpw(*scw2, cv_, grad_f_out, grad_f_out);
    recycleSca(scw2);

    sht_.Filter(grad_f_out, work_arr, shc, grad_f_out);
}

template <typename ScalarContainer, typename VectorContainer> 
void Surface<ScalarContainer, VectorContainer>::div(const VectorContainer 
    &f_in, ScalarContainer &div_f_out) const
{
    if(first_forms_are_stale_)
        updateFirstForms();

    sht_.FirstDerivatives(f_in, work_arr, shc, Xu, Xv);

    Sca* scw(produceSca(x_));
    GeometricDot(Xu,cu_, *scw);
    GeometricDot(Xv,cv_, div_f_out);
    axpy(static_cast<value_type>(1),div_f_out, *scw, div_f_out);
    recycleSca(scw);

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
    
    Sca* scw(produceSca(x_));
    GeometricDot(x_,normal_,*scw);
    axpy(static_cast<value_type>(1)/3, *scw, *scw);

    integrator_(*scw, w_, vol_out);
    recycleSca(scw);
}

template< typename ScalarContainer, typename VectorContainer >
void Surface<ScalarContainer, VectorContainer>::getCenters(Vec &centers) const
{
    if(first_forms_are_stale_)
        updateFirstForms();

    Sca* scw(produceSca(x_));
    GeometricDot(x_, x_, *scw);
    axpy(static_cast<value_type>(.5), *scw, *scw);
    xv(*scw, normal_, Xu);
    recycleSca(scw);

    integrator_(Xu, w_, centers);
    
    scw = produceSca(centers);
    volume(*scw);
    uyInv(centers, *scw, centers);
    recycleSca(scw);
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
    ///@todo change the work_arr require s.t. it is the same size of x_.
    work_arr.resize(2 * x_.getNumSubs(), x_.getShOrder(), x_.getGridDim());
    E.replicate(x_);    
    F.replicate(x_);    
    G.replicate(x_);    
    Xu.replicate(x_);    
    Xv.replicate(x_);

    position_has_changed_outside_ = false;
}

template <typename ScalarContainer, typename VectorContainer>  
ScalarContainer* Surface<ScalarContainer, VectorContainer>::produceSca(
    const VectorContainer &ref) const
{
    Sca* scp;
    
    if(scalar_work_q_.empty())
    {
        scp = new ScalarContainer;
    }
    else
    {
        scp = scalar_work_q_.front();
        scalar_work_q_.pop();
    }
    
    scp->replicate(ref);
    ++checked_out_work_sca_;
    return(scp);
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::recycleSca(
    ScalarContainer* scp) const
{
    scalar_work_q_.push(scp);
    --checked_out_work_sca_;
}
