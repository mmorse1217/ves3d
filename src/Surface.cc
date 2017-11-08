/**
 * @file   Surface.cc
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @date   Tue Feb  2 13:37:21 2010
 *
 * @brief  The implementation of the surface class.
 */
template <typename ScalarContainer, typename VectorContainer>
Surface<ScalarContainer, VectorContainer>::Surface(
    int sh_order, const OperatorsMats<Arr_t> &mats,
    const Vec_t *x_in, int diff_filter_freq,
    int rep_filter_freq, ReparamType rep_type,
    value_type rep_exponent
    ) :
    sh_order_(sh_order),
    diff_filter_freq_((diff_filter_freq == -1) ? 2*sh_order_/3 : diff_filter_freq),
    reparam_filter_freq_((rep_filter_freq == -1) ? sh_order_ : rep_filter_freq),
    reparam_type_(rep_type),
    mats_(&mats),
    sht_(sh_order_, mats.getShMats(sh_order_), diff_filter_freq_), ///@todo make sht_ autonomous
    sht_rep_filter_(sh_order, mats.getShMats(sh_order_), reparam_filter_freq_,rep_exponent),
    sht_resample_(NULL),
    containers_are_stale_(true),
    first_forms_are_stale_(true),
    second_forms_are_stale_(true),
    checked_out_work_sca_(0),
    checked_out_work_vec_(0)
{
    if (sh_order_ == mats.p_ )
	sht_resample_ = new SHMats_t(mats.p_up_, mats.mats_p_up_);
    else
	sht_resample_ = new SHMats_t(mats.p_, mats.mats_p_);

    INFO("Initializing with sh_order="<<sh_order_
	<<", diff_filter_freq="<<diff_filter_freq_
        <<", reparam_filter_freq="<<reparam_filter_freq_
        <<", reparam_type="<<reparam_type_
        <<", sh_order="<<sht_.getShOrder()
        <<", sht_rep_filter="<<sht_rep_filter_.getShOrder()
        <<", sht_resample="<<sht_resample_->getShOrder()
	 );

    x_.set_name("position");
    if (x_in != NULL) setPosition(*x_in);
    
    //contact force
    fc_.replicate(x_);
    fc_.getDevice().Memset(fc_.begin(), 0, sizeof(value_type)*fc_.size());
    fc_.set_name("contact_force");
    velocity_.replicate(x_);
    velocity_.getDevice().Memset(velocity_.begin(), 0, sizeof(value_type)*velocity_.size());
    velocity_.set_name("velocity");
    tension_.replicate(x_);
    tension_.getDevice().Memset(tension_.begin(), 0, sizeof(value_type)*tension_.size());
    tension_.set_name("tension");
}

template <typename ScalarContainer, typename VectorContainer>
Surface<ScalarContainer, VectorContainer>::~Surface()
{
    ASSERT(!checked_out_work_sca_,"All scalar are not returned to work pool");
    ASSERT(!checked_out_work_vec_,"All vectors are not returned to work pool");

    purgeTheWorkSpace();
    delete sht_resample_;
}

template <typename ScalarContainer, typename VectorContainer>
void Surface<ScalarContainer, VectorContainer>::setPosition(const Vec_t& x_in)
{
    PROFILESTART();
    containers_are_stale_ = true;
    first_forms_are_stale_ = true;
    second_forms_are_stale_ = true;

    x_.replicate(x_in);
    axpy(static_cast<value_type>(1), x_in, x_);
    PROFILEEND("",0);
}

template <typename ScalarContainer, typename VectorContainer>
VectorContainer&
Surface<ScalarContainer,VectorContainer>::getPositionModifiable()
{
    containers_are_stale_ = true;
    first_forms_are_stale_ = true;
    second_forms_are_stale_ = true;

    return(x_);
}

template <typename ScalarContainer, typename VectorContainer>
const VectorContainer&
Surface<ScalarContainer, VectorContainer>::getPosition() const
{
    return(x_);
}

template <typename ScalarContainer, typename VectorContainer>
const VectorContainer&
Surface<ScalarContainer, VectorContainer>::getNormal() const
{
    if(first_forms_are_stale_)
        updateFirstForms();
    return(normal_);
}

template <typename ScalarContainer, typename VectorContainer>
const ScalarContainer& Surface<ScalarContainer,VectorContainer>::getAreaElement() const
{
    if(first_forms_are_stale_)
        updateFirstForms();
    return(w_);
}

template <typename ScalarContainer, typename VectorContainer>
const ScalarContainer& Surface<ScalarContainer, VectorContainer>::getMeanCurv() const
{
    if(second_forms_are_stale_)
        updateAll();
    return(h_);
}

template <typename ScalarContainer, typename VectorContainer>
const ScalarContainer& Surface<ScalarContainer,VectorContainer>::getGaussianCurv() const
{
    if(second_forms_are_stale_)
        updateAll();
    return(k_);
}

template <typename ScalarContainer, typename VectorContainer>
void Surface<ScalarContainer, VectorContainer>::
getSmoothedShapePosition(Vec_t &smthd_pos) const
{
    PROFILESTART();
    std::auto_ptr<Vec_t> wrk(checkoutVec());
    std::auto_ptr<Vec_t> shc(checkoutVec());
    sht_rep_filter_.lowPassFilter(x_, *wrk, *shc, smthd_pos);
    recycle(wrk);
    recycle(shc);
    PROFILEEND("",0);
}

template <typename ScalarContainer, typename VectorContainer>
void Surface<ScalarContainer, VectorContainer>::
getSmoothedShapePositionReparam(Vec_t &smthd_pos) const
{
    PROFILESTART();
    std::auto_ptr<Vec_t> wrk(checkoutVec());
    std::auto_ptr<Vec_t> shc(checkoutVec());
    if (reparam_type_ == BoxReparam)
        sht_rep_filter_.lowPassFilter(x_, *wrk, *shc, smthd_pos);
    else
        sht_rep_filter_.lowPassFilterPoly(x_, *wrk, *shc, smthd_pos);
    recycle(wrk);
    recycle(shc);
    PROFILEEND("",0);
}

template <typename ScalarContainer, typename VectorContainer>
void Surface<ScalarContainer, VectorContainer>::
mapToTangentSpace(Vec_t &vec_fld, bool upsample) const
{
    PROFILESTART();
    if (first_forms_are_stale_)
        updateFirstForms();

    std::auto_ptr<Sca_t> scp(checkoutSca());

    if (upsample) {
        /*
         * Upsampling is not preferred and is kept due to backward
         * compatibility. It's better to upsample the surface and do
         * computation there.
         */
        std::auto_ptr<Vec_t> wrk(checkoutVec());
        std::auto_ptr<Vec_t> shc(checkoutVec());
        std::auto_ptr<Vec_t> fld(checkoutVec());

        //up-sampling
        int usf(sht_resample_->getShOrder());

        scp->resize(scp->getNumSubs(), usf);
        wrk->resize(wrk->getNumSubs(), usf);
        shc->resize(shc->getNumSubs(), usf);
        fld->resize(fld->getNumSubs(), usf);
        normal_.resize(normal_.getNumSubs(), usf);

        Resample(vec_fld, sht_, *sht_resample_, *shc, *wrk, *fld);
        Resample(normal_, sht_, *sht_resample_, *shc, *wrk, normal_);

        //re-normalizing
        GeometricDot(normal_, normal_, *scp);
        Sqrt(*scp, *scp);
        uyInv(normal_, *scp, normal_);

        //mapping to tangent
        GeometricDot(*fld, normal_, *scp);
        axpy(static_cast<value_type>(-1.0), *scp, *scp);
        xvpw(*scp, normal_, *fld, *fld);

        //down-sampling
        Resample(*fld   , *sht_resample_, sht_, *shc, *wrk, vec_fld);
        Resample(normal_, *sht_resample_, sht_, *shc, *wrk, normal_);
        normal_.resize(normal_.getNumSubs(), sht_.getShOrder());

        recycle(wrk);
        recycle(shc);
        recycle(fld);
    } else {
        scp->replicate(vec_fld);
        GeometricDot(vec_fld, normal_, *scp);
        axpy(static_cast<value_type>(-1.0), *scp, *scp);
        xvpw(*scp, normal_, vec_fld, vec_fld);
    }

    recycle(scp);
    PROFILEEND("",0);
}

template <typename ScalarContainer, typename VectorContainer>
void Surface<ScalarContainer, VectorContainer>::
updateFirstForms() const
{
    PROFILESTART();
    COUTDEBUG("Updating first fundamental forms");
    if(containers_are_stale_)
        checkContainers();

    std::auto_ptr<Vec_t> wrk(checkoutVec());
    std::auto_ptr<Vec_t> shc(checkoutVec());
    std::auto_ptr<Vec_t> dif(checkoutVec());
    std::auto_ptr<Sca_t> scp(checkoutSca());

    // Spherical harmonic coefficient (dif=du,normal=dv)
    sht_.FirstDerivatives(x_, *wrk, *shc, *dif, normal_);

    // First fundamental coefficients
    GeometricDot(*dif, *dif, E);
    GeometricDot(*dif, normal_, F);
    GeometricDot(normal_, normal_, G);

    // Area element
    xy(E, G, w_);
    xy(F, F, *scp);
    axpy(static_cast<value_type>(-1), *scp, w_, w_);

    // Dividing EFG by W^2
    xyInv(E, w_, E);
    xyInv(F, w_, F);
    xyInv(G, w_, G);
    Sqrt(w_, w_);

    //Div and Grad coefficients
    xv(F, normal_, cu_);
    axpy(static_cast<value_type>(-1), cu_, cu_);
    xvpw(G, *dif, cu_, cu_);

    xv(F, *dif, cv_);
    axpy(static_cast<value_type>(-1), cv_, cv_);
    xvpw(E, normal_, cv_, cv_);
    GeometricCross(*dif, normal_, normal_);
    uyInv(normal_, w_, normal_);

    recycle(wrk);
    recycle(shc);
    recycle(dif);
    recycle(scp);

    first_forms_are_stale_ = false;
    PROFILEEND("",0);
}

template <typename ScalarContainer, typename VectorContainer>
void Surface<ScalarContainer, VectorContainer>::
updateAll() const
{
    PROFILESTART();
    COUTDEBUG("Updating all fundamental forms");
    if(first_forms_are_stale_)
        updateFirstForms();

    std::auto_ptr<Vec_t> wrk(checkoutVec());
    std::auto_ptr<Vec_t> shc(checkoutVec());
    std::auto_ptr<Vec_t> dif(checkoutVec());

    sht_.forward(x_, *wrk, *shc);

    sht_.backward_duv(*shc, *wrk, *dif);
    GeometricDot(*dif, normal_, h_);

    xy(h_, h_, k_);
    axpy(static_cast<value_type>(-1), k_, k_);

    xy(F, h_, h_);
    axpy(static_cast<value_type>(-1), h_, h_);

    std::auto_ptr<Sca_t> L(checkoutSca());
    sht_.backward_d2u(*shc, *wrk, *dif);
    GeometricDot(*dif, normal_, *L);

    std::auto_ptr<Sca_t> N(checkoutSca());
    xy(G, *L, *N);
    axpy(static_cast<value_type>(.5), *N, h_, h_);


    sht_.backward_d2v(*shc, *wrk, *dif);
    GeometricDot(*dif, normal_, *N);

    xy(*L, *N, *L);
    axpy(static_cast<value_type>(1), *L, k_, k_);

    xyInv(k_,w_,k_);
    xyInv(k_,w_,k_);

    xy(E, *N, *N);
    axpy(static_cast<value_type>(.5), *N, h_, h_);

    recycle(L);
    recycle(N);

    sht_.lowPassFilter(k_, *wrk, *shc, k_);
    sht_.lowPassFilter(h_, *wrk, *shc, h_);

    recycle(wrk);
    recycle(shc);
    recycle(dif);

    second_forms_are_stale_ = false;
    PROFILEEND("",0);
}

template <typename ScalarContainer, typename VectorContainer>
void Surface<ScalarContainer, VectorContainer>::
grad(const ScalarContainer &f_in, VectorContainer &grad_f_out) const
{
    PROFILESTART();
    if(first_forms_are_stale_)
        updateFirstForms();

    std::auto_ptr<Sca_t> scw1(checkoutSca());
    std::auto_ptr<Sca_t> scw2(checkoutSca());
    std::auto_ptr<Vec_t> shc(checkoutVec());
    std::auto_ptr<Vec_t> wrk(checkoutVec());

    sht_.FirstDerivatives(f_in, *wrk, *shc, *scw1, *scw2);
    xv(*scw1, cu_, grad_f_out);
    xvpw(*scw2, cv_, grad_f_out, grad_f_out);
    sht_.lowPassFilter(grad_f_out, *wrk, *shc, grad_f_out);

    recycle(scw1);
    recycle(scw2);
    recycle(shc);
    recycle(wrk);
    PROFILEEND("",0);
}

template <typename ScalarContainer, typename VectorContainer>
void Surface<ScalarContainer, VectorContainer>::
div(const VectorContainer &f_in, ScalarContainer &div_f_out) const
{
    PROFILESTART();
    if(first_forms_are_stale_)
        updateFirstForms();

    std::auto_ptr<Vec_t> dif(checkoutVec());
    std::auto_ptr<Vec_t> shc(checkoutVec());
    std::auto_ptr<Vec_t> wrk(checkoutVec());
    std::auto_ptr<Sca_t> scw(checkoutSca());

    sht_.forward(f_in, *wrk, *shc);
    sht_.backward_du(*shc, *wrk, *dif);
    GeometricDot(*dif, cu_, *scw);

    sht_.backward_dv(*shc, *wrk, *dif);
    GeometricDot(*dif, cv_, div_f_out);
    recycle(dif);

    axpy(static_cast<value_type>(1),div_f_out, *scw, div_f_out);
    recycle(scw);

    sht_.lowPassFilter(div_f_out, *wrk, *shc, div_f_out);
    recycle(shc);
    recycle(wrk);

    PROFILEEND("",0);
}

template< typename ScalarContainer, typename VectorContainer>
void Surface<ScalarContainer, VectorContainer>::
area(ScalarContainer &area_out) const
{
    PROFILESTART();
    COUTDEBUG("Computing area");
    if(first_forms_are_stale_)
        updateFirstForms();

    integrator_(w_, area_out);
    PROFILEEND("",0);
}

template< typename ScalarContainer, typename VectorContainer >
void Surface<ScalarContainer, VectorContainer>::
volume(ScalarContainer &vol_out) const
{
    PROFILESTART();
    if(first_forms_are_stale_)
        updateFirstForms();

    COUTDEBUG("Computing volume");
    std::auto_ptr<Sca_t> scw(checkoutSca());
    GeometricDot(x_,normal_,*scw);
    axpy(static_cast<value_type>(1)/3, *scw, *scw);

    integrator_(*scw, w_, vol_out);
    recycle(scw);
    PROFILEEND("",0);
}

template< typename ScalarContainer, typename VectorContainer>
void Surface<ScalarContainer, VectorContainer>::
getCenters(Vec_t &centers) const
{
    PROFILESTART();
    COUTDEBUG("Computing centers of mass");
    if(first_forms_are_stale_)
        updateFirstForms();

    std::auto_ptr<Sca_t> scw(checkoutSca());
    GeometricDot(x_, x_, *scw);
    axpy(static_cast<value_type>(.5), *scw, *scw);

    std::auto_ptr<Vec_t> vcw(checkoutVec());
    xv(*scw, normal_, *vcw);
    recycle(scw);

    integrator_(*vcw, w_, centers);
    recycle(vcw);

    scw = checkoutSca();
    scw->replicate(centers);
    volume(*scw);
    uyInv(centers, *scw, centers);
    recycle(scw);
    PROFILEEND("",0);
}

template <typename ScalarContainer, typename VectorContainer>
void Surface<ScalarContainer, VectorContainer>::
checkContainers() const
{
    normal_.replicate(x_);
    w_.replicate(x_);
    h_.replicate(x_);
    k_.replicate(x_);
    cu_.replicate(x_);
    cv_.replicate(x_);

    E.replicate(x_);
    F.replicate(x_);
    G.replicate(x_);

    containers_are_stale_ = false;
}

template <typename ScalarContainer, typename VectorContainer>
std::auto_ptr<ScalarContainer> Surface<ScalarContainer, VectorContainer>::
checkoutSca() const
{
    std::auto_ptr<Sca_t> scp;

   if(scalar_work_q_.empty())
        scp = static_cast<std::auto_ptr<Sca_t> >(new ScalarContainer);
    else
    {
        scp = static_cast<std::auto_ptr<Sca_t> >(scalar_work_q_.front());
        scalar_work_q_.pop();
    }

    scp->replicate(x_);
    ++checked_out_work_sca_;
    return(scp);
}

template <typename ScalarContainer, typename VectorContainer>
void Surface<ScalarContainer, VectorContainer>::
recycle(std::auto_ptr<Sca_t> scp) const
{
    scalar_work_q_.push(scp.release());
    --checked_out_work_sca_;
}

template <typename ScalarContainer, typename VectorContainer>
std::auto_ptr<VectorContainer> Surface<ScalarContainer, VectorContainer>::
checkoutVec() const
{
    std::auto_ptr<Vec_t> vcp;

    if(vector_work_q_.empty())
        vcp = static_cast<std::auto_ptr<Vec_t> >(new VectorContainer);
    else
    {
        vcp = static_cast<std::auto_ptr<Vec_t> >(vector_work_q_.front());
        vector_work_q_.pop();
    }

    vcp->replicate(x_);
    ++checked_out_work_vec_;

    return(vcp);
}

template <typename ScalarContainer, typename VectorContainer>
void Surface<ScalarContainer, VectorContainer>::recycle(
    std::auto_ptr<Vec_t> vcp) const
{
    vector_work_q_.push(vcp.release());
    --checked_out_work_vec_;
}

template <typename ScalarContainer, typename VectorContainer>
void Surface<ScalarContainer, VectorContainer>::
purgeTheWorkSpace() const
{
    while ( !scalar_work_q_.empty() )
    {
        delete scalar_work_q_.front();
        scalar_work_q_.pop();
    }

    while ( !vector_work_q_.empty() )
    {
        delete vector_work_q_.front();
        vector_work_q_.pop();
    }
}

template <typename ScalarContainer, typename VectorContainer>
void Surface<ScalarContainer, VectorContainer>::
linearizedMeanCurv(const Vec_t &x_new, Sca_t &h_lin) const
{
    PROFILESTART();
    if(first_forms_are_stale_)
        updateFirstForms();

    std::auto_ptr<Vec_t> wrk(checkoutVec());
    std::auto_ptr<Vec_t> shc(checkoutVec());
    std::auto_ptr<Vec_t> dif(checkoutVec());
    std::auto_ptr<Sca_t> scw(checkoutSca());

    sht_.forward(x_new, *wrk, *shc);

    //duu
    sht_.backward_d2u(*shc, *wrk, *dif);
    GeometricDot(*dif, normal_, h_lin);
    xy(G, h_lin, h_lin);

    //duv
    sht_.backward_duv(*shc, *wrk, *dif);
    GeometricDot(*dif, normal_, *scw);
    xy(F, *scw, *scw);
    axpy(static_cast<value_type>(-2), *scw, h_lin, h_lin);

    //dvv
    sht_.backward_d2v(*shc, *wrk, *dif);
    GeometricDot(*dif, normal_, *scw);
    xy(E, *scw, *scw);
    axpy(static_cast<value_type>(1), *scw, h_lin, h_lin);

    axpy(static_cast<value_type>(.5), h_lin, h_lin);

    sht_.lowPassFilter(h_lin, *wrk, *shc, h_lin);

    recycle(wrk);
    recycle(shc);
    recycle(dif);
    recycle(scw);
    PROFILEEND("",0);
}

template <typename S, typename V>
Error_t Surface<S, V>::resample(int new_sh_freq, Surface **new_surf /* delete when you're done */ ) const{
    /* current implementation doesn't permit resampling to arbitrary freq */
    if (new_sh_freq != sht_resample_->getShOrder())
	return ErrorEvent::NotImplementedError;

    PROFILESTART();

    // resample x to the target freq
    std::auto_ptr<Vec_t> wrk(checkoutVec());
    std::auto_ptr<Vec_t> shc(checkoutVec());
    std::auto_ptr<Vec_t> xre(checkoutVec());
    wrk->resize(wrk->getNumSubs(), new_sh_freq);
    shc->resize(shc->getNumSubs(), new_sh_freq);
    xre->resize(xre->getNumSubs(), new_sh_freq);
    Resample(x_, sht_, *sht_resample_, *shc, *wrk, *xre);

    // construct new object
    if (*new_surf == NULL){
        COUTDEBUG("Constructing new surface container");
        *new_surf = new Surface(new_sh_freq, *mats_, xre.get(),
            diff_filter_freq_*new_sh_freq/sh_order_ /* keep it relative */,
            reparam_filter_freq_                    /* not relative     */ );
    } else {
        COUTDEBUG("Reusing the surface container");
        ASSERT((*new_surf)->getShOrder()==new_sh_freq, "Container should have the same order as the argument");
        (*new_surf)->setPosition(*xre);
    }

    recycle(wrk);
    recycle(shc);
    recycle(xre);
    PROFILEEND("",0);

    return ErrorEvent::Success;
}

template <typename S, typename V>
Error_t Surface<S, V>::pack(std::ostream &os, Streamable::Format format) const{

    ASSERT(format==Streamable::ASCII, "BIN is not supported yet");

    os<<"SURFACE\n";
    os<<"version: "<<VERSION<<"\n";
    os<<"name: "<<Streamable::name_<<"\n";
    os<<"sh_order: "<<x_.getShOrder()<<"\n";
    os<<"SHT_order: "<<sht_.getShOrder()<<"\n";
    os<<"SHT_rep_order: "<<sht_rep_filter_.getShOrder()<<"\n";
    os<<"SHT_rep_exponent: "<<sht_rep_filter_.getShFilterExponent()<<"\n";
    os<<"SHT_resample_order: "<<sht_resample_->getShOrder()<<"\n";
    os<<"diff_filter_freq: "<<diff_filter_freq_<<"\n";
    os<<"reparam_filter_freq: "<<reparam_filter_freq_<<"\n";
    os<<"reparam_type: "<<reparam_type_<<"\n";
    x_.pack(os,format);
    fc_.pack(os,format);
    velocity_.pack(os,format);
    tension_.pack(os,format);
    os<<"/SURFACE\n";

    return ErrorEvent::Success;
}
template <typename S, typename V>
Error_t Surface<S, V>::unpack(std::istream &is, Streamable::Format format){

    ASSERT(format==Streamable::ASCII, "BIN is not supported yet");
    std::string s,key;
    int ii,version(0);
    value_type v;

    is>>s;
    ASSERT(s=="SURFACE", "Bad input string (missing header).");

    is>>key;
    if (key=="version:") {
        is>>version>>key;
        if (key=="+") {++version;is>>key;}
    };

    is>>Streamable::name_;
    ASSERT(key=="name:", "bad key name");

    is>>key>>ii;
    ASSERT(key=="sh_order:", "bad key sh_order");

    is>>key>>ii;
    ASSERT(key=="SHT_order:", "bad key SHT_order");
    ASSERT(ii==sht_.getShOrder(), "incompatible data (different sh_order), cannot unpack");

    is>>key>>ii;
    ASSERT(key=="SHT_rep_order:", "bad key SHT_rep_order");
    ASSERT(ii==sht_rep_filter_.getShOrder(), "incompatible data (different sh_order), cannot unpack");

    is>>key>>v;
    ASSERT(key=="SHT_rep_exponent:", "bad key SHT_rep_exponent");
    ASSERT(v==sht_rep_filter_.getShFilterExponent(), "incompatible data (different sh_filter_exponent), cannot unpack");

    is>>key>>ii;
    ASSERT(key=="SHT_resample_order:", "bad key SHT_resample_order");
    ASSERT(ii==sht_resample_->getShOrder(), "incompatible data (different sh_order), cannot unpack");

    is>>key>>ii;
    ASSERT(key=="diff_filter_freq:", "bad key diff_filter_freq");
    ASSERT(ii==diff_filter_freq_, "incompatible data (different filter_freq), cannot unpack");

    is>>key>>ii;
    ASSERT(key=="reparam_filter_freq:", "bad key reparam_filter_freq");
    ASSERT(ii==reparam_filter_freq_, "incompatible data (different rep_filter_freq), cannot unpack");

    is>>key>>s;
    ASSERT(key=="reparam_type:", "bad key reparam_type");
    ReparamType rt(EnumifyReparam(s.c_str()));
    if(rt!=reparam_type_) WARN("Reparametrization type switched from "<<rt<<" to "<<reparam_type_);

    x_.unpack(is, format);
    fc_.unpack(is, format);
    velocity_.unpack(is, format);
    tension_.unpack(is, format);
    is>>s;
    ASSERT(s=="/SURFACE", "Bad input string (missing footer).");

    INFO("Unpacked "<<Streamable::name_<<" data from version "<<version<<" (current version "<<VERSION<<")");

    return ErrorEvent::Success;
}




/////////////////////////////////////////////////////////////////////////////////
template <typename S, typename V>
std::ostream& operator<<(std::ostream& output, const Surface<S, V> &sur)
{
    output<<" =====================================================\n"
          <<"  Number of surfaces     : "<<sur.getNumberOfSurfaces()<<"\n"
          <<"  SH order               : "<<sur.getShOrder()<<"\n"
          <<"  Differentiation freq   : "<<sur.diff_filter_freq_<<"\n"
	  <<"  Reparametrization freq : "<<sur.reparam_filter_freq_<<"\n"
          <<" =====================================================\n";

    return(output);
}
