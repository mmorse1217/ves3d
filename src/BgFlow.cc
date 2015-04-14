#include <Enums.h>
#include <Device.h>

template<typename Vec_t>
Error_t BgFlowFactory(Parameters<typename Vec_t::value_type> &params, BgFlowBase<Vec_t> **vInf){

    ASSERT(*vInf==NULL, "Non-null pointer passed");

    switch (params.bg_flow)
    {
	case ShearFlow:
	    *vInf = new ShearFlowImp<Vec_t>(params.bg_flow_param);
	    break;

	case ExtensionalFlow:
	    *vInf = new ExtensionalFlowImp<Vec_t>(params.bg_flow_param);
	    break;

	case ParabolicFlow:
	    //! @bug hard-coding the extra arguments of parabolic flow; should be fixed
	    *vInf = new ParabolicFlowImp<Vec_t>(10.0, params.bg_flow_param);
	    break;

	case PeriodicFlow:
	    //return ErrorEvent::NotImplementedError;
	    *vInf = new PeriodicFlowImp<Vec_t>(params.bg_flow_param);
	    break;

	case UserDefinedFlow:
	    //don't do anything
	    break;

	default:
	    return ErrorEvent::InvalidParameterError;
    }

    return ErrorEvent::Success;
}



// Shear flow /////////////////////////////////////////////////////////////////
template<typename Vec_t>
ShearFlowImp<Vec_t>::ShearFlowImp(value_type shear_rate) :
    shear_rate_(shear_rate) {}

template<typename Vec_t>
void ShearFlowImp<Vec_t>::operator()(const Vec_t &pos, const value_type time,
        Vec_t &vel_inf) const
{
    int n_surfs = pos.getNumSubs();
    int stride = pos.getStride();
    int idx;

    axpy(0.0, pos, vel_inf);
    for(int ii=0;ii<n_surfs;ii++)
    {
        idx = pos.getTheDim() * ii * stride;
        pos.getDevice().Memcpy(vel_inf.begin() + idx,
            pos.begin() + idx + stride + stride
            ,stride * sizeof(typename Vec_t::value_type),
            pos.getDevice().MemcpyDeviceToDevice);
    }
    axpy(shear_rate_, vel_inf, vel_inf);
}

// Parabolic flow ///////////////////////////////////////////////////////////////
template<typename Vec_t>
ParabolicFlowImp<Vec_t>::ParabolicFlowImp(value_type radius,
    value_type center_vel, value_type flow_dir_x, value_type flow_dir_y,
        value_type flow_dir_z) :
    inv_radius2_( - 1.0 / radius / radius),
    center_vel_(center_vel),
    flow_dir_x_(flow_dir_x),
    flow_dir_y_(flow_dir_y),
    flow_dir_z_(flow_dir_z)
{
    //Normalizing the flow direction unit vector
    value_type norm(flow_dir_x_ * flow_dir_x_ +
        flow_dir_y_ * flow_dir_y_ + flow_dir_z_ * flow_dir_z_);
    norm = sqrt(norm);
    flow_dir_x_ /= norm;
    flow_dir_y_ /= norm;
    flow_dir_z_ /= norm;
}

template<typename Vec_t>
void ParabolicFlowImp<Vec_t>::CheckContainers(
    const Vec_t &ref) const
{
    int ns = flow_direction_.getNumSubs();
    flow_direction_.replicate(ref);
    s_wrk_.replicate(ref);

    if (ns == 0 && flow_direction_.getSubLength() > 0)
    {
        int ll = flow_direction_.getStride();
        value_type* buffer = new value_type[ll * DIM];

        for(int jj=0; jj< ll; ++jj)
        {
            buffer[       jj] = flow_dir_x_;
            buffer[  ll + jj] = flow_dir_y_;
            buffer[2*ll + jj] = flow_dir_z_;
        }

        Vec_t::getDevice().Memcpy(flow_direction_.begin(),
            buffer, DIM *ll * sizeof(value_type),
            Vec_t::getDevice().MemcpyHostToDevice);
        delete[] buffer;
    }

    for(int ii=ns; ii<flow_direction_.getNumSubs(); ++ii)
        Vec_t::getDevice().Memcpy(flow_direction_.getSubN_begin(ii),
            flow_direction_.begin(), flow_direction_.getSubLength() *
            sizeof(value_type),
            Vec_t::getDevice().MemcpyDeviceToDevice);
}

template<typename Vec_t>
void ParabolicFlowImp<Vec_t>::operator()(const
    Vec_t &pos, const value_type time, Vec_t &vel_inf) const
{
    this->CheckContainers(pos);

    GeometricDot(pos, flow_direction_, s_wrk_);
    xv(s_wrk_, flow_direction_, vel_inf);
    axpy(static_cast<value_type>(-1), vel_inf, pos, vel_inf);
    GeometricDot(vel_inf, vel_inf, s_wrk_);
    axpy(inv_radius2_, s_wrk_, s_wrk_);
    xvpw(s_wrk_, flow_direction_, flow_direction_, vel_inf);
    axpy(center_vel_, vel_inf, vel_inf);
};

// Extensional flow /////////////////////////////////////////////////////////////
template<typename Vec_t>
ExtensionalFlowImp<Vec_t>::ExtensionalFlowImp(value_type rate) :
    rate_(rate) {}

template<typename Vec_t>
void ExtensionalFlowImp<Vec_t>::operator()(const Vec_t &pos,
					const value_type time,
					Vec_t &vel_inf) const
{
    //! u=(-x,y/2,z/2)
    if (pos.getStride() != coeffs_.getStride())
	AdjustCoeffs(pos.getShOrder());

    Vec_t::getDevice().ax(coeffs_.begin(), pos.begin(), pos.getSubLength(),
     	pos.getNumSubs(), vel_inf.begin());
}

template<typename Vec_t>
void ExtensionalFlowImp<Vec_t>::AdjustCoeffs(int p) const
{
    COUTDEBUG("Adjusting the coefficient vector size");
    coeffs_.resize(1,p);
    int sz(coeffs_.size()), stride(coeffs_.getStride());

    value_type *buffer = new value_type[sz];
    for(int i(0); i<stride; ++i){
	buffer[i              ]	= -rate_;
	buffer[i+stride       ]	= rate_/2;
	buffer[i+stride+stride]	= rate_/2;
    }

    Vec_t::getDevice().Memcpy(coeffs_.begin(), buffer, sizeof(value_type)*sz, DT::MemcpyHostToDevice);
    delete[] buffer;
}

// Taylor vortex ////////////////////////////////////////////////////////////////
/**
 * @bug The data is manipulated directly that causes segmentation
 * fault on any device other than CPU.
 */
template<typename Vec_t>
TaylorVortexImp<Vec_t>::TaylorVortexImp(value_type strength, value_type x_period,
    value_type y_period) :
    strength_(strength), x_period_(x_period), y_period_(y_period)
{
    assert( Vec_t::getDeviceType() == CPU );
}

template<typename Vec_t>
void TaylorVortexImp<Vec_t>::operator()(const Vec_t &pos, const value_type time,
    Vec_t &vel_inf) const
{
    int n_surfs = pos.getNumSubs();
    int stride = pos.getStride();
    int idx;

    wrk_vec1_.replicate(pos);
    wrk_vec2_.replicate(pos);

    axpy( BGPI / x_period_, pos, wrk_vec1_);
    axpy( BGPI / y_period_, pos, wrk_vec2_);

    for ( size_t ii=0;ii<pos.size(); ++ii )
    {
        wrk_vec1_.begin()[ii] = cos(wrk_vec1_.begin()[ii]);
        wrk_vec2_.begin()[ii] = sin(wrk_vec2_.begin()[ii]);
    }

    axpy(0.0, pos, vel_inf);
    for ( int ss=0; ss<n_surfs; ++ss )
        for ( int ii=0;ii<stride; ++ii)
        {
            idx = ss * DIM * stride + ii;
            vel_inf.begin()[idx          ] = -1 * wrk_vec1_.begin()[idx] * wrk_vec2_.begin()[idx + stride];
            vel_inf.begin()[idx + stride ] =  wrk_vec2_.begin()[idx] * wrk_vec1_.begin()[idx + stride];

        }
    axpy(strength_, vel_inf, vel_inf);
}



// Periodic flow /////////////////////////////////////////////////////////////////
template<typename Vec_t>
PeriodicFlowImp<Vec_t>::PeriodicFlowImp(value_type strength, value_type period) :
    period_(period),
    strength_(strength)
{
    assert( Vec_t::getDevice().type() == CPU );
}

template<typename Vec_t>
void PeriodicFlowImp<Vec_t>::operator()(const Vec_t &pos, const value_type time,
        Vec_t &vel_inf) const
{
    int n_surfs = pos.getNumSubs();
    int stride = pos.getStride();
    int idx;

    axpy(0.0, pos, vel_inf);
    for(int ii=0;ii<n_surfs;ii++){
        idx = pos.getTheDim() * ii * stride;
        const value_type* pos_=pos.begin() + idx;
        value_type* vel_=vel_inf.begin() + idx;
        for(size_t jj=0;jj<stride;jj++){
            vel_[jj+0*stride]=sin(pos_[jj+2*stride]*2.0*M_PI/period_);
        }
    }
    axpy(strength_, vel_inf, vel_inf);
}

