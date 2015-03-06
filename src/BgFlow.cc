#include <Enums.h>
#include <Device.h>

// Shear flow /////////////////////////////////////////////////////////////////
template<typename VecContainer>
ShearFlow<VecContainer>::ShearFlow(value_type shear_rate) :
    shear_rate_(shear_rate) {}


template<typename VecContainer>
void ShearFlow<VecContainer>::operator()(const VecContainer &pos, const value_type time,
        VecContainer &vel_inf) const
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
            ,stride * sizeof(typename VecContainer::value_type),
            pos.getDevice().MemcpyDeviceToDevice);
    }
    axpy(shear_rate_, vel_inf, vel_inf);
}

// Parabolic flow ///////////////////////////////////////////////////////////////
template<typename ScalarContainer, typename VecContainer>
ParabolicFlow<ScalarContainer, VecContainer>::ParabolicFlow(value_type radius,
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

template<typename ScalarContainer, typename VecContainer>
void ParabolicFlow<ScalarContainer, VecContainer>::CheckContainers(
    const VecContainer &ref) const
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

        VecContainer::getDevice().Memcpy(flow_direction_.begin(),
            buffer, DIM *ll * sizeof(value_type),
            VecContainer::getDevice().MemcpyHostToDevice);
        delete[] buffer;
    }

    for(int ii=ns; ii<flow_direction_.getNumSubs(); ++ii)
        VecContainer::getDevice().Memcpy(flow_direction_.getSubN(ii),
            flow_direction_.begin(), flow_direction_.getSubLength() *
            sizeof(value_type),
            VecContainer::getDevice().MemcpyDeviceToDevice);
}

template<typename ScalarContainer, typename VecContainer>
void ParabolicFlow<ScalarContainer, VecContainer>::operator()(const
    VecContainer &pos, const value_type time, VecContainer &vel_inf) const
{
    this->CheckContainers(pos);

    GeometricDot(pos, flow_direction_, s_wrk_);
    xv(s_wrk_, flow_direction_, vel_inf);
    axpy(static_cast<value_type>(-1), vel_inf, pos, vel_inf);
    GeometricDot(vel_inf, vel_inf, s_wrk_);
    axpy(inv_radius2_, s_wrk_, s_wrk_);
    xvpw(s_wrk_, flow_direction_, flow_direction_, vel_inf);
    axpy(center_vel_, vel_inf, vel_inf);
    //  cout<<"Max Vel :"<<MaxAbs(vel_inf)<<" "<<center_vel_<<endl;
};

// Extensional flow /////////////////////////////////////////////////////////////
template<typename Vec_t>
ExtensionalFlow<Vec_t>::ExtensionalFlow(value_type rate) :
    rate_(rate) {}

template<typename Vec_t>
void ExtensionalFlow<Vec_t>::operator()(const Vec_t &pos,
					const value_type time,
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
    axpy(rate_, vel_inf, vel_inf);
}

// Taylor vortex ////////////////////////////////////////////////////////////////
/**
 * @bug The data is manipulated directly that causes segmentation
 * fault on any device other than CPU.
 */
template<typename VecContainer>
TaylorVortex<VecContainer>::TaylorVortex(value_type strength, value_type x_period,
    value_type y_period) :
    strength_(strength), x_period_(x_period), y_period_(y_period)
{
    assert( VecContainer::getDeviceType() == CPU );
}

template<typename VecContainer>
void TaylorVortex<VecContainer>::operator()(const VecContainer &pos, const value_type time,
    VecContainer &vel_inf) const
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
