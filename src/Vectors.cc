/**
 * @file   Vectors.cc
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Fri Jan 29 09:36:58 2010
 * 
 * @brief  The implementation of Vectors.
 */

// Constructors
template<typename T> 
Vectors<T>::Vectors(Device<T> *device_in) :
    Scalars<T>(device_in) , n_vecs_(0){}

template<typename T> 
Vectors<T>::Vectors(Device<T> *device_in, int p_in, int num_vecs_in) :
    Scalars<T>(device_in, p_in, 3*num_vecs_in) , n_vecs_(num_vecs_in){}

template<typename T> 
Vectors<T>::Vectors(Device<T> *device_in, int p_in, int num_vecs_in, const T *vec_data_in) :
    Scalars<T>(device_in, p_in, 3*num_vecs_in, vec_data_in) , n_vecs_(num_vecs_in){}

// Utility functions
template<typename T> 
int Vectors<T>::GetVecLength() const
{
    return(3*Scalars<T>::GetFunLength());
}

template<typename T> 
void Vectors<T>::Resize(int n_vecs_in)
{
    Scalars<T>::Resize(3*n_vecs_in);
    this->n_vecs_ = n_vecs_in;
}

template<typename T> 
void DotProduct(const Vectors<T> &x_in, 
    const Vectors<T> &y_in, Scalars<T> &xDy_out)
{
    assert(x_in.GetDataLength() == y_in.GetDataLength());
    assert(x_in.GetDataLength() <= 3*xDy_out.GetDataLength());
    assert(x_in.GetDevicePtr() == y_in.GetDevicePtr() && y_in.GetDevicePtr() == xDy_out.GetDevicePtr());
    
    x_in.GetDevicePtr()->DotProduct(x_in.begin(), y_in.begin(), x_in.GetFunLength(), x_in.n_vecs_, xDy_out.begin());
}

template<typename T> 
void CrossProduct(const Vectors<T>& x_in, 
    const Vectors<T>& y_in, Vectors<T>& xCy_out)
{
    assert(x_in.GetDataLength() == y_in.GetDataLength() &&
        y_in.GetDataLength() == xCy_out.GetDataLength());
    assert(x_in.GetDevicePtr() == y_in.GetDevicePtr() && y_in.GetDevicePtr() == xCy_out.GetDevicePtr());
    
    x_in.GetDevicePtr()->CrossProduct(x_in.begin(), y_in.begin(), x_in.GetFunLength(), x_in.n_vecs_, xCy_out.begin());
}

template<typename T> 
void xvpw(const Scalars<T>& x_in, const Vectors<T>& v_in, 
    const Vectors<T>& w_in, Vectors<T>& xvpw_out)
{
    assert(3*x_in.GetDataLength() == v_in.GetDataLength());
    assert(v_in.GetDataLength() == w_in.GetDataLength() &&
        w_in.GetDataLength() == xvpw_out.GetDataLength());
    assert(x_in.GetDevicePtr() == v_in.GetDevicePtr() && v_in.GetDevicePtr() == w_in.GetDevicePtr() && w_in.GetDevicePtr() == xvpw_out.GetDevicePtr());

    x_in.GetDevicePtr()->xvpw(x_in.begin(), v_in.begin(), w_in.begin(), v_in.GetFunLength(), v_in.n_vecs_, xvpw_out.begin());
}

template<typename T> 
void xvpb(const Scalars<T>& x_in, const Vectors<T>& v_in, 
    T w_in, Vectors<T>& xvpb_out)
{
    assert(3*x_in.GetDataLength() == v_in.GetDataLength());
    assert(v_in.GetDataLength() == xvpb_out.GetDataLength());
    assert(x_in.GetDevicePtr() == v_in.GetDevicePtr() && v_in.GetDevicePtr() == xvpb_out.GetDevicePtr());

    x_in.GetDevicePtr()->xvpb(x_in.begin(), v_in.begin(), w_in, v_in.GetFunLength(), v_in.n_vecs_, xvpb_out.begin());
}

template<typename T> 
void uyInv(const Vectors<T>& u_in, const Scalars<T>& y_in, 
    Vectors<T>& uyInv_out)
{
    assert(3*y_in.GetDataLength() == u_in.GetDataLength());
    assert(u_in.GetDataLength() == uyInv_out.GetDataLength());
    assert(y_in.GetDevicePtr() == u_in.GetDevicePtr() && u_in.GetDevicePtr() == uyInv_out.GetDevicePtr());

    u_in.GetDevicePtr()->uyInv(u_in.begin(), y_in.begin(), u_in.GetFunLength(), u_in.n_vecs_, uyInv_out.begin());
}

template<typename T> 
void avpw(const T* a_in, const Vectors<T> &v_in, 
    const Scalars<T> &w_in, Vectors<T> &avpw_out) 
{
    ///@todo there is no guarantee  that a_in is on the same device.
    assert(v_in.GetDataLength() == w_in.GetDataLength() &&
        w_in.GetDataLength() == avpw_out.GetDataLength());
    assert(v_in.GetDevicePtr() == w_in.GetDevicePtr() && w_in.GetDevicePtr() == avpw_out.GetDevicePtr());

    v_in.GetDevicePtr()->avpw(a_in, v_in.begin(), w_in.begin(), v_in.GetFunLength(), v_in.n_vecs_, avpw_out.begin());
}
