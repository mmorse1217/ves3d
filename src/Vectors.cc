/**
 * @file   Vectors.cc
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Fri Jan 29 09:36:58 2010
 * 
 * @brief  The implementation of Vectors.
 */

// Constructors
template<typename T> 
Vectors<T>::Vectors(Device<T> &device_in) :
    Scalars<T>(device_in) , n_vecs_(0){}

template<typename T> 
Vectors<T>::Vectors(Device<T> &device_in, int p_in, int num_vecs_in) :
    Scalars<T>(device_in, p_in, 3*num_vecs_in) , n_vecs_(num_vecs_in){}

template<typename T> 
Vectors<T>::Vectors(Device<T> &device_in, int p_in, int num_vecs_in, const T *vec_data_in) :
    Scalars<T>(device_in, p_in, 3*num_vecs_in, vec_data_in) , n_vecs_(num_vecs_in){}

// Utility functions
template<typename T> 
int Vectors<T>::GetVecLength() const
{
    return(3*Scalars<T>::GetFunLength());
}

template<typename T> 
void Vectors<T>::Resize(int p_in, int n_vecs_in)
{
    Scalars<T>::Resize(p_in, 3*n_vecs_in);
    n_vecs_ = n_vecs_in;
}

template<typename T> 
void DotProduct(const Vectors<T> &x_in, 
    const Vectors<T> &y_in, Scalars<T> &xDy_out)
{
    //Vectors<T>  *temp = static_cast<Vectors<T>*> (&xDy_out);
    assert(x_in.GetDataLength() == y_in.GetDataLength());
    assert(x_in.GetDataLength() <= 3*xDy_out.GetDataLength());
    assert(x_in.device_ == y_in.device_ && y_in.device_ == xDy_out.device_);
    
    x_in.device_.DotProduct(x_in.data_, y_in.data_, x_in.GetFunLength(), x_in.n_vecs_, xDy_out.data_);
}

template<typename T> 
void CrossProduct(const Vectors<T>& x_in, 
    const Vectors<T>& y_in, Vectors<T>& xCy_out)
{
    assert(x_in.GetDataLength() == y_in.GetDataLength() &&
        y_in.GetDataLength() == xCy_out.GetDataLength());
    assert(x_in.device_ == y_in.device_ && y_in.device_ == xCy_out.device_);
    
    x_in.device_.CrossProduct(x_in.data_, y_in.data_, x_in.GetFunLength(), x_in.n_vecs_, xCy_out.data_);
}

template<typename T> 
void xvpw(const Scalars<T>& x_in, const Vectors<T>& v_in, 
    const Vectors<T>& w_in, Vectors<T>& xvpw_out)
{
    assert(3*x_in.GetDataLength() == v_in.GetDataLength());
    assert(v_in.GetDataLength() == w_in.GetDataLength() &&
        w_in.GetDataLength() == xvpw_out.GetDataLength());
    assert(x_in.device_ == v_in.device_ && v_in.device_ == w_in.device_ && w_in.device_ == xvpw_out.device_);

    x_in.device_.xvpw(x_in.data_, v_in.data_, w_in.data_, v_in.GetFunLength(), v_in.n_vecs_, xvpw_out.data_);
}

template<typename T> 
void xvpb(const Scalars<T>& x_in, const Vectors<T>& v_in, 
    T w_in, Vectors<T>& xvpb_out)
{
    assert(3*x_in.GetDataLength() == v_in.GetDataLength());
    assert(v_in.GetDataLength() == xvpb_out.GetDataLength());
    assert(x_in.device_ == v_in.device_ && v_in.device_ == xvpb_out.device_);

    x_in.device_.xvpb(x_in.data_, v_in.data_, w_in, v_in.GetFunLength(), v_in.n_vecs_, xvpb_out.data_);
}
