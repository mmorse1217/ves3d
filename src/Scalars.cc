/**
 * @file   Scalars.cc
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Jan 26 13:48:19 2010
 * 
 * @brief  Implementation for the Scalars class. 
 */

#include <iostream>
using namespace std;
//Constructors
template<typename T> 
Scalars<T>::Scalars(Device<T> *device_in, int p_in, int num_funs_in) : 
    device_(device_in), data_(0), p_(p_in), num_funs_(num_funs_in),
    capacity_(0)
{
    Resize(num_funs_);
}

template<typename T> 
Scalars<T>::~Scalars()
{
    device_->Free(data_);
    data_ = 0;
}

//Utility functions
template<typename T> 
void Scalars<T>::Resize(int num_funs_in)
{
    if(num_funs_in > capacity_)
    {
        T *data_old(this->data_);
        this->capacity_ = num_funs_in;
        data_ = device_->Malloc(capacity_ * GetFunLength());
        if(data_old != 0) 
            device_->Memcpy(data_, data_old, GetDataLength(), MemcpyDeviceToDevice);
        device_->Free(data_old);
    }
    this->num_funs_ = num_funs_in;
}

template<typename T> 
size_t Scalars<T>::GetFunLength() const
{
    return(2*p_*(p_ + 1));
}

template<typename T> 
void Scalars<T>::Sqrt()
{
    device_->Sqrt(data_, GetFunLength(), num_funs_, data_);
}

template<typename T> 
T Scalars<T>::Max()
{
    assert(data_ != 0);
    return(device_->Max(data_, GetDataLength()));
}

//Friends
template<typename T> 
void axpy(T a_in, const Scalars<T>& x_in, const Scalars<T>& y_in, Scalars<T>& axpy_out)
{
    assert(x_in.GetDevicePtr() == y_in.GetDevicePtr() && y_in.GetDevicePtr() == axpy_out.GetDevicePtr());
    axpy_out.GetDevicePtr()->axpy(a_in,x_in.begin(), y_in.begin(), x_in.GetFunLength(), x_in.GetNumFuns(), axpy_out.begin());
}

template<typename T> 
void axpb(T a_in, const Scalars<T> &x_in, T y_in, Scalars<T> &axpb_out)
{
    assert(x_in.GetDevicePtr() == axpb_out.GetDevicePtr());
    axpb_out.GetDevicePtr()->axpb(a_in,x_in.begin(), y_in, x_in.GetFunLength(), x_in.GetNumFuns(), axpb_out.begin());
}

template<typename T> 
void xy(const Scalars<T> &x_in, const Scalars<T> &y_in, Scalars<T> &xy_out)
{
    assert(x_in.GetDevicePtr() == y_in.GetDevicePtr() && y_in.GetDevicePtr() == xy_out.GetDevicePtr());
    xy_out.GetDevicePtr()->xy(x_in.begin(), y_in.begin(), x_in.GetFunLength(), x_in.GetNumFuns(),xy_out.begin());
}

template<typename T> 
void xyInv(const Scalars<T> &x_in, const Scalars<T> &y_in, Scalars<T> &xyInv_out)
{
    assert(x_in.GetDevicePtr() == y_in.GetDevicePtr() && y_in.GetDevicePtr() == xyInv_out.GetDevicePtr());
    xyInv_out.GetDevicePtr()->xyInv(x_in.begin(), y_in.begin(), x_in.GetFunLength(), x_in.GetNumFuns(),xyInv_out.begin());
}

template<typename T> 
Scalars<T>::iterator Scalars<T>::begin()
{
    return(data_);
}

template<typename T> 
Scalars<T>::const_iterator Scalars<T>::begin() const
{
    return(data_);
}

template<typename T> 
Scalars<T>::iterator Scalars<T>::end()
{
    return(data_ + stride_ * num_funs_);
}

template<typename T> 
Scalars<T>::const_iterator Scalars<T>::end() const
{
    return(data_ + stride_ * num_funs_);
}

template<typename T> 
const Device<T>* Scalars<T>::GetDevicePtr() const
{
    return(device_);
}

template<typename T> 
int Scalars<T>::GetShOrder() const
{
    return(p_);
}

template<typename T> 
size_t Scalars<T>::GetNumFuns() const
{
    return(num_funs_);
}

template<typename T> 
size_t Scalars<T>::GetDataLength() const
{
    return(num_funs_ * GetFunLength());
}
