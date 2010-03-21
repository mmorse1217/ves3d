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
Scalars<T>::Scalars(Device<T> &device_in) :
    device_(device_in), data_(0), p_(0), n_funs_(0), max_n_funs_(0)
{}

template<typename T> 
Scalars<T>::Scalars(Device<T> &device_in, int p_in, int n_funs_in) : 
    device_(device_in), data_(0), p_(p_in), n_funs_(n_funs_in),
    max_n_funs_(0)
{
    Resize(n_funs_);
}

template<typename T> 
Scalars<T>::Scalars(Device<T> &device_in, int p_in, int n_funs_in, const T *data_in) :
    device_(device_in), data_(0), p_(p_in), n_funs_(n_funs_in),
    max_n_funs_(0)
{
    Resize(n_funs_);
    SetData(data_in);
}

//Destructor
template<typename T> 
Scalars<T>::~Scalars()
{
    device_.Free(data_);
    data_ = 0;
}

//Utility functions
template<typename T> 
void Scalars<T>::Resize(int n_funs_in)
{
    if(n_funs_in > max_n_funs_)
    {
        T *data_old(this->data_);
        this->max_n_funs_ = n_funs_in;
        data_ = device_.Malloc(max_n_funs_ * GetFunLength());
        if(data_old != 0) 
            device_.Memcpy(data_, data_old, GetDataLength(), MemcpyDeviceToDevice);
        device_.Free(data_old);
    }
    this->n_funs_ = n_funs_in;
}

template<typename T> 
int Scalars<T>::GetFunLength() const
{
    return(2*p_*(p_ + 1));
}

template<typename T> 
int Scalars<T>::GetDataLength() const
{
    return(n_funs_ * GetFunLength());
}

template<typename T> 
void Scalars<T>::SetData(const T *data_in)
{
    ///If memory is not allocated return.
    assert(data_ != 0);
    device_.Memcpy(data_, data_in, GetDataLength(), MemcpyHostToDevice);
}

template<typename T> 
const T* Scalars<T>::GetFunctionAt(int fun_idx_in) const
{
    assert(data_ != 0);
    assert(fun_idx_in<n_funs_ && fun_idx_in>=0);
    return(data_ + fun_idx_in*GetFunLength());
}

template<typename T> 
void Scalars<T>::SetFunctionAt(const T *fun_in, int fun_idx_in)
{
    assert(data_ != 0);
    assert(fun_idx_in<n_funs_ && fun_idx_in>=0);

    device_.Memcpy(data_ + fun_idx_in*GetFunLength(), fun_in, GetFunLength(), MemcpyHostToDevice);
}

template<typename T> 
void Scalars<T>::Sqrt()
{
    device_.Sqrt(data_, GetFunLength(), n_funs_, data_);
}

///@todo this may need to be device dependent.
template<typename T> 
T Scalars<T>::Max()
{
    assert(data_ != 0);
    int length = GetDataLength();
    T max_val=*data_;
    
    for(int idx = 0;idx<length;idx++)
        max_val = (max_val > data_[idx]) ? max_val : data_[idx];

    return(max_val);
}

//Friends
template<typename T> 
void axpy(T a_in, const Scalars<T>& x_in, const Scalars<T>& y_in, Scalars<T>& axpy_out)
{
    assert(x_in.device_ == y_in.device_ && y_in.device_ == axpy_out.device_);
    x_in.device_.axpy(a_in,x_in.data_, y_in.data_, x_in.GetFunLength(), x_in.n_funs_, axpy_out.data_);
}

template<typename T> 
void axpb(T a_in, const Scalars<T> &x_in, T y_in, Scalars<T> &axpb_out)
{
    assert(x_in.device_ == axpb_out.device_);
    x_in.device_.axpb(a_in,x_in.data_, y_in, x_in.GetFunLength(), x_in.n_funs_, axpb_out.data_);
}

template<typename T> 
void xy(const Scalars<T> &x_in, const Scalars<T> &y_in, Scalars<T> &xy_out)
{
    assert(x_in.device_ == y_in.device_ && y_in.device_ == xy_out.device_);
    x_in.device_.xy(x_in.data_, y_in.data_, x_in.GetFunLength(), x_in.n_funs_,xy_out.data_);
}

template<typename T> 
void xyInv(const Scalars<T> &x_in, const Scalars<T> &y_in, Scalars<T> &xyInv_out)
{
    assert(x_in.device_ == y_in.device_ && y_in.device_ == xyInv_out.device_);
    x_in.device_.xyInv(x_in.data_, y_in.data_, x_in.GetFunLength(), x_in.n_funs_,xyInv_out.data_);
}
