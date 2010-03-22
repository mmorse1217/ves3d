#ifndef _VESUTIL_H_
#define _VESUTIL_H_

#include "Vectors.h"
#include <iostream>
using namespace std;

// Flow Field /////////////////////////////////////////////////////////////////
template<typename T>
class VelField
{
  public:
    virtual void GetVel(Vectors<T> &x_in, Vectors<T> &vel_out) = 0;
};

template<typename T>
class ShearFlow : public VelField<T>
{
  public:
    T shear_rate_;
    virtual void GetVel(Vectors<T> &x_in, Vectors<T> &vel_out)
    {
#ifndef NDEBUG
    cout<<"ShearFlow::GetVel()"<<endl;
#endif

        assert(x_in.device_ == vel_out.device_);
        assert(x_in.GetDataLength() == vel_out.GetDataLength());

        int n_surfs = x_in.n_vecs_;
        int stride = x_in.GetFunLength();
        int idx;
        ///@bug This may be device dependent.
        for(int ii=0;ii<n_surfs;++ii)
            for(int jj=0;jj<stride;++jj)
            {
                idx = 3*ii*stride +jj;
                vel_out.data_[idx                  ] = shear_rate_ * x_in.data_[idx + stride];
                vel_out.data_[idx + stride         ] = 0.0;
                vel_out.data_[idx + stride + stride] = 0.0;
            }

    }
};

// Interaction /////////////////////////////////////////////////////////////////
template<typename T>
void DirectInteraction(T *x_in, T *density_in, T *vel_out)
{
    cout<<"Interaction"<<endl;
}

#endif //_VESUTIL_H_ 
