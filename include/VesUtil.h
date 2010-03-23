#ifndef _VESUTIL_H_
#define _VESUTIL_H_

#include "Device.h"
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
                vel_out.data_[idx                  ] = shear_rate_ * x_in.data_[idx + stride + stride];
                vel_out.data_[idx + stride         ] = 0.0;
                vel_out.data_[idx + stride + stride] = 0.0;
            }

    }
};

// Interaction /////////////////////////////////////////////////////////////////
template<typename T>
void DirectInteraction(T *x_in, T *density_in, enum CoordinateOrder order, int np, T *vel_out)
{
#ifndef NDEBUG
    cout<<"DirectInteraction()"<<endl;
#endif
    
    assert(order == PointMajor);

    T tx, ty, tz, px, py, pz, dx, dy, dz, invR, cpx, cpy, cpz, cc;

    int t_idx = 0, s_idx;
    //#pragma omp parallel for private(tx, ty, tz, px, py, pz, dx, dy, dz, invR, cpx, cpy, cpz, cc, t_idx, s_idx)
    for(int trg=0; trg<np; trg++)
    {
        tx = x_in[t_idx  ];
        ty = x_in[t_idx+1];
        tz = x_in[t_idx+2];
        
        px = 0;
        py = 0;
        pz = 0;
        
        s_idx = 0;
        for(int src=0; src<np; src++)
        {
            dx  = x_in[s_idx]-tx; 
            cpx = density_in[s_idx++];

            dy=x_in[s_idx]-ty;
            cpy = density_in[s_idx++];

            dz=x_in[s_idx]-tz;
            cpz = density_in[s_idx++];

            if(s_idx>936)
                cerr<<"S "<<s_idx<<endl;

            invR = dx*dx;
            invR+= dy*dy;
            invR+= dz*dz;
                
            if (invR!=0)
                invR = 1.0/sqrt(invR);
            
            cc  = dx*cpx;
            cc += dy*cpy;
            cc += dz*cpz;
            cc *= invR;
            cc *= invR;

            cpx += cc*dx;
            cpy += cc*dy;
            cpz += cc*dz;
                
            px += cpx*invR;
            py += cpy*invR;
            pz += cpz*invR;
        }
        vel_out[t_idx++] = px;
        vel_out[t_idx++] = py;
        vel_out[t_idx++] = pz;
        
        if(t_idx>936)
            cerr<<"T "<<t_idx<<endl;

    }
};

#endif //_VESUTIL_H_ 
