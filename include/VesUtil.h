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
void DirectInteraction(T *x_in, T *density_in, int stride, int n_surfs, T *vel_out)
{
#ifndef NDEBUG
    cout<<"DirectInteraction()"<<endl;
#endif
 
#ifdef PROFILING
    double ss = get_seconds();
#endif

    int trg_idx, src_idx;
    T px, py, pz, tx, ty, tz, dx, dy, dz, invR, cpx, cpy, cpz, cc;
    
#pragma omp parallel for private(trg_idx, src_idx, px, py, pz, tx, ty, tz, dx, dy, dz, invR, cpx, cpy, cpz, cc)
    for(int ii=0;ii<n_surfs;++ii)
        for(int jj=0;jj<stride;++jj) 
        {
            trg_idx = 3*ii*stride + jj;
        
            px = 0;
            py = 0;
            pz = 0;
            
            tx=x_in[trg_idx                   ];
            ty=x_in[trg_idx + stride          ];
            tz=x_in[trg_idx + stride + stride ];
            

            for(int kk=0;kk<n_surfs;++kk)
                for(int ll=0;ll<stride;++ll) 
                {
                    src_idx = 3*kk*stride + ll;
                    
                    dx=x_in[src_idx                  ]-tx;
                    dy=x_in[src_idx + stride         ]-ty;
                    dz=x_in[src_idx + stride + stride]-tz;
                    
                    invR = dx*dx;
                    invR+= dy*dy;
                    invR+= dz*dz;
                
                    if (invR!=0)
                        invR = 1.0/sqrt(invR);
                    
                    cpx = density_in[src_idx                  ];
                    cpy = density_in[src_idx + stride         ];
                    cpz = density_in[src_idx + stride + stride];
                
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
            
            vel_out[trg_idx                   ] = px;
            vel_out[trg_idx + stride          ] = py;
            vel_out[trg_idx + stride + stride ] = pz;
        }

#ifdef PROFILING
    ss = get_seconds()-ss;
    cout<<"DeviceCPU::DirectStokes takes (sec) : "<<ss<<endl;
#endif
    return;
}

///Parser ///////////////////////////////////////////////////////
#include <fstream>
#include <iostream>
#include <cstring>

template<typename T>
T ParserLite(const char *filename)
{
    ifstream in_file(filename);
    
    T class_instance;
    
    
    if ( in_file.is_open() ) {
        
        string line, var_name;
        
        while ( getline ( in_file, line ) ) {
            string::size_type i = line.find_first_not_of ( " \t\n\v" );
            
            if ( i == string::npos || (i !=string::npos && line[i] == '#') )
                continue;
            
            char *cstr = new char [line.size()+1];
            char *val  = new char [line.size()+1];
            
            i = line.find_first_of ( "=" );
            var_name = line.substr(0,i);
            line = line.substr(++i);
            
            strcpy(cstr, var_name.c_str());
            sscanf(cstr,"%s ",cstr);
            
            class_instance.SetMember(cstr,line);

            delete cstr;
        }   
    }
    return class_instance;
}

template <typename T>
bool String2Num(T &num, const string &s)
{
  std::istringstream iss(s);
  return !(iss >> num).fail();
}

#endif //_VESUTIL_H_ 
