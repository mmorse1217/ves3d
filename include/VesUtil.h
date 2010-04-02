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
    virtual void GetVel(Vectors<T> &x_in, Scalars<T> &work, Vectors<T> &vel_out) = 0;
};

template<typename T>
class ShearFlow : public VelField<T>
{
  public:
    T shear_rate_;
    virtual void GetVel(Vectors<T> &x_in, Scalars<T> &work, Vectors<T> &vel_out)
    {
#ifndef NDEBUG
    cout<<"ShearFlow::GetVel()"<<endl;
#endif

        assert(x_in.device_ == vel_out.device_);
        assert(x_in.GetDataLength() == vel_out.GetDataLength());

        int n_surfs = x_in.n_vecs_;
        int stride = x_in.GetFunLength();
        int idx;
        axpb( (T) 0.0, x_in, (T) 0.0, vel_out); //To fill it with zeros
        
        for(int ii=0;ii<n_surfs;ii++)
        {
            idx = 3 * ii * stride;
            x_in.device_.Memcpy(vel_out.data_ + idx, x_in.data_ + idx + stride + stride
                ,stride, MemcpyDeviceToDevice);
        }
        axpb( shear_rate_, vel_out, (T) 0.0, vel_out); 

    }
};

template<typename T>
class ParabolicFlow : public VelField<T>
{
  public:
    T R;
    T U;
    virtual void GetVel(Vectors<T> &x_in, Scalars<T> &work, Vectors<T> &vel_out)
    {
#ifndef NDEBUG
        cout<<"ParabolicFlow::GetVel()"<<endl;
#endif

        assert(x_in.device_ == vel_out.device_);
        assert(x_in.GetDataLength() == vel_out.GetDataLength());

        int n_surfs = x_in.n_vecs_;
        int stride = x_in.GetFunLength();
        int idx;
        axpb( (T) 0.0, x_in, (T) 0.0, vel_out); //To fill it with zeros
        
        for(int ii=0;ii<n_surfs;ii++)
        {
            idx = 3 * ii * stride + stride;
            x_in.device_.Memcpy(vel_out.data_ + idx, x_in.data_ + idx 
                ,2 * stride, MemcpyDeviceToDevice);
        }
        
        DotProduct(vel_out,vel_out,work);
        axpb( (T) -U/R/R, work, (T) U, work);

        axpb( (T) 0.0, x_in, (T) 0.0, vel_out); //To fill it with zeros
        for(int ii=0;ii<n_surfs;ii++)
            x_in.device_.Memcpy(vel_out.data_ + 3 * ii * stride, work.data_ + ii * stride
                , stride, MemcpyDeviceToDevice);
        }
    
};

// Interaction /////////////////////////////////////////////////////////////////
template<typename T>
void DirectInteraction(T *x_in, T *density_in, int stride, int n_surfs, T *vel_out, Device<T> &device, void *user)
{
    T *work = (T*) user;
    size_t np = n_surfs * stride;
    size_t n_surfs_direct = 1;

    //reordering x
    device.ShufflePoints(x_in,  AxisMajor, stride, n_surfs       , work);
    device.ShufflePoints(work, PointMajor, np    , n_surfs_direct, x_in);
    
    //reordering den
    device.ShufflePoints(density_in,  AxisMajor, stride, n_surfs       ,       work);
    device.ShufflePoints(      work, PointMajor, np    , n_surfs_direct, density_in);

    //Direct stokes of Device
    size_t trg_idx_head = 0;
    size_t trg_idx_tail = np;
    device.DirectStokes(np, n_surfs_direct, trg_idx_head, trg_idx_tail, 
        NULL, x_in, x_in, density_in, vel_out);

    //reordering x to the original
    device.ShufflePoints(x_in,  AxisMajor,     np, n_surfs_direct, work);
    device.ShufflePoints(work, PointMajor, stride, n_surfs       , x_in);
    
    //reordering density to the original
    device.ShufflePoints(density_in,  AxisMajor,     np, n_surfs_direct,       work);
    device.ShufflePoints(      work, PointMajor, stride, n_surfs       , density_in);

    //reordering vel to the original
    device.ShufflePoints(vel_out,  AxisMajor,     np, n_surfs_direct, work);
    device.ShufflePoints(   work, PointMajor, stride, n_surfs       , vel_out);
    
}

// #define I_PI 1.0/M_PI/8.0

// template<typename T>
// void DirectInteraction(T *x_in, T *density_in, int stride, int n_surfs, T *vel_out, Device<T> &device, void *user)
// {

// #ifndef NDEBUG
//     cout<<"DirectInteraction()"<<endl;
// #endif
 
// #ifdef PROFILING
//     double ss = get_seconds();
// #endif

//     int trg_idx, src_idx;
//     T px, py, pz, tx, ty, tz, dx, dy, dz, invR, cpx, cpy, cpz, cc;
    
// #pragma omp parallel for private(trg_idx, src_idx, px, py, pz, tx, ty, tz, dx, dy, dz, invR, cpx, cpy, cpz, cc)
//     for(int ii=0;ii<n_surfs;++ii)
//         for(int jj=0;jj<stride;++jj) 
//         {
//             trg_idx = 3*ii*stride + jj;
        
//             px = 0;
//             py = 0;
//             pz = 0;
            
//             tx=x_in[trg_idx                   ];
//             ty=x_in[trg_idx + stride          ];
//             tz=x_in[trg_idx + stride + stride ];
            

//             for(int kk=0;kk<n_surfs;++kk)
//                 for(int ll=0;ll<stride;++ll) 
//                 {
//                     src_idx = 3*kk*stride + ll;
                    
//                     dx=x_in[src_idx                  ]-tx;
//                     dy=x_in[src_idx + stride         ]-ty;
//                     dz=x_in[src_idx + stride + stride]-tz;
                    
//                     invR = dx*dx;
//                     invR+= dy*dy;
//                     invR+= dz*dz;
                
//                     if (invR!=0)
//                         invR = 1.0/sqrt(invR);
                    
//                     cpx = density_in[src_idx                  ];
//                     cpy = density_in[src_idx + stride         ];
//                     cpz = density_in[src_idx + stride + stride];
                
//                     cc  = dx*cpx;
//                     cc += dy*cpy;
//                     cc += dz*cpz;
//                     cc *= invR;
//                     cc *= invR;

//                     cpx += cc*dx;
//                     cpy += cc*dy;
//                     cpz += cc*dz;
                
//                     px += cpx*invR;
//                     py += cpy*invR;
//                     pz += cpz*invR;
//                 }
            
//             vel_out[trg_idx                   ] = px*I_PI;
//             vel_out[trg_idx + stride          ] = py*I_PI;
//             vel_out[trg_idx + stride + stride ] = pz*I_PI;
//         }

// #ifdef PROFILING
//     ss = get_seconds()-ss;
//     cout<<"DeviceCPU::DirectStokes takes (sec) : "<<ss<<endl;
// #endif
//     return;
// }

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
