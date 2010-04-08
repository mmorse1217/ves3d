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

        size_t idx;
        int stride = x_in.GetFunLength();
        int n_surfs= x_in.n_vecs_;

        //setting vel_out to zero
        vel_out.device_.Memset(vel_out.data_, 0, 3 * stride * n_surfs);
        
        //copy y and z to the vel_out
        for(int ii=0;ii<n_surfs;ii++)
        {
            idx = 3 * ii * stride + stride;
            x_in.device_.Memcpy(vel_out.data_ + idx, x_in.data_ + idx 
                ,2 * stride, MemcpyDeviceToDevice);
        }

        //calculating r^2 = y^2 + z^2
        DotProduct(vel_out, vel_out, work);
        //scaling r
        axpb( (T) -U/R/R, work, (T) U, work);

        //setting to zero
        vel_out.device_.Memset(vel_out.data_, 0, 3 * stride * n_surfs);
        
        //setting u_x = r;
        for(int ii=0;ii<n_surfs;ii++)
            vel_out.device_.Memcpy(vel_out.data_ + 3 * ii * stride, 
                work.data_ + ii * stride, stride, MemcpyDeviceToDevice);
    }
    
};

////////////////////////////////////////////////////////////////////////////////
// Interaction /////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<typename T>
void DirectInteraction(T *x_in, T *density_in, int stride, int n_surfs, 
    T *vel_out, Device<T> &device, void *user)
{
#ifndef NDEBUG
    cout<<"DirectInteraction()"<<endl;
#endif
 
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

#ifdef PROFILING
    ss = get_seconds()-ss;
    cout<<"DeviceCPU::DirectInteraction takes (sec) : "<<ss<<endl;
#endif
    return;
}

// Interaction multithreads /////////////////////////////////////////////////
#define DIM 3

template<typename T>
void DirectInteractionMultiThread(T *x_in, T *density_in, int stride, 
    int n_surfs, T *vel_out, Device<T> &device, void *user)
{
#ifndef NDEBUG
    cout<<"DirectInteractionMultiThread()"<<endl;
#endif

    ///To have these visible to all threads
    static size_t **n_surfs_per_thread = new size_t*;
    static T **x_all = new T*;
    static T **d_all = new T*;
    static T **v_all = new T*;
    
    const int nt = omp_get_num_threads();
    size_t NP;    
    
    if(omp_get_thread_num() == 0)
        *n_surfs_per_thread = (size_t*) malloc(nt * sizeof(size_t));
    
    //each thread fills the number of points is has
#pragma omp barrier
    *(*n_surfs_per_thread + omp_get_thread_num()) = n_surfs;

    //Making space for the points on master    
#pragma omp barrier
    if(omp_get_thread_num() == 0)
    {
        size_t tmp1, tmp2;
        tmp1 = **n_surfs_per_thread;
        **n_surfs_per_thread = 0;

        for(int ii=1;ii<nt;++ii)
        {
            tmp2 = *(*n_surfs_per_thread + ii);
            *(*n_surfs_per_thread + ii) = tmp1 +  *(*n_surfs_per_thread + ii-1);
            tmp1 = tmp2;
        }
        
        NP = *(*n_surfs_per_thread + nt -1) + tmp1;
       
        NP *= stride;
        *x_all = (T*) malloc(DIM * NP * sizeof(T));
        *d_all = (T*) malloc(DIM * NP * sizeof(T));
        *v_all = (T*) malloc(DIM * NP * sizeof(T));
    }

    //Shuffling the points and copying to the master location
#pragma omp barrier
    
    T *work = (T*) user;
    size_t np = n_surfs * stride;

    //reordering x and copying
    device.ShufflePoints(x_in,  AxisMajor, stride, n_surfs, work);
    device.Memcpy(*x_all + *(*n_surfs_per_thread + omp_get_thread_num()), work, DIM * np , MemcpyDeviceToHost);
    
    //reordering density and copying
    device.ShufflePoints(density_in,  AxisMajor, stride, n_surfs, work);
    device.Memcpy(*d_all + *(*n_surfs_per_thread + omp_get_thread_num()), work, DIM * np , MemcpyDeviceToHost);

#pragma omp barrier
    //Calculating the potential
    if(omp_get_thread_num() == 0)
    {
        size_t n_surfs_direct = 1;
        
        device.ShufflePoints(*x_all, PointMajor, NP, n_surfs_direct, *v_all);
        device.Memcpy(*x_all, *v_all, DIM * NP, MemcpyHostToHost);

        device.ShufflePoints(*d_all, PointMajor, NP, n_surfs_direct, *v_all);
        device.Memcpy(*d_all, *v_all, DIM * NP, MemcpyHostToHost);

        //Direct stokes of Device
        size_t trg_idx_head = 0;
        size_t trg_idx_tail = NP;
    
        device.DirectStokes(NP, n_surfs_direct, trg_idx_head, trg_idx_tail, 
            NULL, *x_all, *x_all, *d_all, *v_all);

        device.ShufflePoints(*v_all, AxisMajor, NP, n_surfs_direct, *d_all);
        device.Memcpy(*v_all, *d_all, DIM * NP, MemcpyHostToHost);

    }
    
    //Distributing the potential
#pragma omp barrier
    device.Memcpy(work, *v_all + *(*n_surfs_per_thread + omp_get_thread_num()), DIM * np , MemcpyHostToDevice);
    device.ShufflePoints(work, PointMajor, stride, n_surfs, vel_out);

    //Freeing memory
#pragma omp barrier
    if(omp_get_thread_num() == 0)
    {
        free(*n_surfs_per_thread);
        device.Free(*x_all);
        device.Free(*d_all);
        device.Free(*v_all);
    }

#ifdef PROFILING
    ss = get_seconds()-ss;
    cout<<"DeviceCPU::DirectInteractionMultiThread takes (sec) : "<<ss<<endl;
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
