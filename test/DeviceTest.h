/**
 * @test
 * @file   DeviceTest.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Sun Feb 28 12:54:20 2010
 * 
 * @brief The tester class for the device class.
 */

#include "Device.h"
#include <iostream>
#include <math.h>
#include <time.h>
#include <typeinfo>

using namespace std;

template<enum DeviceType DT, typename T>
class DeviceTest
{
  private:
    T eps;
    Device<DT> *device;
    
  public:
    DeviceTest(Device<DT> *device_in);

    bool PerformAll();
    bool TestMalloc();
    bool TestCalloc(); 
    bool TestMemcpy(); 
    bool TestDotProduct(); 
    bool TestCrossProduct(); 
    bool TestSqrt();
    bool Testxy(); 
    bool TestxyInv(); 
    bool TestuyInv(); 
    bool Testaxpy();
    bool Testavpw();
    bool Testxvpw();
    bool TestReduce();
    bool TestTranspose();
    bool TestMax();
};

template<enum DeviceType DT, typename T>
DeviceTest<DT, T>::DeviceTest(Device<DT> *device_in) :
    device(device_in)
{
    eps = 1.0;
    
    do { eps /= 2.0; }
    while ((T)(1.0 + (eps/2.0)) != 1.0);
    
    eps *= 1e2;
}

template<enum DeviceType DT, typename T>
bool DeviceTest<DT,T>::PerformAll()
{
    bool test_result;
    test_result = TestMalloc()
        && TestCalloc() 
        && TestMemcpy() 
        && TestDotProduct() 
        && TestCrossProduct() 
        && TestSqrt()
        && Testxy() 
        && TestxyInv() 
        && TestuyInv() 
        && Testaxpy()
        && Testavpw()
        && Testxvpw()
        && TestReduce()
        && TestTranspose()
        && TestMax();
    
    string res_print = (test_result) ? "Passed" : "Failed";
    cout<<"\n *** Device Class tests with DT="<<DT
        <<" and T="<< typeid(T).name()<< ": " + res_print + " ***\n"<<endl;
    
    return test_result;
}

template<enum DeviceType DT, typename T>
bool DeviceTest<DT,T>::TestMalloc()
{
    size_t arr_size(1e6);
    bool res(false);
    T* a = (T*) device->Malloc(arr_size * sizeof(T));
    if(a != NULL)
        res = true;
    
    device->Free(a);
    
    string res_print = (res) ? "Passed" : "Failed";
    cout<<" * Device::Malloc: " + res_print + " *"<<endl;
    return res;
}

template<enum DeviceType DT, typename T>
bool DeviceTest<DT,T>::TestCalloc()
{
    size_t arr_size = int(1e3); 
    T* a = (T*) device->Calloc(arr_size, sizeof(T));
    T* b = (T*) malloc(arr_size * sizeof(T));
    
    for(int ii=0;ii<arr_size;++ii)
        b[ii]=1.0;
    
    device->Memcpy(b, a, arr_size * sizeof(T), MemcpyDeviceToHost);

    bool res;
    for(size_t idx=0;idx<arr_size;++idx)
        res = (b[idx] == 0) ? true : false;

    string res_print = (res) ? "Passed" : "Failed";
    device->Free(a);
    free(b);

    cout<<" * Device::Calloc: " + res_print + " *"<<endl;
    return (res);
}

template<enum DeviceType DT, typename T>
bool DeviceTest<DT,T>::TestMemcpy()
{
    int arr_size = int(1e6);
    T* a = (T*) device->Malloc(arr_size * sizeof(T));
    T* b = (T*) malloc(arr_size * sizeof(T));

    for(int idx=0;idx<arr_size;++idx)
        b[idx] = 1;
    device->Memcpy(a, b, arr_size * sizeof(T), MemcpyHostToDevice);
        
    for(int idx=0;idx<arr_size;++idx)
        b[idx] = 2;
    device->Memcpy(b, a, arr_size * sizeof(T), MemcpyDeviceToHost);

    bool res = true;
    for(int idx=0;idx<arr_size;++idx)
        res = res && (b[idx] == (T) 1) ? true : false;

    device->Free(a);
    free(b);
        
    string res_print = (res) ? "Passed" : "Failed";
    cout<<" * Device::Memcpy: " + res_print + " *"<<endl;
    return res;
}

template<enum DeviceType DT, typename T>
bool DeviceTest<DT,T>::TestDotProduct()
{
    bool res = true;
    srand48(time(0));
    
    {//Orthogonality (random single)
        int stride = 312, num_vecs = 100, idx;
        int sc_length = stride*num_vecs;
        int arr_length = DIM*num_vecs*stride;
        T* a = (T*) device->Malloc(arr_length * sizeof(T));
        T* b = (T*) device->Malloc(arr_length * sizeof(T));
        T* c = (T*) device->Malloc( sc_length * sizeof(T));

        T* a_host = (T*) malloc(arr_length * sizeof(T));
        T* b_host = (T*) malloc(arr_length * sizeof(T));
        T* c_host = (T*) malloc( sc_length * sizeof(T));
            
        // a dot b should be zero since they are geometrically orthogonal 
        T rnum;
        for(int ii=0;ii<num_vecs;++ii)
            for(int jj=0;jj<stride;++jj)
            {
                idx = DIM*ii*stride+jj;
                rnum = (T) drand48();
                    
                a_host[idx] = rnum; a_host[idx+stride] = -1.0;;
                b_host[idx] = 1.0 ; b_host[idx+stride] = rnum;

                if(DIM==3)
                {
                    a_host[idx+2*stride] = rnum;
                    b_host[idx+2*stride] = 0.0;
                }
                
                c_host[ii*stride+jj] = 1;
            }
            
        device->Memcpy(a, a_host, arr_length * sizeof(T), MemcpyHostToDevice);
        device->Memcpy(b, b_host, arr_length * sizeof(T), MemcpyHostToDevice);
        device->Memcpy(c, c_host,  sc_length * sizeof(T), MemcpyHostToDevice);
            
        device->DotProduct(a,b,stride,num_vecs,c);
        device->Memcpy(c_host, c, sc_length * sizeof(T), MemcpyDeviceToHost);
            
        T err = 0;
        for(int idx=0;idx<sc_length;idx++)
            err = (c_host[idx]>err) ? c_host[idx] : err ;

        device->Free(a);
        device->Free(b);
        device->Free(c);
            
        free(a_host);
        free(b_host);
        free(c_host);
        res = res && ((err<eps) ? true : false);

        string res_print = (res) ? "Passed" : "Failed";
        cout<<" * Device::DotProduct (Orthogonality): " + res_print + " *"<<endl;
    }

    {//Normalization (random, single)
        int stride = 401, num_vecs = 1;
        int sc_length = stride*num_vecs;
        int arr_length = DIM*num_vecs*stride;
        T* a = (T*) device->Malloc(arr_length * sizeof(T));
        T* c = (T*) device->Malloc(sc_length * sizeof(T));

        T* a_host = (T*) malloc(arr_length * sizeof(T));
        T* c_host = (T*) malloc( sc_length * sizeof(T));

        for(int idx=0;idx<stride;idx++)
        {                
            a_host[idx       ] = drand48();
            a_host[idx+stride] = drand48();

            if(DIM==3)
                a_host[idx+stride+stride] = drand48();
            c_host[idx] = 1.0;
        }
            
        device->Memcpy(a, a_host, arr_length * sizeof(T), MemcpyHostToDevice);
        device->Memcpy(c, c_host,  sc_length * sizeof(T), MemcpyHostToDevice);
        device->DotProduct(a,a,stride,num_vecs,c);
        device->Memcpy(c_host, c, sc_length * sizeof(T), MemcpyDeviceToHost);

        T nn;
        for(int idx=0;idx<sc_length;idx++)
        {
            nn = sqrt(c_host[idx]);
            a_host[idx       ] /= nn;
            a_host[idx+stride] /= nn;

            if(DIM==3)
                a_host[idx+stride+stride] /= nn;
        }
                        
        device->Memcpy(a, a_host, arr_length * sizeof(T), MemcpyHostToDevice);
        device->DotProduct(a,a,stride,num_vecs,c);
        device->Memcpy(c_host, c, sc_length * sizeof(T), MemcpyDeviceToHost);

        T err = 0;
        for(int idx=0;idx<sc_length;idx++)
            err = (c_host[idx]-1>err) ? c_host[idx]-1 : err ;
      
        device->Free(a);
        device->Free(c);
        free(a_host);
        free(c_host);

        res = res && (err<eps) ? true : false;
            
        string res_print = (res) ? "Passed" : "Failed";
        cout<<" * Device::DotProduct (Normalization): " + res_print + " *"<<endl;
    }
    return res;
}

template<enum DeviceType DT, typename T>
bool DeviceTest<DT,T>::TestCrossProduct()
{
    bool res = true;
    {//Self-product 
        int stride = 312, num_vecs = 100, idx;
        int arr_length = DIM*num_vecs*stride;
        T* a = (T*) device->Malloc(arr_length * sizeof(T));
        T* c = (T*) device->Malloc(arr_length * sizeof(T));
        T* a_host = (T*) malloc(arr_length * sizeof(T));
        T* c_host = (T*) malloc(arr_length * sizeof(T));
        
        for(int ii=0;ii<num_vecs;++ii)
            for(int idx=0;idx<stride;idx++)
            {                
                a_host[DIM*ii*stride + idx              ] = drand48();
                a_host[DIM*ii*stride + idx+stride       ] = drand48();
                a_host[DIM*ii*stride + idx+stride+stride] = drand48();
            }
        
        device->Memcpy(a, a_host, arr_length * sizeof(T), MemcpyHostToDevice);
        device->CrossProduct(a,a,stride,num_vecs,c);
        device->Memcpy(c_host, c, arr_length * sizeof(T), MemcpyDeviceToHost);
        
        T err = 0;
        for(int idx=0;idx<stride*num_vecs;idx++)
            err = (c_host[idx]>err) ? c_host[idx] : err ;
        
        device->Free(a);
        device->Free(c);
        free(a_host);
        free(c_host);
        
        res = res && (err<eps) ? true : false;
        
        string res_print = (res) ? "Passed" : "Failed";
        cout<<" * Device::CrossProduct (Self-product): " + res_print + " *"<<endl;
    }
    
    {//triple product
        int stride = 723, num_vecs = 1;
        int arr_length = DIM*num_vecs*stride;
        int sc_length = num_vecs*stride;
            
        T* a = (T*) device->Malloc(arr_length * sizeof(T));
        T* b = (T*) device->Malloc(arr_length * sizeof(T));
        T* c = (T*) device->Malloc(arr_length * sizeof(T));
        T* d = (T*) device->Malloc(arr_length * sizeof(T));
            
        T* e = (T*) device->Malloc(sc_length  * sizeof(T));
        T* f = (T*) device->Malloc(sc_length  * sizeof(T));            

        T* a_host = (T*) malloc(arr_length * sizeof(T));
        T* b_host = (T*) malloc(arr_length * sizeof(T));
        T* c_host = (T*) malloc(arr_length * sizeof(T));
        T* d_host = (T*) malloc(arr_length * sizeof(T));

        T* e_host = (T*) malloc(sc_length * sizeof(T));
        T* f_host = (T*) malloc(sc_length * sizeof(T));

        for(int idx=0;idx<stride;idx++)
        {                
            a_host[idx              ] = drand48();
            a_host[idx+stride       ] = drand48();
            a_host[idx+stride+stride] = drand48();

            b_host[idx              ] = drand48();
            b_host[idx+stride       ] = drand48();
            b_host[idx+stride+stride] = drand48();
                
            c_host[idx              ] = drand48();
            c_host[idx+stride       ] = drand48();
            c_host[idx+stride+stride] = drand48();
        }

        device->Memcpy(a, a_host, arr_length * sizeof(T), MemcpyHostToDevice);
        device->Memcpy(b, c_host, arr_length * sizeof(T), MemcpyHostToDevice);
        device->Memcpy(c, c_host, arr_length * sizeof(T), MemcpyHostToDevice);

        // (a x b).c
        device->CrossProduct(a,b,stride,num_vecs,d);
        device->DotProduct(d,c,stride,num_vecs,e);
        device->Memcpy(e_host, e, sc_length * sizeof(T), MemcpyDeviceToHost);

        // (b x a).c
        device->CrossProduct(b,a,stride,num_vecs,d);
        device->DotProduct(d,c,stride,num_vecs,f);
        device->Memcpy(f_host, f, sc_length  * sizeof(T), MemcpyDeviceToHost);

        T err = 0;
        for(int idx=0;idx<sc_length;idx++)
            err = (e_host[idx]+f_host[idx]>err) ? e_host[idx]+f_host[idx] : err ;

        // (c x a).b
        device->CrossProduct(c,a,stride,num_vecs,d);
        device->DotProduct(d,b,stride,num_vecs,f);
        device->Memcpy(f_host, f, sc_length  * sizeof(T), MemcpyDeviceToHost);
            
        for(int idx=0;idx<sc_length;idx++)
            err = ((e_host[idx]-f_host[idx])>err) ? e_host[idx]-f_host[idx] : err ;

        device->Free(a);
        device->Free(b);
        device->Free(c);
        device->Free(d);
        device->Free(e);
        device->Free(f);

        free(a_host);
        free(b_host);
        free(c_host);
        free(d_host);
        free(e_host);
        free(f_host);

        res = res && (err<eps) ? true : false;
            
        string res_print = (res) ? "Passed" : "Failed";
        cout<<" * Device::CrossProduct (Triple product): " + res_print + " *"<<endl;
    }
    return res;
}

template<enum DeviceType DT, typename T>
bool DeviceTest<DT,T>::TestSqrt()
{
    bool res = true;
    {
        int stride = 413, num_vecs = 7;
        int sc_length = stride*num_vecs;
        T* x = (T*) device->Malloc(sc_length * sizeof(T));
        T* y = (T*) device->Malloc(sc_length * sizeof(T));
            
        T* x_host = (T*) malloc(sc_length * sizeof(T));
        T* y_host = (T*) malloc(sc_length * sizeof(T));
        for(int idx=0;idx<sc_length;idx++)
            x_host[idx] = (T) drand48();

        device->Memcpy(x, x_host, sc_length * sizeof(T), MemcpyHostToDevice);
        device->Sqrt(x,sc_length,y);
        device->Memcpy(y_host, y, sc_length * sizeof(T), MemcpyDeviceToHost);
        T err = 0;
        for(int idx=0;idx<sc_length;idx++)
        {
            T d = x_host[idx]-y_host[idx]*y_host[idx];
            err = (d>err) ? d : err ;
        }
        device->Free(x);
        device->Free(y);
        free(x_host);
        free(y_host);
            
        res = res && (err<eps) ? true : false;
            
        string res_print = (res) ? "Passed" : "Failed";
        cout<<" * Device::Sqrt : " + res_print + " *"<<endl;
    }
    return res;
}

template<enum DeviceType DT, typename T>
bool DeviceTest<DT,T>::Testxy()
{
    bool res = true;
    {
        int stride = 65, num_vecs = 4;
        int sc_length = stride*num_vecs;
        T* x = (T*) device->Malloc(sc_length * sizeof(T));
        T* y = (T*) device->Malloc(sc_length * sizeof(T));
        T* z = (T*) device->Malloc(sc_length * sizeof(T));
        
        T* x_host = (T*) malloc(sc_length * sizeof(T));
        T* y_host = (T*) malloc(sc_length * sizeof(T));
        T* z_host = (T*) malloc(sc_length * sizeof(T));
        
        
        for(int idx=0;idx<sc_length;idx++)
        {
            x_host[idx] = (T) drand48();
            y_host[idx] = (T) drand48();
        }
        
        device->Memcpy(x, x_host, sc_length * sizeof(T), MemcpyHostToDevice);
        device->Memcpy(y, y_host, sc_length * sizeof(T), MemcpyHostToDevice);
        device->xy(x,y,sc_length,z);
        device->Memcpy(z_host, z, sc_length * sizeof(T), MemcpyDeviceToHost);
            
        T err = 0, diff;
        for(int idx=0;idx<sc_length;idx++)
        {
            diff = fabs(x_host[idx]*y_host[idx]-z_host[idx]);
            err = (diff>err) ? diff : err ;
        }
        device->Free(x);
        device->Free(y);
        device->Free(z);
        free(x_host);
        free(y_host);
        free(z_host);
            
            
        res = res && (err<eps) ? true : false;
            
        string res_print = (res) ? "Passed" : "Failed";
        cout<<" * Device::xy : " + res_print + " *"<<endl;
    }
    return res;
}

template<enum DeviceType DT, typename T>
bool DeviceTest<DT,T>::TestxyInv()
{
    bool res = true;
    {
        int stride = 413, num_vecs = 1;
        int sc_length = stride*num_vecs;
        T* x = (T*) device->Malloc(sc_length * sizeof(T));
        T* y = (T*) device->Malloc(sc_length * sizeof(T));
            
        T* x_host = (T*) malloc(sc_length * sizeof(T));
        T* y_host = (T*) malloc(sc_length * sizeof(T));
        for(int idx=0;idx<sc_length;idx++)
            x_host[idx] = (T) drand48();
            
        device->Memcpy(x, x_host, sc_length * sizeof(T), MemcpyHostToDevice);
        device->xyInv((T*) NULL, x, sc_length,y);
        device->xyInv((T*) NULL, y, sc_length,y);
        device->Memcpy(y_host, y, sc_length * sizeof(T), MemcpyDeviceToHost);
        T err = 0;
        for(int idx=0;idx<sc_length;idx++)
            err = (x_host[idx]-y_host[idx]>err) ? x_host[idx]-y_host[idx] : err ;
            
        device->Free(x);
        device->Free(y);
        free(x_host);
        free(y_host);
            
        res = res && (err<eps) ? true : false;
            
        string res_print = (res) ? "Passed" : "Failed";
        cout<<" * Device::xInv : " + res_print + " *"<<endl;
    }
        
    {
        int stride = 413, num_vecs = 2;
        int sc_length = stride*num_vecs;
        T* x = (T*) device->Malloc(sc_length * sizeof(T));
        T* y = (T*) device->Malloc(sc_length * sizeof(T));
        T* z = (T*) device->Malloc(sc_length * sizeof(T));

        T* x_host = (T*) malloc(sc_length * sizeof(T));
        T* y_host = (T*) malloc(sc_length * sizeof(T));
        T* z_host = (T*) malloc(sc_length * sizeof(T));
            
        for(int idx=0;idx<sc_length;idx++)
            x_host[idx] = (T) drand48();
                      
        device->Memcpy(x, x_host, sc_length * sizeof(T), MemcpyHostToDevice);
        device->xyInv((T*) NULL, x,sc_length,y); 
        device->xyInv(x,y,sc_length,z);
        device->Memcpy(z_host, z, sc_length * sizeof(T), MemcpyDeviceToHost);
            
        T err = 0;
        for(int idx=0;idx<sc_length;idx++)
            err = (z_host[idx]-1.0>err) ? z_host[idx]-1 : err ;
            
        device->Free(x);
        device->Free(y);
        device->Free(z);
            
        free(x_host);
        free(y_host);
        free(z_host);

        res = res && (err<eps) ? true : false;
            
        string res_print = (res) ? "Passed" : "Failed";
        cout<<" * Device::xyInv : " + res_print + " *"<<endl;
    }
    return res;
}

template<enum DeviceType DT, typename T>
bool DeviceTest<DT,T>::TestuyInv()
{
    bool res = true;
    {
        int stride = 413, num_vecs = 2;
        int vec_length = DIM*stride*num_vecs;
        int sc_length = stride*num_vecs;
        T* u = (T*) device->Malloc(vec_length * sizeof(T));
        T* y = (T*) device->Malloc(sc_length * sizeof(T));

        T* u_host = (T*) malloc(vec_length * sizeof(T));
        T* y_host = (T*) malloc(sc_length * sizeof(T));
            
        for(int idx=0;idx<vec_length;idx++)
            u_host[idx] = (T) drand48();
            
        device->Memcpy(u, u_host, vec_length * sizeof(T), MemcpyHostToDevice);
        device->DotProduct(u, u, stride, num_vecs, y);
        device->Sqrt(y,sc_length,y);
        device->uyInv(u,y,stride,num_vecs,u);
        device->DotProduct(u,u,stride, num_vecs,y);
        device->Memcpy(y_host, y, sc_length * sizeof(T), MemcpyDeviceToHost);
            
        T err = 0;
        for(int idx=0;idx<sc_length;idx++)
            err = (y_host[idx]-1.0>err) ? y_host[idx]-1 : err ;
            
        device->Free(u);
        device->Free(y);
            
        free(u_host);
        free(y_host);
            
        res = res && (err<eps) ? true : false;
            
        string res_print = (res) ? "Passed" : "Failed";
        cout<<" * Device::uyInv : " + res_print + " *"<<endl;
    }
    return res;
}

template<enum DeviceType DT, typename T>
bool DeviceTest<DT,T>::Testaxpy()
{
    bool res = true;
    {
        int stride = 531, num_vecs = 3;
        int sc_length = stride*num_vecs;
        T* x = (T*) device->Malloc(sc_length * sizeof(T));
        T* y = (T*) device->Malloc(sc_length * sizeof(T));

        T* x_host = (T*) malloc(sc_length * sizeof(T));
        T* y_host = (T*) malloc(sc_length * sizeof(T));

            
        for(int idx=0;idx<sc_length;idx++)
            x_host[idx] = (T) drand48();
                        
        device->Memcpy(x, x_host, sc_length * sizeof(T), MemcpyHostToDevice);
        device->axpy((T) -1.0,x, (T*) NULL, sc_length,y);
        device->Memcpy(y_host, y, sc_length * sizeof(T), MemcpyDeviceToHost);
            
        T err = 0;
        for(int idx=0;idx<sc_length;idx++)
            err = (x_host[idx]+y_host[idx]>err) ? (x_host[idx]+y_host[idx]) : err ;
            
        device->Free(x);
        device->Free(y);
        free(x_host);
        free(y_host);

        res = res && (err<eps) ? true : false;
            
        string res_print = (res) ? "Passed" : "Failed";
        cout<<" * Device::axpb : " + res_print + " *"<<endl;
    }

    {
        int stride = 531, num_vecs = 3;
        int sc_length = stride*num_vecs;
        T* x = (T*) device->Malloc(sc_length * sizeof(T));
        T* y = (T*) device->Malloc(sc_length * sizeof(T));
            
        T* x_host = (T*) malloc(sc_length * sizeof(T));
        T* y_host = (T*) malloc(sc_length * sizeof(T));

        for(int idx=0;idx<sc_length;idx++)
            x_host[idx] = (T) drand48();
                        
        device->Memcpy(x, x_host, sc_length * sizeof(T), MemcpyHostToDevice);
        device->axpy((T) -1.0,x,x,sc_length,y);
        device->Memcpy(y_host, y, sc_length * sizeof(T), MemcpyDeviceToHost);

        T err = 0;
        for(int idx=0;idx<sc_length;idx++)
            err = (y_host[idx]>err) ? y_host[idx] : err ;
            
        device->Free(x);
        device->Free(y);
        free(x_host);
        free(y_host);
            
        res = res && (err<eps) ? true : false;
            
        string res_print = (res) ? "Passed" : "Failed";
        cout<<" * Device::axpy : " + res_print + " *"<<endl;
    }
    return res;
}
    
template<enum DeviceType DT, typename T>
bool DeviceTest<DT,T>::Testavpw()
{
    bool res = true;
    {
        
        int stride = 531, num_vecs = 3;
        int sc_length = stride*num_vecs;
        int vec_length = DIM*sc_length;
            
        T* v = (T*) device->Malloc(vec_length * sizeof(T));
        T* a = (T*) device->Malloc(num_vecs * sizeof(T));

        T* v_host = (T*) malloc(vec_length * sizeof(T));
        T* w_host = (T*) malloc(vec_length * sizeof(T));
        T* a_host = (T*) malloc(num_vecs * sizeof(T));

        for(int idx=0;idx<vec_length;idx++)
            v_host[idx] = (T) drand48();
        for(int idx=0;idx<num_vecs;++idx)
            a_host[idx] = idx;
                        
        device->Memcpy(v,v_host, vec_length * sizeof(T), MemcpyHostToDevice);
        device->Memcpy(a,a_host, num_vecs * sizeof(T), MemcpyHostToDevice);
        device->avpw(a,v,v,stride,num_vecs,v);
        device->Memcpy(w_host, v, vec_length * sizeof(T), MemcpyDeviceToHost);

        T err = 0;
        for(int ii=0;ii<num_vecs;++ii)
            for(int jj=0;jj<DIM*stride;jj++)
            {
                size_t idx=ii*DIM*stride + jj;
                T dd = (ii+1)*v_host[idx]-w_host[idx];
                err = (dd>err) ? dd : err ;
            }
        
        device->Free(v);
        device->Free(a);
        free(v_host);
        free(w_host);
        free(a_host);
            
        res = res && (err<eps) ? true : false;

        string res_print = (res) ? "Passed" : "Failed";
        cout<<" * Device::avpw : " + res_print + " *"<<endl;
        
    }
    return res;
}

template<enum DeviceType DT, typename T>
bool DeviceTest<DT,T>::Testxvpw()
{
    bool res = true;
    {
        int stride = 1001, num_vecs = 7;
        int sc_length = stride*num_vecs;
        int vec_length = DIM*sc_length;

        T* x = (T*) device->Malloc(vec_length * sizeof(T));
        T* y = (T*) device->Malloc(vec_length * sizeof(T));
        T* z = (T*) device->Malloc(vec_length * sizeof(T));

        T* a = (T*) device->Malloc(sc_length * sizeof(T) );
        T* b = (T*) device->Malloc(sc_length * sizeof(T) );
        T* c = (T*) device->Malloc(sc_length * sizeof(T) );
        T* d = (T*) device->Malloc(sc_length * sizeof(T) );

        T* x_host = (T*) malloc(vec_length * sizeof(T) );
        T* y_host = (T*) malloc(vec_length * sizeof(T) );
        T* z_host = (T*) malloc(vec_length * sizeof(T) );

        T* a_host = (T*) malloc(sc_length  * sizeof(T) );
        T* b_host = (T*) malloc(sc_length  * sizeof(T) );
        T* c_host = (T*) malloc(sc_length  * sizeof(T) );
        T* d_host = (T*) malloc(sc_length  * sizeof(T) );

            
        for(int idx=0;idx<vec_length;idx++)
        {
            x_host[idx] = (T) drand48();
            y_host[idx] = (T) drand48();
        }

        for(int idx=0;idx<sc_length;idx++)
            a_host[idx] = (T) drand48();

        device->Memcpy(a, a_host, sc_length * sizeof(T), MemcpyHostToDevice);
        device->Memcpy(x, x_host, vec_length * sizeof(T), MemcpyHostToDevice);
        device->Memcpy(y, y_host, vec_length * sizeof(T), MemcpyHostToDevice);

        device->xvpw(a,x, (T*) NULL,stride,num_vecs,z);
            
        device->DotProduct(x,y,stride,num_vecs,b);
        device->DotProduct(z,y,stride,num_vecs,d);

        device->Memcpy(b_host, b, sc_length * sizeof(T), MemcpyDeviceToHost);
        device->Memcpy(d_host, d, sc_length * sizeof(T), MemcpyDeviceToHost);

        T err = 0, diff;
        for(int idx=0;idx<sc_length;idx++)
        {
            diff = fabs(a_host[idx]*b_host[idx]-d_host[idx]);
            err = (diff>err) ? diff : err ;
        }
            
        device->xvpw(a,x,y,stride,num_vecs,z);
        device->DotProduct(x,y,stride,num_vecs,b);
        device->DotProduct(y,y,stride,num_vecs,c);
        device->DotProduct(z,y,stride,num_vecs,d);
            
        device->Memcpy(a_host, a, sc_length * sizeof(T), MemcpyDeviceToHost);
        device->Memcpy(b_host, b, sc_length * sizeof(T), MemcpyDeviceToHost);
        device->Memcpy(c_host, c, sc_length * sizeof(T), MemcpyDeviceToHost);
        device->Memcpy(d_host, d, sc_length * sizeof(T), MemcpyDeviceToHost);

        for(int idx=0;idx<sc_length;idx++)
        {
            diff = fabs(a_host[idx]*b_host[idx]+c_host[idx]-d_host[idx]);
            err = (diff>err) ? diff : err ;
        }

        device->Free(x);
        device->Free(y);
        device->Free(z);
        device->Free(a);
        device->Free(b);
        device->Free(c);
        device->Free(d);

        free(x_host);
        free(y_host);
        free(z_host);
        free(a_host);
        free(b_host);
        free(c_host);
        free(d_host);
            
        res = res && (err<eps) ? true : false;
            
        string res_print = (res) ? "Passed" : "Failed";
        cout<<" * Device::xvpw : " + res_print + " *"<<endl;
    }
    return res;
}    

template<enum DeviceType DT, typename T>
bool DeviceTest<DT,T>::TestReduce()
{
    bool res = true;
    
    int stride = 312;
    int ns = 5;
    int length = ns*stride;
        
    T *x = (T*) device->Malloc(length * sizeof(T));
    T *w = (T*) device->Malloc(length * sizeof(T));
    T *q = (T*) device->Malloc(stride * sizeof(T));
    T *I = (T*) device->Malloc(ns * sizeof(T));        
        
    T *x_host = (T*) malloc(length * sizeof(T));
    T *w_host = (T*) malloc(length * sizeof(T));
    T *q_host = (T*) malloc(length * sizeof(T));
    T *I_host = (T*) malloc(ns * sizeof(T));

    for(int ii=0;ii<ns;++ii)
        for(int jj=0;jj<stride;++jj)
        {
            x_host[ii*stride+jj] = jj;
            w_host[ii*stride+jj] = .5;
            q_host[ii*stride+jj] = .25;
        }

    device->Memcpy(x,x_host,length * sizeof(T),MemcpyHostToDevice);
    device->Memcpy(w,w_host,length * sizeof(T),MemcpyHostToDevice);
    device->Memcpy(q,q_host,stride * sizeof(T),MemcpyHostToDevice);
         
    device->Reduce(x, 1, w, q, stride, ns, I);
    device->Memcpy(I_host,I,ns * sizeof(T),MemcpyDeviceToHost);

    T err = 0;
    T II = (stride-1)*stride/16.0;
    for(int ii=0;ii<ns;++ii)
    {
        T diff = fabs(I_host[ii]-II);
        err = (err>diff) ? err : diff;
    }
    res = res && (err<eps) ? true : false;


    device->Reduce((T*) NULL, 1, w, q, stride, ns, I);
    device->Memcpy(I_host,I,ns * sizeof(T),MemcpyDeviceToHost);

    II = stride/8.0;
    for(int ii=0;ii<ns;++ii)
    {
        T diff = fabs(I_host[ii]-II);
        err = (err>diff) ? err : diff;
    }
    res = res && (err<eps) ? true : false;

    device->Reduce(w, 1, x, q, stride, ns, I);
    device->Memcpy(I_host,I,ns * sizeof(T),MemcpyDeviceToHost);

    II = (stride-1)*stride/16.0;
    for(int ii=0;ii<ns;++ii)
    {
        T diff = fabs(I_host[ii]-II);
        err = (err>diff) ? err : diff;
    }
    res = res && (err<eps) ? true : false;

    free(x_host);
    free(w_host);
    free(q_host);
    free(I_host);

    device->Free(x);
    device->Free(w);
    device->Free(q);
    device->Free(I);

    string res_print = (res) ? "Passed" : "Failed";
    cout<<" * Device::Reduce : " + res_print + " *"<<endl;
    return res;
}

template<enum DeviceType DT, typename T>
bool DeviceTest<DT,T>::TestTranspose()
{
    bool res(true);
    {
        int n = 31, m = 11;
        int length = n*m;
        
        T* x = (T*) device->Malloc(length * sizeof(T));
        T* y = (T*) device->Malloc(length * sizeof(T));
        
        T* x_host = (T*) malloc(length * sizeof(T) );
        
        for(int ii=0;ii<length;++ii)
            x_host[ii] = ii;
        
        device->Memcpy(x, x_host, length * sizeof(T), MemcpyHostToDevice);
        device->Transpose(x, n, m, y);
        device->Transpose(y, m, n, x);
        device->Memcpy(x_host, x, length * sizeof(T), MemcpyDeviceToHost);
        
        T err=0, diff;
        for(int jj=0;jj<length;jj++)
        {     
            diff = fabs(x_host[jj] - jj);
            err = (diff>err) ? diff : err ;
        }
        
        res = res && (err<eps) ? true : false;
        string res_print = (res) ? "Passed" : "Failed";
        cout<<" * Device::Transpose : " + res_print + " *"<<endl;
        
        device->Free(x);
        device->Free(y);
        free(x_host);
    }

    {
        int n = 2, m = 3;
        int length = n*m;
        
        T* x = (T*) device->Malloc(length * sizeof(T));
        T* y = (T*) device->Malloc(length * sizeof(T));
        
        T x_host[] = {0,2,4,1,3,5};
        
        device->Memcpy(x, x_host, length * sizeof(T), MemcpyHostToDevice);
        device->Transpose(x, n, m, y);
        device->Memcpy(x_host, y, length * sizeof(T), MemcpyDeviceToHost);
        
        T err=0, diff;
        for(int jj=0;jj<length;jj++)
        {     
            diff = fabs(x_host[jj] - jj);
            err = (diff>err) ? diff : err ;
        }

        res = res && (err<eps) ? true : false;
        string res_print = (res) ? "Passed" : "Failed";
        cout<<" * Device::Transpose : " + res_print + " *"<<endl;
        
        device->Free(x);
        device->Free(y);
    }

    return res;
}


template<enum DeviceType DT, typename T>
bool DeviceTest<DT,T>::TestMax()
{
    bool res = true;
    {
        int length = 10012;
        T* x = (T*) device->Malloc(length * sizeof(T));
        T* x_host = (T*) malloc(length * sizeof(T));
            
        T max = abs(x_host[0]);
        for(int idx=0;idx<length;idx++)
        {
            x_host[idx] = (T) drand48() * 10 - 5;
            max = (max > abs(x_host[idx])) ? max : abs(x_host[idx]);
        }
            
        device->Memcpy(x, x_host, length * sizeof(T), MemcpyHostToDevice);
        T mx = device->MaxAbs(x,length);

        device->Free(x);
        free(x_host);
                        
        T err = fabs(mx-max);
        res = res && (err<eps) ? true : false;
            
        string res_print = (res) ? "Passed" : "Failed";
        cout<<" * Device::Max : " + res_print + " *"<<endl;
    }
    return res;
}

