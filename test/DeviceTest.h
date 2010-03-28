/**
 * @file   DeviceTest.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Sun Feb 28 12:54:20 2010
 * 
 * @brief The tester class for the device class and its subclasses.
 */

#include <iostream>
#include <math.h>
#include <time.h>
#include "Device.h"

using namespace std;

template<typename T>
class DeviceTest
{
  private:
    T eps_;
    Device<T>& device_;
    
  public:

    DeviceTest(Device<T>& device_in) :
        device_(device_in)
    {
        eps_ = 1.0;
        
        do { eps_ /= 2.0; }
        while ((T)(1.0 + (eps_/2.0)) != 1.0);
        
        eps_ *= 1e2;
    };

    bool PerformAll()
    {
        bool test_result;
        test_result = TestMalloc() && TestCalloc() && TestMemcpy()&& TestDotProduct() && TestCrossProduct()
            && TestxInv() && Testxy() && TestxyInv() && Testaxpy() && Testaxpb() && Testxvpw()
            && TestShufflePoints();
        
        
        string res_print = (test_result) ? "Passed" : "Failed";
        cout<<"\n*** Device Class tests: " + res_print + " ***\n"<<endl;
        
        return test_result;
    };

    bool TestMalloc()
    {
        int arr_size = (int) 1e6;
        T* a = device_.Malloc(arr_size);
        T* b = (T*) malloc(arr_size * sizeof(T));

        for(int idx=0;idx<arr_size;++idx)
            b[idx] = idx;
        device_.Memcpy(a, b, arr_size, MemcpyHostToDevice);
        device_.Free(a);
        
        cout<<"* Device::Malloc: Passed *"<<endl;
        return true;
    };

    bool TestCalloc()
    {
			int arr_size = int(1e6); 
        T* a = device_.Calloc(arr_size);
        T* b = (T*) malloc(arr_size * sizeof(T));

        for(int idx=0;idx<arr_size;++idx)
            b[idx] = 1;
        device_.Memcpy(b, a, arr_size, MemcpyDeviceToHost);

        bool res;
        for(int idx=0;idx<arr_size;++idx)
            res = (b[idx] == 0) ? true : false;

        string res_print = (res) ? "Passed" : "Failed";
        device_.Free(a);
        
        cout<<"* Device::Calloc: " + res_print + " *"<<endl;
        return (res);
    };

    bool TestMemcpy()
    {
			int arr_size = int(1e6);
        T* a = device_.Malloc(arr_size);
        T* b = (T*) malloc(arr_size * sizeof(T));

        for(int idx=0;idx<arr_size;++idx)
            b[idx] = 1;
        device_.Memcpy(a, b, arr_size, MemcpyHostToDevice);
        
        for(int idx=0;idx<arr_size;++idx)
            b[idx] = 2;
        device_.Memcpy(b, a, arr_size, MemcpyDeviceToHost);

        bool res = true;
        for(int idx=0;idx<arr_size;++idx)
            res = res && (b[idx] == (T) 1) ? true : false;

        device_.Free(a);
        device_.Free(b);
        
        string res_print = (res) ? "Passed" : "Failed";
        cout<<"* Device::Memcpy: " + res_print + " *"<<endl;
        return res;
    };

    bool TestDotProduct()
    {
        bool res = true;
        srand48(time(0));

        {//Orthogonality (random single)
            int stride = 312, num_vecs = 100, idx;
            int sc_length = stride*num_vecs;
            int arr_length = 3*num_vecs*stride;
            T* a = device_.Malloc(arr_length);
            T* b = device_.Malloc(arr_length);
            T* c = device_.Malloc(sc_length);

            T* a_host = (T*) malloc(arr_length * sizeof(T));
            T* b_host = (T*) malloc(arr_length * sizeof(T));
            T* c_host = (T*) malloc(sc_length * sizeof(T));
            
            // a dot b should be zero since they are geometrically orthogonal 
            T rnum;
            for(int ii=0;ii<num_vecs;++ii)
                for(int jj=0;jj<stride;++jj)
                {
                    idx = 3*ii*stride+jj;
                    rnum = (T) drand48();
                    
                    a_host[idx] = rnum; a_host[idx+stride] = 1   ; a_host[idx+2*stride] = rnum;
                    b_host[idx] = 1   ; b_host[idx+stride] = rnum; b_host[idx+2*stride] = -2;
                    c_host[ii*stride+jj] = 1;
                }
            
            device_.Memcpy(a, a_host, arr_length, MemcpyHostToDevice);
            device_.Memcpy(b, b_host, arr_length, MemcpyHostToDevice);
            device_.Memcpy(c, c_host, sc_length , MemcpyHostToDevice);
            
            device_.DotProduct(a,b,stride,num_vecs,c);
            device_.Memcpy(c_host, c, sc_length, MemcpyDeviceToHost);
            
            T err = 0;
            for(int idx=0;idx<sc_length;idx++)
                err = (c_host[idx]>err) ? c_host[idx] : err ;

            device_.Free(a);
            device_.Free(b);
            device_.Free(c);
            
            free(a_host);
            free(b_host);
            free(c_host);
            res = res && ((err<eps_) ? true : false);

            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* Device::DotProduct (Orthogonality): " + res_print + " *"<<endl;
        }

        {//Normalization (random, single)
            int stride = 401, num_vecs = 1;
            int sc_length = stride*num_vecs;
            int arr_length = 3*num_vecs*stride;
            T* a = device_.Malloc(arr_length);
            T* c = device_.Malloc(sc_length);

            T* a_host = (T*) malloc(arr_length * sizeof(T));
            T* c_host = (T*) malloc(sc_length * sizeof(T));

            for(int idx=0;idx<stride;idx++)
            {                
                a_host[idx              ] = drand48();
                a_host[idx+stride       ] = drand48();
                a_host[idx+stride+stride] = drand48();
                c_host[idx] = 1.0;
            }
            
            device_.Memcpy(a, a_host, arr_length, MemcpyHostToDevice);
            device_.Memcpy(c, c_host, sc_length , MemcpyHostToDevice);
            device_.DotProduct(a,a,stride,num_vecs,c);
            device_.Memcpy(c_host, c, sc_length, MemcpyDeviceToHost);

            T nn;
            for(int idx=0;idx<sc_length;idx++)
            {
                nn = sqrt(c_host[idx]);
                a_host[idx              ] /= nn;
                a_host[idx+stride       ] /= nn;
                a_host[idx+stride+stride] /= nn;
            }
                        
            device_.Memcpy(a, a_host, arr_length, MemcpyHostToDevice);
            device_.DotProduct(a,a,stride,num_vecs,c);
            device_.Memcpy(c_host, c, sc_length, MemcpyDeviceToHost);

            T err = 0;
            for(int idx=0;idx<sc_length;idx++)
                err = (c_host[idx]-1>err) ? c_host[idx]-1 : err ;
      
            device_.Free(a);
            device_.Free(c);
            free(a_host);
            free(c_host);

            res = res && (err<eps_) ? true : false;
            
            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* Device::DotProduct (Normalization): " + res_print + " *"<<endl;
        }
        return res;
    };

    bool TestCrossProduct()
    {
        bool res = true;
        {//Self-product 
            int stride = 312, num_vecs = 100, idx;
            int arr_length = 3*num_vecs*stride;
            T* a = device_.Malloc(arr_length);
            T* c = device_.Malloc(arr_length);
            
            for(int ii=0;ii<num_vecs;++ii)
                for(int jj=0;idx<stride;jj++)
                {                
                    idx = ii*3*stride + jj;
                    a[idx              ] = drand48();
                    a[idx+stride       ] = drand48();
                    a[idx+stride+stride] = drand48();
                }

            device_.CrossProduct(a,a,stride,num_vecs,c);
            T err = 0;
            for(int idx=0;idx<stride*num_vecs;idx++)
                err = (c[idx]>err) ? c[idx] : err ;
            
            device_.Free(a);
            device_.Free(c);
            
            res = res && (err<eps_) ? true : false;
            
            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* Device::CrossProduct (Self-product): " + res_print + " *"<<endl;
        }
        
        {//triple product
            int stride = 723, num_vecs = 1;
            int arr_length = 3*num_vecs*stride;
            int sc_length = num_vecs*stride;
            
            T* a = device_.Malloc(arr_length);
            T* b = device_.Malloc(arr_length);
            T* c = device_.Malloc(arr_length);
            T* d = device_.Malloc(arr_length);
            
            T* e = device_.Malloc(sc_length );
            T* f = device_.Malloc(sc_length );            

            for(int idx=0;idx<stride;idx++)
            {                
                a[idx              ] = drand48();
                a[idx+stride       ] = drand48();
                a[idx+stride+stride] = drand48();

                b[idx              ] = drand48();
                b[idx+stride       ] = drand48();
                b[idx+stride+stride] = drand48();
                
                c[idx              ] = drand48();
                c[idx+stride       ] = drand48();
                c[idx+stride+stride] = drand48();
            }

            // (a x b).c
            device_.CrossProduct(a,b,stride,num_vecs,d);
            device_.DotProduct(d,c,stride,num_vecs,e);

            // (b x a).c
            device_.CrossProduct(b,a,stride,num_vecs,d);
            device_.DotProduct(d,c,stride,num_vecs,f);

            T err = 0;
            for(int idx=0;idx<sc_length;idx++)
                err = (e[idx]+f[idx]>err) ? e[idx]+f[idx] : err ;

            // (c x a).b
            device_.CrossProduct(c,a,stride,num_vecs,d);
            device_.DotProduct(d,b,stride,num_vecs,f);
            
            for(int idx=0;idx<sc_length;idx++)
                err = ((e[idx]-f[idx])>err) ? e[idx]-f[idx] : err ;

            device_.Free(a);
            device_.Free(b);
            device_.Free(c);
            device_.Free(d);
            device_.Free(e);
            device_.Free(f);

            res = res && (err<eps_) ? true : false;
            
            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* Device::CrossProduct (Triple product): " + res_print + " *"<<endl;
        }
        return res;
    };

    bool TestxInv()
    {
        bool res = true;
        {
            int stride = 413, num_vecs = 1;
            int sc_length = stride*num_vecs;
            T* x = device_.Malloc(sc_length);
            T* y = device_.Malloc(sc_length);
            
            for(int idx=0;idx<sc_length;idx++)
            x[idx] = (T) drand48();
                        
            device_.xInv(x,stride,num_vecs,y);
            device_.xInv(y,stride,num_vecs,y);
            
            T err = 0;
            for(int idx=0;idx<sc_length;idx++)
                err = (x[idx]-y[idx]>err) ? x[idx]-y[idx] : err ;
            
            device_.Free(x);
            device_.Free(y);
            
            res = res && (err<eps_) ? true : false;
            
            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* Device::xInv : " + res_print + " *"<<endl;
        }
        return res;
    };

    bool Testxy()
    {
        bool res = true;
        {
            int stride = 65, num_vecs = 4;
            int sc_length = stride*num_vecs;
            T* x = device_.Malloc(sc_length);
            T* y = device_.Malloc(sc_length);
            T* z = device_.Malloc(sc_length);
            
            for(int idx=0;idx<sc_length;idx++)
            {
                x[idx] = (T) drand48();
                y[idx] = (T) drand48();
            }
                        
            device_.xy(x,y,stride,num_vecs,z);
            
            T err = 0, diff;
            for(int idx=0;idx<sc_length;idx++)
            {
                diff = fabs(x[idx]*y[idx]-z[idx]);
                err = (diff>err) ? diff : err ;
            }
            device_.Free(x);
            device_.Free(y);
            device_.Free(z);
            
            
            res = res && (err<eps_) ? true : false;
            
            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* Device::xy : " + res_print + " *"<<endl;
        }
        return res;
    }

    bool TestxyInv()
    {
        bool res = true;
        {
            int stride = 413, num_vecs = 2;
            int sc_length = stride*num_vecs;
            T* x = device_.Malloc(sc_length);
            T* y = device_.Malloc(sc_length);
            T* z = device_.Malloc(sc_length);
            
            for(int idx=0;idx<sc_length;idx++)
            x[idx] = (T) drand48();
                        
            device_.xInv(x,stride,num_vecs,y);
            device_.xyInv(x,y,stride,num_vecs,z);
            
            T err = 0;
            for(int idx=0;idx<sc_length;idx++)
                err = (z[idx]-1.0>err) ? z[idx]-1 : err ;
            
            device_.Free(x);
            device_.Free(y);
            device_.Free(z);
            
            res = res && (err<eps_) ? true : false;
            
            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* Device::xyInv : " + res_print + " *"<<endl;
        }
        return res;
    };

    bool Testaxpy()
    {
        bool res = true;
        {
            int stride = 531, num_vecs = 3;
            int sc_length = stride*num_vecs;
            T* x = device_.Malloc(sc_length);
            T* y = device_.Malloc(sc_length);
            
            for(int idx=0;idx<sc_length;idx++)
                x[idx] = (T) drand48();
                        
            device_.axpy(-1.0,x,x,stride,num_vecs,y);
            
            T err = 0;
            for(int idx=0;idx<sc_length;idx++)
                err = (y[idx]>err) ? y[idx] : err ;
            
            device_.Free(x);
            device_.Free(y);
            
            res = res && (err<eps_) ? true : false;
            
            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* Device::axpy: " + res_print + " *"<<endl;
        }
        return res;
    }

    bool Testaxpb()
    {
        bool res = true;
        {
            int stride = 531, num_vecs = 3;
            int sc_length = stride*num_vecs;
            T* x = device_.Malloc(sc_length);
            T* y = device_.Malloc(sc_length);
            
            for(int idx=0;idx<sc_length;idx++)
                x[idx] = (T) drand48();
                        
            device_.axpb(-1.0,x,0.0,stride,num_vecs,y);
            
            T err = 0;
            for(int idx=0;idx<sc_length;idx++)
                err = (x[idx]+y[idx]>err) ? (x[idx]+y[idx]) : err ;
            
            device_.Free(x);
            device_.Free(y);
            
            res = res && (err<eps_) ? true : false;
            
            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* Device::axpb : " + res_print + " *"<<endl;
        }
        return res;
    }
    
    bool Testxvpw()
    {
        bool res = true;
        {
            int stride = 1001, num_vecs = 7;
            int sc_length = stride*num_vecs;
            int vec_length = 3*sc_length;

            T* x = device_.Malloc(vec_length);
            T* y = device_.Malloc(vec_length);
            T* z = device_.Malloc(vec_length);

            T* a = device_.Malloc(sc_length );
            T* b = device_.Malloc(sc_length );
            T* c = device_.Malloc(sc_length );
            T* d = device_.Malloc(sc_length );
            
            for(int idx=0;idx<vec_length;idx++)
            {
                x[idx] = (T) drand48();
                y[idx] = (T) drand48();
            }

            for(int idx=0;idx<sc_length;idx++)
                a[idx] = (T) drand48();

            device_.xvpb(a,x,0.0,stride,num_vecs,z);
            
            device_.DotProduct(x,y,stride,num_vecs,b);
            device_.DotProduct(z,y,stride,num_vecs,d);
            
            T err = 0, diff;
            for(int idx=0;idx<sc_length;idx++)
            {
                diff = fabs(a[idx]*b[idx]-d[idx]);
                err = (diff>err) ? diff : err ;
            }
            
            device_.xvpw(a,x,y,stride,num_vecs,z);
            
            device_.DotProduct(x,y,stride,num_vecs,b);
            device_.DotProduct(y,y,stride,num_vecs,c);
            device_.DotProduct(z,y,stride,num_vecs,d);
            
            for(int idx=0;idx<sc_length;idx++)
            {
                diff = fabs(a[idx]*b[idx]+c[idx]-d[idx]);
                err = (diff>err) ? diff : err ;
            }

            device_.Free(x);
            device_.Free(y);
            device_.Free(z);
            device_.Free(a);
            device_.Free(b);
            device_.Free(c);
            device_.Free(d);
            
            res = res && (err<eps_) ? true : false;
            
            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* Device::xvpw : " + res_print + " *"<<endl;
        }
        return res;
    }    

    bool TestShufflePoints()
    {
        int stride = 11, num_vecs = 1;
        int vec_length = 3*stride*num_vecs;

        T* x = device_.Malloc(vec_length);
        T* y = device_.Malloc(vec_length);

        int idx;
        for(int ii=0;ii<num_vecs;++ii)
            for(int jj=0;jj<stride;jj++)
            {                
                idx = ii*3*stride + jj;
                x[idx              ] = jj;
                x[idx+stride       ] = jj;
                x[idx+stride+stride] = jj;
            }

        device_.ShufflePoints(x, AxisMajor, stride, num_vecs, y);

        T diff, err=0;
        for(int ii=0;ii<num_vecs;++ii)
            for(int jj=0;jj<stride;jj++)
            {                
                idx = ii*3*stride + 3*jj;
                diff = fabs(y[idx] + y[++idx] + y[++idx] - 3*jj);
                err = (diff>err) ? diff : err ;
            }
        bool res = (err<eps_) ? true : false;

        device_.ShufflePoints(y, PointMajor, stride, num_vecs, x);

        err=0;
        for(int ii=0;ii<num_vecs;++ii)
            for(int jj=0;jj<stride;jj++)
            {                
                idx = ii*3*stride + jj;
                diff = fabs(x[idx] + x[idx+stride] + x[idx+stride+stride] - 3*jj);
                err = (diff>err) ? diff : err ;
            }

        res = res && (err<eps_) ? true : false;
        string res_print = (res) ? "Passed" : "Failed";
        cout<<"* Device::ShufflePoints : " + res_print + " *"<<endl;
                
        return res;
    }
};

