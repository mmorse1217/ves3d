/**
 * @file   DeviceTest.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Sun Feb 28 12:54:20 2010
 * 
 * @brief The tester class for the device class and its subclasses.
 */

#include <iostream>
#include <math.h>
#include "Device.h"

using namespace std;

template<typename T>
class DeviceTest
{
  private:
    T eps_;
    Device<T>* device_;

  public:
    void SetDevice(Device<T>* device_in)
    {
        device_ = device_in;
    };


    DeviceTest()
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
            && TestxInv() && Testxy() && TestxyInv() && Testaxpy() && Testaxpb() && Testxvpw();
        
        
        string res_print = (test_result) ? "Passed" : "Failed";
        cout<<"\n*** Device Class tests: " + res_print + " ***\n"<<endl;
        
        return test_result;
    };

    bool TestMalloc()
    {
        int arr_size = (int) 1e6;
        T* a = device_->Malloc(arr_size);
        
        int idx;
        for(idx=0;idx<arr_size;++idx)
            a[idx] = idx;
        device_->Free(a);
        
        cout<<"* DeviceCPU::Malloc: Passed *"<<endl;
        return true;
    };

    bool TestCalloc()
    {
        int arr_size = 1e6;
        T* a = device_->Calloc(arr_size);
     
        bool res;
        for(int idx=0;idx<arr_size;++idx)
            res = (a[idx] == 0) ? true : false;

        string res_print = (res) ? "Passed" : "Failed";
        device_->Free(a);
        
        cout<<"* DeviceCPU::Calloc: " + res_print + " *"<<endl;
        return (res);
    };

    bool TestMemcpy()
    {
        int arr_size = 1e6;
        T* a = device_->Malloc(arr_size);
        
        int idx;
        for(idx=0;idx<arr_size;++idx)
            a[idx] = idx;
        
        T* b = device_->Calloc(arr_size);

        device_->Memcpy(b, a, arr_size, MemcpyHostToHost);

        bool res = true;
        for(idx=0;idx<arr_size;++idx)
            res = res && (a[idx] == b[idx]) ? true : false;

        device_->Free(a);
        device_->Free(b);
        
        string res_print = (res) ? "Passed" : "Failed";
        cout<<"* DeviceCPU::Memcpy: " + res_print + " *"<<endl;
        return res;
    };

    bool TestDotProduct()
    {
        bool res = true;
        srand48(time(0));

        {//Orthogonality (random single)
            int stride = 312, num_vecs = 1;
            int sc_length = stride*num_vecs;
            int arr_length = 3*num_vecs*stride;
            T* a = device_->Malloc(arr_length);
            T* b = device_->Malloc(arr_length);
            T* c = device_->Malloc(sc_length);
            
            // a dot b should be zero since they are geometrically orthogonal 
            T rnum;
            for(int idx=0;idx<stride;idx++)
            {
                rnum = (T) drand48();
                
                a[idx] = rnum; a[idx+stride] = 1   ; a[idx+2*stride] = rnum;
                b[idx] = 1   ; b[idx+stride] = rnum; b[idx+2*stride] = -2;
                c[idx] = 1;
            }
            
            device_->DotProduct(a,b,stride,num_vecs,c);
            T err = 0;
            for(int idx=0;idx<sc_length;idx++)
                err = (c[idx]>err) ? c[idx] : err ;

            device_->Free(a);
            device_->Free(b);
            device_->Free(c);
            
            res = res && ((err<eps_) ? true : false);

            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* DeviceCPU::DotProduct (Orthogonality): " + res_print + " *"<<endl;
        }

        {//Normalization (random, single)
            int stride = 401, num_vecs = 1;
            int sc_length = stride*num_vecs;
            int arr_length = 3*num_vecs*stride;
            T* a = device_->Malloc(arr_length);
            T* c = device_->Malloc(sc_length);

            for(int idx=0;idx<stride;idx++)
            {                
                a[idx              ] = drand48();
                a[idx+stride       ] = drand48();
                a[idx+stride+stride] = drand48();
                c[idx] = 1.0;
            }
            
            device_->DotProduct(a,a,stride,num_vecs,c);

            T nn;
            for(int idx=0;idx<sc_length;idx++)
            {
                nn = sqrt(c[idx]);
                a[idx              ] /= nn;
                a[idx+stride       ] /= nn;
                a[idx+stride+stride] /= nn;
            }
                        
            device_->DotProduct(a,a,stride,num_vecs,c);

            T err = 0;
            for(int idx=0;idx<sc_length;idx++)
                err = (c[idx]-1>err) ? c[idx]-1 : err ;
      
            device_->Free(a);
            device_->Free(c);
            
            res = res && (err<eps_) ? true : false;
            
            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* DeviceCPU::DotProduct (Normalization): " + res_print + " *"<<endl;
        }
        return res;
    };

    bool TestCrossProduct()
    {
        bool res = true;
        {//Self-product 
            int stride = 312, num_vecs = 1;
            int arr_length = 3*num_vecs*stride;
            T* a = device_->Malloc(arr_length);
            T* c = device_->Malloc(arr_length);
            
            
            for(int idx=0;idx<stride;idx++)
            {                
                a[idx              ] = drand48();
                a[idx+stride       ] = drand48();
                a[idx+stride+stride] = drand48();
            }

            device_->CrossProduct(a,a,stride,num_vecs,c);
            T err = 0;
            for(int idx=0;idx<stride*num_vecs;idx++)
                err = (c[idx]>err) ? c[idx] : err ;
            
            device_->Free(a);
            device_->Free(c);
            
            res = res && (err<eps_) ? true : false;
            
            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* DeviceCPU::CrossProduct (Self-product): " + res_print + " *"<<endl;
        }
        
        {//triple product
            int stride = 723, num_vecs = 1;
            int arr_length = 3*num_vecs*stride;
            int sc_length = num_vecs*stride;
            
            T* a = device_->Malloc(arr_length);
            T* b = device_->Malloc(arr_length);
            T* c = device_->Malloc(arr_length);
            T* d = device_->Malloc(arr_length);
            
            T* e = device_->Malloc(sc_length );
            T* f = device_->Malloc(sc_length );            

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
            device_->CrossProduct(a,b,stride,num_vecs,d);
            device_->DotProduct(d,c,stride,num_vecs,e);

            // (b x a).c
            device_->CrossProduct(b,a,stride,num_vecs,d);
            device_->DotProduct(d,c,stride,num_vecs,f);

            T err = 0;
            for(int idx=0;idx<sc_length;idx++)
                err = (e[idx]+f[idx]>err) ? e[idx]+f[idx] : err ;

            // (c x a).b
            device_->CrossProduct(c,a,stride,num_vecs,d);
            device_->DotProduct(d,b,stride,num_vecs,f);
            
            for(int idx=0;idx<sc_length;idx++)
                err = ((e[idx]-f[idx])>err) ? e[idx]-f[idx] : err ;

            device_->Free(a);
            device_->Free(b);
            device_->Free(c);
            device_->Free(d);
            device_->Free(e);
            device_->Free(f);

            res = res && (err<eps_) ? true : false;
            
            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* DeviceCPU::CrossProduct (Triple product): " + res_print + " *"<<endl;
        }
        return res;
    };

    bool TestxInv()
    {
        bool res = true;
        {
            int stride = 413, num_vecs = 1;
            int sc_length = stride*num_vecs;
            T* x = device_->Malloc(sc_length);
            T* y = device_->Malloc(sc_length);
            
            for(int idx=0;idx<sc_length;idx++)
            x[idx] = (T) drand48();
                        
            device_->xInv(x,stride,num_vecs,y);
            device_->xInv(y,stride,num_vecs,y);
            
            T err = 0;
            for(int idx=0;idx<sc_length;idx++)
                err = (x[idx]-y[idx]>err) ? x[idx]-y[idx] : err ;
            
            device_->Free(x);
            device_->Free(y);
            
            res = res && (err<eps_) ? true : false;
            
            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* DeviceCPU::xInv : " + res_print + " *"<<endl;
        }
        return res;
    };

    bool Testxy()
    {
        bool res = true;
        {
            int stride = 65, num_vecs = 4;
            int sc_length = stride*num_vecs;
            T* x = device_->Malloc(sc_length);
            T* y = device_->Malloc(sc_length);
            T* z = device_->Malloc(sc_length);
            
            for(int idx=0;idx<sc_length;idx++)
            {
                x[idx] = (T) drand48();
                y[idx] = (T) drand48();
            }
                        
            device_->xy(x,y,stride,num_vecs,z);
            
            T err = 0, diff;
            for(int idx=0;idx<sc_length;idx++)
            {
                diff = fabs(x[idx]*y[idx]-z[idx]);
                err = (diff>err) ? diff : err ;
            }
            device_->Free(x);
            device_->Free(y);
            device_->Free(z);
            
            
            res = res && (err<eps_) ? true : false;
            
            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* DeviceCPU::xy : " + res_print + " *"<<endl;
        }
        return res;
    }

    bool TestxyInv()
    {
        bool res = true;
        {
            int stride = 413, num_vecs = 2;
            int sc_length = stride*num_vecs;
            T* x = device_->Malloc(sc_length);
            T* y = device_->Malloc(sc_length);
            T* z = device_->Malloc(sc_length);
            
            for(int idx=0;idx<sc_length;idx++)
            x[idx] = (T) drand48();
                        
            device_->xInv(x,stride,num_vecs,y);
            device_->xyInv(x,y,stride,num_vecs,z);
            
            T err = 0;
            for(int idx=0;idx<sc_length;idx++)
                err = (z[idx]-1.0>err) ? z[idx]-1 : err ;
            
            device_->Free(x);
            device_->Free(y);
            device_->Free(z);
            
            res = res && (err<eps_) ? true : false;
            
            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* DeviceCPU::xyInv : " + res_print + " *"<<endl;
        }
        return res;
    };

    bool Testaxpy()
    {
        bool res = true;
        {
            int stride = 531, num_vecs = 3;
            int sc_length = stride*num_vecs;
            T* x = device_->Malloc(sc_length);
            T* y = device_->Malloc(sc_length);
            
            for(int idx=0;idx<sc_length;idx++)
                x[idx] = (T) drand48();
                        
            device_->axpy(-1.0,x,x,stride,num_vecs,y);
            
            T err = 0;
            for(int idx=0;idx<sc_length;idx++)
                err = (y[idx]>err) ? y[idx] : err ;
            
            device_->Free(x);
            device_->Free(y);
            
            res = res && (err<eps_) ? true : false;
            
            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* DeviceCPU::axpy: " + res_print + " *"<<endl;
        }
        return res;
    }

    bool Testaxpb()
    {
        bool res = true;
        {
            int stride = 531, num_vecs = 3;
            int sc_length = stride*num_vecs;
            T* x = device_->Malloc(sc_length);
            T* y = device_->Malloc(sc_length);
            
            for(int idx=0;idx<sc_length;idx++)
                x[idx] = (T) drand48();
                        
            device_->axpb(-1.0,x,0.0,stride,num_vecs,y);
            
            T err = 0;
            for(int idx=0;idx<sc_length;idx++)
                err = (x[idx]+y[idx]>err) ? (x[idx]+y[idx]) : err ;
            
            device_->Free(x);
            device_->Free(y);
            
            res = res && (err<eps_) ? true : false;
            
            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* DeviceCPU::axpb : " + res_print + " *"<<endl;
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

            T* x = device_->Malloc(vec_length);
            T* y = device_->Malloc(vec_length);
            T* z = device_->Malloc(vec_length);

            T* a = device_->Malloc(sc_length );
            T* b = device_->Malloc(sc_length );
            T* c = device_->Malloc(sc_length );
            T* d = device_->Malloc(sc_length );
            
            for(int idx=0;idx<vec_length;idx++)
            {
                x[idx] = (T) drand48();
                y[idx] = (T) drand48();
            }

            for(int idx=0;idx<sc_length;idx++)
                a[idx] = (T) drand48();

            device_->xvpb(a,x,0.0,stride,num_vecs,z);
            
            device_->DotProduct(x,y,stride,num_vecs,b);
            device_->DotProduct(z,y,stride,num_vecs,d);
            
            T err = 0, diff;
            for(int idx=0;idx<sc_length;idx++)
            {
                diff = fabs(a[idx]*b[idx]-d[idx]);
                err = (diff>err) ? diff : err ;
            }
            
            device_->xvpw(a,x,y,stride,num_vecs,z);
            
            device_->DotProduct(x,y,stride,num_vecs,b);
            device_->DotProduct(y,y,stride,num_vecs,c);
            device_->DotProduct(z,y,stride,num_vecs,d);
            
            for(int idx=0;idx<sc_length;idx++)
            {
                diff = fabs(a[idx]*b[idx]+c[idx]-d[idx]);
                err = (diff>err) ? diff : err ;
            }

            device_->Free(x);
            device_->Free(y);
            device_->Free(z);
            device_->Free(a);
            device_->Free(b);
            device_->Free(c);
            device_->Free(d);
            
            res = res && (err<eps_) ? true : false;
            
            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* DeviceCPU::xvpw : " + res_print + " *"<<endl;
        }
        return res;
    }    
};

