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
        test_result = 
            TestMalloc() && 
            TestCalloc() && 
            TestMemcpy() && 
            TestDotProduct() && 
            TestCrossProduct() && 
            TestxInv() &&
            Testxy() && 
            TestxyInv() && 
            TestuyInv() && 
            Testaxpy() &&
            Testaxpb() &&
            Testxvpw() &&
            TestShufflePoints() &&
            TestMax() &&
            TestReduce() &&
            TestResample();
        
        
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
            T* a_host = (T*) malloc(arr_length * sizeof(T));
            T* c_host = (T*) malloc(arr_length * sizeof(T));

            for(int ii=0;ii<num_vecs;++ii)
                for(int idx=0;idx<stride;idx++)
                {                
                    a_host[3*ii*stride + idx              ] = drand48();
                    a_host[3*ii*stride + idx+stride       ] = drand48();
                    a_host[3*ii*stride + idx+stride+stride] = drand48();
                }
            
            device_.Memcpy(a, a_host, arr_length, MemcpyHostToDevice);
            device_.CrossProduct(a,a,stride,num_vecs,c);
            device_.Memcpy(c_host, c, arr_length , MemcpyDeviceToHost);

            T err = 0;
            for(int idx=0;idx<stride*num_vecs;idx++)
                err = (c_host[idx]>err) ? c_host[idx] : err ;
            
            device_.Free(a);
            device_.Free(c);
            free(a_host);
            free(c_host);
            
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

            device_.Memcpy(a, a_host, arr_length, MemcpyHostToDevice);
            device_.Memcpy(b, c_host, arr_length, MemcpyHostToDevice);
            device_.Memcpy(c, c_host, arr_length, MemcpyHostToDevice);

            // (a x b).c
            device_.CrossProduct(a,b,stride,num_vecs,d);
            device_.DotProduct(d,c,stride,num_vecs,e);
            device_.Memcpy(e_host, e, sc_length , MemcpyDeviceToHost);

            // (b x a).c
            device_.CrossProduct(b,a,stride,num_vecs,d);
            device_.DotProduct(d,c,stride,num_vecs,f);
            device_.Memcpy(f_host, f, sc_length , MemcpyDeviceToHost);

            T err = 0;
            for(int idx=0;idx<sc_length;idx++)
                err = (e_host[idx]+f_host[idx]>err) ? e_host[idx]+f_host[idx] : err ;

            // (c x a).b
            device_.CrossProduct(c,a,stride,num_vecs,d);
            device_.DotProduct(d,b,stride,num_vecs,f);
            device_.Memcpy(f_host, f, sc_length , MemcpyDeviceToHost);
            
            for(int idx=0;idx<sc_length;idx++)
                err = ((e_host[idx]-f_host[idx])>err) ? e_host[idx]-f_host[idx] : err ;

            device_.Free(a);
            device_.Free(b);
            device_.Free(c);
            device_.Free(d);
            device_.Free(e);
            device_.Free(f);

            free(a_host);
            free(b_host);
            free(c_host);
            free(d_host);
            free(e_host);
            free(f_host);

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
            
            T* x_host = (T*) malloc(sc_length * sizeof(T));
            T* y_host = (T*) malloc(sc_length * sizeof(T));
            for(int idx=0;idx<sc_length;idx++)
                x_host[idx] = (T) drand48();
            
            device_.Memcpy(x, x_host, sc_length, MemcpyHostToDevice);
            device_.xInv(x,stride,num_vecs,y);
            device_.xInv(y,stride,num_vecs,y);
            device_.Memcpy(y_host, y, sc_length, MemcpyDeviceToHost);
            T err = 0;
            for(int idx=0;idx<sc_length;idx++)
                err = (x_host[idx]-y_host[idx]>err) ? x_host[idx]-y_host[idx] : err ;
            
            device_.Free(x);
            device_.Free(y);
            free(x_host);
            free(y_host);
            
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

            T* x_host = (T*) malloc(sc_length * sizeof(T));
            T* y_host = (T*) malloc(sc_length * sizeof(T));
            T* z_host = (T*) malloc(sc_length * sizeof(T));

            
            for(int idx=0;idx<sc_length;idx++)
            {
                x_host[idx] = (T) drand48();
                y_host[idx] = (T) drand48();
            }
            
            device_.Memcpy(x, x_host, sc_length, MemcpyHostToDevice);
            device_.Memcpy(y, y_host, sc_length, MemcpyHostToDevice);
            device_.xy(x,y,stride,num_vecs,z);
            device_.Memcpy(z_host, z, sc_length, MemcpyDeviceToHost);
            
            T err = 0, diff;
            for(int idx=0;idx<sc_length;idx++)
            {
                diff = fabs(x_host[idx]*y_host[idx]-z_host[idx]);
                err = (diff>err) ? diff : err ;
            }
            device_.Free(x);
            device_.Free(y);
            device_.Free(z);
            free(x_host);
            free(y_host);
            free(z_host);
            
            
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

            T* x_host = (T*) malloc(sc_length * sizeof(T));
            T* y_host = (T*) malloc(sc_length * sizeof(T));
            T* z_host = (T*) malloc(sc_length * sizeof(T));
            
            for(int idx=0;idx<sc_length;idx++)
                x_host[idx] = (T) drand48();
                      
            device_.Memcpy(x, x_host, sc_length, MemcpyHostToDevice);
            device_.xInv(x,stride,num_vecs,y); 
            device_.xyInv(x,y,stride,num_vecs,z);
            device_.Memcpy(z_host, z, sc_length, MemcpyDeviceToHost);
            
            T err = 0;
            for(int idx=0;idx<sc_length;idx++)
                err = (z_host[idx]-1.0>err) ? z_host[idx]-1 : err ;
            
            device_.Free(x);
            device_.Free(y);
            device_.Free(z);
            
            free(x_host);
            free(y_host);
            free(z_host);

            res = res && (err<eps_) ? true : false;
            
            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* Device::xyInv : " + res_print + " *"<<endl;
        }
        return res;
     };

    bool TestuyInv()
    {
        bool res = true;
        {
            int stride = 413, num_vecs = 2;
            int vec_length = 3*stride*num_vecs;
            int sc_length = stride*num_vecs;
            T* u = device_.Malloc(vec_length);
            T* y = device_.Malloc(sc_length);

            T* u_host = (T*) malloc(vec_length * sizeof(T));
            T* y_host = (T*) malloc(sc_length * sizeof(T));
            
            for(int idx=0;idx<vec_length;idx++)
                u_host[idx] = (T) drand48();
            
            device_.Memcpy(u, u_host, vec_length, MemcpyHostToDevice);
            device_.DotProduct(u, u, stride, num_vecs, y);
            device_.Memcpy(y_host, y, sc_length, MemcpyDeviceToHost);
            
            for(int idx=0;idx<sc_length;idx++)
                y_host[idx] = sqrt(y_host[idx]);

            device_.Memcpy(y, y_host, sc_length, MemcpyHostToDevice);
            device_.uyInv(u,y,stride,num_vecs,u);
            device_.DotProduct(u,u,stride, num_vecs,y);
            device_.Memcpy(y_host, y, sc_length, MemcpyDeviceToHost);
            
            T err = 0;
            for(int idx=0;idx<sc_length;idx++)
                err = (y_host[idx]-1.0>err) ? y_host[idx]-1 : err ;
            
            device_.Free(u);
            device_.Free(y);
            
            free(u_host);
            free(y_host);
            
            res = res && (err<eps_) ? true : false;
            
            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* Device::uyInv : " + res_print + " *"<<endl;
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
            
            T* x_host = (T*) malloc(sc_length * sizeof(T));
            T* y_host = (T*) malloc(sc_length * sizeof(T));

            for(int idx=0;idx<sc_length;idx++)
                x_host[idx] = (T) drand48();
                        
            device_.Memcpy(x, x_host, sc_length, MemcpyHostToDevice);
            device_.axpy(-1.0,x,x,stride,num_vecs,y);
            device_.Memcpy(y_host, y, sc_length, MemcpyDeviceToHost);

            T err = 0;
            for(int idx=0;idx<sc_length;idx++)
                err = (y_host[idx]>err) ? y_host[idx] : err ;
            
            device_.Free(x);
            device_.Free(y);
            free(x_host);
            free(y_host);
            
            res = res && (err<eps_) ? true : false;
            
            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* Device::axpy: " + res_print + " *"<<endl;
        }
        return res;
         return true;
     }

     bool Testaxpb()
     {
        bool res = true;
        {
            int stride = 531, num_vecs = 3;
            int sc_length = stride*num_vecs;
            T* x = device_.Malloc(sc_length);
            T* y = device_.Malloc(sc_length);

            T* x_host = (T*) malloc(sc_length * sizeof(T));
            T* y_host = (T*) malloc(sc_length * sizeof(T));

            
            for(int idx=0;idx<sc_length;idx++)
                x_host[idx] = (T) drand48();
                        
            device_.Memcpy(x, x_host, sc_length, MemcpyHostToDevice);
            device_.axpb(-1.0,x,0.0,stride,num_vecs,y);
            device_.Memcpy(y_host, y, sc_length, MemcpyDeviceToHost);
            
            T err = 0;
            for(int idx=0;idx<sc_length;idx++)
                err = (x_host[idx]+y_host[idx]>err) ? (x_host[idx]+y_host[idx]) : err ;
            
            device_.Free(x);
            device_.Free(y);
            free(x_host);
            free(y_host);

            res = res && (err<eps_) ? true : false;
            
            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* Device::axpb : " + res_print + " *"<<endl;
        }
        return res;
        return true;
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

            device_.Memcpy(a, a_host, sc_length, MemcpyHostToDevice);
            device_.Memcpy(x, x_host, vec_length, MemcpyHostToDevice);
            device_.Memcpy(y, y_host, vec_length, MemcpyHostToDevice);

            device_.xvpb(a,x,0.0,stride,num_vecs,z);
            
            device_.DotProduct(x,y,stride,num_vecs,b);
            device_.DotProduct(z,y,stride,num_vecs,d);

            device_.Memcpy(b_host, b, sc_length, MemcpyDeviceToHost);
            device_.Memcpy(d_host, d, sc_length, MemcpyDeviceToHost);

            T err = 0, diff;
            for(int idx=0;idx<sc_length;idx++)
            {
                diff = fabs(a_host[idx]*b_host[idx]-d_host[idx]);
                err = (diff>err) ? diff : err ;
            }
            
            device_.xvpw(a,x,y,stride,num_vecs,z);
            device_.DotProduct(x,y,stride,num_vecs,b);
            device_.DotProduct(y,y,stride,num_vecs,c);
            device_.DotProduct(z,y,stride,num_vecs,d);
            
            device_.Memcpy(a_host, a, sc_length, MemcpyDeviceToHost);
            device_.Memcpy(b_host, b, sc_length, MemcpyDeviceToHost);
            device_.Memcpy(c_host, c, sc_length, MemcpyDeviceToHost);
            device_.Memcpy(d_host, d, sc_length, MemcpyDeviceToHost);

            for(int idx=0;idx<sc_length;idx++)
            {
                diff = fabs(a_host[idx]*b_host[idx]+c_host[idx]-d_host[idx]);
                err = (diff>err) ? diff : err ;
            }

            device_.Free(x);
            device_.Free(y);
            device_.Free(z);
            device_.Free(a);
            device_.Free(b);
            device_.Free(c);
            device_.Free(d);

            free(x_host);
            free(y_host);
            free(z_host);
            free(a_host);
            free(b_host);
            free(c_host);
            free(d_host);
            
            res = res && (err<eps_) ? true : false;
            
            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* Device::xvpw : " + res_print + " *"<<endl;
        }
        return res;
    }    

     bool TestShufflePoints()
     {
        int stride = 312, num_vecs = 8;
        int vec_length = 3*stride*num_vecs;

        T* x = device_.Malloc(vec_length);
        T* y = device_.Malloc(vec_length);
        
        T* x_host = (T*) malloc(vec_length * sizeof(T) );
        T* y_host = (T*) malloc(vec_length * sizeof(T) );

        int idx;
        for(int ii=0;ii<vec_length;++ii)
            x_host[ii] = ii;

        device_.Memcpy(x, x_host, vec_length, MemcpyHostToDevice);
        device_.ShufflePoints(x, AxisMajor, stride, num_vecs, y);
        device_.Memcpy(y_host, y, vec_length, MemcpyDeviceToHost);

        T diff, err=0;
        for(int ii=0;ii<num_vecs;++ii)
            for(int jj=0;jj<stride;++jj)
            {                
                idx = 3*ii*stride+3*jj;
                T val1 = ii*3*stride+jj;
                T val2 = ii*3*stride+jj+stride;
                T val3 = ii*3*stride+jj+2*stride;
                diff =  fabs(y_host[idx] -val1) + fabs(y_host[idx+1] -val2) + fabs(y_host[idx+2] -val3);
                err = (diff>err) ? diff : err ;
            }

        bool res = (err<eps_) ? true : false;

        device_.ShufflePoints(y, PointMajor, stride, num_vecs, x);
        device_.Memcpy(x_host, x, vec_length, MemcpyDeviceToHost);

        err=0;
        for(int jj=0;jj<3*num_vecs * stride;jj++)
            {     
                diff = fabs(x_host[jj] - jj);
                err = (diff>err) ? diff : err ;
            }
        
        res = res && (err<eps_) ? true : false;
        string res_print = (res) ? "Passed" : "Failed";
        cout<<"* Device::ShufflePoints : " + res_print + " *"<<endl;

        device_.Free(x);
        device_.Free(y);
        free(x_host);
        free(y_host);

        return res;
     }


    bool TestMax()
    {
        bool res = true;
        {
            int length = 10012;
            T* x = device_.Malloc(length);
            T* y = device_.Malloc(length);

            T* x_host = (T*) malloc(length * sizeof(T));
            T* y_host = (T*) malloc(length * sizeof(T));
            
            T max = 0;
            for(int idx=0;idx<length;idx++)
            {
                x_host[idx] = (T) drand48();
                y_host[idx] = 0.5*x_host[idx] - 1.0;
                max = (max > x_host[idx]) ? max : x_host[idx];
            }
            
            device_.Memcpy(x, x_host, length, MemcpyHostToDevice);
            device_.Memcpy(y, y_host, length, MemcpyHostToDevice);
            T mx = device_.Max(x,length);
            T my = device_.Max(y,length);
            
            device_.Free(x);
            device_.Free(y);
            free(x_host);
            free(y_host);
                        
            T err = fabs(mx-max) + fabs(.5*max - 1.0-my);
            res = res && (err<eps_) ? true : false;
            
            string res_print = (res) ? "Passed" : "Failed";
            cout<<"* Device::Max : " + res_print + " *"<<endl;
        }
        return res;
    }

    bool TestReduce()
    {
        bool res = true;
        
        int stride = 312;
        int ns = 5;
        int length = ns*stride;
        
        T *x = device_.Malloc(length);
        T *w = device_.Malloc(length);
        T *q = device_.Malloc(stride);
        T *I = device_.Malloc(ns);        
        
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

         device_.Memcpy(x,x_host,length,MemcpyHostToDevice);
         device_.Memcpy(w,w_host,length,MemcpyHostToDevice);
         device_.Memcpy(q,q_host,stride,MemcpyHostToDevice);
         
         device_.Reduce(x, w, q, stride, ns, I);
         device_.Memcpy(I_host,I,ns,MemcpyDeviceToHost);

         T err = 0;
         T II = (stride-1)*stride/16.0;
         for(int ii=0;ii<ns;++ii)
         {
             T diff = fabs(I_host[ii]-II);
             err = (err>diff) ? err : diff;
         }
         res = res && (err<eps_) ? true : false;


         device_.Reduce(NULL, w, q, stride, ns, I);
         device_.Memcpy(I_host,I,ns,MemcpyDeviceToHost);

         II = stride/8.0;
         for(int ii=0;ii<ns;++ii)
         {
             T diff = fabs(I_host[ii]-II);
             err = (err>diff) ? err : diff;
         }
         res = res && (err<eps_) ? true : false;

         device_.Reduce(w, x, q, stride, ns, I);
         device_.Memcpy(I_host,I,ns,MemcpyDeviceToHost);

         II = (stride-1)*stride/16.0;
         for(int ii=0;ii<ns;++ii)
         {
             T diff = fabs(I_host[ii]-II);
             err = (err>diff) ? err : diff;
         }
         res = res && (err<eps_) ? true : false;

        free(x_host);
        free(w_host);
        free(q_host);
        free(I_host);

        device_.Free(x);
        device_.Free(w);
        device_.Free(q);
        device_.Free(I);
        
        string res_print = (res) ? "Passed" : "Failed";
        cout<<"* Device::Reduce : " + res_print + " *"<<endl;
        return res;
    }

    bool TestResample()
    {
//         int p = 12;
//         int q = 16;
//         int n_fun = 3;
//         int lp = p * (p + 2) * n_fun;
//         int lq = q * (q + 2) * n_fun;
//         int idx = 0;
    
//         T *sp_host = (T*) malloc(lp * sizeof(T));
//         T *sq_host = (T*) malloc(lq * sizeof(T));
//         T *al_host = (T*) malloc(p * (p + 2) *sizeof(T));

//         T *sp = device_.Malloc(lp);
//         T *sq = device_.Malloc(lq);
//         T *al = device_.Malloc(p * (p + 2));

//         for(int ff=0; ff< n_fun; ++ff)
//             for(int ii=0; ii< 2 * p; ++ii)
//             {
//                 int len = p + 1 - (ii+1)/2;
//                 for(int jj=0; jj < len; ++jj)
//                     sp_host[idx++] = jj;
//             }

//         idx = 0;
//         for(int ii=0; ii< 2 * p; ++ii)
//         {
//             int len = p + 1 - (ii+1)/2;
//             for(int jj=0; jj < len; ++jj)
//                 al_host[idx++] = (len-jj)<=(p-2*p/3) ? 0 : 1;
//         }

//         device_.Memcpy(sp, sp_host, lp, MemcpyHostToDevice);
//         device_.Memcpy(al, al_host, p * (p + 2), MemcpyHostToDevice);
//         device_.ScaleFreqs(p, n_fun, sp, al, sp);
//         device_.Resample(p, n_fun, q, sp, sq);
//         device_.Memcpy(sq_host, sq, lq, MemcpyDeviceToHost);
    
//         idx = 0;
//         for(int ff=0; ff< n_fun; ++ff)
//             for(int ii=0; ii< 2 * q; ++ii)
//             {
//                 int len = q + 1 - (ii+1)/2;
//                 for(int jj=0; jj < len; ++jj)
//                     cout<<sq_host[idx++]<<" ";
//                 cout<<endl;
//             }
        
//         device_.Free(al);
//         device_.Free(sp);
//         device_.Free(sq);

//         free(sp_host);
//         free(sq_host);
//         free(al_host);

        return false;
    }
};

