/**
 * @file   DeviceTest.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Sun Feb 28 12:54:20 2010
 * 
 * @brief The tester class for the device class and its subclasses.
 */

#include <iostream>
#include "Device.h"
#include <math.h>

using namespace std;

class DeviceTest
{
  private:
    DeviceCPU<float> cpu_device_;
    
  public:
    bool PerformAll()
    {
        bool test_result;
        test_result = TestMalloc();// && TestCalloc() && TestMemcpy();
        // && TestDotProduct() && TestCrossProduct()
//             && TestxInv() && TestAxPyScA() && TestAxPyScAScY() && TestAxPy() && TestxTy();
        
        
        string res_print = (test_result) ? "Passed" : "Failed";
        cout<<"\n*** Device Class tests: " + res_print + " ***\n"<<endl;

        return test_result;
    };

    bool TestMalloc()
    {
        int arr_size = 1e6;
        float* a = (float*) cpu_device_.Malloc(arr_size * sizeof(float));
        
        int idx;
        for(idx=0;idx<arr_size;++idx)
            a[idx] = idx;
        //cpu_device_.Free(a);
        
        cout<<"* DeviceCPU::Malloc: Passed *"<<endl;
        return true;
    };

//     bool TestCalloc()
//     {
//         int arr_size = 1e6;
//         float* a = (float*) cpu_device_.Calloc(arr_size, sizeof(float));
     
//         bool res;
//         for(int idx=0;idx<arr_size;++idx)
//             res = (a[idx] == 0) ? true : false;

//         string res_print = (res) ? "Passed" : "Failed";
//         cpu_device_.Free(a);
        
//         cout<<"* DeviceCPU::Calloc: " + res_print + " *"<<endl;
//         return (res);
//     };

//     bool TestMemcpy()
//     {
//         int arr_size = 1e6;
//         float* a = (float*) cpu_device_.Malloc(arr_size * sizeof(float));
        
//         int idx;
//         for(idx=0;idx<arr_size;++idx)
//             a[idx] = idx;
        
//         float* b = (float*) cpu_device_.Calloc(arr_size, sizeof(float));

//         cpu_device_.Memcpy(b, a, arr_size * sizeof(float));

//         bool res;
//         for(idx=0;idx<arr_size;++idx)
//             res = (a[idx] == b[idx]) ? true : false;

//         cpu_device_.Free(a);
//         cpu_device_.Free(b);
        
//         string res_print = (res) ? "Passed" : "Failed";
//         cout<<"* DeviceCPU::Memcpy: " + res_print + " *"<<endl;
//         return res;
//     };

//     bool TestDotProduct()
//     {

//         bool res = true;
//         srand48(time(0));

//         {//Orthogonality (random single)
//             int stride = 312, num_vecs = 1;
//             int sc_length = stride*num_vecs;
//             int arr_length = 3*num_vecs*stride;
//             float* a = (float*) cpu_device_.Malloc(arr_length * sizeof(float));
//             float* b = (float*) cpu_device_.Malloc(arr_length * sizeof(float));
//             float* c = (float*) cpu_device_.Malloc(sc_length  * sizeof(float));
            
//             // a dot b should be zero since they are geometrically orthogonal 
//             float rnum;
//             for(int idx=0;idx<stride;idx++)
//             {
//                 rnum = (float) drand48();
                
//                 a[idx] = rnum; a[idx+stride] = 1   ; a[idx+2*stride] = rnum;
//                 b[idx] = 1   ; b[idx+stride] = rnum; b[idx+2*stride] = -2;
//                 c[idx] = 1;
//             }
            
//             cpu_device_.DotProduct(a,b,stride,num_vecs,c);
//             float err = 0;
//             for(int idx=0;idx<sc_length;idx++)
//                 err = (c[idx]>err) ? c[idx] : err ;
            
//             cpu_device_.Free(a);
//             cpu_device_.Free(b);
//             cpu_device_.Free(c);
            
//             res = res && (err<1e-7) ? true : false;
            
//             string res_print = (res) ? "Passed" : "Failed";
//             cout<<"* DeviceCPU::DotProduct (Orthogonality, single): " + res_print + " *"<<endl;
//         }

// //         {//Orthogonality (random, double)
// //             int stride = 312, num_vecs = 1;
// //             int sc_length = stride*num_vecs;
// //             int arr_length = 3*num_vecs*stride;
// //             double* a = (double*) cpu_device_.Malloc(arr_length * sizeof(double));
// //             double* b = (double*) cpu_device_.Malloc(arr_length * sizeof(double));
// //             double* c = (double*) cpu_device_.Malloc(sc_length  * sizeof(double));
            
// //             // a dot b should be zero since they are geometrically orthogonal 
// //             double rnum;
// //             for(int idx=0;idx<stride;idx++)
// //             {
// //                 rnum = drand48();
                
// //                 a[idx] = rnum; a[idx+stride] = 1.0 ; a[idx+2*stride] = rnum;
// //                 b[idx] = 1.0 ; b[idx+stride] = rnum; b[idx+2*stride] = -2.0;
// //                 c[idx] = 1.0;
// //             }
            
// //             cpu_device_.DotProduct(a,b,stride,num_vecs,c);
// //             double err = 0;
// //             for(int idx=0;idx<sc_length;idx++)
// //                 err = (c[idx]>err) ? c[idx] : err ;
                        
// //             cpu_device_.Free(a);
// //             cpu_device_.Free(b);
// //             cpu_device_.Free(c);
            
// //             res = res && (err<1e-14) ? true : false;
            
// //             string res_print = (res) ? "Passed" : "Failed";
// //             cout<<"* DeviceCPU::DotProduct (Orthogonality, double): " + res_print + " *"<<endl;
// //         }

//         {//Normalization (random, single)
//             int stride = 401, num_vecs = 1;
//             int sc_length = stride*num_vecs;
//             int arr_length = 3*num_vecs*stride;
//             float* a = (float*) cpu_device_.Malloc(arr_length * sizeof(float));
//             float* c = (float*) cpu_device_.Malloc(sc_length  * sizeof(float));

//             for(int idx=0;idx<stride;idx++)
//             {                
//                 a[idx              ] = drand48();
//                 a[idx+stride       ] = drand48();
//                 a[idx+stride+stride] = drand48();
//                 c[idx] = 1.0;
//             }
            
//             cpu_device_.DotProduct(a,a,stride,num_vecs,c);

//             float nn;
//             for(int idx=0;idx<sc_length;idx++)
//             {
//                 nn = sqrt(c[idx]);
//                 a[idx              ] /= nn;
//                 a[idx+stride       ] /= nn;
//                 a[idx+stride+stride] /= nn;
//             }
                        
//             cpu_device_.DotProduct(a,a,stride,num_vecs,c);

//             float err = 0;
//             for(int idx=0;idx<sc_length;idx++)
//                 err = (c[idx]-1>err) ? c[idx]-1 : err ;
      
//             cpu_device_.Free(a);
//             cpu_device_.Free(c);
            
//             res = res && (err<1e-6) ? true : false;
            
//             string res_print = (res) ? "Passed" : "Failed";
//             cout<<"* DeviceCPU::DotProduct (Normalization, single): " + res_print + " *"<<endl;
//         }

// //         {//Normalization (random, double)
// //             int stride = 401, num_vecs = 1;
// //             int sc_length = stride*num_vecs;
// //             int arr_length = 3*num_vecs*stride;
// //             double* a = (double*) cpu_device_.Malloc(arr_length * sizeof(double));
// //             double* c = (double*) cpu_device_.Malloc(sc_length  * sizeof(double));

// //             for(int idx=0;idx<stride;idx++)
// //             {                
// //                 a[idx              ] = drand48();
// //                 a[idx+stride       ] = drand48();
// //                 a[idx+stride+stride] = drand48();
// //                 c[idx] = 1.0;
// //             }
            
// //             cpu_device_.DotProduct(a,a,stride,num_vecs,c);

// //             double nn;
// //             for(int idx=0;idx<sc_length;idx++)
// //             {
// //                 nn = sqrt(c[idx]);
// //                 a[idx              ] /= nn;
// //                 a[idx+stride       ] /= nn;
// //                 a[idx+stride+stride] /= nn;
// //             }
                        
// //             cpu_device_.DotProduct(a,a,stride,num_vecs,c);

// //             double err = 0;
// //             for(int idx=0;idx<sc_length;idx++)
// //                 err = (c[idx]-1>err) ? c[idx]-1 : err ;
      
// //             cpu_device_.Free(a);
// //             cpu_device_.Free(c);
            
// //             res = res && (err<1e-14) ? true : false;
            
// //             string res_print = (res) ? "Passed" : "Failed";
// //             cout<<"* DeviceCPU::DotProduct (Normalization, double): " + res_print + " *"<<endl;
// //         }
        
//         return res;
//     };

//     bool TestCrossProduct()
//     {
//         bool res = true;
//         {//Self-product (single)
//             int stride = 312, num_vecs = 1;
//             int arr_length = 3*num_vecs*stride;
//             float* a = (float*) cpu_device_.Malloc(arr_length * sizeof(float));
//             float* c = (float*) cpu_device_.Malloc(arr_length * sizeof(float));
            
            
//             for(int idx=0;idx<stride;idx++)
//             {                
//                 a[idx              ] = drand48();
//                 a[idx+stride       ] = drand48();
//                 a[idx+stride+stride] = drand48();
//             }

//             cpu_device_.CrossProduct(a,a,stride,num_vecs,c);
//             float err = 0;
//             for(int idx=0;idx<stride*num_vecs;idx++)
//                 err = (c[idx]>err) ? c[idx] : err ;
            
//             cpu_device_.Free(a);
//             cpu_device_.Free(c);
            
//             res = res && (err<1e-14) ? true : false;
            
//             string res_print = (res) ? "Passed" : "Failed";
//             cout<<"* DeviceCPU::CrossProduct (Self-product, float): " + res_print + " *"<<endl;
//         }

// //         {//Self-product (double)
// //             int stride = 312, num_vecs = 1;
// //             int arr_length = 3*num_vecs*stride;
// //             double* a = (double*) cpu_device_.Malloc(arr_length * sizeof(double));
// //             double* c = (double*) cpu_device_.Malloc(arr_length * sizeof(double));
            
            
// //             for(int idx=0;idx<stride;idx++)
// //             {                
// //                 a[idx              ] = drand48();
// //                 a[idx+stride       ] = drand48();
// //                 a[idx+stride+stride] = drand48();
// //             }

// //             cpu_device_.CrossProduct(a,a,stride,num_vecs,c);
// //             double err = 0;
// //             for(int idx=0;idx<stride*num_vecs;idx++)
// //                 err = (c[idx]>err) ? c[idx] : err ;
            
// //             cpu_device_.Free(a);
// //             cpu_device_.Free(c);
            
// //             res = res && (err<1e-14) ? true : false;
            
// //             string res_print = (res) ? "Passed" : "Failed";
// //             cout<<"* DeviceCPU::CrossProduct (Self-product, double): " + res_print + " *"<<endl;
// //         }

        
//         {//triple product (single)
//             int stride = 723, num_vecs = 1;
//             int arr_length = 3*num_vecs*stride;
//             int sc_length = num_vecs*stride;
            
//             float* a = (float*) cpu_device_.Malloc(arr_length * sizeof(float));
//             float* b = (float*) cpu_device_.Malloc(arr_length * sizeof(float));
//             float* c = (float*) cpu_device_.Malloc(arr_length * sizeof(float));
//             float* d = (float*) cpu_device_.Malloc(arr_length * sizeof(float));
            
//             float* e = (float*) cpu_device_.Malloc(sc_length  * sizeof(float));
//             float* f = (float*) cpu_device_.Malloc(sc_length  * sizeof(float));            

//             for(int idx=0;idx<stride;idx++)
//             {                
//                 a[idx              ] = drand48();
//                 a[idx+stride       ] = drand48();
//                 a[idx+stride+stride] = drand48();

//                 b[idx              ] = drand48();
//                 b[idx+stride       ] = drand48();
//                 b[idx+stride+stride] = drand48();
                
//                 c[idx              ] = drand48();
//                 c[idx+stride       ] = drand48();
//                 c[idx+stride+stride] = drand48();
//             }

//             // (a x b).c
//             cpu_device_.CrossProduct(a,b,stride,num_vecs,d);
//             cpu_device_.DotProduct(d,c,stride,num_vecs,e);

//             // (b x a).c
//             cpu_device_.CrossProduct(b,a,stride,num_vecs,d);
//             cpu_device_.DotProduct(d,c,stride,num_vecs,f);

//             float err = 0;
//             for(int idx=0;idx<sc_length;idx++)
//                 err = (e[idx]+f[idx]>err) ? e[idx]+f[idx] : err ;

//             // (c x a).b
//             cpu_device_.CrossProduct(c,a,stride,num_vecs,d);
//             cpu_device_.DotProduct(d,b,stride,num_vecs,f);
            
//             for(int idx=0;idx<sc_length;idx++)
//                 err = ((e[idx]-f[idx])>err) ? e[idx]-f[idx] : err ;

//             cpu_device_.Free(a);
//             cpu_device_.Free(b);
//             cpu_device_.Free(c);
//             cpu_device_.Free(d);
//             cpu_device_.Free(e);
//             cpu_device_.Free(f);

//             res = res && (err<1e-6) ? true : false;
            
//             string res_print = (res) ? "Passed" : "Failed";
//             cout<<"* DeviceCPU::CrossProduct (Triple product, float): " + res_print + " *"<<endl;
//         }

// //         {//triple product (double)
// //             int stride = 723, num_vecs = 1;
// //             int arr_length = 3*num_vecs*stride;
// //             int sc_length = num_vecs*stride;
            
// //             double* a = (double*) cpu_device_.Malloc(arr_length * sizeof(double));
// //             double* b = (double*) cpu_device_.Malloc(arr_length * sizeof(double));
// //             double* c = (double*) cpu_device_.Malloc(arr_length * sizeof(double));
// //             double* d = (double*) cpu_device_.Malloc(arr_length * sizeof(double));
            
// //             double* e = (double*) cpu_device_.Malloc(sc_length  * sizeof(double));
// //             double* f = (double*) cpu_device_.Malloc(sc_length  * sizeof(double));            

// //             for(int idx=0;idx<stride;idx++)
// //             {                
// //                 a[idx              ] = drand48();
// //                 a[idx+stride       ] = drand48();
// //                 a[idx+stride+stride] = drand48();

// //                 b[idx              ] = drand48();
// //                 b[idx+stride       ] = drand48();
// //                 b[idx+stride+stride] = drand48();
                
// //                 c[idx              ] = drand48();
// //                 c[idx+stride       ] = drand48();
// //                 c[idx+stride+stride] = drand48();
// //             }

// //             // (a x b).c
// //             cpu_device_.CrossProduct(a,b,stride,num_vecs,d);
// //             cpu_device_.DotProduct(d,c,stride,num_vecs,e);

// //             // (b x a).c
// //             cpu_device_.CrossProduct(b,a,stride,num_vecs,d);
// //             cpu_device_.DotProduct(d,c,stride,num_vecs,f);

// //             double err = 0;
// //             for(int idx=0;idx<sc_length;idx++)
// //                 err = (e[idx]+f[idx]>err) ? e[idx]+f[idx] : err ;

// //             // (c x a).b
// //             cpu_device_.CrossProduct(c,a,stride,num_vecs,d);
// //             cpu_device_.DotProduct(d,b,stride,num_vecs,f);
            
// //             for(int idx=0;idx<sc_length;idx++)
// //                 err = ((e[idx]-f[idx])>err) ? e[idx]-f[idx] : err ;

// //             cpu_device_.Free(a);
// //             cpu_device_.Free(b);
// //             cpu_device_.Free(c);
// //             cpu_device_.Free(d);
// //             cpu_device_.Free(e);
// //             cpu_device_.Free(f);

// //             res = res && (err<1e-14) ? true : false;
            
// //             string res_print = (res) ? "Passed" : "Failed";
// //             cout<<"* DeviceCPU::CrossProduct (Triple product, double): " + res_print + " *"<<endl;
// //         }

//         return res;
//     };

//     bool TestxInv()
//     {
//         bool res = true;
//         {
//             int stride = 413, num_vecs = 1;
//             int sc_length = stride*num_vecs;
//             float* x = (float*) cpu_device_.Malloc(sc_length  * sizeof(float));
//             float* y = (float*) cpu_device_.Malloc(sc_length  * sizeof(float));
            
//             for(int idx=0;idx<sc_length;idx++)
//             x[idx] = (float) drand48();
                        
//             cpu_device_.xInv(x,stride,num_vecs,y);
//             cpu_device_.xInv(y,stride,num_vecs,y);
            
//             float err = 0;
//             for(int idx=0;idx<sc_length;idx++)
//                 err = (x[idx]-y[idx]>err) ? x[idx]-y[idx] : err ;
            
//             cpu_device_.Free(x);
//             cpu_device_.Free(y);
            
//             res = res && (err<1e-6) ? true : false;
            
//             string res_print = (res) ? "Passed" : "Failed";
//             cout<<"* DeviceCPU::xInv : " + res_print + " *"<<endl;
//         }
//         return res;
//     }

//     bool TestAxPyScA()
//     {
//         bool res = true;
//         {
//             int stride = 531, num_vecs = 3;
//             int sc_length = stride*num_vecs;
//             float* x = (float*) cpu_device_.Malloc(sc_length  * sizeof(float));
//             float* y = (float*) cpu_device_.Malloc(sc_length  * sizeof(float));
            
//             for(int idx=0;idx<sc_length;idx++)
//                 x[idx] = (float) drand48();
                        
//             cpu_device_.AxPy(-1.0,x,x,stride,num_vecs,y);
            
//             float err = 0;
//             for(int idx=0;idx<sc_length;idx++)
//                 err = (y[idx]>err) ? y[idx] : err ;
            
//             cpu_device_.Free(x);
//             cpu_device_.Free(y);
            
//             res = res && (err<1e-6) ? true : false;
            
//             string res_print = (res) ? "Passed" : "Failed";
//             cout<<"* DeviceCPU::AxPy (scalar field) : " + res_print + " *"<<endl;
//         }
//         return res;
//     }

//     bool TestAxPyScAScY()
//     {
//         bool res = true;
//         {
//             int stride = 531, num_vecs = 3;
//             int sc_length = stride*num_vecs;
//             float* x = (float*) cpu_device_.Malloc(sc_length  * sizeof(float));
//             float* y = (float*) cpu_device_.Malloc(sc_length  * sizeof(float));
            
//             for(int idx=0;idx<sc_length;idx++)
//                 x[idx] = (float) drand48();
                        
//             cpu_device_.AxPy(-1.0,x,0.0,stride,num_vecs,y);
            
//             float err = 0;
//             for(int idx=0;idx<sc_length;idx++)
//                 err = (x[idx]+y[idx]>err) ? (x[idx]+y[idx]) : err ;
            
//             cpu_device_.Free(x);
//             cpu_device_.Free(y);
            
//             res = res && (err<1e-6) ? true : false;
            
//             string res_print = (res) ? "Passed" : "Failed";
//             cout<<"* DeviceCPU::AxPy (scalar field, fixed y) : " + res_print + " *"<<endl;
//         }
//         return res;
//     }

//     bool TestAxPy()
//     {
//         bool res = true;
//         {
//             int stride = 1001, num_vecs = 7;
//             int sc_length = stride*num_vecs;
//             int vec_length = 3*sc_length;

//             float* x = (float*) cpu_device_.Malloc(vec_length  * sizeof(float));
//             float* y = (float*) cpu_device_.Malloc(vec_length  * sizeof(float));
//             float* z = (float*) cpu_device_.Malloc(vec_length  * sizeof(float));

//             float* a = (float*) cpu_device_.Malloc(sc_length   * sizeof(float));
//             float* b = (float*) cpu_device_.Malloc(sc_length   * sizeof(float));
//             float* c = (float*) cpu_device_.Malloc(sc_length   * sizeof(float));
//             float* d = (float*) cpu_device_.Malloc(sc_length   * sizeof(float));
            
//             for(int idx=0;idx<vec_length;idx++)
//             {
//                 x[idx] = (float) drand48();
//                 y[idx] = (float) drand48();
//             }

//             for(int idx=0;idx<sc_length;idx++)
//                 a[idx] = (float) drand48();

//             cpu_device_.AxPy(a,x,0.0,stride,num_vecs,z);
            
//             cpu_device_.DotProduct(x,y,stride,num_vecs,b);
//             cpu_device_.DotProduct(z,y,stride,num_vecs,d);
            
//             float err = 0, diff;
//             for(int idx=0;idx<sc_length;idx++)
//             {
//                 diff = fabs(a[idx]*b[idx]-d[idx]);
//                 err = (diff>err) ? diff : err ;
//             }
            
//             cpu_device_.AxPy(a,x,y,stride,num_vecs,z);
            
//             cpu_device_.DotProduct(x,y,stride,num_vecs,b);
//             cpu_device_.DotProduct(y,y,stride,num_vecs,c);
//             cpu_device_.DotProduct(z,y,stride,num_vecs,d);
            
//             for(int idx=0;idx<sc_length;idx++)
//             {
//                 diff = fabs(a[idx]*b[idx]+c[idx]-d[idx]);
//                 err = (diff>err) ? diff : err ;
//             }

//             cpu_device_.Free(x);
//             cpu_device_.Free(y);
//             cpu_device_.Free(z);
//             cpu_device_.Free(a);
//             cpu_device_.Free(b);
//             cpu_device_.Free(c);
//             cpu_device_.Free(d);
            
//             res = res && (err<1e-6) ? true : false;
            
//             string res_print = (res) ? "Passed" : "Failed";
//             cout<<"* DeviceCPU::AxPy : " + res_print + " *"<<endl;
//         }
//         return res;
//     }
    
//     bool TestxTy()
//     {
//         bool res = true;
//         {
//             int stride = 65, num_vecs = 4;
//             int sc_length = stride*num_vecs;
//             float* x = (float*) cpu_device_.Malloc(sc_length  * sizeof(float));
//             float* y = (float*) cpu_device_.Malloc(sc_length  * sizeof(float));
//             float* z = (float*) cpu_device_.Malloc(sc_length  * sizeof(float));
            
//             for(int idx=0;idx<sc_length;idx++)
//             {
//                 x[idx] = (float) drand48();
//                 y[idx] = (float) drand48();
//             }
                        
//             cpu_device_.xTy(x,y,stride,num_vecs,z);
            
//             float err = 0, diff;
//             for(int idx=0;idx<sc_length;idx++)
//             {
//                 diff = fabs(x[idx]*y[idx]-z[idx]);
//                 err = (diff>err) ? diff : err ;
//             }
//             cpu_device_.Free(x);
//             cpu_device_.Free(y);
//             cpu_device_.Free(z);
            
            
//             res = res && (err<1e-6) ? true : false;
            
//             string res_print = (res) ? "Passed" : "Failed";
//             cout<<"* DeviceCPU::xTy : " + res_print + " *"<<endl;
//         }
//         return res;
//     }
    
};

