#include <iostream>
#include "DeviceCPU.h"
#include "Scalars.h"
#include "Vectors.h"
#include "SHTrans.h"
#include "Surface.h"
#include "DataIO.h"
#include <math.h>

using namespace std;
typedef float T;

int main(int argc, char ** argv)
{
    DeviceCPU<T> cpu;
    DataIO<T> myIO;

    int p(12), nVec(1);
    
    int fLen(2*p*(p+1));
    int dLen(6*p*(p+1)*nVec);

    // reading data
    T* pos_arr = new T[dLen];
    myIO.ReadData("../data/dumbbell_cart12_single.txt",dLen,pos_arr);
    Vectors<T> pos_vec(cpu,p,nVec,pos_arr);
    delete []pos_arr;

    Surface<T> S(cpu,p,nVec,pos_vec);
    
    // Checking the grad and div operator
    Scalars<T> div_n(cpu,p,nVec);
    S.SurfDiv(S.normal_, div_n);
    axpy((T) .5, div_n, S.h_,div_n);
        
    T err = 0, err2;
    for(int ii=0;ii<div_n.GetDataLength(); ii++)
    {
        err2 = fabs(div_n.data_[ii]);
        err = (err>err2) ? err : err2;
    }
    cout<<"\n The error in the surface divergence (For the dumbbell .09313 expected)= "<<err<<endl;
    
    Scalars<T> X(cpu,p,nVec,S.x_.GetFunctionAt(0));
    Scalars<T> Y(cpu,p,nVec,S.x_.GetFunctionAt(1));
    Scalars<T> Z(cpu,p,nVec,S.x_.GetFunctionAt(2));

    Vectors<T> grad(cpu,p,nVec), lap(cpu,p,nVec);
    Scalars<T> lap_x(cpu,p,nVec), lap_y(cpu,p,nVec), lap_z(cpu,p,nVec);

    S.SurfGrad(X,grad);
    S.SurfDiv(grad,lap_x);
    lap.SetFunctionAt(lap_x.data_,0);
    
    S.SurfGrad(Y,grad);
    S.SurfDiv(grad,lap_y);
    lap.SetFunctionAt(lap_y.data_,1);
    
    S.SurfGrad(Z,grad);
    S.SurfDiv(grad,lap_z);
    lap.SetFunctionAt(lap_z.data_,2);

    Scalars<T> hh(cpu, p, nVec);
    DotProduct(lap,S.normal_,hh);
    axpy((T) -.5, hh, S.h_,hh);
        
    err = 0;
    for(int ii=0;ii<hh.GetDataLength(); ii++)
    {
        err2 = fabs(hh.data_[ii]);
        err = (err>err2) ? err : err2;
    }
    cout<<"\n The error in the surface grad (For the dumbbell .13703 expected)= "<<err<<endl;

    return 0;
}


//     }

    

//     T ts=.01;
//     int m=20;

//     if(argc>1)
//         ts = (T) atof(argv[1]);
 
//     if(argc>2)
//         m = atoi(argv[2]);
    
//     cout<<"ts = "<<ts<<", m = "<<m<<endl;
//     // reading data
//     T* pos_arr = new T[vec_len];
//     myIO.ReadData("../data/dumbbell_cart12",vec_len,pos_arr);
    
//     // Checking the grad and div operator
//     Scalars<T> div_n(p,nVec);
//     S.SurfDiv(S.normal_, div_n);
//     AxPy((T) -.5, div_n, (T) 0.0,div_n);
    
//   //   err = 0;
//     for(int ii=0;ii<div_n.GetDataLength(); ii+=2*p)
//         cout<<div_n.data_[ii]<<endl;
// //     {
// //         err2 = fabs(div_n.data_[ii]);
// //         err = (err>err2) ? err : err2;
// //     }
    

    
//     //

//     for(int ii=1;ii<nVec;++ii)
//         for(int jj=0;jj<np;++jj)
//         {
//             pos_arr[ii*3*np+jj     ] = 1.0 + pos_arr[jj     ];
//             pos_arr[ii*3*np+jj+  np] = 2.0 + pos_arr[jj+  np];
//             pos_arr[ii*3*np+jj+2*np] = 3.0 + pos_arr[jj+2*np];
//         }
    
//     Vectors<T> pos_vec(p,nVec,pos_arr);
//     delete []pos_arr;
    
//     // initializing Vesicle
//     Surface<T> S(p,nVec,pos_vec);
//     Vectors<T> bend_force(p,nVec);

//     for(int tt=0;tt<m;++tt)
//     {
//         AxPy(S.h_, S.normal_, (float) 0.0, bend_force);
//         //Add pressure to preserve volume

//         AxPy(ts,bend_force, S.x_, pos_vec);
//         S.SetX(pos_vec);
//         cout<<tt<<" ";
//     }
//     cout<<endl;
    
//     T err = 0, err2;
//     for(int ii=1;ii<nVec;++ii)
//         for(int jj=0;jj<np;++jj)
//         {
//             err2 = fabs(S.x_.data_[ii*3*np+jj     ] - S.x_.data_[jj     ] - 1.0);
//             err2+= fabs(S.x_.data_[ii*3*np+jj+  np] - S.x_.data_[jj+  np] - 2.0);
//             err2+= fabs(S.x_.data_[ii*3*np+jj+2*np] - S.x_.data_[jj+2*np] - 3.0);
            
//             err = (err>err2) ? err : err2;
//         }

//     cout<<err<<endl;
//     myIO.WriteData("./cart12_final",3*np,S.x_.data_);

