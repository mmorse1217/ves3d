#include <iostream>
#include "DeviceCPU.h"
#include "Scalars.h"
#include "Vectors.h"
#include "SHTrans.h"
#include "Surface.h"
//#include "DataIO.h"
#include <math.h>

using namespace std;
typedef float T;

int main(int argc, char ** argv)
{
//     int p=12;
//     int n_surfs = 1;
    
//     DeviceCPU<float> cpu1, cpu2;
//     Scalars<float> sc(cpu1, p, n_surfs);

//     SHTrans<float> diff(cpu1, p, n_surfs);
    
//     diff.AllDerivatives(sc,sc,sc,sc,sc,sc);
//     diff.FirstDerivatives(sc,sc,sc);


    DeviceCPU<T> cpu;
    
    int p(12), nVec(10000);
    
    int fLen(2*p*(p+1));
    int dLen(6*p*(p+1)*nVec);
    
    T *data_in = new T[dLen];
    for(int idx=0;idx<dLen;++idx)
        data_in[idx] = idx;
    
    Vectors<T> pos_vec(cpu,p,nVec,data_in);
    delete[] data_in;

    Surface<T> S(cpu,p,nVec,pos_vec);

    return 0;
}

     

//     {
//         int p = 12;
//         int np = 2*p*(p+1);
//         int num_vecs = 1;
//         int vec_len = 3*num_vecs*np;
//         DataIO<scType> myIO;
        
//         // reading data
//         scType* pos_arr = new scType[vec_len];
//         myIO.ReadData("../data/dumbbell_cart12",vec_len,pos_arr);
    
//         Vectors<scType> pos_vec(p,num_vecs,pos_arr);
//         delete []pos_arr;
    
//         // initializing Vesicle
//         Surface<scType> S(p,num_vecs,pos_vec);
        
//         // Checking the grad and div operator
//         Scalars<scType> div_n(p,num_vecs);
//         S.SurfDiv(S.normal_, div_n);
//         AxPy((scType) .5, div_n, S.h_,div_n);
        
//         scType err = 0, err2;
//         for(int ii=0;ii<div_n.GetDataLength(); ii++)
//         {
//             err2 = fabs(div_n.data_[ii]);
//             err = (err>err2) ? err : err2;
//         }
//         cout<<err<<endl;

//         Scalars<scType> X(p,num_vecs,S.x_.GetFunctionAt(0));
//         Scalars<scType> Y(p,num_vecs,S.x_.GetFunctionAt(1));
//         Scalars<scType> Z(p,num_vecs,S.x_.GetFunctionAt(2));

//         Vectors<scType> grad(p,num_vecs), lap(p,num_vecs);
//         Scalars<scType> lap_x(p,num_vecs), lap_y(p,num_vecs), lap_z(p,num_vecs);

//         S.SurfGrad(X,grad);
//         S.SurfDiv(grad,lap_x);
//         lap.SetFunctionAt(lap_x.data_,0);

//         S.SurfGrad(Y,grad);
//         S.SurfDiv(grad,lap_y);
//         lap.SetFunctionAt(lap_y.data_,1);

//         S.SurfGrad(Z,grad);
//         S.SurfDiv(grad,lap_z);
//         lap.SetFunctionAt(lap_z.data_,2);

//         Scalars<scType> hh(p, num_vecs);
//         DotProduct(lap,S.normal_,hh);
//         AxPy((scType) -.5, hh, S.h_,hh);
        
//         err = 0;
//         for(int ii=0;ii<hh.GetDataLength(); ii++)
//         {
//             err2 = fabs(hh.data_[ii]);
//             err = (err>err2) ? err : err2;
//         }
//         cout<<err<<endl;

//     }

    

//     scType ts=.01;
//     int m=20;

//     if(argc>1)
//         ts = (scType) atof(argv[1]);
 
//     if(argc>2)
//         m = atoi(argv[2]);
    
//     cout<<"ts = "<<ts<<", m = "<<m<<endl;
//     // reading data
//     scType* pos_arr = new scType[vec_len];
//     myIO.ReadData("../data/dumbbell_cart12",vec_len,pos_arr);
    
//     // Checking the grad and div operator
//     Scalars<scType> div_n(p,num_vecs);
//     S.SurfDiv(S.normal_, div_n);
//     AxPy((scType) -.5, div_n, (scType) 0.0,div_n);
    
//   //   err = 0;
//     for(int ii=0;ii<div_n.GetDataLength(); ii+=2*p)
//         cout<<div_n.data_[ii]<<endl;
// //     {
// //         err2 = fabs(div_n.data_[ii]);
// //         err = (err>err2) ? err : err2;
// //     }
    

    
//     //

//     for(int ii=1;ii<num_vecs;++ii)
//         for(int jj=0;jj<np;++jj)
//         {
//             pos_arr[ii*3*np+jj     ] = 1.0 + pos_arr[jj     ];
//             pos_arr[ii*3*np+jj+  np] = 2.0 + pos_arr[jj+  np];
//             pos_arr[ii*3*np+jj+2*np] = 3.0 + pos_arr[jj+2*np];
//         }
    
//     Vectors<scType> pos_vec(p,num_vecs,pos_arr);
//     delete []pos_arr;
    
//     // initializing Vesicle
//     Surface<scType> S(p,num_vecs,pos_vec);
//     Vectors<scType> bend_force(p,num_vecs);

//     for(int tt=0;tt<m;++tt)
//     {
//         AxPy(S.h_, S.normal_, (float) 0.0, bend_force);
//         //Add pressure to preserve volume

//         AxPy(ts,bend_force, S.x_, pos_vec);
//         S.SetX(pos_vec);
//         cout<<tt<<" ";
//     }
//     cout<<endl;
    
//     scType err = 0, err2;
//     for(int ii=1;ii<num_vecs;++ii)
//         for(int jj=0;jj<np;++jj)
//         {
//             err2 = fabs(S.x_.data_[ii*3*np+jj     ] - S.x_.data_[jj     ] - 1.0);
//             err2+= fabs(S.x_.data_[ii*3*np+jj+  np] - S.x_.data_[jj+  np] - 2.0);
//             err2+= fabs(S.x_.data_[ii*3*np+jj+2*np] - S.x_.data_[jj+2*np] - 3.0);
            
//             err = (err>err2) ? err : err2;
//         }

//     cout<<err<<endl;
//     myIO.WriteData("./cart12_final",3*np,S.x_.data_);

