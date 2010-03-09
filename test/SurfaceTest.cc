#include <iostream>
#include "DeviceCPU.h"
#include "Scalars.h"
#include "Vectors.h"
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
    cpu.InitializeSHT(p, "../data/legTrans12_single.txt",
            "../data/legTransInv12_single.txt",
            "../data/d1legTrans12_single.txt",
            "../data/d2legTrans12_single.txt");

    int fLen(2*p*(p+1));
    int dLen(6*p*(p+1)*nVec);

    // reading data
    Surface<T> S(cpu,p,nVec);
    myIO.ReadData("../data/dumbbell_cart12_single.txt",dLen,S.x_.data_);
    S.UpdateProps();
    
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
