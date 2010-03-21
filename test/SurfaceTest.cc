 #include <iostream>
#include "DeviceCPU.h"
#include "DataIO.h"
#include "Scalars.h"
#include "Vectors.h"
#include "Surface.h"
#include <math.h>

using namespace std;
typedef float T;

int main(int argc, char ** argv)
{
    DeviceCPU<T> cpu;
    DataIO<T> myIO;

    int p(12), nVec(5);
    cpu.InitializeSHT(p,2*p);

    int fLen(2*p*(p+1));
    int dLen(6*p*(p+1));

    //Setting surface parameters
    SurfaceParams<T> par;
    
    par.p_ = p;
    par.n_surfs_ = 1;
    par.kappa_ = 1e-2;
    par.filter_freq_ = 8;
    par.rep_ts_ = 1e-1;
    par.rep_max_vel_ = 1e-1;
    par.rep_iter_max_ = 100;
    par.rep_up_freq_ = 24;
    par.rep_filter_freq_ = 4;

    // memory allocation
    Surface<T> S(cpu,par);
    Scalars<T> X(cpu,p,nVec);
    Scalars<T> Y(cpu,p,nVec);
    Scalars<T> Z(cpu,p,nVec);
    
    // initializing vesicle positions from text file
    myIO.ReadData("../data/dumbbell_cart12_single.txt",dLen,S.x_.data_);

    S.Resize(nVec);
    for(int ii=1;ii<nVec;ii++)
        for(int idx=0;idx<dLen;idx++)
            S.x_.data_[ ii*dLen + idx] = S.x_.data_[idx];
    S.UpdateAll();
    
    for(int ii=0;ii<nVec;ii++)
        for(int idx=0;idx<fLen;idx++)
        {
            X.data_[ ii*fLen + idx] = S.x_.data_[idx];
            Y.data_[ ii*fLen + idx] = S.x_.data_[idx + fLen];
            Z.data_[ ii*fLen + idx] = S.x_.data_[idx + fLen + fLen];
        }
    
    // Checking the grad and div operator
    S.Resize(3); 
    Scalars<T> div_n(cpu,p,3);
    S.SurfDiv(S.normal_, div_n);
    axpy((T) .5, div_n, S.h_,div_n);
        
    T err = 0, err2;
    for(int ii=0;ii<div_n.GetDataLength(); ii++)
    {
        err2 = fabs(div_n.data_[ii]);
        err = (err>err2) ? err : err2;
    }
    cout<<"\n The error in the surface divergence (For the dumbbell .02964 expected)= "<<err<<endl;
    
    S.Resize(nVec);
    for(int ii=1;ii<nVec;ii++)
        for(int idx=0;idx<dLen;idx++)
            S.x_.data_[ ii*dLen + idx] = S.x_.data_[idx];

    dLen *=nVec;
    S.UpdateAll();
    Vectors<T> grad(cpu,p,nVec), lap(cpu,p,nVec);

    S.SurfGrad(X,grad);
    S.SurfDiv(grad,X);
    
    S.SurfGrad(Y,grad);
    S.SurfDiv(grad,Y);
    
    S.SurfGrad(Z,grad);
    S.SurfDiv(grad,Z);
    
    for(int ii=0;ii<nVec;ii++)
        for(int idx=0;idx<fLen;idx++)
        {
            lap.data_[ 3*ii*fLen + idx              ] = X.data_[ii*fLen + idx];
            lap.data_[ 3*ii*fLen + idx + fLen       ] = Y.data_[ii*fLen + idx];
            lap.data_[ 3*ii*fLen + idx + fLen + fLen] = Z.data_[ii*fLen + idx];
        }

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

//     S.Area();
//     S.Volume();

//     S.Reparam();
//     char fname[300];
//     sprintf(fname,"X.txt");
//     myIO.WriteData(fname, 3*fLen, S.x_.data_);

//     Vectors<T> den(cpu,p,nVec),vel(cpu,p,nVec);
//     xvpb(S.h_,S.normal_, (T) 0,den);
//     S.StokesMatVec(den,vel);
    
//     char fname[300];
//     sprintf(fname,"X.txt");
//     myIO.WriteData(fname, 3*fLen, S.x_.data_);
    
//     sprintf(fname,"V.txt");
//     myIO.WriteData(fname, 3*fLen, vel.data_);

    return 0;
}
