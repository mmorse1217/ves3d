#include <iostream>
#include "DeviceGPU.h"
#include "DataIO.h"
#include "Scalars.h"
#include "Vectors.h"
#include "Surface.h"
#include <math.h>

using namespace std;
typedef float T;

int main(int argc, char ** argv)
{

    int p(12), nVec(5);

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
    par.rep_up_freq_ = 12;
    par.rep_filter_freq_ = 4;

    //Device
    DeviceGPU<T> gpu;

    //IO
    DataIO<T> myIO(gpu,"",0);

    //Reading data
    bool readFromFile = true;
    OperatorsMats<T> mats(myIO, par.p_, 2*par.p_, readFromFile);
    
    //initialing device
    gpu.InitializeSHT(mats);
    
    // memory allocation
    Surface<T> S(gpu,par,mats);
    Scalars<T> X(gpu,p,nVec);
    Scalars<T> Y(gpu,p,nVec);
    Scalars<T> Z(gpu,p,nVec);
    
    // initializing vesicle positions from text file
    myIO.ReadData("precomputed/dumbbell_cart12_single.txt",dLen,S.x_.data_);
    //myIO.ReadData("precomputed/biconcave_ra95_6",dLen,S.x_.data_);
    
    S.Resize(nVec);
    for(int ii=1;ii<nVec;ii++)
        S.device_.Memcpy(S.x_.data_ + ii*dLen, S.x_.data_, dLen, MemcpyDeviceToDevice);
    S.UpdateAll();
    
    for(int ii=0;ii<nVec;ii++)
    {
        X.device_.Memcpy(X.data_ + ii*fLen, S.x_.data_          , fLen, MemcpyDeviceToDevice);
        X.device_.Memcpy(Y.data_ + ii*fLen, S.x_.data_+fLen     , fLen, MemcpyDeviceToDevice);
        X.device_.Memcpy(Z.data_ + ii*fLen, S.x_.data_+fLen+fLen, fLen, MemcpyDeviceToDevice);
    }
    
    // Checking the grad and div operator
    S.Resize(3); 
    Scalars<T> div_n(gpu,p,3);
    S.SurfDiv(S.normal_, div_n);
    axpy((T) .5, div_n, S.h_,div_n);
    
    T *dat = (T*) malloc(div_n.GetDataLength() * sizeof(T));
    div_n.device_.Memcpy(dat,div_n.data_, div_n.GetDataLength(), MemcpyDeviceToHost);
    
    T err = 0, err2;
    for(int ii=0;ii<div_n.GetDataLength(); ii++)
    {
        err2 = fabs(dat[ii]);
        err = (err>err2) ? err : err2;
    }
    cout<<"\n The error in the surface divergence (For the dumbbell .02964 expected)= "<<err<<endl;
    
    S.Resize(nVec);
    for(int ii=1;ii<nVec;ii++)
        S.device_.Memcpy(S.x_.data_ + ii*dLen, S.x_.data_, dLen, MemcpyDeviceToDevice);

    dLen *=nVec;
    S.UpdateAll();
    Vectors<T> grad(gpu,p,nVec), lap(gpu,p,nVec);

    S.SurfGrad(X,grad);
    S.SurfDiv(grad,X);
    
    S.SurfGrad(Y,grad);
    S.SurfDiv(grad,Y);
    
    S.SurfGrad(Z,grad);
    S.SurfDiv(grad,Z);
    
    for(int ii=0;ii<nVec;ii++)
    {
        lap.device_.Memcpy(lap.data_ + 3*ii*fLen          , X.data_+ii*fLen, fLen, MemcpyDeviceToDevice);
        lap.device_.Memcpy(lap.data_ + 3*ii*fLen+fLen     , Y.data_+ii*fLen, fLen, MemcpyDeviceToDevice);
        lap.device_.Memcpy(lap.data_ + 3*ii*fLen+fLen+fLen, Z.data_+ii*fLen, fLen, MemcpyDeviceToDevice);
    }

    Scalars<T> hh(gpu, p, nVec);
    DotProduct(lap,S.normal_,hh);
    axpy((T) -.5, hh, S.h_,hh);

    free(dat);
    dat = (T*) malloc(hh.GetDataLength() * sizeof(T));
    hh.device_.Memcpy(dat,hh.data_, hh.GetDataLength(), MemcpyDeviceToHost);
    err = 0;
    for(int ii=0;ii<hh.GetDataLength(); ii++)
    {
        err2 = fabs(dat[ii]);
        err = (err>err2) ? err : err2;
    }
    cout<<"\n The error in the surface grad (For the dumbbell .13703 expected)= "<<err<<endl;

    free(dat);
//     S.Area();
//     S.Volume();

//     S.Reparam();
//     char fname[300];
//     sprintf(fname,"X.txt");
//     myIO.WriteData(fname, 3*fLen, S.x_.data_);

//     Vectors<T> den(gpu,p,nVec),vel(gpu,p,nVec);
//     xvpb(S.h_,S.normal_, (T) 0,den);
//     S.StokesMatVec(den,vel);
    
//     char fname[300];
//     sprintf(fname,"X.txt");
//     myIO.WriteData(fname, 3*fLen, S.x_.data_);
    
//     sprintf(fname,"V.txt");
//     myIO.WriteData(fname, 3*fLen, vel.data_);

    return 0;
}
