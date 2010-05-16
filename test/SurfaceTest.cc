#include <iostream>
#include "Device.h"
#include "DataIO.h"
#include "Scalars.h"
#include "Vectors.h"
#include "Surface.h"
#include <math.h>

using namespace std;
typedef float T;

int main(int argc, char ** argv)
{
    {
        int p(12), nVec(5);

        int fLen(2*p*(p+1));
        int dLen(6*p*(p+1));
        
        //Setting surface parameters
        SurfaceParams<T> par; 
        
        par.p_ = p;
        par.n_surfs_ = 1;
        par.kappa_ = 1e-2;
        par.filter_freq_ = 2*p/3;
        par.rep_ts_ = 1e-1;
        par.rep_max_vel_ = 1e-1;
        par.rep_iter_max_ = 100;
        par.rep_up_freq_ = 2*p;
        par.rep_filter_freq_ = p/3;
        
        //Device
        Device<CPU> cpu;
        
        //IO
        DataIO<T,CPU> myIO(&cpu,"",0);
        
        //Reading data
        bool readFromFile = true;
        OperatorsMats<T> mats(myIO, par.p_, 2*par.p_, readFromFile);
        
        // memory allocation
        Surface<T,CPU> S(&cpu,par,mats);
        Scalars<T,CPU> X(&cpu,p,nVec);
        Scalars<T,CPU> Y(&cpu,p,nVec);
        Scalars<T,CPU> Z(&cpu,p,nVec);
        
        // initializing vesicle positions from text file
        char fname[400];
        sprintf(fname,"%s/precomputed/dumbbell_cart12_single.txt",getenv("VES3D_DIR"));
        myIO.ReadData(fname,dLen,S.x_.begin());
        //myIO.ReadData("precomputed/biconcave_ra95_12",dLen,S.x_.begin());
        
        S.Resize(nVec);
        for(int ii=1;ii<nVec;ii++)
            for(int idx=0;idx<dLen;idx++)
                *(S.x_.begin() + ii*dLen + idx) = *(S.x_.begin() + idx);
        S.UpdateAll();
    
        for(int ii=0;ii<nVec;ii++)
            for(int idx=0;idx<fLen;idx++)
            {
                *(X.begin() +  ii*fLen + idx) = *(S.x_.begin() + idx);
                *(Y.begin() +  ii*fLen + idx) = *(S.x_.begin() + idx + fLen);
                *(Z.begin() +  ii*fLen + idx) = *(S.x_.begin() + idx + fLen + fLen);
            }
    
        // Checking the grad and div operator
        S.Resize(3); 
        Scalars<T,CPU> div_n(&cpu,p,3);
        S.SurfDiv(S.normal_, div_n);
        axpy((T) .5, div_n, S.h_,div_n);
        
        T err = 0, err2;
        for(int ii=0;ii<div_n.Size(); ii++)
        {
            err2 = fabs(*(div_n.begin() + ii));
            err = (err>err2) ? err : err2;
        }
        cout<<"\n The error in the surface divergence (For the dumbbell .02964 expected)= "<<err<<endl;
    
        S.Resize(nVec);
        for(int ii=1;ii<nVec;ii++)
            for(int idx=0;idx<dLen;idx++)
                *(S.x_.begin() +  ii*dLen + idx) = *(S.x_.begin() + idx);

        dLen *=nVec;
        S.UpdateAll();
        Vectors<T,CPU> grad(&cpu,p,nVec), lap(&cpu,p,nVec);

        S.SurfGrad(X,grad);
        S.SurfDiv(grad,X);
    
        S.SurfGrad(Y,grad);
        S.SurfDiv(grad,Y);
    
        S.SurfGrad(Z,grad);
        S.SurfDiv(grad,Z);
    
        for(int ii=0;ii<nVec;ii++)
            for(int idx=0;idx<fLen;idx++)
            {
                *(lap.begin() +  3*ii*fLen + idx              ) = *(X.begin() + ii*fLen + idx);
                *(lap.begin() +  3*ii*fLen + idx + fLen       ) = *(Y.begin() + ii*fLen + idx);
                *(lap.begin() +  3*ii*fLen + idx + fLen + fLen) = *(Z.begin() + ii*fLen + idx);
            }

        Scalars<T,CPU> hh(&cpu, p, nVec);
        DotProduct(lap,S.normal_,hh);
        axpy((T) -.5, hh, S.h_,hh);
        
        err = 0;
        for(int ii=0;ii<hh.Size(); ii++)
        {
            err2 = fabs(hh.begin()[ii]);
            err = (err>err2) ? err : err2;
        }
        cout<<"\n The error in the surface grad (For the dumbbell .13703 expected)= "<<err<<endl;

        //     S.Area();
        //     S.Volume();

        //     S.Reparam(); char fname[300]; sprintf(fname,"X.txt");
        //     myIO.WriteData(fname, 3*fLen, S.x_.data_);

        //     Vectors<T,CPU> den(cpu,p,nVec),vel(cpu,p,nVec);
        //     xvpb(S.h_,S.normal_, (T) 0,den);
        //     S.StokesMatVec(den,vel);
    
        //     char fname[300];
        //     sprintf(fname,"X.txt");
        //     myIO.WriteData(fname, 3*fLen, S.x_.data_);
    
        //     sprintf(fname,"V.txt");
        //     myIO.WriteData(fname, 3*fLen, vel.data_);
    }

    PROFILEREPORT(SortTime);
    return 0;
}
