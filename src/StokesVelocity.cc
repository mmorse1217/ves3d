#include <omp.h>
#include <iostream>
#include "Surface.h"
#include <profile.hpp>

#define __ENABLE_PVFMM_PROFILER__
#define __USE_NEW_SINGULAR_INTEG__
//#define __CHECK_SINGULAR_INTEG__




#include <legendre_rule.hpp>
#include <matrix.hpp>
#define SHMAXDEG 256

template <class Real>
class SphericalHarmonics{

  public:

    static void SHC2Grid(const pvfmm::Vector<Real>& S, long p0, long p1, pvfmm::Vector<Real>& X, pvfmm::Vector<Real>* X_theta=NULL, pvfmm::Vector<Real>* X_phi=NULL){
      pvfmm::Matrix<Real>& Mf =SphericalHarmonics<Real>::MatFourier(p0,p1);
      pvfmm::Matrix<Real>& Mdf=SphericalHarmonics<Real>::MatFourierGrad(p0,p1);
      std::vector<pvfmm::Matrix<Real> >& Ml =SphericalHarmonics<Real>::MatLegendre(p0,p1);
      std::vector<pvfmm::Matrix<Real> >& Mdl=SphericalHarmonics<Real>::MatLegendreGrad(p0,p1);
      assert(p0==Ml.size()-1);
      assert(p0==Mf.Dim(0)/2);
      assert(p1==Mf.Dim(1)/2);

      long N=S.Dim()/(p0*(p0+2));
      assert(N*p0*(p0+2)==S.Dim());

      if(X.Dim()!=N*2*p1*(p1+1)) X.ReInit(N*2*p1*(p1+1));
      if(X_phi   && X_phi  ->Dim()!=N*2*p1*(p1+1)) X_phi  ->ReInit(N*2*p1*(p1+1));
      if(X_theta && X_theta->Dim()!=N*2*p1*(p1+1)) X_theta->ReInit(N*2*p1*(p1+1));

      static pvfmm::Vector<Real> B0, B1;
      B0.ReInit(N*  p0*(p0+2));
      B1.ReInit(N*2*p0*(p1+1));

      #pragma omp parallel
      { // B0 <-- Rearrange(S)
        long tid=omp_get_thread_num();
        long omp_p=omp_get_num_threads();

        long a=(tid+0)*N/omp_p;
        long b=(tid+1)*N/omp_p;
        for(long i=a;i<b;i++){
          long offset=0;
          for(long j=0;j<2*p0;j++){
            long len=p0+1-(j+1)/2;
            Real* B_=&B0[i*len+N*offset];
            Real* S_=&S[i*p0*(p0+2)+offset];
            for(long k=0;k<len;k++) B_[k]=S_[k];
            offset+=len;
          }
        }
      }

      #pragma omp parallel
      { // Evaluate Legendre polynomial
        long tid=omp_get_thread_num();
        long omp_p=omp_get_num_threads();

        long offset0=0;
        long offset1=0;
        for(long i=0;i<p0+1;i++){
          long N0=2*N;
          if(i==0 || i==p0) N0=N;
          pvfmm::Matrix<Real> Min (N0, p0+1-i,&B0[0]+offset0,false);
          pvfmm::Matrix<Real> Mout(N0, p1+1  ,&B1[0]+offset1,false);
          { // Mout = Min * Ml[i]  // split between threads
            long a=(tid+0)*N0/omp_p;
            long b=(tid+1)*N0/omp_p;
            if(a<b){
              pvfmm::Matrix<Real> Min_ (b-a, Min .Dim(1), &Min [a][0],false);
              pvfmm::Matrix<Real> Mout_(b-a, Mout.Dim(1), &Mout[a][0],false);
              pvfmm::Matrix<Real>::GEMM(Mout_,Min_,Ml[i]);
            }
          }
          offset0+=Min .Dim(0)*Min .Dim(1);
          offset1+=Mout.Dim(0)*Mout.Dim(1);
        }
      }

      #pragma omp parallel
      { // Transpose and evaluate Fourier
        long tid=omp_get_thread_num();
        long omp_p=omp_get_num_threads();

        long a=(tid+0)*N*(p1+1)/omp_p;
        long b=(tid+1)*N*(p1+1)/omp_p;

        const long block_size=16;
        pvfmm::Matrix<Real> B2(block_size,2*p0);
        for(long i0=a;i0<b;i0+=block_size){
          long i1=std::min(b,i0+block_size);
          for(long i=i0;i<i1;i++){
            for(long j=0;j<2*p0;j++){
              B2[i-i0][j]=B1[j*N*(p1+1)+i];
            }
          }

          pvfmm::Matrix<Real> Min (i1-i0,2*p0,&B2[0][0]  , false);
          pvfmm::Matrix<Real> Mout(i1-i0,2*p1,&X[i0*2*p1], false);
          pvfmm::Matrix<Real>::GEMM(Mout, Min, Mf);

          if(X_theta){ // Evaluate Fourier gradient
            pvfmm::Matrix<Real> Mout(i1-i0,2*p1,&(*X_theta)[i0*2*p1], false);
            pvfmm::Matrix<Real>::GEMM(Mout, Min, Mdf);
          }
        }
      }

      if(X_phi){
        #pragma omp parallel
        { // Evaluate Legendre gradient
          long tid=omp_get_thread_num();
          long omp_p=omp_get_num_threads();

          long offset0=0;
          long offset1=0;
          for(long i=0;i<p0+1;i++){
            long N0=2*N;
            if(i==0 || i==p0) N0=N;
            pvfmm::Matrix<Real> Min (N0, p0+1-i,&B0[0]+offset0,false);
            pvfmm::Matrix<Real> Mout(N0, p1+1  ,&B1[0]+offset1,false);
            { // Mout = Min * Mdl[i]  // split between threads
              long a=(tid+0)*N0/omp_p;
              long b=(tid+1)*N0/omp_p;
              if(a<b){
                pvfmm::Matrix<Real> Min_ (b-a, Min .Dim(1), &Min [a][0],false);
                pvfmm::Matrix<Real> Mout_(b-a, Mout.Dim(1), &Mout[a][0],false);
                pvfmm::Matrix<Real>::GEMM(Mout_,Min_,Mdl[i]);
              }
            }
            offset0+=Min .Dim(0)*Min .Dim(1);
            offset1+=Mout.Dim(0)*Mout.Dim(1);
          }
        }

        #pragma omp parallel
        { // Transpose and evaluate Fourier
          long tid=omp_get_thread_num();
          long omp_p=omp_get_num_threads();

          long a=(tid+0)*N*(p1+1)/omp_p;
          long b=(tid+1)*N*(p1+1)/omp_p;

          const long block_size=16;
          pvfmm::Matrix<Real> B2(block_size,2*p0);
          for(long i0=a;i0<b;i0+=block_size){
            long i1=std::min(b,i0+block_size);
            for(long i=i0;i<i1;i++){
              for(long j=0;j<2*p0;j++){
                B2[i-i0][j]=B1[j*N*(p1+1)+i];
              }
            }

            pvfmm::Matrix<Real> Min (i1-i0,2*p0,&B2[0][0]         , false);
            pvfmm::Matrix<Real> Mout(i1-i0,2*p1,&(*X_phi)[i0*2*p1], false);
            pvfmm::Matrix<Real>::GEMM(Mout, Min, Mf);
          }
        }
      }
    }

    static void Grid2SHC(const pvfmm::Vector<Real>& X, long p0, long p1, pvfmm::Vector<Real>& S){
      pvfmm::Matrix<Real> Mf =SphericalHarmonics<Real>::MatFourier(p1,p0).pinv();
      std::vector<pvfmm::Matrix<Real> > Ml =SphericalHarmonics<Real>::MatLegendre(p1,p0);
      for(long i=0;i<Ml.size();i++) Ml[i]=Ml[i].pinv();
      assert(p1==Ml.size()-1);
      assert(p0==Mf.Dim(0)/2);
      assert(p1==Mf.Dim(1)/2);

      long N=X.Dim()/(2*p0*(p0+1));
      assert(N*2*p0*(p0+1)==X.Dim());
      if(S.Dim()!=N*(p1*(p1+2))) S.ReInit(N*(p1*(p1+2)));

      static pvfmm::Vector<Real> B0, B1;
      B0.ReInit(N*  p1*(p1+2));
      B1.ReInit(N*2*p1*(p0+1));

      #pragma omp parallel
      { // Evaluate Fourier and transpose
        long tid=omp_get_thread_num();
        long omp_p=omp_get_num_threads();

        long a=(tid+0)*N*(p0+1)/omp_p;
        long b=(tid+1)*N*(p0+1)/omp_p;

        const long block_size=16;
        pvfmm::Matrix<Real> B2(block_size,2*p1);
        for(long i0=a;i0<b;i0+=block_size){
          long i1=std::min(b,i0+block_size);
          pvfmm::Matrix<Real> Min (i1-i0,2*p0,&X[i0*2*p0], false);
          pvfmm::Matrix<Real> Mout(i1-i0,2*p1,&B2[0][0]  , false);
          pvfmm::Matrix<Real>::GEMM(Mout, Min, Mf);

          for(long i=i0;i<i1;i++){
            for(long j=0;j<2*p1;j++){
              B1[j*N*(p0+1)+i]=B2[i-i0][j];
            }
          }
        }
      }

      #pragma omp parallel
      { // Evaluate Legendre polynomial
        long tid=omp_get_thread_num();
        long omp_p=omp_get_num_threads();

        long offset0=0;
        long offset1=0;
        for(long i=0;i<p1+1;i++){
          long N0=2*N;
          if(i==0 || i==p1) N0=N;
          pvfmm::Matrix<Real> Min (N0, p0+1  ,&B1[0]+offset0,false);
          pvfmm::Matrix<Real> Mout(N0, p1+1-i,&B0[0]+offset1,false);
          { // Mout = Min * Ml[i]  // split between threads
            long a=(tid+0)*N0/omp_p;
            long b=(tid+1)*N0/omp_p;
            if(a<b){
              pvfmm::Matrix<Real> Min_ (b-a, Min .Dim(1), &Min [a][0],false);
              pvfmm::Matrix<Real> Mout_(b-a, Mout.Dim(1), &Mout[a][0],false);
              pvfmm::Matrix<Real>::GEMM(Mout_,Min_,Ml[i]);
            }
          }
          offset0+=Min .Dim(0)*Min .Dim(1);
          offset1+=Mout.Dim(0)*Mout.Dim(1);
        }
      }

      #pragma omp parallel
      { // S <-- Rearrange(B0)
        long tid=omp_get_thread_num();
        long omp_p=omp_get_num_threads();

        long a=(tid+0)*N/omp_p;
        long b=(tid+1)*N/omp_p;
        for(long i=a;i<b;i++){
          long offset=0;
          for(long j=0;j<2*p1;j++){
            long len=p1+1-(j+1)/2;
            Real* B_=&B0[i*len+N*offset];
            Real* S_=&S[i*p1*(p1+2)+offset];
            for(long k=0;k<len;k++) S_[k]=B_[k];
            offset+=len;
          }
        }
      }
    }

    static void SHC2GridTranspose(const pvfmm::Vector<Real>& X, long p0, long p1, pvfmm::Vector<Real>& S){
      pvfmm::Matrix<Real> Mf =SphericalHarmonics<Real>::MatFourier(p1,p0).Transpose();
      std::vector<pvfmm::Matrix<Real> > Ml =SphericalHarmonics<Real>::MatLegendre(p1,p0);
      for(long i=0;i<Ml.size();i++) Ml[i]=Ml[i].Transpose();
      assert(p1==Ml.size()-1);
      assert(p0==Mf.Dim(0)/2);
      assert(p1==Mf.Dim(1)/2);

      long N=X.Dim()/(2*p0*(p0+1));
      assert(N*2*p0*(p0+1)==X.Dim());
      if(S.Dim()!=N*(p1*(p1+2))) S.ReInit(N*(p1*(p1+2)));

      static pvfmm::Vector<Real> B0, B1;
      B0.ReInit(N*  p1*(p1+2));
      B1.ReInit(N*2*p1*(p0+1));

      #pragma omp parallel
      { // Evaluate Fourier and transpose
        long tid=omp_get_thread_num();
        long omp_p=omp_get_num_threads();

        long a=(tid+0)*N*(p0+1)/omp_p;
        long b=(tid+1)*N*(p0+1)/omp_p;

        const long block_size=16;
        pvfmm::Matrix<Real> B2(block_size,2*p1);
        for(long i0=a;i0<b;i0+=block_size){
          long i1=std::min(b,i0+block_size);
          pvfmm::Matrix<Real> Min (i1-i0,2*p0,&X[i0*2*p0], false);
          pvfmm::Matrix<Real> Mout(i1-i0,2*p1,&B2[0][0]  , false);
          pvfmm::Matrix<Real>::GEMM(Mout, Min, Mf);

          for(long i=i0;i<i1;i++){
            for(long j=0;j<2*p1;j++){
              B1[j*N*(p0+1)+i]=B2[i-i0][j];
            }
          }
        }
      }

      #pragma omp parallel
      { // Evaluate Legendre polynomial
        long tid=omp_get_thread_num();
        long omp_p=omp_get_num_threads();

        long offset0=0;
        long offset1=0;
        for(long i=0;i<p1+1;i++){
          long N0=2*N;
          if(i==0 || i==p1) N0=N;
          pvfmm::Matrix<Real> Min (N0, p0+1  ,&B1[0]+offset0,false);
          pvfmm::Matrix<Real> Mout(N0, p1+1-i,&B0[0]+offset1,false);
          { // Mout = Min * Ml[i]  // split between threads
            long a=(tid+0)*N0/omp_p;
            long b=(tid+1)*N0/omp_p;
            if(a<b){
              pvfmm::Matrix<Real> Min_ (b-a, Min .Dim(1), &Min [a][0],false);
              pvfmm::Matrix<Real> Mout_(b-a, Mout.Dim(1), &Mout[a][0],false);
              pvfmm::Matrix<Real>::GEMM(Mout_,Min_,Ml[i]);
            }
          }
          offset0+=Min .Dim(0)*Min .Dim(1);
          offset1+=Mout.Dim(0)*Mout.Dim(1);
        }
      }

      #pragma omp parallel
      { // S <-- Rearrange(B0)
        long tid=omp_get_thread_num();
        long omp_p=omp_get_num_threads();

        long a=(tid+0)*N/omp_p;
        long b=(tid+1)*N/omp_p;
        for(long i=a;i<b;i++){
          long offset=0;
          for(long j=0;j<2*p1;j++){
            long len=p1+1-(j+1)/2;
            Real* B_=&B0[i*len+N*offset];
            Real* S_=&S[i*p1*(p1+2)+offset];
            for(long k=0;k<len;k++) S_[k]=B_[k];
            offset+=len;
          }
        }
      }
    }

    static void SHC2Pole(const pvfmm::Vector<Real>& S, long p0, pvfmm::Vector<Real>& P){
      pvfmm::Vector<Real> QP[2];
      { // Set QP
        Real x[2]={-1,1};
        std::vector<Real> alp((p0+1)*(p0+2)/2);
        const Real SQRT2PI=sqrt(2*M_PI);
        for(long i=0;i<2;i++){
          LegPoly(&alp[0], &x[i], 1, p0);
          QP[i].ReInit(p0+1,&alp[0]);
          for(long j=0;j<p0+1;j++) QP[i][j]*=SQRT2PI;
        }
      }

      long N=S.Dim()/(p0*(p0+2));
      assert(N*p0*(p0+2)==S.Dim());
      if(P.Dim()!=N*2) P.ReInit(N*2);

      #pragma omp parallel
      { // Compute pole
        long tid=omp_get_thread_num();
        long omp_p=omp_get_num_threads();

        long a=(tid+0)*N/omp_p;
        long b=(tid+1)*N/omp_p;

        for(long i=a;i<b;i++){
          Real P_[2]={0,0};
          for(long j=0;j<p0+1;j++){
            P_[0]+=S[i*p0*(p0+2)+j]*QP[0][j];
            P_[1]+=S[i*p0*(p0+2)+j]*QP[1][j];
          }
          P[2*i+0]=P_[0];
          P[2*i+1]=P_[1];
        }
      }
    }

    static void RotateAll(const pvfmm::Vector<Real>& S, long p0, long dof, pvfmm::Vector<Real>& S_){
      std::vector<pvfmm::Matrix<Real> >& Mr=MatRotate(p0);
      long Ncoef=p0*(p0+2);

      long N=S.Dim()/Ncoef/dof;
      assert(N*Ncoef*dof==S.Dim());
      if(S_.Dim()!=N*dof*Ncoef*p0*(p0+1)) S_.ReInit(N*dof*Ncoef*p0*(p0+1));
      pvfmm::Matrix<Real> S0(N*dof          ,Ncoef, &S [0], false);
      pvfmm::Matrix<Real> S1(N*dof*p0*(p0+1),Ncoef, &S_[0], false);

      #pragma omp parallel
      { // Construct all p0*(p0+1) rotations
        long tid=omp_get_thread_num();
        long omp_p=omp_get_num_threads();
        pvfmm::Matrix<Real> B0(dof*p0,Ncoef); // memory buffer

        long a=(tid+0)*N/omp_p;
        long b=(tid+1)*N/omp_p;
        for(long i=a;i<b;i++){
          for(long d=0;d<dof;d++){
            for(long j=0;j<p0;j++){
              long offset=0;
              for(long k=0;k<p0+1;k++){
                Real r[2]={cos(k*j*M_PI/p0),-sin(k*j*M_PI/p0)}; // exp(i*k*theta)
                long len=p0+1-k;
                if(k!=0 && k!=p0){
                  for(long l=0;l<len;l++){
                    Real x[2];
                    x[0]=S0[i*dof+d][offset+len*0+l];
                    x[1]=S0[i*dof+d][offset+len*1+l];
                    B0[j*dof+d][offset+len*0+l]=x[0]*r[0]-x[1]*r[1];
                    B0[j*dof+d][offset+len*1+l]=x[0]*r[1]+x[1]*r[0];
                  }
                  offset+=2*len;
                }else{
                  for(long l=0;l<len;l++){
                    B0[j*dof+d][offset+l]=S0[i*dof+d][offset+l];
                  }
                  offset+=len;
                }
              }
              assert(offset==Ncoef);
            }
          }
          for(long t=0;t<p0+1;t++){
            pvfmm::Matrix<Real> Mout(dof*p0,Ncoef,&S1[(i*(p0+1)+t)*dof*p0][0],false);
            pvfmm::Matrix<Real>::GEMM(Mout,B0,Mr[t]);
          }
        }
      }
    }

    static void RotateTranspose(const pvfmm::Vector<Real>& S_, long p0, long dof, pvfmm::Vector<Real>& S){
      std::vector<pvfmm::Matrix<Real> > Mr=MatRotate(p0);
      for(long i=0;i<p0+1;i++) Mr[i]=Mr[i].Transpose();
      long Ncoef=p0*(p0+2);

      long N=S_.Dim()/Ncoef/dof/(p0*(p0+1));
      assert(N*Ncoef*dof*(p0*(p0+1))==S_.Dim());
      if(S.Dim()!=N*dof*Ncoef*p0*(p0+1)) S.ReInit(N*dof*Ncoef*p0*(p0+1));
      pvfmm::Matrix<Real> S0(N*dof*p0*(p0+1),Ncoef, &S [0], false);
      pvfmm::Matrix<Real> S1(N*dof*p0*(p0+1),Ncoef, &S_[0], false);

      #pragma omp parallel
      { // Transpose all p0*(p0+1) rotations
        long tid=omp_get_thread_num();
        long omp_p=omp_get_num_threads();
        pvfmm::Matrix<Real> B0(dof*p0,Ncoef); // memory buffer

        long a=(tid+0)*N/omp_p;
        long b=(tid+1)*N/omp_p;
        for(long i=a;i<b;i++){
          for(long t=0;t<p0+1;t++){
            long idx0=(i*(p0+1)+t)*p0*dof;
            pvfmm::Matrix<Real> Min(p0*dof,Ncoef, &S1[idx0][0],false);
            pvfmm::Matrix<Real>::GEMM(B0,Min,Mr[t]);

            for(long j=0;j<p0;j++){
              for(long d=0;d<dof;d++){
                long idx1=idx0+j*dof+d;
                long offset=0;
                for(long k=0;k<p0+1;k++){
                  Real r[2]={cos(k*j*M_PI/p0),sin(k*j*M_PI/p0)}; // exp(i*k*theta)
                  long len=p0+1-k;
                  if(k!=0 && k!=p0){
                    for(long l=0;l<len;l++){
                      Real x[2];
                      x[0]=B0[j*dof+d][offset+len*0+l];
                      x[1]=B0[j*dof+d][offset+len*1+l];
                      S0[idx1][offset+len*0+l]=x[0]*r[0]-x[1]*r[1];
                      S0[idx1][offset+len*1+l]=x[0]*r[1]+x[1]*r[0];
                    }
                    offset+=2*len;
                  }else{
                    for(long l=0;l<len;l++){
                      S0[idx1][offset+l]=B0[j*dof+d][offset+l];
                    }
                    offset+=len;
                  }
                }
                assert(offset==Ncoef);
              }
            }
          }
        }
      }
    }

    static pvfmm::Vector<Real>& LegendreNodes(long p1){
      assert(p1<SHMAXDEG);
      matrix.Qx_.resize(SHMAXDEG);
      pvfmm::Vector<Real>& Qx=matrix.Qx_[p1];
      if(!Qx.Dim()){
        std::vector<Real> qx1(p1+1);
        std::vector<Real> qw1(p1+1);
        cgqf(p1+1, 1, 0.0, 0.0, -1.0, 1.0, &qx1[0], &qw1[0]);
        Qx=qx1;
      }
      return Qx;
    }

    static pvfmm::Vector<Real>& LegendreWeights(long p1){
      assert(p1<SHMAXDEG);
      matrix.Qw_.resize(SHMAXDEG);
      pvfmm::Vector<Real>& Qw=matrix.Qw_[p1];
      if(!Qw.Dim()){
        std::vector<Real> qx1(p1+1);
        std::vector<Real> qw1(p1+1);
        cgqf(p1+1, 1, 0.0, 0.0, -1.0, 1.0, &qx1[0], &qw1[0]);
        for(long i=0;i<qw1.size();i++) qw1[i]*=M_PI/p1/sqrt(1-qx1[i]*qx1[i]);
        Qw=qw1;
      }
      return Qw;
    }

    static pvfmm::Vector<Real>& SingularWeights(long p1){
      assert(p1<SHMAXDEG);
      matrix.Sw_.resize(SHMAXDEG);
      pvfmm::Vector<Real>& Sw=matrix.Sw_[p1];
      if(!Sw.Dim()){
        std::vector<Real> qx1(p1+1);
        std::vector<Real> qw1(p1+1);
        cgqf(p1+1, 1, 0.0, 0.0, -1.0, 1.0, &qx1[0], &qw1[0]);

        std::vector<Real> Yf(p1+1,0);
        { // Set Yf
          Real x0=1.0;
          std::vector<Real> alp0((p1+1)*(p1+2)/2);
          LegPoly(&alp0[0], &x0, 1, p1);

          std::vector<Real> alp((p1+1) * (p1+1)*(p1+2)/2);
          LegPoly(&alp[0], &qx1[0], p1+1, p1);

          for(long j=0;j<p1+1;j++){
            for(long i=0;i<p1+1;i++){
              Yf[i]+=4*M_PI/(2*j+1) * alp0[j] * alp[j*(p1+1)+i];
            }
          }
        }

        Sw.ReInit(p1+1);
        for(long i=0;i<p1+1;i++){
          Sw[i]=(qw1[i]*M_PI/p1)*Yf[i]/cos(acos(qx1[i])/2);
        }
      }
      return Sw;
    }

    static pvfmm::Matrix<Real>& MatFourier(long p0, long p1){
      assert(p0<SHMAXDEG && p1<SHMAXDEG);
      matrix.Mf_ .resize(SHMAXDEG*SHMAXDEG);
      pvfmm::Matrix<Real>& Mf =matrix.Mf_ [p0*SHMAXDEG+p1];
      if(!Mf.Dim(0)){
        const Real SQRT2PI=sqrt(2*M_PI);
        { // Set Mf
          pvfmm::Matrix<Real> M(2*p0,2*p1);
          for(long j=0;j<2*p1;j++){
            M[0][j]=SQRT2PI*1.0;
            for(long k=1;k<p0;k++){
              M[2*k-1][j]=SQRT2PI*cos(j*k*M_PI/p1);
              M[2*k-0][j]=SQRT2PI*sin(j*k*M_PI/p1);
            }
            M[2*p0-1][j]=SQRT2PI*cos(j*p0*M_PI/p1);
          }
          Mf=M;
        }
      }
      return Mf;
    }

    static pvfmm::Matrix<Real>& MatFourierGrad(long p0, long p1){
      assert(p0<SHMAXDEG && p1<SHMAXDEG);
      matrix.Mdf_.resize(SHMAXDEG*SHMAXDEG);
      pvfmm::Matrix<Real>& Mdf=matrix.Mdf_[p0*SHMAXDEG+p1];
      if(!Mdf.Dim(0)){
        const Real SQRT2PI=sqrt(2*M_PI);
        { // Set Mdf_
          pvfmm::Matrix<Real> M(2*p0,2*p1);
          for(long j=0;j<2*p1;j++){
            M[0][j]=SQRT2PI*0.0;
            for(long k=1;k<p0;k++){
              M[2*k-1][j]=-SQRT2PI*k*sin(j*k*M_PI/p1);
              M[2*k-0][j]= SQRT2PI*k*cos(j*k*M_PI/p1);
            }
            M[2*p0-1][j]=-SQRT2PI*p0*sin(j*p0*M_PI/p1);
          }
          Mdf=M;
        }
      }
      return Mdf;
    }

    static std::vector<pvfmm::Matrix<Real> >& MatLegendre(long p0, long p1){
      assert(p0<SHMAXDEG && p1<SHMAXDEG);
      matrix.Ml_ .resize(SHMAXDEG*SHMAXDEG);
      std::vector<pvfmm::Matrix<Real> >& Ml =matrix.Ml_ [p0*SHMAXDEG+p1];
      if(!Ml.size()){
        std::vector<Real> qx1(p1+1);
        std::vector<Real> qw1(p1+1);
        cgqf(p1+1, 1, 0.0, 0.0, -1.0, 1.0, &qx1[0], &qw1[0]);

        { // Set Ml
          std::vector<Real> alp(qx1.size()*(p0+1)*(p0+2)/2);
          LegPoly(&alp[0], &qx1[0], qx1.size(), p0);

          Ml.resize(p0+1);
          Real* ptr=&alp[0];
          for(long i=0;i<=p0;i++){
            Ml[i].ReInit(p0+1-i, qx1.size(), ptr);
            ptr+=Ml[i].Dim(0)*Ml[i].Dim(1);
          }
        }
      }
      return Ml;
    }

    static std::vector<pvfmm::Matrix<Real> >& MatLegendreGrad(long p0, long p1){
      assert(p0<SHMAXDEG && p1<SHMAXDEG);
      matrix.Mdl_.resize(SHMAXDEG*SHMAXDEG);
      std::vector<pvfmm::Matrix<Real> >& Mdl=matrix.Mdl_[p0*SHMAXDEG+p1];
      if(!Mdl.size()){
        std::vector<Real> qx1(p1+1);
        std::vector<Real> qw1(p1+1);
        cgqf(p1+1, 1, 0.0, 0.0, -1.0, 1.0, &qx1[0], &qw1[0]);

        { // Set Mdl
          std::vector<Real> alp(qx1.size()*(p0+1)*(p0+2)/2);
          LegPolyDeriv(&alp[0], &qx1[0], qx1.size(), p0);

          Mdl.resize(p0+1);
          Real* ptr=&alp[0];
          for(long i=0;i<=p0;i++){
            Mdl[i].ReInit(p0+1-i, qx1.size(), ptr);
            ptr+=Mdl[i].Dim(0)*Mdl[i].Dim(1);
          }
        }
      }
      return Mdl;
    }

    static std::vector<pvfmm::Matrix<Real> >& MatRotate(long p0){
      assert(p0<SHMAXDEG);
      matrix.Mr_.resize(SHMAXDEG);
      std::vector<pvfmm::Matrix<Real> >& Mr=matrix.Mr_[p0];
      if(!Mr.size()){
        const Real SQRT2PI=sqrt(2*M_PI);
        long Ncoef=p0*(p0+2);
        long Ngrid=2*p0*(p0+1);
        long Naleg=(p0+1)*(p0+2)/2;

        pvfmm::Matrix<Real> Mcoord0(3,Ngrid);
        pvfmm::Vector<Real>& x=LegendreNodes(p0);
        for(long i=0;i<p0+1;i++){ // Set Mcoord0
          for(long j=0;j<2*p0;j++){
            Mcoord0[0][i*2*p0+j]=x[i];
            Mcoord0[1][i*2*p0+j]=sqrt(1-x[i]*x[i])*sin(M_PI*j/p0);
            Mcoord0[2][i*2*p0+j]=sqrt(1-x[i]*x[i])*cos(M_PI*j/p0);
          }
        }

        for(long l=0;l<p0+1;l++){ // For each rotation angle
          pvfmm::Matrix<Real> Mcoord1;
          { // Rotate coordinates
            pvfmm::Matrix<Real> M(COORD_DIM, COORD_DIM);
            Real cos_=-x[l];
            Real sin_=-sqrt(1.0-x[l]*x[l]);
            M[0][0]= cos_; M[0][1]=0; M[0][2]=-sin_;
            M[1][0]=    0; M[1][1]=1; M[1][2]=    0;
            M[2][0]= sin_; M[2][1]=0; M[2][2]= cos_;
            Mcoord1=M*Mcoord0;
          }

          pvfmm::Matrix<Real> Mleg(Naleg, Ngrid);
          { // Set Mleg
            LegPoly(&Mleg[0][0], &Mcoord1[0][0], Ngrid, p0);
          }

          pvfmm::Vector<Real> theta(Ngrid);
          for(long i=0;i<theta.Dim();i++){ // Set theta
            theta[i]=atan2(Mcoord1[1][i],Mcoord1[2][i]);
          }

          pvfmm::Matrix<Real> Mcoef2grid(Ncoef, Ngrid);
          { // Build Mcoef2grid
            long offset0=0;
            long offset1=0;
            for(long i=0;i<p0+1;i++){
              long len=p0+1-i;
              { // P * cos
                for(long j=0;j<len;j++){
                  for(long k=0;k<Ngrid;k++){
                    Mcoef2grid[offset1+j][k]=SQRT2PI*Mleg[offset0+j][k]*cos(i*theta[k]);
                  }
                }
                offset1+=len;
              }
              if(i!=0 && i!=p0){ // P * sin
                for(long j=0;j<len;j++){
                  for(long k=0;k<Ngrid;k++){
                    Mcoef2grid[offset1+j][k]=SQRT2PI*Mleg[offset0+j][k]*sin(i*theta[k]);
                  }
                }
                offset1+=len;
              }
              offset0+=len;
            }
            assert(offset0==Naleg);
            assert(offset1==Ncoef);
          }

          pvfmm::Vector<Real> Vcoef2coef(Ncoef*Ncoef);
          pvfmm::Vector<Real> Vcoef2grid(Ncoef*Ngrid, &Mcoef2grid[0][0],false);
          Grid2SHC(Vcoef2grid, p0, p0, Vcoef2coef);

          pvfmm::Matrix<Real> Mcoef2coef(Ncoef, Ncoef, &Vcoef2coef[0],false);
          Mr.push_back(Mcoef2coef);
        }
      }
      return Mr;
    }

    static void StokesSingularInteg(const pvfmm::Vector<Real>& S, long p0, long p1, pvfmm::Vector<Real>* SLMatrix=NULL, pvfmm::Vector<Real>* DLMatrix=NULL){
      long Ngrid=2*p0*(p0+1);
      long Ncoef=  p0*(p0+2);
      long Nves=S.Dim()/(Ngrid*COORD_DIM);
      if(SLMatrix) SLMatrix->ReInit(Nves*(Ncoef*COORD_DIM)*(Ncoef*COORD_DIM));
      if(DLMatrix) DLMatrix->ReInit(Nves*(Ncoef*COORD_DIM)*(Ncoef*COORD_DIM));

      long BLOCK_SIZE=omp_get_max_threads();
      for(long a=0;a<Nves;a+=BLOCK_SIZE){
        long b=std::min(a+BLOCK_SIZE, Nves);

        pvfmm::Vector<Real> _SLMatrix, _DLMatrix, _S;
        if(SLMatrix) _SLMatrix.ReInit((b-a)*(Ncoef*COORD_DIM)*(Ncoef*COORD_DIM),&SLMatrix[0][a*(Ncoef*COORD_DIM)*(Ncoef*COORD_DIM)],false);
        if(DLMatrix) _DLMatrix.ReInit((b-a)*(Ncoef*COORD_DIM)*(Ncoef*COORD_DIM),&DLMatrix[0][a*(Ncoef*COORD_DIM)*(Ncoef*COORD_DIM)],false);
        _S                    .ReInit((b-a)*(Ngrid*COORD_DIM)                  ,&S          [a*(Ngrid*COORD_DIM)                  ],false);

        if(SLMatrix && DLMatrix) StokesSingularInteg_< true,  true>(_S, p0, p1, _SLMatrix, _DLMatrix);
        else        if(SLMatrix) StokesSingularInteg_< true, false>(_S, p0, p1, _SLMatrix, _DLMatrix);
        else        if(DLMatrix) StokesSingularInteg_<false,  true>(_S, p0, p1, _SLMatrix, _DLMatrix);
      }
    }

  private:

    /**
     * \brief Computes all the Associated Legendre Polynomials (normalized) upto the specified degree.
     * \param[in] degree The degree upto which the legendre polynomials have to be computed.
     * \param[in] X The input values for which the polynomials have to be computed.
     * \param[in] N The number of input points.
     * \param[out] poly_val The output array of size (degree+1)*(degree+2)*N/2 containing the computed polynomial values.
     * The output values are in the order:
     * P(n,m)[i] => {P(0,0)[0], P(0,0)[1], ..., P(0,0)[N-1], P(1,0)[0], ..., P(1,0)[N-1],
     * P(2,0)[0], ..., P(degree,0)[N-1], P(1,1)[0], ...,P(2,1)[0], ..., P(degree,degree)[N-1]}
     */
    static void LegPoly(Real* poly_val, const Real* X, long N, long degree){
      Real* p_val=poly_val;
      Real fact=1.0/(Real)sqrt(4*M_PI);

      std::vector<Real> u(N);
      for(int n=0;n<N;n++){
        u[n]=sqrt(1-X[n]*X[n]);
        if(X[n]*X[n]>1.0) u[n]=0;
        p_val[n]=fact;
      }

      Real* p_val_nxt=poly_val;
      for(int i=1;i<=degree;i++){
        p_val_nxt=&p_val_nxt[N*(degree-i+2)];
        Real c=(i==1?sqrt(3.0/2.0):1);
        if(i>1)c*=sqrt((Real)(2*i+1)/(2*i));
        for(int n=0;n<N;n++){
          p_val_nxt[n]=-p_val[n]*u[n]*c;
        }
        p_val=p_val_nxt;
      }

      p_val=poly_val;
      for(int m=0;m<degree;m++){
        for(int n=0;n<N;n++){
          Real pmm=0;
          Real pmmp1=p_val[n];
          Real pll;
          for(int ll=m+1;ll<=degree;ll++){
            Real a=sqrt(((Real)(2*ll-1)*(2*ll+1))/((ll-m)*(ll+m)));
            Real b=sqrt(((Real)(2*ll+1)*(ll+m-1)*(ll-m-1))/((ll-m)*(ll+m)*(2*ll-3)));
            pll=X[n]*a*pmmp1-b*pmm;
            pmm=pmmp1;
            pmmp1=pll;
            p_val[N*(ll-m)+n]=pll;
          }
        }
        p_val=&p_val[N*(degree-m+1)];
      }
    }

    static void LegPolyDeriv(Real* poly_val, const Real* X, long N, long degree){
      std::vector<Real> leg_poly((degree+1)*(degree+2)*N/2);
      LegPoly(&leg_poly[0], X, N, degree);

      for(long m=0;m<=degree;m++){
        for(long n=0;n<=degree;n++) if(m<=n){
          const Real* Pn =&leg_poly[0];
          const Real* Pn_=&leg_poly[0];
          if((m+0)<=(n+0)) Pn =&leg_poly[N*(((degree*2-abs(m+0)+1)*abs(m+0))/2+(n+0))];
          if((m+1)<=(n+0)) Pn_=&leg_poly[N*(((degree*2-abs(m+1)+1)*abs(m+1))/2+(n+0))];
          Real*            Hn =&poly_val[N*(((degree*2-abs(m+0)+1)*abs(m+0))/2+(n+0))];

          Real c1=(abs(m+0)<=(n+0)?1.0:0)*m;
          Real c2=(abs(m+1)<=(n+0)?1.0:0)*sqrt(n+m+1)*sqrt(n>m?n-m:1);
          for(int i=0;i<N;i++){
            Hn[i]=-(c1*X[i]*Pn[i]+c2*sqrt(1-X[i]*X[i])*Pn_[i])/sqrt(1-X[i]*X[i]);
          }
        }
      }
    }

    template <bool SLayer, bool DLayer>
    static void StokesSingularInteg_(const pvfmm::Vector<Real>& X0, long p0, long p1, pvfmm::Vector<Real>& SL, pvfmm::Vector<Real>& DL){

      pvfmm::Profile::Tic("Rotate");
      static pvfmm::Vector<Real> S0, S;
      SphericalHarmonics<Real>::Grid2SHC(X0, p0, p0, S0);
      SphericalHarmonics<Real>::RotateAll(S0, p0, COORD_DIM, S);
      pvfmm::Profile::Toc();


      pvfmm::Profile::Tic("Upsample");
      static pvfmm::Vector<Real> X, X_phi, X_theta, trg;
      SphericalHarmonics<Real>::SHC2Grid(S, p0, p1, X, &X_theta, &X_phi);
      SphericalHarmonics<Real>::SHC2Pole(S, p0, trg);
      pvfmm::Profile::Toc();


      pvfmm::Profile::Tic("Stokes");
      static pvfmm::Vector<Real> SL0, DL0;
      { // Stokes kernel
        long M0=2*p0*(p0+1);
        long M1=2*p1*(p1+1);
        long N=trg.Dim()/(2*COORD_DIM);
        assert(X.Dim()==M1*COORD_DIM*N);
        if(SLayer && SL0.Dim()!=N*2*COORD_DIM*COORD_DIM*M1) SL0.ReInit(2*N*COORD_DIM*COORD_DIM*M1);
        if(DLayer && DL0.Dim()!=N*2*COORD_DIM*COORD_DIM*M1) DL0.ReInit(2*N*COORD_DIM*COORD_DIM*M1);
        pvfmm::Vector<Real>& qw=SphericalHarmonics<Real>::SingularWeights(p1);

        const Real scal_const_dl = 3.0/(4.0*M_PI);
        const Real scal_const_sl = 1.0/(8.0*M_PI);
        Real eps=-1;
        if(eps<0){
          eps=1;
          while(eps*(Real)0.5+(Real)1.0>1.0) eps*=0.5;
        }

        #pragma omp parallel
        {
          long tid=omp_get_thread_num();
          long omp_p=omp_get_num_threads();

          long a=(tid+0)*N/omp_p;
          long b=(tid+1)*N/omp_p;
          for(long i=a;i<b;i++){
            for(long t=0;t<2;t++){
              Real tx, ty, tz;
              { // Read target coordinates
                tx=trg[i*2*COORD_DIM+0*2+t];
                ty=trg[i*2*COORD_DIM+1*2+t];
                tz=trg[i*2*COORD_DIM+2*2+t];
              }

              for(long j0=0;j0<p1+1;j0++){
                for(long j1=0;j1<2*p1;j1++){
                  long s=2*p1*j0+j1;

                  Real dx, dy, dz;
                  { // Compute dx, dy, dz
                    dx=tx-X[(i*COORD_DIM+0)*M1+s];
                    dy=ty-X[(i*COORD_DIM+1)*M1+s];
                    dz=tz-X[(i*COORD_DIM+2)*M1+s];
                  }

                  Real nx, ny, nz;
                  { // Compute source normal
                    Real x_theta=X_theta[(i*COORD_DIM+0)*M1+s];
                    Real y_theta=X_theta[(i*COORD_DIM+1)*M1+s];
                    Real z_theta=X_theta[(i*COORD_DIM+2)*M1+s];

                    Real x_phi=X_phi[(i*COORD_DIM+0)*M1+s];
                    Real y_phi=X_phi[(i*COORD_DIM+1)*M1+s];
                    Real z_phi=X_phi[(i*COORD_DIM+2)*M1+s];

                    nx=(y_theta*z_phi-z_theta*y_phi);
                    ny=(z_theta*x_phi-x_theta*z_phi);
                    nz=(x_theta*y_phi-y_theta*x_phi);
                  }

                  Real area_elem=1.0;
                  if(SLayer){ // Compute area_elem
                    area_elem=sqrt(nx*nx+ny*ny+nz*nz);
                  }

                  Real rinv, rinv2;
                  { // Compute rinv, rinv2
                    Real r2=dx*dx+dy*dy+dz*dz;
                    rinv=1.0/sqrt(r2);
                    if(r2<=eps) rinv=0;
                    rinv2=rinv*rinv;
                  }

                  if(DLayer){
                    Real rinv5=rinv2*rinv2*rinv;
                    Real r_dot_n_rinv5=scal_const_dl*qw[j0*t+(p1-j0)*(1-t)] * (nx*dx+ny*dy+nz*dz)*rinv5;
                    DL0[((i*2+t)*COORD_DIM*COORD_DIM+0)*M1+s]=dx*dx*r_dot_n_rinv5;
                    DL0[((i*2+t)*COORD_DIM*COORD_DIM+1)*M1+s]=dx*dy*r_dot_n_rinv5;
                    DL0[((i*2+t)*COORD_DIM*COORD_DIM+2)*M1+s]=dx*dz*r_dot_n_rinv5;
                    DL0[((i*2+t)*COORD_DIM*COORD_DIM+3)*M1+s]=dy*dx*r_dot_n_rinv5;
                    DL0[((i*2+t)*COORD_DIM*COORD_DIM+4)*M1+s]=dy*dy*r_dot_n_rinv5;
                    DL0[((i*2+t)*COORD_DIM*COORD_DIM+5)*M1+s]=dy*dz*r_dot_n_rinv5;
                    DL0[((i*2+t)*COORD_DIM*COORD_DIM+6)*M1+s]=dz*dx*r_dot_n_rinv5;
                    DL0[((i*2+t)*COORD_DIM*COORD_DIM+7)*M1+s]=dz*dy*r_dot_n_rinv5;
                    DL0[((i*2+t)*COORD_DIM*COORD_DIM+8)*M1+s]=dz*dz*r_dot_n_rinv5;
                  }
                  if(SLayer){
                    Real area_rinv =scal_const_sl*qw[j0*t+(p1-j0)*(1-t)] * area_elem*rinv;
                    Real area_rinv2=area_rinv*rinv2;
                    SL0[((i*2+t)*COORD_DIM*COORD_DIM+0)*M1+s]=area_rinv+dx*dx*area_rinv2;
                    SL0[((i*2+t)*COORD_DIM*COORD_DIM+1)*M1+s]=          dx*dy*area_rinv2;
                    SL0[((i*2+t)*COORD_DIM*COORD_DIM+2)*M1+s]=          dx*dz*area_rinv2;
                    SL0[((i*2+t)*COORD_DIM*COORD_DIM+3)*M1+s]=          dy*dx*area_rinv2;
                    SL0[((i*2+t)*COORD_DIM*COORD_DIM+4)*M1+s]=area_rinv+dy*dy*area_rinv2;
                    SL0[((i*2+t)*COORD_DIM*COORD_DIM+5)*M1+s]=          dy*dz*area_rinv2;
                    SL0[((i*2+t)*COORD_DIM*COORD_DIM+6)*M1+s]=          dz*dx*area_rinv2;
                    SL0[((i*2+t)*COORD_DIM*COORD_DIM+7)*M1+s]=          dz*dy*area_rinv2;
                    SL0[((i*2+t)*COORD_DIM*COORD_DIM+8)*M1+s]=area_rinv+dz*dz*area_rinv2;
                  }
                }
              }
            }
          }
        }
      }
      pvfmm::Profile::Toc();


      pvfmm::Profile::Tic("UpsampleTranspose");
      static pvfmm::Vector<Real> SL1, DL1;
      SphericalHarmonics<Real>::SHC2GridTranspose(SL0, p1, p0, SL1);
      SphericalHarmonics<Real>::SHC2GridTranspose(DL0, p1, p0, DL1);
      pvfmm::Profile::Toc();


      pvfmm::Profile::Tic("RotateTranspose");
      static pvfmm::Vector<Real> SL2, DL2;
      SphericalHarmonics<Real>::RotateTranspose(SL1, p0, 2*COORD_DIM*COORD_DIM, SL2);
      SphericalHarmonics<Real>::RotateTranspose(DL1, p0, 2*COORD_DIM*COORD_DIM, DL2);
      pvfmm::Profile::Toc();


      pvfmm::Profile::Tic("Rearrange");
      { // Transpose
        long Ncoef=p0*(p0+2);
        long Ngrid=2*p0*(p0+1);
        { // Transpose SL2
          long N=SL2.Dim()/(COORD_DIM*COORD_DIM*Ncoef*Ngrid);
          #pragma omp parallel
          {
            long tid=omp_get_thread_num();
            long omp_p=omp_get_num_threads();
            pvfmm::Matrix<Real> B(COORD_DIM*Ncoef,Ngrid*COORD_DIM);

            long a=(tid+0)*N/omp_p;
            long b=(tid+1)*N/omp_p;
            for(long i=a;i<b;i++){
              pvfmm::Matrix<Real> M0(Ngrid*COORD_DIM, COORD_DIM*Ncoef, &SL2[i*COORD_DIM*Ngrid*COORD_DIM*Ncoef], false);
              for(long k=0;k<B.Dim(0);k++){ // Transpose
                for(long j=0;j<B.Dim(1);j++){ // TODO: needs blocking
                  B[k][j]=M0[j][k];
                }
              }
              pvfmm::Matrix<Real> M1(Ncoef*COORD_DIM, COORD_DIM*Ngrid, &SL2[i*COORD_DIM*Ncoef*COORD_DIM*Ngrid], false);
              for(long k=0;k<B.Dim(0);k++){ // Rearrange
                for(long j0=0;j0<COORD_DIM;j0++){
                  for(long j1=0;j1<p0+1;j1++){
                    for(long j2=0;j2<p0;j2++) M1[k][((j0*(p0+1)+   j1)*2+0)*p0+j2]=B[k][((j1*p0+j2)*2+0)*COORD_DIM+j0];
                    for(long j2=0;j2<p0;j2++) M1[k][((j0*(p0+1)+p0-j1)*2+1)*p0+j2]=B[k][((j1*p0+j2)*2+1)*COORD_DIM+j0];
                  }
                }
              }
            }
          }
        }
        { // Transpose DL2
          long N=DL2.Dim()/(COORD_DIM*COORD_DIM*Ncoef*Ngrid);
          #pragma omp parallel
          {
            long tid=omp_get_thread_num();
            long omp_p=omp_get_num_threads();
            pvfmm::Matrix<Real> B(COORD_DIM*Ncoef,Ngrid*COORD_DIM);

            long a=(tid+0)*N/omp_p;
            long b=(tid+1)*N/omp_p;
            for(long i=a;i<b;i++){
              pvfmm::Matrix<Real> M0(Ngrid*COORD_DIM, COORD_DIM*Ncoef, &DL2[i*COORD_DIM*Ngrid*COORD_DIM*Ncoef], false);
              for(long k=0;k<B.Dim(0);k++){ // Transpose
                for(long j=0;j<B.Dim(1);j++){ // TODO: needs blocking
                  B[k][j]=M0[j][k];
                }
              }
              pvfmm::Matrix<Real> M1(Ncoef*COORD_DIM, COORD_DIM*Ngrid, &DL2[i*COORD_DIM*Ncoef*COORD_DIM*Ngrid], false);
              for(long k=0;k<B.Dim(0);k++){ // Rearrange
                for(long j0=0;j0<COORD_DIM;j0++){
                  for(long j1=0;j1<p0+1;j1++){
                    for(long j2=0;j2<p0;j2++) M1[k][((j0*(p0+1)+   j1)*2+0)*p0+j2]=B[k][((j1*p0+j2)*2+0)*COORD_DIM+j0];
                    for(long j2=0;j2<p0;j2++) M1[k][((j0*(p0+1)+p0-j1)*2+1)*p0+j2]=B[k][((j1*p0+j2)*2+1)*COORD_DIM+j0];
                  }
                }
              }
            }
          }
        }
      }
      pvfmm::Profile::Toc();


      pvfmm::Profile::Tic("Grid2SHC");
      SphericalHarmonics<Real>::Grid2SHC(SL2, p0, p0, SL);
      SphericalHarmonics<Real>::Grid2SHC(DL2, p0, p0, DL);
      pvfmm::Profile::Toc();

    }

    static struct MatrixStorage{
      MatrixStorage(int size){
        Qx_ .resize(size);
        Qw_ .resize(size);
        Sw_ .resize(size);
        Mf_ .resize(size*size);
        Mdf_.resize(size*size);
        Ml_ .resize(size*size);
        Mdl_.resize(size*size);
        Mr_ .resize(size);
      }
      std::vector<pvfmm::Vector<Real> > Qx_;
      std::vector<pvfmm::Vector<Real> > Qw_;
      std::vector<pvfmm::Vector<Real> > Sw_;
      std::vector<pvfmm::Matrix<Real> > Mf_ ;
      std::vector<pvfmm::Matrix<Real> > Mdf_;
      std::vector<std::vector<pvfmm::Matrix<Real> > > Ml_ ;
      std::vector<std::vector<pvfmm::Matrix<Real> > > Mdl_;
      std::vector<std::vector<pvfmm::Matrix<Real> > > Mr_;
    } matrix;

};

SphericalHarmonics<double>::MatrixStorage SphericalHarmonics<double>::matrix(SHMAXDEG);
SphericalHarmonics<float >::MatrixStorage SphericalHarmonics<float >::matrix(SHMAXDEG);




template<typename Real_t>
static void LegPoly(const Real_t* x, size_t n, size_t q, Real_t* y){
  if(q>0) for(size_t i=0;i<n;i++) y[i]=1.0;
  if(q>1) for(size_t i=0;i<n;i++) y[n+i]=x[i];
  for(size_t j=2;j<q;j++){
    Real_t inv_j=1.0/j;
    for(size_t i=0;i<n;i++){
      y[j*n+i]=((2.0*j-1.0)*x[i]*y[(j-1)*n+i] - (j-1.0)*y[(j-2)*n+i])*inv_j;
    }
  }

  // Normalize
  for(size_t j=0;j<q;j++){
    for(size_t i=0;i<n;i++){
      y[j*n+i]*=sqrt(2.0*j+1.0);
    }
  }
}

template<typename Surf_t>
StokesVelocity<Surf_t>::StokesVelocity(
    const OperatorsMats<Arr_t> &mats,
    const Parameters<Real_t> &sim_par_,
    MPI_Comm c):
  near_singular(sht_up_.getShOrder(),sim_par_.periodic_length,c),
  sht_   (mats.p_   , mats.mats_p_   ),
  sht_up_(mats.p_up_, mats.mats_p_up_),
  box_size_(sim_par_.periodic_length),
  sim_par(sim_par_),
  move_pole(mats),
  S_self(NULL),
  S_far(NULL),
  comm(c)
{
  S=NULL;
  S_vel_ptr=NULL;
  force_single=NULL;
  force_double=NULL;

  fmm_flag=fmm_flag | StokesVelocity::UpdateSrcCoord;
  fmm_flag=fmm_flag | StokesVelocity::UpdateDensitySL;
  fmm_flag=fmm_flag | StokesVelocity::UpdateDensityDL;
  fmm_flag=fmm_flag | StokesVelocity::UpdateTrgCoord;
  self_flag=self_flag | StokesVelocity::UpdateSrcCoord;
  self_flag=self_flag | StokesVelocity::UpdateDensitySL;
  self_flag=self_flag | StokesVelocity::UpdateDensityDL;
  self_flag=self_flag | StokesVelocity::UpdateSurfaceVel;

  sh_order     =sht_   .getShOrder();
  sh_order_far =sht_up_.getShOrder();
  #ifdef __USE_NEW_SINGULAR_INTEG__
  sh_order_self=sht_up_.getShOrder();
  #else
  sh_order_self=sht_   .getShOrder();
  #endif
  assert(sim_par.sh_order     ==sht_   .getShOrder());
  assert(sim_par.upsample_freq==sht_up_.getShOrder());
  { // Set quad_weights_ (sh_order_far)
    quad_weights_.resize(1,sh_order_far);
    quad_weights_.getDevice().Memcpy(quad_weights_.begin(),
        mats.quad_weights_p_up_,
        quad_weights_.size() * sizeof(Real_t),
        Surf_t::device_type::MemcpyDeviceToDevice);
  }

  { // Set w_sph_inv_ (sh_order)
    Sca_t w_sph_;
    w_sph_    .resize(1, sh_order);
    w_sph_inv_.resize(1, sh_order);
    int np = w_sph_.getStride();

    w_sph_.getDevice().Memcpy(w_sph_.begin(), mats.w_sph_,
        np * sizeof(Real_t), device_type::MemcpyDeviceToDevice);
    xInv(w_sph_,w_sph_inv_);
  }
  { // Set sing_quad_weights_, w_sph_sing_quad_weights_ (sh_order_self)
    Sca_t w_sph_;
    w_sph_.resize(1, sh_order_self);
    int np = w_sph_.getStride();

    if(sh_order_self==sht_up_.getShOrder()){
      w_sph_.getDevice().Memcpy(w_sph_.begin(), mats.w_sph_up_,
          np * sizeof(Real_t), device_type::MemcpyDeviceToDevice);
    }else{
      w_sph_.getDevice().Memcpy(w_sph_.begin(), mats.w_sph_,
          np * sizeof(Real_t), device_type::MemcpyDeviceToDevice);
    }

    //Singular quadrature weights
    sing_quad_weights_.resize(1,sh_order_self);
    if(sh_order_self==sht_up_.getShOrder()){
      sing_quad_weights_.getDevice().Memcpy(sing_quad_weights_.begin(),
          mats.sing_quad_weights_up_, sing_quad_weights_.size() *
          sizeof(Real_t),
          device_type::MemcpyDeviceToDevice);
    }else{
      sing_quad_weights_.getDevice().Memcpy(sing_quad_weights_.begin(),
          mats.sing_quad_weights_, sing_quad_weights_.size() *
          sizeof(Real_t),
          device_type::MemcpyDeviceToDevice);
    }

    // Precompute w_sph_sing_quad_
    w_sph_sing_quad_weights_.resize(1, sh_order_self);
    ax(sing_quad_weights_, w_sph_, w_sph_sing_quad_weights_);
  }

  { // compute quadrature to find pole
    size_t p=sh_order;

    //Gauss-Legendre quadrature nodes and weights
    std::vector<double> x(p+1),w(p+1);
    cgqf(p+1, 1, 0.0, 0.0, -1.0, 1.0, &x[0], &w[0]);

    std::vector<double> leg((p+1)*(p+1));
    LegPoly(&x[0],x.size(),p+1,&leg[0]);
    pole_quad.ReInit(p+1);
    pole_quad.SetZero();
    for(size_t j=0;j<p+1;j++){
      for(size_t i=0;i<p+1;i++){
        pole_quad[i]+=leg[j*(p+1)+i]*sqrt(2.0*j+1.0);
      }
    }
    for(size_t i=0;i<p+1;i++){
      pole_quad[i]*=w[i]*0.25/p;
    }
  }

  pvfmm_ctx=PVFMMCreateContext<Real_t>(box_size_);
}

template<typename Surf_t>
StokesVelocity<Surf_t>::~StokesVelocity(){
  PVFMMDestroyContext<Real_t>(&pvfmm_ctx);
  if(sh_order_far!=sh_order && S_far) delete S_far;
  if(S_self) delete S_self;
}



template<typename Surf_t>
void StokesVelocity<Surf_t>::SetSrcCoord(const Surf_t& S_){
  self_flag=self_flag | StokesVelocity::UpdateSrcCoord;
  fmm_flag=fmm_flag | StokesVelocity::UpdateSrcCoord;
  near_singular.SetSrcCoord(src_coord_up);
  S=&S_;

  //#ifndef NDEBUG
  MonitorError();
  //#endif
}

template<typename Surf_t>
void StokesVelocity<Surf_t>::SetSurfaceVel(const Vec_t& S_vel_){
  self_flag=self_flag | StokesVelocity::UpdateSurfaceVel;
  S_vel_ptr=&S_vel_;
}

template<typename Surf_t>
void StokesVelocity<Surf_t>::SetDensitySL(const Vec_t* force_single_){
  self_flag=self_flag | StokesVelocity::UpdateDensitySL;
  fmm_flag=fmm_flag | StokesVelocity::UpdateDensitySL;
  force_single=force_single_;
}

template<typename Surf_t>
void StokesVelocity<Surf_t>::SetDensityDL(const Vec_t* force_double_){
  self_flag=self_flag | StokesVelocity::UpdateDensityDL;
  fmm_flag=fmm_flag | StokesVelocity::UpdateDensityDL;
  force_double=force_double_;
}



template<typename Surf_t>
void StokesVelocity<Surf_t>::SetTrgCoord(const Surf_t& T_){
  fmm_flag=fmm_flag | StokesVelocity::UpdateTrgCoord;
  size_t omp_p=omp_get_max_threads();

  trg_is_surf=true;
  const Vec_t& x=T_.getPosition();
  size_t N_ves = x.getNumSubs(); // Number of vesicles
  size_t M_ves = x.getStride(); // Points per vesicle
  assert(M_ves==x.getGridDim().first*x.getGridDim().second);
  trg_coord.ReInit(N_ves*M_ves*COORD_DIM);

  #pragma omp parallel for
  for(size_t tid=0;tid<omp_p;tid++){ // Set tree pt data
    size_t a=((tid+0)*N_ves)/omp_p;
    size_t b=((tid+1)*N_ves)/omp_p;
    for(size_t i=a;i<b;i++){
      // read each component of x
      const Real_t* xk=x.getSubN_begin(i)+0*M_ves;
      const Real_t* yk=x.getSubN_begin(i)+1*M_ves;
      const Real_t* zk=x.getSubN_begin(i)+2*M_ves;

      for(size_t j=0;j<M_ves;j++){
        trg_coord[(i*M_ves+j)*COORD_DIM+0]=xk[j];
        trg_coord[(i*M_ves+j)*COORD_DIM+1]=yk[j];
        trg_coord[(i*M_ves+j)*COORD_DIM+2]=zk[j];
      }
    }
  }
  near_singular.SetTrgCoord(&trg_coord[0],trg_coord.Dim()/COORD_DIM,trg_is_surf);
}

template<typename Surf_t>
void StokesVelocity<Surf_t>::SetTrgCoord(Real_t* trg_coord_, size_t N){ // TODO: change to const Real_t*
  fmm_flag=fmm_flag | StokesVelocity::UpdateTrgCoord;
  trg_is_surf=false;
  trg_coord.ReInit(N*COORD_DIM, trg_coord_);
  near_singular.SetTrgCoord(&trg_coord[0],trg_coord.Dim()/COORD_DIM,trg_is_surf);
}



template<typename Surf_t>
void StokesVelocity<Surf_t>::GetPole(const Vec_t& v,  PVFMMVec_t& pvfmm_v){
  assert(v.getShOrder()==sh_order);
  size_t omp_p=omp_get_max_threads();
  size_t N_ves   =v.getNumSubs(); // Number of vesicles
  size_t M_ves   =(1+sh_order    )*(2*sh_order    );
  size_t M_ves_up=(1+sh_order_far)*(2*sh_order_far);
  if(pvfmm_v.Dim()!=N_ves*(2+M_ves_up)*COORD_DIM){
    pvfmm_v.ReInit(N_ves*(2+M_ves_up)*COORD_DIM);
  }

  #pragma omp parallel for
  for(size_t tid=0;tid<omp_p;tid++){
    size_t a=((tid+0)*N_ves)/omp_p;
    size_t b=((tid+1)*N_ves)/omp_p;
    for(size_t i=a;i<b;i++){
      const Real_t* x=v.getSubN_begin(i)+0*M_ves;
      const Real_t* y=v.getSubN_begin(i)+1*M_ves;
      const Real_t* z=v.getSubN_begin(i)+2*M_ves;

      Real_t pole[2*COORD_DIM];
      for(size_t j=0;j<2*COORD_DIM;j++) pole[j]=0;
      for(size_t k0=0;k0<(1+sh_order);k0++){
        for(size_t k1=0;k1<(2*sh_order);k1++){
          size_t k=k1+k0*(2*sh_order);
          pole[0*COORD_DIM+0]+=pole_quad[sh_order-k0]*x[k];
          pole[0*COORD_DIM+1]+=pole_quad[sh_order-k0]*y[k];
          pole[0*COORD_DIM+2]+=pole_quad[sh_order-k0]*z[k];
          pole[1*COORD_DIM+0]+=pole_quad[         k0]*x[k];
          pole[1*COORD_DIM+1]+=pole_quad[         k0]*y[k];
          pole[1*COORD_DIM+2]+=pole_quad[         k0]*z[k];
        }
      }

      pvfmm_v[(i*(2+M_ves_up)+0)*COORD_DIM+0]=pole[0*COORD_DIM+0];
      pvfmm_v[(i*(2+M_ves_up)+0)*COORD_DIM+1]=pole[0*COORD_DIM+1];
      pvfmm_v[(i*(2+M_ves_up)+0)*COORD_DIM+2]=pole[0*COORD_DIM+2];
      pvfmm_v[(i*(2+M_ves_up)+1)*COORD_DIM+0]=pole[1*COORD_DIM+0];
      pvfmm_v[(i*(2+M_ves_up)+1)*COORD_DIM+1]=pole[1*COORD_DIM+1];
      pvfmm_v[(i*(2+M_ves_up)+1)*COORD_DIM+2]=pole[1*COORD_DIM+2];
    }
  }
}

template<typename Surf_t>
void StokesVelocity<Surf_t>::Vec2PVFMMVec(const Vec_t& v,  PVFMMVec_t& pvfmm_v){
  size_t omp_p=omp_get_max_threads();
  size_t N_ves = v.getNumSubs(); // Number of vesicles
  size_t M_ves = v.getStride(); // Points per vesicle
  if(pvfmm_v.Dim()!=N_ves*(2+M_ves)*COORD_DIM){
    pvfmm_v.ReInit(N_ves*(2+M_ves)*COORD_DIM);
  }

  #pragma omp parallel for
  for(size_t tid=0;tid<omp_p;tid++){
    size_t a=((tid+0)*N_ves)/omp_p;
    size_t b=((tid+1)*N_ves)/omp_p;
    for(size_t i=a;i<b;i++){
      const Real_t* x=v.getSubN_begin(i)+0*M_ves;
      const Real_t* y=v.getSubN_begin(i)+1*M_ves;
      const Real_t* z=v.getSubN_begin(i)+2*M_ves;
      for(size_t j=0;j<M_ves;j++){
        pvfmm_v[(i*(2+M_ves)+j+2)*COORD_DIM+0]=x[j];
        pvfmm_v[(i*(2+M_ves)+j+2)*COORD_DIM+1]=y[j];
        pvfmm_v[(i*(2+M_ves)+j+2)*COORD_DIM+2]=z[j];
      }
    }
  }
}

template<typename Surf_t>
void StokesVelocity<Surf_t>::Upsample(const Vec_t& v, Vec_t* v_out,  PVFMMVec_t* pvfmm_v){
  assert(v.getShOrder()==sh_order);

  Vec_t wrk[3]; // TODO: Pre-allocate
  if(!v_out) v_out=&wrk[2];
  wrk[0].resize(v.getNumSubs(), std::max(sh_order_far,sh_order));
  wrk[1].resize(v.getNumSubs(), std::max(sh_order_far,sh_order));
  v_out->resize(v.getNumSubs(), std::max(sh_order_far,sh_order));
  Resample(v, sht_, sht_up_, wrk[0], wrk[1], *v_out);

  if(pvfmm_v){
    Vec2PVFMMVec(*v_out, *pvfmm_v);
    GetPole     ( v    , *pvfmm_v);
  }
}

template<typename Surf_t>
void StokesVelocity<Surf_t>::Setup(){

  if(force_single && force_single!=&F_repl){ // Add repulsion
    const PVFMMVec_t& f_repl=near_singular.ForceRepul();
    F_repl.resize(force_single->getNumSubs(), sh_order);
    size_t omp_p=omp_get_max_threads();

    Vec_t& x=F_repl;
    size_t N_ves = x.getNumSubs(); // Number of vesicles
    size_t M_ves = x.getStride(); // Points per vesicle
    assert(M_ves==x.getGridDim().first*x.getGridDim().second);
    assert(N_ves*M_ves*COORD_DIM==f_repl.Dim());

    #pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){ // Set tree pt data
      size_t a=((tid+0)*N_ves)/omp_p;
      size_t b=((tid+1)*N_ves)/omp_p;
      for(size_t i=a;i<b;i++){
        // read each component of x
        Real_t* xk=x.getSubN_begin(i)+0*M_ves;
        Real_t* yk=x.getSubN_begin(i)+1*M_ves;
        Real_t* zk=x.getSubN_begin(i)+2*M_ves;

        for(size_t j=0;j<M_ves;j++){
          xk[j]=f_repl[(i*M_ves+j)*COORD_DIM+0];
          yk[j]=f_repl[(i*M_ves+j)*COORD_DIM+1];
          zk[j]=f_repl[(i*M_ves+j)*COORD_DIM+2];
        }
      }
    }

    axpy(1.0,*force_single,F_repl,F_repl);
    force_single=&F_repl;
  }

  assert(S);
  size_t omp_p=omp_get_max_threads();
  if(fmm_flag & StokesVelocity::UpdateSrcCoord){ // Compute src_coord_up
    fmm_flag=fmm_flag & ~StokesVelocity::UpdateSrcCoord;

    if(sh_order_far==sh_order) S_far=S;
    else CHK(S->resample(sh_order_far, (Surf_t**)&S_far));
    const Vec_t& S_coord   =S    ->getPosition();
    const Vec_t& S_coord_up=S_far->getPosition();
    Vec2PVFMMVec(S_coord_up, src_coord_up);
    GetPole     (S_coord   , src_coord_up);
  }

  Vec_t qforce; // TODO: Pre-allocate
  const int force_dim=COORD_DIM;
  if(fmm_flag & StokesVelocity::UpdateDensitySL){ // Compute qforce_single_up
    fmm_flag=fmm_flag & ~StokesVelocity::UpdateDensitySL;
    fmm_flag=fmm_flag | StokesVelocity::UpdateqDensitySL;
    if(force_single){
      Upsample(*force_single, &qforce);
      size_t N_ves = qforce.getNumSubs(); // Number of vesicles
      size_t M_ves = qforce.getStride(); // Points per vesicle

      xv(S_far->getAreaElement(), qforce, qforce);
      ax<Sca_t>(   quad_weights_, qforce, qforce);

      qforce_single_up.ReInit(N_ves*(M_ves+2)*(force_dim));
      #pragma omp parallel for
      for(size_t tid=0;tid<omp_p;tid++){
        size_t a=((tid+0)*N_ves)/omp_p;
        size_t b=((tid+1)*N_ves)/omp_p;
        for(size_t i=a;i<b;i++){ // loop over all vesicles
          // read each component of force_single
          const Real_t* fxk=qforce.getSubN_begin(i)+0*M_ves;
          const Real_t* fyk=qforce.getSubN_begin(i)+1*M_ves;
          const Real_t* fzk=qforce.getSubN_begin(i)+2*M_ves;

          for(size_t j=0;j<M_ves;j++){
            qforce_single_up[(i*(M_ves+2)+j+2)*force_dim+0]=fxk[j];
            qforce_single_up[(i*(M_ves+2)+j+2)*force_dim+1]=fyk[j];
            qforce_single_up[(i*(M_ves+2)+j+2)*force_dim+2]=fzk[j];
          }
          for(size_t j=0;j<force_dim;j++){
            qforce_single_up[(i*(M_ves+2)+0)*force_dim+j]=0;
            qforce_single_up[(i*(M_ves+2)+1)*force_dim+j]=0;
          }
        }
      }
      near_singular.SetDensitySL(&qforce_single_up);
    }else{
      qforce_single_up.ReInit(0);
      near_singular.SetDensitySL(NULL);
    }
  }
  if(fmm_flag & StokesVelocity::UpdateDensityDL){ // Compute qforce_double_up
    fmm_flag=fmm_flag & ~StokesVelocity::UpdateDensityDL;
    fmm_flag=fmm_flag | StokesVelocity::UpdateqDensityDL;
    if(force_double){
      Upsample(*force_double, &qforce, &force_double_up);
      size_t N_ves = qforce.getNumSubs(); // Number of vesicles
      size_t M_ves = qforce.getStride(); // Points per vesicle

      xv(S_far->getAreaElement(), qforce, qforce);
      ax<Sca_t>(   quad_weights_, qforce, qforce);

      qforce_double_up.ReInit(N_ves*(M_ves+2)*(force_dim+COORD_DIM));
      const Vec_t* normal=&S_far->getNormal();
      #pragma omp parallel for
      for(size_t tid=0;tid<omp_p;tid++){
        size_t a=((tid+0)*N_ves)/omp_p;
        size_t b=((tid+1)*N_ves)/omp_p;
        for(size_t i=a;i<b;i++){ // loop over all vesicles
          // read each component of force_double
          const Real_t* fxk=qforce.getSubN_begin(i)+0*M_ves;
          const Real_t* fyk=qforce.getSubN_begin(i)+1*M_ves;
          const Real_t* fzk=qforce.getSubN_begin(i)+2*M_ves;

          // read each component of normal
          const Real_t* nxk=normal->getSubN_begin(i)+0*M_ves;
          const Real_t* nyk=normal->getSubN_begin(i)+1*M_ves;
          const Real_t* nzk=normal->getSubN_begin(i)+2*M_ves;

          for(size_t j=0;j<M_ves;j++){
            qforce_double_up[(i*(M_ves+2)+j+2)*(force_dim+COORD_DIM)+0]=fxk[j];
            qforce_double_up[(i*(M_ves+2)+j+2)*(force_dim+COORD_DIM)+1]=fyk[j];
            qforce_double_up[(i*(M_ves+2)+j+2)*(force_dim+COORD_DIM)+2]=fzk[j];
            qforce_double_up[(i*(M_ves+2)+j+2)*(force_dim+COORD_DIM)+3]=nxk[j];
            qforce_double_up[(i*(M_ves+2)+j+2)*(force_dim+COORD_DIM)+4]=nyk[j];
            qforce_double_up[(i*(M_ves+2)+j+2)*(force_dim+COORD_DIM)+5]=nzk[j];
          }
          for(size_t j=0;j<force_dim+COORD_DIM;j++){
            qforce_double_up[(i*(M_ves+2)+0)*(force_dim+COORD_DIM)+j]=0;
            qforce_double_up[(i*(M_ves+2)+1)*(force_dim+COORD_DIM)+j]=0;
          }
        }
      }
      near_singular.SetDensityDL(&qforce_double_up, &force_double_up);
    }else{
      qforce_double_up.ReInit(0);
      near_singular.SetDensityDL(NULL, NULL);
    }
  }
}


template<typename Surf_t>
template<bool SL, bool DL>
void StokesVelocity<Surf_t>::SingularInteg(){ // (deprecated)
  if(!SL && !DL) return;

  assert(S);
  CHK(S->resample(sh_order_self, (Surf_t**)&S_self));

  int imax(S->getPosition().getGridDim().first);
  int jmax(S->getPosition().getGridDim().second);
  int nv     = S     ->getPosition().getNumSubs();
  int np     = S     ->getPosition().getStride();
  int np_src = S_self->getPosition().getStride();

  Vec_t  coord_up,  force_single_up,  force_double_up;
  Vec_t        coord_out(nv, sh_order);
  Vec_t force_single_out(nv, sh_order);
  Vec_t force_double_out(nv, sh_order);
  Vec_t force_single_area_elem(nv, sh_order_self);
  Vec_t force_double_area_elem(nv, sh_order_self);

  std::vector<const Sca_t*>  inputs;
  std::vector<      Sca_t*> outputs;
  { // set Position
    inputs .push_back(&S->getPosition());
    outputs.push_back(&       coord_out);
  }
  if(SL){ // set SL force
    inputs .push_back( force_single    );
    outputs.push_back(&force_single_out);
  }
  if(DL){ // set DL force
    inputs .push_back( force_double    );
    outputs.push_back(&force_double_out);
  }
  move_pole.setOperands(&inputs[0], inputs.size(), sim_par.singular_stokes);

  for(int ii=0;ii < imax; ++ii)
  for(int jj=0;jj < jmax; ++jj){
    move_pole(ii, jj, &outputs[0]);

    Vec_t *force_single, *force_double;
    if(sh_order==sh_order_self){
      coord_out.getDevice().Memcpy(S_self->getPositionModifiable().begin(),
          coord_out.begin(), coord_out.size()*sizeof(Real_t), device_type::MemcpyDeviceToDevice);
      force_single=&force_single_out;
      force_double=&force_double_out;
    }else{
      Upsample(coord_out, &S_self->getPositionModifiable());
      if(SL) Upsample(force_single_out, &force_single_up);
      if(DL) Upsample(force_double_out, &force_double_up);
      force_single=&force_single_up;
      force_double=&force_double_up;
    }

    if(SL){
      xv(S_self->getAreaElement(), *force_single, force_single_area_elem);
      S->getPosition().getDevice().DirectStokes           (
          S_self->getPosition().begin(),                              force_single_area_elem.begin(),
          sing_quad_weights_.begin(), np_src, np, nv, S->getPosition().begin(),
          ii * jmax + jj, ii * jmax + jj + 1, SL_vel.begin());
    }

    if(DL){
      xv(S_self->getAreaElement(), *force_double, force_double_area_elem);
      S->getPosition().getDevice().DirectStokesDoubleLayer(
          S_self->getPosition().begin(), S_self->getNormal().begin(), force_double_area_elem.begin(),
          sing_quad_weights_.begin(), np_src, np, nv, S->getPosition().begin(),
          ii * jmax + jj, ii * jmax + jj + 1, DL_vel.begin());
    }
  }
}

template<typename Surf_t>
const StokesVelocity<Surf_t>::Vec_t& StokesVelocity<Surf_t>::SelfInteraction(bool update_self){
  if(self_flag & StokesVelocity::UpdateSurfaceVel){
    if(S_vel_ptr){ // Nothing to do; it has already been computed for us.
      Upsample(*S_vel_ptr, NULL,  &surf_vel_up);
      near_singular.SetSurfaceVel(&surf_vel_up);
      self_flag=StokesVelocity::UpdateNone;
    }else{ // Compute from SL and DL density
      self_flag=self_flag | StokesVelocity::UpdateDensitySL;
      self_flag=self_flag | StokesVelocity::UpdateDensityDL;
    }
  }
  if(update_self && self_flag){ // Compute self interaction
    { // Original self-interaction code (deprecated)
      #ifdef __CHECK_SINGULAR_INTEG__
      pvfmm::Profile::Tic("SelfInteraction",&comm);
      bool prof_state=pvfmm::Profile::Enable(true);//false);
      assert(S);
      int imax(S->getPosition().getGridDim().first);
      int jmax(S->getPosition().getGridDim().second);
      int np = S->getPosition().getStride();
      int nv = S->getPosition().getNumSubs();
      if((self_flag & (StokesVelocity::UpdateDensitySL | StokesVelocity::UpdateSrcCoord)) &&
         (self_flag & (StokesVelocity::UpdateDensityDL | StokesVelocity::UpdateSrcCoord)) &&
         force_single && force_double){ // Single and Double layer
        SL_vel.resize(nv, sh_order);
        DL_vel.resize(nv, sh_order);
        {
          #ifdef __USE_NEW_SINGULAR_INTEG__
          SingularInteg<true,true>();
          #else
          assert(sh_order==sh_order_self);
          Vec_t coord_out(nv, sh_order);
          Vec_t force_single_in (nv, sh_order);
          Vec_t force_single_out(nv, sh_order);

          Vec_t normal_out(nv, sh_order);
          Vec_t force_double_in (nv, sh_order);
          Vec_t force_double_out(nv, sh_order);

          Sca_t t1(nv, sh_order);
          ax(w_sph_inv_, S->getAreaElement(), t1);
          xv(        t1,       *force_single, force_single_in);
          xv(        t1,       *force_double, force_double_in);

          int numinputs = 4;
          const Sca_t* inputs[] = {&S->getPosition(), &force_single_in , &S->getNormal(), &force_double_in };
          Sca_t*      outputs[] = {&       coord_out, &force_single_out, &    normal_out, &force_double_out};
          move_pole.setOperands(inputs, numinputs, sim_par.singular_stokes);

          for(int ii=0;ii < imax; ++ii)
          for(int jj=0;jj < jmax; ++jj){
            move_pole(ii, jj, outputs);

            S->getPosition().getDevice().DirectStokes(coord_out.begin(), force_single_out.begin(),
                w_sph_sing_quad_weights_.begin(), np, np, nv, S->getPosition().begin(),
                ii * jmax + jj, ii * jmax + jj + 1, SL_vel.begin());

            S->getPosition().getDevice().DirectStokesDoubleLayer(coord_out.begin(), normal_out.begin(), force_double_out.begin(),
                w_sph_sing_quad_weights_.begin(), np, np, nv, S->getPosition().begin(),
                ii * jmax + jj, ii * jmax + jj + 1, DL_vel.begin());
          }
          #endif
        }
      }else{
        if(self_flag & (StokesVelocity::UpdateDensitySL | StokesVelocity::UpdateSrcCoord)){ // Single layer
          SL_vel.resize(nv, sh_order);
          if(force_single){
            #ifdef __USE_NEW_SINGULAR_INTEG__
            SingularInteg<true,false>();
            #else
            assert(sh_order==sh_order_self);
            Vec_t coord_out(nv, sh_order);
            Vec_t force_single_in (nv, sh_order);
            Vec_t force_single_out(nv, sh_order);

            Sca_t t1(nv, sh_order);
            ax(w_sph_inv_, S->getAreaElement(), t1);
            xv(        t1,       *force_single, force_single_in);

            int numinputs = 2;
            const Sca_t* inputs[] = {&S->getPosition(), &force_single_in };
            Sca_t*      outputs[] = {&       coord_out, &force_single_out};
            move_pole.setOperands(inputs, numinputs, sim_par.singular_stokes);

            for(int ii=0;ii < imax; ++ii)
            for(int jj=0;jj < jmax; ++jj){
              move_pole(ii, jj, outputs);

              S->getPosition().getDevice().DirectStokes(coord_out.begin(), force_single_out.begin(),
                  w_sph_sing_quad_weights_.begin(), np, np, nv, S->getPosition().begin(),
                  ii * jmax + jj, ii * jmax + jj + 1, SL_vel.begin());
            }
            #endif
          }else{
            Vec_t::getDevice().Memset(SL_vel.begin(),0,SL_vel.size()*sizeof(Real_t));
          }
        }
        if(self_flag & (StokesVelocity::UpdateDensityDL | StokesVelocity::UpdateSrcCoord)){ // Double layer
          DL_vel.resize(nv, sh_order);
          if(force_double){
            #ifdef __USE_NEW_SINGULAR_INTEG__
            SingularInteg<false,true>();
            #else
            assert(sh_order==sh_order_self);
            Vec_t coord_out(nv, sh_order);
            Vec_t normal_out(nv, sh_order);
            Vec_t force_double_in (nv, sh_order);
            Vec_t force_double_out(nv, sh_order);

            Sca_t t1(nv, sh_order);
            ax(w_sph_inv_, S->getAreaElement(), t1);
            xv(        t1,       *force_double, force_double_in);

            int numinputs = 3;
            const Sca_t* inputs[] = {&S->getPosition(), &S->getNormal(), &force_double_in };
            Sca_t*      outputs[] = {&       coord_out, &    normal_out, &force_double_out};
            move_pole.setOperands(inputs, numinputs, sim_par.singular_stokes);

            for(int ii=0;ii < imax; ++ii)
            for(int jj=0;jj < jmax; ++jj){
              move_pole(ii, jj, outputs);

              S->getPosition().getDevice().DirectStokesDoubleLayer(coord_out.begin(), normal_out.begin(), force_double_out.begin(),
                  w_sph_sing_quad_weights_.begin(), np, np, nv, S->getPosition().begin(),
                  ii * jmax + jj, ii * jmax + jj + 1, DL_vel.begin());
            }
            #endif
          }else{
            Vec_t::getDevice().Memset(DL_vel.begin(),0,DL_vel.size()*sizeof(Real_t));
          }
        }
      }
      S_vel.resize(nv, sh_order);
      axpy(1.0,SL_vel,DL_vel,S_vel);
      pvfmm::Profile::Enable(prof_state);
      pvfmm::Profile::Toc();
      #endif
    }
    { // Compute S_vel
      pvfmm::Profile::Tic("SelfInteraction",&comm);
      bool prof_state=pvfmm::Profile::Enable(false);

      assert(S);
      assert(S->getShOrder().first==sh_order+1);
      long Ngrid = 2*sh_order*(sh_order+1);
      long Ncoef = sh_order*(sh_order+2);
      const Vec_t& x=S->getPosition();
      long nv = x.getNumSubs();


      { // Compute self-interaction matrices
        if(self_flag & StokesVelocity::UpdateSrcCoord){
          if(!force_single) SLMatrix.ReInit(0);
          if(!force_double) DLMatrix.ReInit(0);
        }
        pvfmm::Vector<Real_t> *SLMatrix_=NULL, *DLMatrix_=NULL;
        if(force_single && (!SLMatrix.Dim() || (self_flag & StokesVelocity::UpdateSrcCoord))) SLMatrix_=&SLMatrix;
        if(force_double && (!DLMatrix.Dim() || (self_flag & StokesVelocity::UpdateSrcCoord))) DLMatrix_=&DLMatrix;
        pvfmm::Vector<Real_t> x_vec(nv*COORD_DIM*Ngrid,(Real_t*)x.begin(),false);
        SphericalHarmonics<Real_t>::StokesSingularInteg(x_vec, sh_order, sh_order_self, SLMatrix_, DLMatrix_);
      }

      SL_vel.resize(nv, sh_order);
      if(force_single){
        pvfmm::Vector<Real_t> F, F_;
        F.ReInit(nv*COORD_DIM*Ngrid,(Real_t*)force_single->begin(),false);
        SphericalHarmonics<Real_t>::Grid2SHC(F,sh_order,sh_order,F_);

        pvfmm::Vector<Real_t> V_(nv*COORD_DIM*Ncoef);
        #pragma omp parallel
        { // mat-vec
          long tid=omp_get_thread_num();
          long omp_p=omp_get_num_threads();

          long a=(tid+0)*nv/omp_p;
          long b=(tid+1)*nv/omp_p;
          for(long i=a;i<b;i++){
            pvfmm::Matrix<Real_t> Mv(1,COORD_DIM*Ncoef,&V_[i*COORD_DIM*Ncoef],false);
            pvfmm::Matrix<Real_t> Mf(1,COORD_DIM*Ncoef,&F_[i*COORD_DIM*Ncoef],false);
            pvfmm::Matrix<Real_t> M(COORD_DIM*Ncoef,COORD_DIM*Ncoef,&SLMatrix[i*COORD_DIM*Ncoef*COORD_DIM*Ncoef],false);
            pvfmm::Matrix<Real_t>::GEMM(Mv,Mf,M);
          }
        }
        { // print error
          #ifdef __CHECK_SINGULAR_INTEG__
          Real_t err=0;
          pvfmm::Vector<Real_t> Vgrid(nv*COORD_DIM*Ngrid,(Real_t*)SL_vel.begin(),false), Vshc;
          SphericalHarmonics<Real_t>::Grid2SHC(Vgrid,sh_order,sh_order,Vshc);
          for(long i=0;i<nv*COORD_DIM*Ncoef;i++) err=std::max(err,fabs(Vshc[i]-V_[i]));
          INFO("StokesVelocity: SL-Error: "<<err);
          #endif
        }
        pvfmm::Vector<Real_t> V(nv*COORD_DIM*Ngrid,(Real_t*)SL_vel.begin(),false);
        SphericalHarmonics<Real_t>::SHC2Grid(V_,sh_order,sh_order,V);
      }else{
        Vec_t::getDevice().Memset(SL_vel.begin(),0,SL_vel.size()*sizeof(Real_t));
      }

      DL_vel.resize(nv, sh_order);
      if(force_double){
        pvfmm::Vector<Real_t> F, F_;
        F.ReInit(nv*COORD_DIM*Ngrid,(Real_t*)force_double->begin(),false);
        SphericalHarmonics<Real_t>::Grid2SHC(F,sh_order,sh_order,F_);

        pvfmm::Vector<Real_t> V_(nv*COORD_DIM*Ncoef);
        #pragma omp parallel
        { // mat-vec
          long tid=omp_get_thread_num();
          long omp_p=omp_get_num_threads();

          long a=(tid+0)*nv/omp_p;
          long b=(tid+1)*nv/omp_p;
          for(long i=a;i<b;i++){
            pvfmm::Matrix<Real_t> Mv(1,COORD_DIM*Ncoef,&V_[i*COORD_DIM*Ncoef],false);
            pvfmm::Matrix<Real_t> Mf(1,COORD_DIM*Ncoef,&F_[i*COORD_DIM*Ncoef],false);
            pvfmm::Matrix<Real_t> M(COORD_DIM*Ncoef,COORD_DIM*Ncoef,&DLMatrix[i*COORD_DIM*Ncoef*COORD_DIM*Ncoef],false);
            pvfmm::Matrix<Real_t>::GEMM(Mv,Mf,M);
          }
        }
        { // print error
          #ifdef __CHECK_SINGULAR_INTEG__
          Real_t err=0;
          pvfmm::Vector<Real_t> Vgrid(nv*COORD_DIM*Ngrid,(Real_t*)DL_vel.begin(),false), Vshc;
          SphericalHarmonics<Real_t>::Grid2SHC(Vgrid,sh_order,sh_order,Vshc);
          for(long i=0;i<nv*COORD_DIM*Ncoef;i++) err=std::max(err,fabs(Vshc[i]-V_[i]));
          INFO("StokesVelocity: DL-Error: "<<err);
          #endif
        }
        pvfmm::Vector<Real_t> V(nv*COORD_DIM*Ngrid,(Real_t*)DL_vel.begin(),false);
        SphericalHarmonics<Real_t>::SHC2Grid(V_,sh_order,sh_order,V);
      }else{
        Vec_t::getDevice().Memset(DL_vel.begin(),0,DL_vel.size()*sizeof(Real_t));
      }

      S_vel.resize(nv, sh_order);
      axpy(1.0,SL_vel,DL_vel,S_vel);
      pvfmm::Profile::Enable(prof_state);
      pvfmm::Profile::Toc();
    }
    S_vel_ptr=&S_vel;
    Upsample(*S_vel_ptr, NULL,  &surf_vel_up);
    near_singular.SetSurfaceVel(&surf_vel_up);
    self_flag=StokesVelocity::UpdateNone;
  }
  return *S_vel_ptr;
}

template<typename Surf_t>
const StokesVelocity<Surf_t>::PVFMMVec_t& StokesVelocity<Surf_t>::NearInteraction(bool update_near){
  if(update_near) Setup();
  return near_singular(update_near);
}

template<typename Surf_t>
const StokesVelocity<Surf_t>::PVFMMVec_t& StokesVelocity<Surf_t>::FarInteraction(bool update_far){
  if(update_far && fmm_flag){ // Compute velocity with FMM
    pvfmm::Profile::Tic("FarInteraction",&comm);
    bool prof_state=pvfmm::Profile::Enable(false);
    Setup();
    fmm_vel.ReInit(trg_coord.Dim());
    PVFMMEval(&src_coord_up[0],
              (qforce_single_up.Dim()?&qforce_single_up[0]:NULL),
              (qforce_double_up.Dim()?&qforce_double_up[0]:NULL),
              src_coord_up.Dim()/COORD_DIM,
              &trg_coord[0], &fmm_vel[0], trg_coord.Dim()/COORD_DIM, &pvfmm_ctx);
    near_singular.SubtractDirect(fmm_vel);
    fmm_flag=StokesVelocity::UpdateNone;
    pvfmm::Profile::Enable(prof_state);
    pvfmm::Profile::Toc();
  }
  return fmm_vel;
}



template<typename Surf_t>
StokesVelocity<Surf_t>::Real_t* StokesVelocity<Surf_t>::operator()(unsigned int flag){
  bool update_self=(flag & StokesVelocity::UpdateSelf);
  bool update_near=(flag & StokesVelocity::UpdateNear);
  bool update_far =(flag & StokesVelocity::UpdateFar );
  size_t omp_p=omp_get_max_threads();

#ifdef __ENABLE_PVFMM_PROFILER__
  bool prof_state=pvfmm::Profile::Enable(true);
#endif
  SelfInteraction(update_self);
  const PVFMMVec_t& near_vel=NearInteraction(update_near);
  const PVFMMVec_t& far_vel=FarInteraction(update_far);
  { // Compute trg_vel = far_vel + near_vel
    trg_vel.ReInit(trg_coord.Dim());
    assert(trg_vel.Dim()==far_vel.Dim());
    assert(trg_vel.Dim()==near_vel.Dim());
    #pragma omp parallel for
    for(size_t i=0;i<trg_vel.Dim();i++){
      trg_vel[i]=far_vel[i]+near_vel[i];
    }
  }
  if(trg_is_surf){
    Vec_t& x=S_vel;
    size_t N_ves = x.getNumSubs(); // Number of vesicles
    size_t M_ves = x.getStride(); // Points per vesicle
    assert(M_ves==x.getGridDim().first*x.getGridDim().second);
    assert(N_ves*M_ves*COORD_DIM==trg_vel.Dim());

    #pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){ // Set tree pt data
      size_t a=((tid+0)*N_ves)/omp_p;
      size_t b=((tid+1)*N_ves)/omp_p;
      for(size_t i=a;i<b;i++){
        // read each component of x
        Real_t* xk=x.getSubN_begin(i)+0*M_ves;
        Real_t* yk=x.getSubN_begin(i)+1*M_ves;
        Real_t* zk=x.getSubN_begin(i)+2*M_ves;

        for(size_t j=0;j<M_ves;j++){
          trg_vel[(i*M_ves+j)*COORD_DIM+0]+=xk[j];
          trg_vel[(i*M_ves+j)*COORD_DIM+1]+=yk[j];
          trg_vel[(i*M_ves+j)*COORD_DIM+2]+=zk[j];
        }
      }
    }
  }
#ifdef __ENABLE_PVFMM_PROFILER__
  pvfmm::Profile::Enable(prof_state);
  pvfmm::Profile::print(&comm);
#endif
  return &trg_vel[0];
}

template<typename Surf_t>
void StokesVelocity<Surf_t>::operator()(Vec_t& T_vel, unsigned int flag){
  size_t omp_p=omp_get_max_threads();
  (*this)(flag);

  Vec_t& x=T_vel;
  size_t N_ves = x.getNumSubs(); // Number of vesicles
  size_t M_ves = x.getStride(); // Points per vesicle
  assert(M_ves==x.getGridDim().first*x.getGridDim().second);
  assert(N_ves*M_ves*COORD_DIM==trg_vel.Dim());

  #pragma omp parallel for
  for(size_t tid=0;tid<omp_p;tid++){ // Set tree pt data
    size_t a=((tid+0)*N_ves)/omp_p;
    size_t b=((tid+1)*N_ves)/omp_p;
    for(size_t i=a;i<b;i++){
      // read each component of x
      Real_t* xk=x.getSubN_begin(i)+0*M_ves;
      Real_t* yk=x.getSubN_begin(i)+1*M_ves;
      Real_t* zk=x.getSubN_begin(i)+2*M_ves;

      for(size_t j=0;j<M_ves;j++){
        xk[j]=trg_vel[(i*M_ves+j)*COORD_DIM+0];
        yk[j]=trg_vel[(i*M_ves+j)*COORD_DIM+1];
        zk[j]=trg_vel[(i*M_ves+j)*COORD_DIM+2];
      }
    }
  }
}

template<typename Surf_t>
StokesVelocity<Surf_t>::Real_t StokesVelocity<Surf_t>::MonitorError(){
  if(!S) return -1;

  // Save original state
  const Vec_t* force_single_orig=force_single;
  const Vec_t* force_double_orig=force_double;
  PVFMMVec_t trg_coord_orig=trg_coord;
  bool trg_is_surf_orig=trg_is_surf;

  Vec_t force;
  force.replicate(S->getPosition());
  size_t N_ves = force.getNumSubs(); // Number of vesicles
  size_t M_ves = force.getStride(); // Points per vesicle
  for(size_t i=0;i<N_ves;i++){ // Set force
    Real_t* xk=force.getSubN_begin(i)+0*M_ves;
    Real_t* yk=force.getSubN_begin(i)+1*M_ves;
    Real_t* zk=force.getSubN_begin(i)+2*M_ves;
    for(size_t j=0;j<M_ves;j++){
      xk[j]=1.0;
      yk[j]=1.0;
      zk[j]=1.0;
    }
  }

  SetDensitySL(NULL);
  SetDensityDL(&force);
  SetTrgCoord(*S);

  Real_t* velocity=this->operator()();
  Real_t* velocity_near=&(NearInteraction(true)[0]);
  Real_t* velocity_self=S_vel.begin();

#if HAVE_PVFMM
  if(1){ // Write VTK file
    static unsigned long iter=0;
    unsigned long skip=1;
    if(iter%skip==0){
      char fname[1000];
      sprintf(fname, "vis/test1_%05d", (int)(iter/skip));
      WriteVTK(*S, fname, MPI_COMM_WORLD, &S_vel, S->getShOrder()*2);
    }
    iter++;
  }
#endif // HAVE_PVFMM

  double norm_glb[3]={0,0,0};
  { // Compute error norm
    double norm_local[3]={0,0,0};
    for(size_t i=0;i<N_ves*M_ves*COORD_DIM;i++){
      norm_local[0]=std::max(norm_local[0],fabs(velocity_self[i]+0.5));
      norm_local[1]=std::max(norm_local[1],fabs(velocity_near[i]    ));
      norm_local[2]=std::max(norm_local[2],fabs(velocity     [i]+0.5));
    }
    MPI_Allreduce(&norm_local, &norm_glb, 3, pvfmm::par::Mpi_datatype<Real_t>::value(), pvfmm::par::Mpi_datatype<Real_t>::max(), comm);
  }

  { // Restore original state
    SetDensitySL(force_single_orig);
    SetDensityDL(force_double_orig);

    fmm_flag=fmm_flag | StokesVelocity::UpdateTrgCoord;
    trg_is_surf=trg_is_surf_orig;
    trg_coord.ReInit(trg_coord_orig.Dim(), &trg_coord_orig[0]);
    near_singular.SetTrgCoord(&trg_coord[0],trg_coord.Dim()/COORD_DIM,trg_is_surf);
  }

  INFO("StokesVelocity: Double-layer integration error: "<<norm_glb[0]<<' '<<norm_glb[1]<<' '<<norm_glb[2]);
  return norm_glb[0];
}



template<typename Surf_t>
void StokesVelocity<Surf_t>::u_ref(const Real_t* coord, int n, Real_t* out){ //Analytical velocity for sphere
  Real_t R0=1.0;
  for(int i=0;i<n;i++){
    const Real_t* c=&coord[i*COORD_DIM];
    Real_t r_2=0;
    r_2+=c[0]*c[0];
    r_2+=c[1]*c[1];
    r_2+=c[2]*c[2];
    Real_t r=sqrt(r_2);

    Real_t cos_t=c[0]/r;
    Real_t sin_t=sqrt(1-cos_t*cos_t);
    Real_t R0_r=R0/r;
    Real_t ur= cos_t*(1.0-1.50*R0_r+0.50*R0_r*R0_r*R0_r);
    Real_t ut=-sin_t*(1.0-0.75*R0_r-0.25*R0_r*R0_r*R0_r);
    out[i*COORD_DIM+0]=cos_t*ur-sin_t*ut-1.0;

    Real_t r_yz=sqrt(c[1]*c[1]+c[2]*c[2]);
    out[i*COORD_DIM+1]=(sin_t*ur+cos_t*ut)*c[1]/r_yz;
    out[i*COORD_DIM+2]=(sin_t*ur+cos_t*ut)*c[2]/r_yz;
  }
}

template<typename Surf_t>
void StokesVelocity<Surf_t>::force(const Real_t* coord, int n, Real_t* out){ // Force on sphere
  Real_t R0=1.0;
  for(int i=0;i<n;i++){
    const Real_t* c=&coord[i*COORD_DIM];
    Real_t r_2=0;
    r_2+=c[0]*c[0];
    r_2+=c[1]*c[1];
    r_2+=c[2]*c[2];
    Real_t r=sqrt(r_2);

    Real_t cos_t=c[0]/r;
    Real_t sin_t=sqrt(1-cos_t*cos_t);
    Real_t R0_r=R0/r;
    Real_t ur=-1.5*cos_t;
    Real_t ut=+1.5*sin_t;
    out[i*COORD_DIM+0]=cos_t*ur-sin_t*ut;

    Real_t r_yz=sqrt(c[1]*c[1]+c[2]*c[2]);
    out[i*COORD_DIM+1]=(sin_t*ur+cos_t*ut)*c[1]/r_yz;
    out[i*COORD_DIM+2]=(sin_t*ur+cos_t*ut)*c[2]/r_yz;
  }
}

template<typename Surf_t>
void StokesVelocity<Surf_t>::Test(){
  int rank, size;
  MPI_Comm comm=MPI_COMM_WORLD;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);
  srand48(rank);

  // Set parameters
  Parameters<Real_t> sim_par;
  sim_par.sh_order = 6;
  sim_par.upsample_freq = 32;
  sim_par.periodic_length = -5;

  // Reading operators from file
  bool readFromFile = true;
  Mats_t mats(readFromFile, sim_par);

  //Setting the background flow
  //ShearFlow<Vec_t> vInf(sim_par.bg_flow_param);

  //StokesVelocity<Surf_t> stokes_vel(mats, sim_par, vInf, comm);
  StokesVelocity<Surf_t> stokes_vel(mats, sim_par, comm);

  //=================================================================//

  // Create vectors
  size_t nVec=1;
  Vec_t x0(nVec, sim_par.sh_order); // coordinates
  Vec_t f0(nVec, sim_par.sh_order); // force single layer
  Vec_t f1(nVec, sim_par.sh_order); // force double layer
  { // Set coordinate values
    int imax(x0.getGridDim().first);
    int jmax(x0.getGridDim().second);

    std::vector<Real_t> qx;
    { // compute legendre node points
      size_t p=x0.getShOrder();
      qx.resize(p+1);
      std::vector<Real_t> qw(p+1);
      cgqf(p+1, 1, 0.0, 0.0, -1.0, 1.0, &qx[0], &qw[0]);
    }

    int fLen = x0.getStride();
    for(size_t k=0;k<nVec;k++){
      Real_t* x_k=x0.getSubN_begin(k)+0*fLen;
      Real_t* y_k=x0.getSubN_begin(k)+1*fLen;
      Real_t* z_k=x0.getSubN_begin(k)+2*fLen;

      Real_t* f0x_k=f0.getSubN_begin(k)+0*fLen;
      Real_t* f0y_k=f0.getSubN_begin(k)+1*fLen;
      Real_t* f0z_k=f0.getSubN_begin(k)+2*fLen;

      Real_t* f1x_k=f1.getSubN_begin(k)+0*fLen;
      Real_t* f1y_k=f1.getSubN_begin(k)+1*fLen;
      Real_t* f1z_k=f1.getSubN_begin(k)+2*fLen;
      for(size_t i=0;i<imax;i++){
        Real_t cos_t=qx[i];
        //Real_t cos_t=cos((i+1)*M_PI/(imax+1));
        Real_t sin_t=sqrt(1.0-cos_t*cos_t);
        for(size_t j=0;j<jmax;j++){
          { // Set x, y, z
            x_k[j+i*jmax]=-cos_t;
            y_k[j+i*jmax]=sin_t*sin(j*2*M_PI/jmax);
            z_k[j+i*jmax]=sin_t*cos(j*2*M_PI/jmax);
          }
          if(!rank && !k){ // Set fx, fy, fz
            Real_t c[3]={x_k[j+i*jmax],
                         y_k[j+i*jmax],
                         z_k[j+i*jmax]};

            Real_t f0[3]={0,0,0};
            force(c,1,f0);
            f0x_k[j+i*jmax]=f0[0];
            f0y_k[j+i*jmax]=f0[1];
            f0z_k[j+i*jmax]=f0[2];

            Real_t f1[3]={1,1,1};
            f1x_k[j+i*jmax]=f1[0];
            f1y_k[j+i*jmax]=f1[1];
            f1z_k[j+i*jmax]=f1[2];
          }else{
            f0x_k[j+i*jmax]=0;
            f0y_k[j+i*jmax]=0;
            f0z_k[j+i*jmax]=0;

            f1x_k[j+i*jmax]=0;
            f1y_k[j+i*jmax]=0;
            f1z_k[j+i*jmax]=0;
          }
        }
      }
    }
  }

  // Creating surface objects
  Surf_t S(sim_par.sh_order,mats,&x0);

  Vec_t vel_surf(nVec, sim_par.sh_order);
  { // Set analytical surface velocity
    int imax(x0.getGridDim().first);
    int jmax(x0.getGridDim().second);
    int fLen = x0.getStride();
    for(size_t k=0;k<nVec;k++){
      Real_t* x_k=x0.getSubN_begin(k)+0*fLen;
      Real_t* y_k=x0.getSubN_begin(k)+1*fLen;
      Real_t* z_k=x0.getSubN_begin(k)+2*fLen;

      Real_t* vx_k=vel_surf.getSubN_begin(k)+0*fLen;
      Real_t* vy_k=vel_surf.getSubN_begin(k)+1*fLen;
      Real_t* vz_k=vel_surf.getSubN_begin(k)+2*fLen;
      for(size_t i=0;i<imax;i++){
        for(size_t j=0;j<jmax;j++){
          if(!rank && !k){
            Real_t c[3]={x_k[j+i*jmax],
                         y_k[j+i*jmax],
                         z_k[j+i*jmax]};
            Real_t vel[3]={0,0,0};
            u_ref(c,1,vel);
            vx_k[j+i*jmax]=vel[0];
            vy_k[j+i*jmax]=vel[1];
            vz_k[j+i*jmax]=vel[2];
          }else{
            vx_k[j+i*jmax]=0;
            vy_k[j+i*jmax]=0;
            vz_k[j+i*jmax]=0;
          }
        }
      }
    }
  }

  // Set target points
  PVFMMVec_t trg_coord; Real_t R0=1.0, R1=1.6;
  if(1){ // Set trg coord (random points)
    size_t N=10000;
    trg_coord.ReInit(N*COORD_DIM);
    trg_coord.SetZero();
    for(size_t i=0;i<N;i++){
      Real_t r=0;
      Real_t* c=&trg_coord[i*COORD_DIM+0];
      while(r<0.5 || r>1.0){
        c[0]=(drand48()-0.5)*2.0;
        c[1]=(drand48()-0.5)*2.0;
        c[2]=(drand48()-0.5)*2.0;
        r=sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);
      }
      Real_t R=R0+i*(R1-R0)/N;
      c[0]=c[0]*R/r;
      c[1]=c[1]*R/r;
      c[2]=c[2]*R/r;
    }
  }else{ // Surface mesh
    int imax(x0.getGridDim().first);
    int jmax(x0.getGridDim().second);
    size_t N=imax*jmax;
    trg_coord.ReInit(N*COORD_DIM);

    std::vector<Real_t> qx;
    { // compute legendre node points
      size_t p=x0.getShOrder();
      qx.resize(p+1);
      std::vector<Real_t> qw(p+1);
      cgqf(p+1, 1, 0.0, 0.0, -1.0, 1.0, &qx[0], &qw[0]);
    }

    for(size_t i=0;i<imax;i++){
      Real_t cos_t=qx[i];
      Real_t sin_t=sqrt(1.0-cos_t*cos_t);
      for(size_t j=0;j<jmax;j++){
        Real_t R=R0+(j+i*jmax)*(R1-R0)/N;
        if((j+i*jmax)*(R1-R0)/N<0.05) R=R0;
        trg_coord[(j+i*jmax)*COORD_DIM+0]=R*cos_t;
        trg_coord[(j+i*jmax)*COORD_DIM+1]=R*sin_t*sin(j*2*M_PI/jmax);
        trg_coord[(j+i*jmax)*COORD_DIM+2]=R*sin_t*cos(j*2*M_PI/jmax);
      }
    }
  }

  stokes_vel.SetSrcCoord(S);
  //stokes_vel.SetSurfaceVel(vel_surf); // optional
  stokes_vel.SetDensitySL(&f0);
  stokes_vel.SetDensityDL(&f1);

  stokes_vel.SetTrgCoord(&trg_coord[0],trg_coord.Dim()/COORD_DIM);
  //stokes_vel.SetTrgCoord(S);

  pvfmm::Profile::Tic("StokesVelocity",&comm);
  Real_t* v=stokes_vel();
  { // Compute error
    size_t N=trg_coord.Dim()/COORD_DIM;
    pvfmm::Matrix<Real_t> trg_vel(N,COORD_DIM,v,false);
    pvfmm::Matrix<Real_t> ref_vel(N,COORD_DIM);
    u_ref(&trg_coord[0], N, ref_vel[0]);

    size_t M=15;
    PVFMMVec_t max_err(M); max_err.SetZero();
    PVFMMVec_t max_vel(M); max_vel.SetZero();
    pvfmm::Matrix<Real_t> err=trg_vel-ref_vel;
    for(size_t j=0;j<M;j++){
      size_t a=((j+0)*err.Dim(0)*err.Dim(1))/M;
      size_t b=((j+1)*err.Dim(0)*err.Dim(1))/M;
      for(size_t i=a;i<b;i++){
        if(fabs(err    [0][i])>max_err[j]) max_err[j]=fabs(err    [0][i]);
        if(fabs(ref_vel[0][i])>max_vel[j]) max_vel[j]=fabs(ref_vel[0][i]);
      }
      Real_t R=R0+(a+b)*(R1-R0)/N/6.0;
      max_err[j]=log(max_err[j]/max_vel[j])/log(10.0);
      max_vel[j]=R;
    }
    std::cout<<max_err;
    std::cout<<max_vel;
  }
  pvfmm::Profile::Toc();
}

template<typename Surf_t>
void WriteVTK(const Surf_t& S, const char* fname, MPI_Comm comm=VES3D_COMM_WORLD, const typename Surf_t::Vec_t* v_ptr=NULL, int order=-1){
  typedef typename Surf_t::value_type Real_t;
  typedef typename Surf_t::Vec_t Vec_t;
  typedef double VTKReal_t;
  int data__dof=COORD_DIM;
  size_t p0=S.getShOrder();
  size_t p1=(order>0?order:p0); // upsample

  pvfmm::Vector<Real_t> X, Xp, V, Vp;
  { // Upsample X
    const Vec_t& x=S.getPosition();
    pvfmm::Vector<Real_t> X0(x.size(),(Real_t*)x.begin(),false);
    pvfmm::Vector<Real_t> X1;
    SphericalHarmonics<Real_t>::Grid2SHC(X0,p0,p0,X1);
    SphericalHarmonics<Real_t>::SHC2Grid(X1,p0,p1,X);
    SphericalHarmonics<Real_t>::SHC2Pole(X1, p0, Xp);
  }
  if(v_ptr){ // Upsample V
    pvfmm::Vector<Real_t> X0(v_ptr->size(),(Real_t*)v_ptr->begin(),false);
    pvfmm::Vector<Real_t> X1;
    SphericalHarmonics<Real_t>::Grid2SHC(X0,p0,p0,X1);
    SphericalHarmonics<Real_t>::SHC2Grid(X1,p0,p1,V);
    SphericalHarmonics<Real_t>::SHC2Pole(X1, p0, Vp);
  }

  std::vector<VTKReal_t> point_coord;
  std::vector<VTKReal_t> point_value;
  std::vector< int32_t> poly_connect;
  std::vector< int32_t> poly_offset;
  { // Set point_coord, point_value, poly_connect
    size_t N_ves = X.Dim()/(2*p1*(p1+1)*COORD_DIM); // Number of vesicles
    assert(Xp.Dim() == N_ves*(2*p1*(p1+1)*COORD_DIM));
    for(size_t k=0;k<N_ves;k++){ // Set point_coord
      for(size_t i=0;i<p1+1;i++){
        for(size_t j=0;j<2*p1;j++){
          for(size_t l=0;l<COORD_DIM;l++){
            point_coord.push_back(X[j+2*p1*(i+(p1+1)*(l+k*COORD_DIM))]);
          }
        }
      }
      for(size_t l=0;l<COORD_DIM;l++) point_coord.push_back(Xp[0+2*(l+k*COORD_DIM)]);
      for(size_t l=0;l<COORD_DIM;l++) point_coord.push_back(Xp[1+2*(l+k*COORD_DIM)]);
    }

    if(v_ptr)
    for(size_t k=0;k<N_ves;k++){ // Set point_value
      for(size_t i=0;i<p1+1;i++){
        for(size_t j=0;j<2*p1;j++){
          for(size_t l=0;l<data__dof;l++){
            point_value.push_back(V[j+2*p1*(i+(p1+1)*(l+k*data__dof))]);
          }
        }
      }
      for(size_t l=0;l<data__dof;l++) point_value.push_back(Vp[0+2*(l+k*data__dof)]);
      for(size_t l=0;l<data__dof;l++) point_value.push_back(Vp[1+2*(l+k*data__dof)]);
    }

    for(size_t k=0;k<N_ves;k++){
      for(size_t j=0;j<2*p1;j++){
        size_t i0= 0;
        size_t i1=p1;
        size_t j0=((j+0)       );
        size_t j1=((j+1)%(2*p1));

        poly_connect.push_back((2*p1*(p1+1)+2)*k + 2*p1*(p1+1)+0);
        poly_connect.push_back((2*p1*(p1+1)+2)*k + 2*p1*i0+j0);
        poly_connect.push_back((2*p1*(p1+1)+2)*k + 2*p1*i0+j1);
        poly_offset.push_back(poly_connect.size());

        poly_connect.push_back((2*p1*(p1+1)+2)*k + 2*p1*(p1+1)+1);
        poly_connect.push_back((2*p1*(p1+1)+2)*k + 2*p1*i1+j0);
        poly_connect.push_back((2*p1*(p1+1)+2)*k + 2*p1*i1+j1);
        poly_offset.push_back(poly_connect.size());
      }
      for(size_t i=0;i<p1;i++){
        for(size_t j=0;j<2*p1;j++){
          size_t i0=((i+0)       );
          size_t i1=((i+1)       );
          size_t j0=((j+0)       );
          size_t j1=((j+1)%(2*p1));
          poly_connect.push_back((2*p1*(p1+1)+2)*k + 2*p1*i0+j0);
          poly_connect.push_back((2*p1*(p1+1)+2)*k + 2*p1*i1+j0);
          poly_connect.push_back((2*p1*(p1+1)+2)*k + 2*p1*i1+j1);
          poly_connect.push_back((2*p1*(p1+1)+2)*k + 2*p1*i0+j1);
          poly_offset.push_back(poly_connect.size());
        }
      }
    }
  }

  int myrank, np;
  MPI_Comm_size(comm,&np);
  MPI_Comm_rank(comm,&myrank);

  std::vector<VTKReal_t>& coord=point_coord;
  std::vector<VTKReal_t>& value=point_value;
  std::vector<int32_t> connect=poly_connect;
  std::vector<int32_t> offset=poly_offset;

  int pt_cnt=coord.size()/COORD_DIM;
  int poly_cnt=poly_offset.size();

  //Open file for writing.
  std::stringstream vtufname;
  vtufname<<fname<<"_"<<std::setfill('0')<<std::setw(6)<<myrank<<".vtp";
  std::ofstream vtufile;
  vtufile.open(vtufname.str().c_str());
  if(vtufile.fail()) return;

  bool isLittleEndian;
  {
    uint16_t number = 0x1;
    uint8_t *numPtr = (uint8_t*)&number;
    isLittleEndian=(numPtr[0] == 1);
  }

  //Proceed to write to file.
  size_t data_size=0;
  vtufile<<"<?xml version=\"1.0\"?>\n";
  if(isLittleEndian) vtufile<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  else               vtufile<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"BigEndian\">\n";
  //===========================================================================
  vtufile<<"  <PolyData>\n";
  vtufile<<"    <Piece NumberOfPoints=\""<<pt_cnt<<"\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\""<<poly_cnt<<"\">\n";

  //---------------------------------------------------------------------------
  vtufile<<"      <Points>\n";
  vtufile<<"        <DataArray type=\"Float"<<sizeof(VTKReal_t)*8<<"\" NumberOfComponents=\""<<COORD_DIM<<"\" Name=\"Position\" format=\"appended\" offset=\""<<data_size<<"\" />\n";
  data_size+=sizeof(uint32_t)+coord.size()*sizeof(VTKReal_t);
  vtufile<<"      </Points>\n";
  //---------------------------------------------------------------------------
  if(value.size()){ // value
    vtufile<<"      <PointData>\n";
    vtufile<<"        <DataArray type=\"Float"<<sizeof(VTKReal_t)*8<<"\" NumberOfComponents=\""<<value.size()/pt_cnt<<"\" Name=\""<<"value"<<"\" format=\"appended\" offset=\""<<data_size<<"\" />\n";
    data_size+=sizeof(uint32_t)+value.size()*sizeof(VTKReal_t);
    vtufile<<"      </PointData>\n";
  }
  //---------------------------------------------------------------------------
  vtufile<<"      <Polys>\n";
  vtufile<<"        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<<data_size<<"\" />\n";
  data_size+=sizeof(uint32_t)+connect.size()*sizeof(int32_t);
  vtufile<<"        <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<<data_size<<"\" />\n";
  data_size+=sizeof(uint32_t)+offset.size() *sizeof(int32_t);
  vtufile<<"      </Polys>\n";
  //---------------------------------------------------------------------------

  vtufile<<"    </Piece>\n";
  vtufile<<"  </PolyData>\n";
  //===========================================================================
  vtufile<<"  <AppendedData encoding=\"raw\">\n";
  vtufile<<"    _";

  int32_t block_size;
  block_size=coord.size()*sizeof(VTKReal_t); vtufile.write((char*)&block_size, sizeof(int32_t)); vtufile.write((char*)&coord  [0], coord.size()*sizeof(VTKReal_t));
  if(value.size()){ // value
    block_size=value.size()*sizeof(VTKReal_t); vtufile.write((char*)&block_size, sizeof(int32_t)); vtufile.write((char*)&value  [0], value.size()*sizeof(VTKReal_t));
  }
  block_size=connect.size()*sizeof(int32_t); vtufile.write((char*)&block_size, sizeof(int32_t)); vtufile.write((char*)&connect[0], connect.size()*sizeof(int32_t));
  block_size=offset .size()*sizeof(int32_t); vtufile.write((char*)&block_size, sizeof(int32_t)); vtufile.write((char*)&offset [0], offset .size()*sizeof(int32_t));

  vtufile<<"\n";
  vtufile<<"  </AppendedData>\n";
  //===========================================================================
  vtufile<<"</VTKFile>\n";
  vtufile.close();


  if(myrank) return;
  std::stringstream pvtufname;
  pvtufname<<fname<<".pvtp";
  std::ofstream pvtufile;
  pvtufile.open(pvtufname.str().c_str());
  if(pvtufile.fail()) return;
  pvtufile<<"<?xml version=\"1.0\"?>\n";
  pvtufile<<"<VTKFile type=\"PPolyData\">\n";
  pvtufile<<"  <PPolyData GhostLevel=\"0\">\n";
  pvtufile<<"      <PPoints>\n";
  pvtufile<<"        <PDataArray type=\"Float"<<sizeof(VTKReal_t)*8<<"\" NumberOfComponents=\""<<COORD_DIM<<"\" Name=\"Position\"/>\n";
  pvtufile<<"      </PPoints>\n";
  if(value.size()){ // value
    pvtufile<<"      <PPointData>\n";
    pvtufile<<"        <PDataArray type=\"Float"<<sizeof(VTKReal_t)*8<<"\" NumberOfComponents=\""<<value.size()/pt_cnt<<"\" Name=\""<<"value"<<"\"/>\n";
    pvtufile<<"      </PPointData>\n";
  }
  {
    // Extract filename from path.
    std::stringstream vtupath;
    vtupath<<'/'<<fname;
    std::string pathname = vtupath.str();
    unsigned found = pathname.find_last_of("/\\");
    std::string fname_ = pathname.substr(found+1);
    for(int i=0;i<np;i++) pvtufile<<"      <Piece Source=\""<<fname_<<"_"<<std::setfill('0')<<std::setw(6)<<i<<".vtp\"/>\n";
  }
  pvtufile<<"  </PPolyData>\n";
  pvtufile<<"</VTKFile>\n";
  pvtufile.close();
}
