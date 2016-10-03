#include <legendre_rule.hpp>

template <class Real>
void SphericalHarmonics<Real>::SHC2Grid(const pvfmm::Vector<Real>& S, long p0, long p1, pvfmm::Vector<Real>& X, pvfmm::Vector<Real>* X_theta, pvfmm::Vector<Real>* X_phi){
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

template <class Real>
void SphericalHarmonics<Real>::Grid2SHC(const pvfmm::Vector<Real>& X, long p0, long p1, pvfmm::Vector<Real>& S){
  pvfmm::Matrix<Real> Mf =SphericalHarmonics<Real>::MatFourierInv(p0,p1);
  std::vector<pvfmm::Matrix<Real> > Ml =SphericalHarmonics<Real>::MatLegendreInv(p0,p1);
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

template <class Real>
void SphericalHarmonics<Real>::SHC2GridTranspose(const pvfmm::Vector<Real>& X, long p0, long p1, pvfmm::Vector<Real>& S){
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

template <class Real>
void SphericalHarmonics<Real>::SHC2Pole(const pvfmm::Vector<Real>& S, long p0, pvfmm::Vector<Real>& P){
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

template <class Real>
void SphericalHarmonics<Real>::RotateAll(const pvfmm::Vector<Real>& S, long p0, long dof, pvfmm::Vector<Real>& S_){
  std::vector<pvfmm::Matrix<Real> >& Mr=MatRotate(p0);
  std::vector<std::vector<long> > coeff_perm(p0+1);
  { // Set coeff_perm
    for(long n=0;n<=p0;n++) coeff_perm[n].resize(std::min(2*n+1,2*p0));
    long itr=0;
    for(long i=0;i<2*p0;i++){
      long m=(i+1)/2;
      for(long n=m;n<=p0;n++){
        coeff_perm[n][i]=itr;
        itr++;
      }
    }
  }
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

    std::vector<pvfmm::Matrix<Real> > Bi(p0+1), Bo(p0+1); // memory buffers
    for(long i=0;i<=p0;i++){ // initialize Bi, Bo
      Bi[i].ReInit(dof*p0,coeff_perm[i].size());
      Bo[i].ReInit(dof*p0,coeff_perm[i].size());
    }

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
      { // Fast rotation
        for(long k=0;k<dof*p0;k++){ // forward permutation
          for(long l=0;l<=p0;l++){
            for(long j=0;j<coeff_perm[l].size();j++){
              Bi[l][k][j]=B0[k][coeff_perm[l][j]];
            }
          }
        }
        for(long t=0;t<=p0;t++){
          for(long l=0;l<=p0;l++){ // mat-vec
            pvfmm::Matrix<Real>::GEMM(Bo[l],Bi[l],Mr[t*(p0+1)+l]);
          }
          pvfmm::Matrix<Real> Mout(dof*p0,Ncoef,&S1[(i*(p0+1)+t)*dof*p0][0],false);
          for(long k=0;k<dof*p0;k++){ // reverse permutation
            for(long l=0;l<=p0;l++){
              for(long j=0;j<coeff_perm[l].size();j++){
                Mout[k][coeff_perm[l][j]]=Bo[l][k][j];
              }
            }
          }
        }
      }
    }
  }
}

template <class Real>
void SphericalHarmonics<Real>::RotateTranspose(const pvfmm::Vector<Real>& S_, long p0, long dof, pvfmm::Vector<Real>& S){
  std::vector<pvfmm::Matrix<Real> > Mr=MatRotate(p0);
  for(long i=0;i<Mr.size();i++) Mr[i]=Mr[i].Transpose();
  std::vector<std::vector<long> > coeff_perm(p0+1);
  { // Set coeff_perm
    for(long n=0;n<=p0;n++) coeff_perm[n].resize(std::min(2*n+1,2*p0));
    long itr=0;
    for(long i=0;i<2*p0;i++){
      long m=(i+1)/2;
      for(long n=m;n<=p0;n++){
        coeff_perm[n][i]=itr;
        itr++;
      }
    }
  }
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

    std::vector<pvfmm::Matrix<Real> > Bi(p0+1), Bo(p0+1); // memory buffers
    for(long i=0;i<=p0;i++){ // initialize Bi, Bo
      Bi[i].ReInit(dof*p0,coeff_perm[i].size());
      Bo[i].ReInit(dof*p0,coeff_perm[i].size());
    }

    long a=(tid+0)*N/omp_p;
    long b=(tid+1)*N/omp_p;
    for(long i=a;i<b;i++){
      for(long t=0;t<p0+1;t++){
        long idx0=(i*(p0+1)+t)*p0*dof;
        { // Fast rotation
          pvfmm::Matrix<Real> Min(p0*dof,Ncoef, &S1[idx0][0],false);
          for(long k=0;k<dof*p0;k++){ // forward permutation
            for(long l=0;l<=p0;l++){
              for(long j=0;j<coeff_perm[l].size();j++){
                Bi[l][k][j]=Min[k][coeff_perm[l][j]];
              }
            }
          }
          for(long l=0;l<=p0;l++){ // mat-vec
            pvfmm::Matrix<Real>::GEMM(Bo[l],Bi[l],Mr[t*(p0+1)+l]);
          }
          for(long k=0;k<dof*p0;k++){ // reverse permutation
            for(long l=0;l<=p0;l++){
              for(long j=0;j<coeff_perm[l].size();j++){
                B0[k][coeff_perm[l][j]]=Bo[l][k][j];
              }
            }
          }
        }
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

template <class Real>
pvfmm::Vector<Real>& SphericalHarmonics<Real>::LegendreNodes(long p1){
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

template <class Real>
pvfmm::Vector<Real>& SphericalHarmonics<Real>::LegendreWeights(long p1){
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

template <class Real>
pvfmm::Vector<Real>& SphericalHarmonics<Real>::SingularWeights(long p1){
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

template <class Real>
pvfmm::Matrix<Real>& SphericalHarmonics<Real>::MatFourier(long p0, long p1){
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

template <class Real>
pvfmm::Matrix<Real>& SphericalHarmonics<Real>::MatFourierInv(long p0, long p1){
  assert(p0<SHMAXDEG && p1<SHMAXDEG);
  matrix.Mfinv_ .resize(SHMAXDEG*SHMAXDEG);
  pvfmm::Matrix<Real>& Mf =matrix.Mfinv_ [p0*SHMAXDEG+p1];
  if(!Mf.Dim(0)){
    const Real INVSQRT2PI=1.0/sqrt(2*M_PI)/p0;
    { // Set Mf
      pvfmm::Matrix<Real> M(2*p0,2*p1);
      M.SetZero();
      if(p1>p0) p1=p0;
      for(long j=0;j<2*p0;j++){
        M[j][0]=INVSQRT2PI*0.5;
        for(long k=1;k<p1;k++){
          M[j][2*k-1]=INVSQRT2PI*cos(j*k*M_PI/p0);
          M[j][2*k-0]=INVSQRT2PI*sin(j*k*M_PI/p0);
        }
        M[j][2*p1-1]=INVSQRT2PI*cos(j*p1*M_PI/p0);
      }
      if(p1==p0) for(long j=0;j<2*p0;j++) M[j][2*p1-1]*=0.5;
      Mf=M;
    }
  }
  return Mf;
}

template <class Real>
pvfmm::Matrix<Real>& SphericalHarmonics<Real>::MatFourierGrad(long p0, long p1){
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

template <class Real>
std::vector<pvfmm::Matrix<Real> >& SphericalHarmonics<Real>::MatLegendre(long p0, long p1){
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

template <class Real>
std::vector<pvfmm::Matrix<Real> >& SphericalHarmonics<Real>::MatLegendreInv(long p0, long p1){
  assert(p0<SHMAXDEG && p1<SHMAXDEG);
  matrix.Mlinv_ .resize(SHMAXDEG*SHMAXDEG);
  std::vector<pvfmm::Matrix<Real> >& Ml =matrix.Mlinv_ [p0*SHMAXDEG+p1];
  if(!Ml.size()){
    std::vector<Real> qx1(p0+1);
    std::vector<Real> qw1(p0+1);
    cgqf(p0+1, 1, 0.0, 0.0, -1.0, 1.0, &qx1[0], &qw1[0]);

    { // Set Ml
      std::vector<Real> alp(qx1.size()*(p1+1)*(p1+2)/2);
      LegPoly(&alp[0], &qx1[0], qx1.size(), p1);

      Ml.resize(p1+1);
      Real* ptr=&alp[0];
      for(long i=0;i<=p1;i++){
        Ml[i].ReInit(qx1.size(), p1+1-i);
        pvfmm::Matrix<Real> M(p1+1-i, qx1.size(), ptr, false);
        for(long j=0;j<p1+1-i;j++){ // Transpose and weights
          for(long k=0;k<qx1.size();k++){
            Ml[i][k][j]=M[j][k]*qw1[k]*2*M_PI;
          }
        }
        ptr+=Ml[i].Dim(0)*Ml[i].Dim(1);
      }
    }
  }
  return Ml;
}

template <class Real>
std::vector<pvfmm::Matrix<Real> >& SphericalHarmonics<Real>::MatLegendreGrad(long p0, long p1){
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

template <class Real>
std::vector<pvfmm::Matrix<Real> >& SphericalHarmonics<Real>::MatRotate(long p0){
  std::vector<std::vector<long> > coeff_perm(p0+1);
  { // Set coeff_perm
    for(long n=0;n<=p0;n++) coeff_perm[n].resize(std::min(2*n+1,2*p0));
    long itr=0;
    for(long i=0;i<2*p0;i++){
      long m=(i+1)/2;
      for(long n=m;n<=p0;n++){
        coeff_perm[n][i]=itr;
        itr++;
      }
    }
  }

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
      for(long n=0;n<=p0;n++){ // Create matrices for fast rotation
        pvfmm::Matrix<Real> M(coeff_perm[n].size(),coeff_perm[n].size());
        for(long i=0;i<coeff_perm[n].size();i++){
          for(long j=0;j<coeff_perm[n].size();j++){
            M[i][j]=Mcoef2coef[coeff_perm[n][i]][coeff_perm[n][j]];
          }
        }
        Mr.push_back(M);
      }
    }
  }
  return Mr;
}

template <class Real>
void SphericalHarmonics<Real>::StokesSingularInteg(const pvfmm::Vector<Real>& S, long p0, long p1, pvfmm::Vector<Real>* SLMatrix, pvfmm::Vector<Real>* DLMatrix){
  long Ngrid=2*p0*(p0+1);
  long Ncoef=  p0*(p0+2);
  long Nves=S.Dim()/(Ngrid*COORD_DIM);
  if(SLMatrix) SLMatrix->ReInit(Nves*(Ncoef*COORD_DIM)*(Ncoef*COORD_DIM));
  if(DLMatrix) DLMatrix->ReInit(Nves*(Ncoef*COORD_DIM)*(Ncoef*COORD_DIM));

  long BLOCK_SIZE=6e9/((3*2*p1*(p1+1))*(3*2*p0*(p0+1))*2*8); // Limit memory usage to 6GB
  BLOCK_SIZE=std::min<long>(BLOCK_SIZE,omp_get_max_threads());
  BLOCK_SIZE=std::max<long>(BLOCK_SIZE,1);

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

template <class Real>
void SphericalHarmonics<Real>::LegPoly(Real* poly_val, const Real* X, long N, long degree){
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

template <class Real>
void SphericalHarmonics<Real>::LegPolyDeriv(Real* poly_val, const Real* X, long N, long degree){
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

template <class Real>
template <bool SLayer, bool DLayer>
void SphericalHarmonics<Real>::StokesSingularInteg_(const pvfmm::Vector<Real>& X0, long p0, long p1, pvfmm::Vector<Real>& SL, pvfmm::Vector<Real>& DL){

  pvfmm::Profile::Tic("Rotate");
  static pvfmm::Vector<Real> S0, S;
  SphericalHarmonics<Real>::Grid2SHC(X0, p0, p0, S0);
  SphericalHarmonics<Real>::RotateAll(S0, p0, COORD_DIM, S);
  pvfmm::Profile::Toc();


  pvfmm::Profile::Tic("Upsample");
  pvfmm::Vector<Real> X, X_phi, X_theta, trg;
  SphericalHarmonics<Real>::SHC2Grid(S, p0, p1, X, &X_theta, &X_phi);
  SphericalHarmonics<Real>::SHC2Pole(S, p0, trg);
  pvfmm::Profile::Toc();


  pvfmm::Profile::Tic("Stokes");
  pvfmm::Vector<Real> SL0, DL0;
  { // Stokes kernel
    long M0=2*p0*(p0+1);
    long M1=2*p1*(p1+1);
    long N=trg.Dim()/(2*COORD_DIM);
    assert(X.Dim()==M1*COORD_DIM*N);
    if(SLayer && SL0.Dim()!=N*2*6*M1) SL0.ReInit(2*N*6*M1);
    if(DLayer && DL0.Dim()!=N*2*6*M1) DL0.ReInit(2*N*6*M1);
    pvfmm::Vector<Real>& qw=SphericalHarmonics<Real>::SingularWeights(p1);

    const Real scal_const_dl = 3.0/(4.0*M_PI);
    const Real scal_const_sl = 1.0/(8.0*M_PI);
    static Real eps=-1;
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
                DL0[((i*2+t)*6+0)*M1+s]=dx*dx*r_dot_n_rinv5;
                DL0[((i*2+t)*6+1)*M1+s]=dx*dy*r_dot_n_rinv5;
                DL0[((i*2+t)*6+2)*M1+s]=dx*dz*r_dot_n_rinv5;
                DL0[((i*2+t)*6+3)*M1+s]=dy*dy*r_dot_n_rinv5;
                DL0[((i*2+t)*6+4)*M1+s]=dy*dz*r_dot_n_rinv5;
                DL0[((i*2+t)*6+5)*M1+s]=dz*dz*r_dot_n_rinv5;
              }
              if(SLayer){
                Real area_rinv =scal_const_sl*qw[j0*t+(p1-j0)*(1-t)] * area_elem*rinv;
                Real area_rinv2=area_rinv*rinv2;
                SL0[((i*2+t)*6+0)*M1+s]=area_rinv+dx*dx*area_rinv2;
                SL0[((i*2+t)*6+1)*M1+s]=          dx*dy*area_rinv2;
                SL0[((i*2+t)*6+2)*M1+s]=          dx*dz*area_rinv2;
                SL0[((i*2+t)*6+3)*M1+s]=area_rinv+dy*dy*area_rinv2;
                SL0[((i*2+t)*6+4)*M1+s]=          dy*dz*area_rinv2;
                SL0[((i*2+t)*6+5)*M1+s]=area_rinv+dz*dz*area_rinv2;
              }
            }
          }
        }
      }
    }
    pvfmm::Profile::Add_FLOP(20*(2*p1)*(p1+1)*2*N);
    if(SLayer) pvfmm::Profile::Add_FLOP((19+6)*(2*p1)*(p1+1)*2*N);
    if(DLayer) pvfmm::Profile::Add_FLOP( 22   *(2*p1)*(p1+1)*2*N);
  }
  pvfmm::Profile::Toc();


  pvfmm::Profile::Tic("UpsampleTranspose");
  static pvfmm::Vector<Real> SL1, DL1;
  SphericalHarmonics<Real>::SHC2GridTranspose(SL0, p1, p0, SL1);
  SphericalHarmonics<Real>::SHC2GridTranspose(DL0, p1, p0, DL1);
  pvfmm::Profile::Toc();


  pvfmm::Profile::Tic("RotateTranspose");
  static pvfmm::Vector<Real> SL2, DL2;
  SphericalHarmonics<Real>::RotateTranspose(SL1, p0, 2*6, SL2);
  SphericalHarmonics<Real>::RotateTranspose(DL1, p0, 2*6, DL2);
  pvfmm::Profile::Toc();


  pvfmm::Profile::Tic("Rearrange");
  static pvfmm::Vector<Real> SL3, DL3;
  { // Transpose
    long Ncoef=p0*(p0+2);
    long Ngrid=2*p0*(p0+1);
    { // Transpose SL2
      long N=SL2.Dim()/(6*Ncoef*Ngrid);
      SL3.ReInit(N*COORD_DIM*Ncoef*COORD_DIM*Ngrid);
      #pragma omp parallel
      {
        long tid=omp_get_thread_num();
        long omp_p=omp_get_num_threads();
        pvfmm::Matrix<Real> B(COORD_DIM*Ncoef,Ngrid*COORD_DIM);

        long a=(tid+0)*N/omp_p;
        long b=(tid+1)*N/omp_p;
        for(long i=a;i<b;i++){
          pvfmm::Matrix<Real> M0(Ngrid*6, Ncoef, &SL2[i*Ngrid*6*Ncoef], false);
          for(long k=0;k<Ncoef;k++){ // Transpose
            for(long j=0;j<Ngrid;j++){ // TODO: needs blocking
              B[k+Ncoef*0][j*COORD_DIM+0]=M0[j*6+0][k];
              B[k+Ncoef*1][j*COORD_DIM+0]=M0[j*6+1][k];
              B[k+Ncoef*2][j*COORD_DIM+0]=M0[j*6+2][k];
              B[k+Ncoef*0][j*COORD_DIM+1]=M0[j*6+1][k];
              B[k+Ncoef*1][j*COORD_DIM+1]=M0[j*6+3][k];
              B[k+Ncoef*2][j*COORD_DIM+1]=M0[j*6+4][k];
              B[k+Ncoef*0][j*COORD_DIM+2]=M0[j*6+2][k];
              B[k+Ncoef*1][j*COORD_DIM+2]=M0[j*6+4][k];
              B[k+Ncoef*2][j*COORD_DIM+2]=M0[j*6+5][k];
            }
          }
          pvfmm::Matrix<Real> M1(Ncoef*COORD_DIM, COORD_DIM*Ngrid, &SL3[i*COORD_DIM*Ncoef*COORD_DIM*Ngrid], false);
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
      long N=DL2.Dim()/(6*Ncoef*Ngrid);
      DL3.ReInit(N*COORD_DIM*Ncoef*COORD_DIM*Ngrid);
      #pragma omp parallel
      {
        long tid=omp_get_thread_num();
        long omp_p=omp_get_num_threads();
        pvfmm::Matrix<Real> B(COORD_DIM*Ncoef,Ngrid*COORD_DIM);

        long a=(tid+0)*N/omp_p;
        long b=(tid+1)*N/omp_p;
        for(long i=a;i<b;i++){
          pvfmm::Matrix<Real> M0(Ngrid*6, Ncoef, &DL2[i*Ngrid*6*Ncoef], false);
          for(long k=0;k<Ncoef;k++){ // Transpose
            for(long j=0;j<Ngrid;j++){ // TODO: needs blocking
              B[k+Ncoef*0][j*COORD_DIM+0]=M0[j*6+0][k];
              B[k+Ncoef*1][j*COORD_DIM+0]=M0[j*6+1][k];
              B[k+Ncoef*2][j*COORD_DIM+0]=M0[j*6+2][k];
              B[k+Ncoef*0][j*COORD_DIM+1]=M0[j*6+1][k];
              B[k+Ncoef*1][j*COORD_DIM+1]=M0[j*6+3][k];
              B[k+Ncoef*2][j*COORD_DIM+1]=M0[j*6+4][k];
              B[k+Ncoef*0][j*COORD_DIM+2]=M0[j*6+2][k];
              B[k+Ncoef*1][j*COORD_DIM+2]=M0[j*6+4][k];
              B[k+Ncoef*2][j*COORD_DIM+2]=M0[j*6+5][k];
            }
          }
          pvfmm::Matrix<Real> M1(Ncoef*COORD_DIM, COORD_DIM*Ngrid, &DL3[i*COORD_DIM*Ncoef*COORD_DIM*Ngrid], false);
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
  SphericalHarmonics<Real>::Grid2SHC(SL3, p0, p0, SL);
  SphericalHarmonics<Real>::Grid2SHC(DL3, p0, p0, DL);
  pvfmm::Profile::Toc();

}

