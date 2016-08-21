#include <omp.h>
#include <iostream>
#include <profile.hpp>
#include <SphericalHarmonics.h>
#include <legendre_rule.h>

#define __ENABLE_PVFMM_PROFILER__
//#define __SH_FILTER__

template<class Real>
StokesVelocity<Real>::StokesVelocity(int sh_order_, int sh_order_up_, Real box_size_, Real repul_dist_, MPI_Comm comm_):
  sh_order(sh_order_), sh_order_up_self(sh_order_up_), sh_order_up(sh_order_up_), box_size(box_size_), comm(comm_), trg_is_surf(true), near_singular0(box_size_, repul_dist_, comm_), near_singular1(box_size_, repul_dist_, comm_)
{
  pvfmm_ctx=PVFMMCreateContext<Real>(box_size_);
  fmm_setup=true;
}

template <class Real>
StokesVelocity<Real>::~StokesVelocity(){
  PVFMMDestroyContext<Real>(&pvfmm_ctx);
}


template <class Real>
void StokesVelocity<Real>::SetSrcCoord(const PVFMMVec& S){
  scoord.ReInit(S.Dim(), (Real*)&S[0], true);
  { // filter
    #ifdef __SH_FILTER__
    pvfmm::Vector<Real>& V=scoord;
    pvfmm::Vector<Real> tmp;
    SphericalHarmonics<Real>::Grid2SHC(V,sh_order,sh_order,tmp);
    SphericalHarmonics<Real>::SHC2Grid(tmp,sh_order,sh_order,V);
    #endif
  }
  fmm_setup=true;

  SLMatrix.ReInit(0);
  DLMatrix.ReInit(0);

  scoord_far.ReInit(0);
  tcoord_repl.ReInit(0);
  scoord_norm.ReInit(0);
  scoord_area.ReInit(0);

  rforce_single.ReInit(0);
  qforce_single.ReInit(0);
  uforce_double.ReInit(0);
  qforce_double.ReInit(0);

  S_vel.ReInit(0);
  S_vel_up.ReInit(0);
  fmm_vel.ReInit(0);
  trg_vel.ReInit(0);
}

template <class Real>
template <class Vec>
void StokesVelocity<Real>::SetSrcCoord(const Vec& S){
  PVFMMVec tmp(S.size(), (Real*)S.begin(), false);
  SetSrcCoord(tmp);
}

template <class Real>
void StokesVelocity<Real>::SetDensitySL(const PVFMMVec* f, bool add_repul_){
  if(f){
    if(force_single.Dim()!=f->Dim()) fmm_setup=true;
    force_single.ReInit(f->Dim(), (Real*)&f[0][0], true);
  }else if(force_single.Dim()!=0){
    fmm_setup=true;
    force_single.ReInit(0);
    near_singular0.SetDensitySL(NULL);
    near_singular1.SetDensitySL(NULL);
  }

  if(force_single.Dim()){ // filter
    #ifdef __SH_FILTER__
    pvfmm::Vector<Real>& V=force_single;
    pvfmm::Vector<Real> tmp;
    SphericalHarmonics<Real>::Grid2SHC(V,sh_order,sh_order,tmp);
    SphericalHarmonics<Real>::SHC2Grid(tmp,sh_order,sh_order,V);
    #endif
  }

  rforce_single.ReInit(0);
  qforce_single.ReInit(0);
  add_repul=add_repul_;

  S_vel.ReInit(0);
  S_vel_up.ReInit(0);
  fmm_vel.ReInit(0);
  trg_vel.ReInit(0);
}

template <class Real>
template <class Vec>
void StokesVelocity<Real>::SetDensitySL(const Vec* f, bool add_repul_){
  if(f){
    PVFMMVec tmp(f->size(), (Real*)f->begin(), false);
    SetDensitySL(&tmp, add_repul_);
  }else SetDensitySL((const PVFMMVec*)NULL, add_repul_);
}

template <class Real>
void StokesVelocity<Real>::SetDensityDL(const PVFMMVec* f){
  if(f){
    if(force_double.Dim()!=f->Dim()) fmm_setup=true;
    force_double.ReInit(f->Dim(), (Real*)&f[0][0], true);
  }else if(force_double.Dim()!=0){
    fmm_setup=true;
    force_double.ReInit(0);
    near_singular0.SetDensityDL(NULL, NULL);
    near_singular1.SetDensityDL(NULL, NULL);
  }

  if(force_double.Dim()){ // filter
    #ifdef __SH_FILTER__
    pvfmm::Vector<Real>& V=force_double;
    pvfmm::Vector<Real> tmp;
    SphericalHarmonics<Real>::Grid2SHC(V,sh_order,sh_order,tmp);
    SphericalHarmonics<Real>::SHC2Grid(tmp,sh_order,sh_order,V);
    #endif
  }

  uforce_double.ReInit(0);
  qforce_double.ReInit(0);

  S_vel.ReInit(0);
  S_vel_up.ReInit(0);
  fmm_vel.ReInit(0);
  trg_vel.ReInit(0);
}

template <class Real>
template <class Vec>
void StokesVelocity<Real>::SetDensityDL(const Vec* f){
  if(f){
    PVFMMVec tmp(f->size(), (Real*)f->begin(), false);
    SetDensityDL(&tmp);
  }else SetDensityDL(NULL);
}

template <class Real>
void StokesVelocity<Real>::SetTrgCoord(const PVFMMVec* T){
  if(T){
    trg_is_surf=false;
    tcoord.ReInit(T->Dim(),&T[0][0]);
    near_singular1.SetTrgCoord(&tcoord[0],tcoord.Dim()/COORD_DIM,false);
  }else{
    trg_is_surf=true;
    tcoord.ReInit(0);
  }

  fmm_setup=true;
  S_vel.ReInit(0);
  S_vel_up.ReInit(0);
  fmm_vel.ReInit(0);
  trg_vel.ReInit(0);
}

template <class Real>
const StokesVelocity<Real>::PVFMMVec& StokesVelocity<Real>::operator()(){
#ifdef __ENABLE_PVFMM_PROFILER__
  bool prof_state=pvfmm::Profile::Enable(true);
  pvfmm::Profile::Tic("StokesVelocity",&comm, true);
#endif

  { // Setup
    pvfmm::Profile::Tic("Setup",&comm, true);
    bool prof_state=pvfmm::Profile::Enable(false);

    if(!SLMatrix.Dim() || !DLMatrix.Dim()){
      pvfmm::Profile::Tic("SelfMatrix",&comm, true);
      if(!SLMatrix.Dim() && !DLMatrix.Dim() && force_single.Dim() && force_double.Dim()){
        SphericalHarmonics<Real>::StokesSingularInteg(scoord, sh_order, sh_order_up_self, &SLMatrix, &DLMatrix);
      }else if(!SLMatrix.Dim() && force_single.Dim()){
        SphericalHarmonics<Real>::StokesSingularInteg(scoord, sh_order, sh_order_up_self, &SLMatrix, NULL);
      }else if(!DLMatrix.Dim() && force_double.Dim()){
        SphericalHarmonics<Real>::StokesSingularInteg(scoord, sh_order, sh_order_up_self, NULL, &DLMatrix);
      }
      pvfmm::Profile::Toc();
    }

    if(!scoord_far.Dim()){
      assert(!scoord_norm.Dim());
      assert(!scoord_area.Dim());
      assert(!tcoord_repl.Dim());

      pvfmm::Profile::Tic("SCoordFar",&comm, true);
      static PVFMMVec scoord_shc, scoord_up, X_theta, X_phi, scoord_pole;
      SphericalHarmonics<Real>::Grid2SHC(scoord    , sh_order, sh_order   , scoord_shc);
      SphericalHarmonics<Real>::SHC2Grid(scoord_shc, sh_order, sh_order_up, scoord_up, &X_theta, &X_phi);
      SphericalHarmonics<Real>::SHC2Pole(scoord_shc, sh_order, scoord_pole);
      { // Set scoord_far
        long Nves=scoord_pole.Dim()/COORD_DIM/2;
        long Mves=2*sh_order_up*(1+sh_order_up);
        scoord_far.ReInit(Nves*(Mves+2)*COORD_DIM);
        #pragma omp parallel for
        for(long i=0;i<Nves;i++){
          for(long k=0;k<COORD_DIM;k++){
            scoord_far[(i*(Mves+2)+0)*COORD_DIM+k]=scoord_pole[(i*COORD_DIM+k)*2+0];
            scoord_far[(i*(Mves+2)+1)*COORD_DIM+k]=scoord_pole[(i*COORD_DIM+k)*2+1];
          }
          for(long j=0;j<Mves;j++){
            for(long k=0;k<COORD_DIM;k++){
              scoord_far[(i*(Mves+2)+(j+2))*COORD_DIM+k]=scoord_up[(i*COORD_DIM+k)*Mves+j];
            }
          }
        }
      }
      pvfmm::Profile::Toc();

      pvfmm::Profile::Tic("SCoordNear",&comm, true);
      { // Set tcoord_repl
        long Nves=scoord_pole.Dim()/COORD_DIM/2;
        long Mves=2*sh_order*(1+sh_order);
        tcoord_repl.ReInit(Nves*Mves*COORD_DIM);
        #pragma omp parallel for
        for(long i=0;i<Nves;i++){
          for(long j=0;j<Mves;j++){
            for(long k=0;k<COORD_DIM;k++){
              tcoord_repl[(i*Mves+j)*COORD_DIM+k]=scoord[(i*COORD_DIM+k)*Mves+j];
            }
          }
        }
      }
      near_singular0.SetTrgCoord(&tcoord_repl[0],tcoord_repl.Dim()/COORD_DIM,true);
      near_singular0.SetSrcCoord(scoord_far,sh_order_up);
      near_singular1.SetSrcCoord(scoord_far,sh_order_up);
      pvfmm::Profile::Toc();

      pvfmm::Profile::Tic("AreaNormal",&comm, true);
      { // Set scoord_norm, scoord_area
        long Mves=2*sh_order_up*(sh_order_up+1);
        long N=X_theta.Dim()/Mves/COORD_DIM;
        scoord_norm.ReInit(N*COORD_DIM*Mves);
        scoord_area.ReInit(N*Mves);
        #pragma omp parallel for
        for(long i=0;i<N;i++){
          for(long j=0;j<Mves;j++){
            Real nx, ny, nz;
            { // Compute source normal
              Real x_theta=X_theta[(i*COORD_DIM+0)*Mves+j];
              Real y_theta=X_theta[(i*COORD_DIM+1)*Mves+j];
              Real z_theta=X_theta[(i*COORD_DIM+2)*Mves+j];

              Real x_phi=X_phi[(i*COORD_DIM+0)*Mves+j];
              Real y_phi=X_phi[(i*COORD_DIM+1)*Mves+j];
              Real z_phi=X_phi[(i*COORD_DIM+2)*Mves+j];

              nx=(y_theta*z_phi-z_theta*y_phi);
              ny=(z_theta*x_phi-x_theta*z_phi);
              nz=(x_theta*y_phi-y_theta*x_phi);
            }
            Real area=sqrt(nx*nx+ny*ny+nz*nz);
            scoord_area[i*Mves+j]=area;
            Real inv_area=1.0/area;
            scoord_norm[(i*COORD_DIM+0)*Mves+j]=nx*inv_area;
            scoord_norm[(i*COORD_DIM+1)*Mves+j]=ny*inv_area;
            scoord_norm[(i*COORD_DIM+2)*Mves+j]=nz*inv_area;
          }
        }
      }
      pvfmm::Profile::Toc();
    }

    if(!rforce_single.Dim() && add_repul){ // Add repulsion
      pvfmm::Profile::Tic("Repulsion",&comm);
      const PVFMMVec& f_repl=near_singular0.ForceRepul();
      long Mves=2*sh_order*(sh_order+1);
      long Nves=f_repl.Dim()/Mves/COORD_DIM;
      assert(f_repl.Dim()==Nves*Mves*COORD_DIM);
      rforce_single.ReInit(Nves*Mves*COORD_DIM);
      #pragma omp parallel for
      for(long i=0;i<Nves;i++){
        for(long j=0;j<COORD_DIM;j++){
          for(long k=0;k<Mves;k++){
            rforce_single[(i*COORD_DIM+j)*Mves+k]=force_single[(i*COORD_DIM+j)*Mves+k]+f_repl[(i*Mves+k)*COORD_DIM+j];
          }
        }
      }
      pvfmm::Profile::Toc();
      if(force_single.Dim()){
        assert(force_single.Dim()==Nves*Mves*COORD_DIM);
        #pragma omp parallel for
        for(long i=0;i<Nves*COORD_DIM*Mves;i++){
          rforce_single[i]+=force_single[i];
        }
      }
    }else if(!rforce_single.Dim() && force_single.Dim()){
      rforce_single.ReInit(force_single.Dim(),&force_single[0],false);
    }

    { // Compute qforce_single, qforce_double
      pvfmm::Vector<Real>& qw=SphericalHarmonics<Real>::LegendreWeights(sh_order_up);
      if(!qforce_single.Dim() && rforce_single.Dim()){ // Compute qforce_single
        static PVFMMVec shc, grid;
        SphericalHarmonics<Real>::Grid2SHC(rforce_single, sh_order, sh_order, shc);
        SphericalHarmonics<Real>::SHC2Grid(shc, sh_order, sh_order_up, grid);

        long Mves=2*sh_order_up*(sh_order_up+1);
        long Nves=grid.Dim()/Mves/COORD_DIM;
        assert(scoord_area.Dim()==Nves*Mves);

        qforce_single.ReInit(Nves*(Mves+2)*COORD_DIM);
        #pragma omp parallel for
        for(long i=0;i<Nves;i++){
          for(long k=0;k<COORD_DIM;k++){
            qforce_single[(i*(Mves+2)+0)*COORD_DIM+k]=0;
            qforce_single[(i*(Mves+2)+1)*COORD_DIM+k]=0;
          }
          for(long j0=0;j0<sh_order_up+1;j0++){
            for(long j1=0;j1<sh_order_up*2;j1++){
              long j=j0*sh_order_up*2+j1;
              Real w=scoord_area[i*Mves+j]*qw[j0];
              for(long k=0;k<COORD_DIM;k++){
                qforce_single[(i*(Mves+2)+(j+2))*COORD_DIM+k]=grid[(i*COORD_DIM+k)*Mves+j]*w;
              }
            }
          }
        }
        near_singular0.SetDensitySL(&qforce_single);
        near_singular1.SetDensitySL(&qforce_single);
      }
      if(!qforce_double.Dim() &&  force_double.Dim()){ // Compute qforce_double
        assert(!uforce_double.Dim());

        static PVFMMVec shc, grid, pole;
        SphericalHarmonics<Real>::Grid2SHC(force_double, sh_order, sh_order, shc);
        SphericalHarmonics<Real>::SHC2Grid(shc, sh_order, sh_order_up, grid);
        SphericalHarmonics<Real>::SHC2Pole(shc, sh_order, pole);

        long Mves=2*sh_order_up*(sh_order_up+1);
        long Nves=grid.Dim()/Mves/COORD_DIM;
        assert(scoord_norm.Dim()==Nves*Mves*COORD_DIM);
        assert(scoord_area.Dim()==Nves*Mves);

        uforce_double.ReInit(Nves*(Mves+2)*(1*COORD_DIM));
        qforce_double.ReInit(Nves*(Mves+2)*(2*COORD_DIM));
        #pragma omp parallel for
        for(long i=0;i<Nves;i++){
          for(long k=0;k<COORD_DIM;k++){
            uforce_double[(i*(Mves+2)+0)*COORD_DIM+k]=pole[(i*COORD_DIM+k)*2+0];
            uforce_double[(i*(Mves+2)+1)*COORD_DIM+k]=pole[(i*COORD_DIM+k)*2+1];

            qforce_double[(i*(Mves+2)+0)*2*COORD_DIM+0*COORD_DIM+k]=0;
            qforce_double[(i*(Mves+2)+0)*2*COORD_DIM+1*COORD_DIM+k]=0;
            qforce_double[(i*(Mves+2)+1)*2*COORD_DIM+0*COORD_DIM+k]=0;
            qforce_double[(i*(Mves+2)+1)*2*COORD_DIM+1*COORD_DIM+k]=0;
          }
          for(long j0=0;j0<sh_order_up+1;j0++){
            for(long j1=0;j1<sh_order_up*2;j1++){
              long j=j0*sh_order_up*2+j1;
              Real w=scoord_area[i*Mves+j]*qw[j0];
              for(long k=0;k<COORD_DIM;k++){
                uforce_double[(i*(Mves+2)+(j+2))*1*COORD_DIM+0*COORD_DIM+k]=       grid[(i*COORD_DIM+k)*Mves+j];
                qforce_double[(i*(Mves+2)+(j+2))*2*COORD_DIM+0*COORD_DIM+k]=       grid[(i*COORD_DIM+k)*Mves+j]*w;
                qforce_double[(i*(Mves+2)+(j+2))*2*COORD_DIM+1*COORD_DIM+k]=scoord_norm[(i*COORD_DIM+k)*Mves+j];
              }
            }
          }
        }
        near_singular0.SetDensityDL(&qforce_double, &uforce_double);
        near_singular1.SetDensityDL(&qforce_double, &uforce_double);
      }
    }

    pvfmm::Profile::Enable(prof_state);
    pvfmm::Profile::Toc();
  }

  if(!S_vel.Dim()){ // Compute self interaction
    pvfmm::Profile::Tic("SelfInteraction",&comm);
    bool prof_state=pvfmm::Profile::Enable(false);
    assert(!S_vel_up.Dim());
    static PVFMMVec vel_up, vel_pole, Vcoef;
    { // Compute Vcoeff
      long Ncoef =   sh_order*(sh_order+2);
      long Ngrid = 2*sh_order*(sh_order+1);
      static PVFMMVec SL_vel, DL_vel;
      SL_vel.ReInit(0);
      DL_vel.ReInit(0);

      if(rforce_single.Dim()){ // Set SL_vel
        static pvfmm::Vector<Real> F;
        SphericalHarmonics<Real>::Grid2SHC(rforce_single,sh_order,sh_order,F);

        long nv = rforce_single.Dim()/Ngrid/COORD_DIM;
        SL_vel.ReInit(nv*COORD_DIM*Ncoef);
        #pragma omp parallel
        { // mat-vec
          long tid=omp_get_thread_num();
          long omp_p=omp_get_num_threads();

          long a=(tid+0)*nv/omp_p;
          long b=(tid+1)*nv/omp_p;
          for(long i=a;i<b;i++){
            pvfmm::Matrix<Real> Mv(1,COORD_DIM*Ncoef,&SL_vel[i*COORD_DIM*Ncoef],false);
            pvfmm::Matrix<Real> Mf(1,COORD_DIM*Ncoef,&F     [i*COORD_DIM*Ncoef],false);
            pvfmm::Matrix<Real> M(COORD_DIM*Ncoef,COORD_DIM*Ncoef,&SLMatrix[i*COORD_DIM*Ncoef*COORD_DIM*Ncoef],false);
            pvfmm::Matrix<Real>::GEMM(Mv,Mf,M);
          }
        }
      }
      if(force_double.Dim()){ // Set DL_vel
        static pvfmm::Vector<Real> F;
        SphericalHarmonics<Real>::Grid2SHC(force_double,sh_order,sh_order,F);

        long nv = force_double.Dim()/Ngrid/COORD_DIM;
        DL_vel.ReInit(nv*COORD_DIM*Ncoef);
        #pragma omp parallel
        { // mat-vec
          long tid=omp_get_thread_num();
          long omp_p=omp_get_num_threads();

          long a=(tid+0)*nv/omp_p;
          long b=(tid+1)*nv/omp_p;
          for(long i=a;i<b;i++){
            pvfmm::Matrix<Real> Mv(1,COORD_DIM*Ncoef,&DL_vel[i*COORD_DIM*Ncoef],false);
            pvfmm::Matrix<Real> Mf(1,COORD_DIM*Ncoef,&F     [i*COORD_DIM*Ncoef],false);
            pvfmm::Matrix<Real> M(COORD_DIM*Ncoef,COORD_DIM*Ncoef,&DLMatrix[i*COORD_DIM*Ncoef*COORD_DIM*Ncoef],false);
            pvfmm::Matrix<Real>::GEMM(Mv,Mf,M);
          }
        }
      }
      if(SL_vel.Dim() && DL_vel.Dim()){ // Vcoef=SL_vel+DL_vel
        Vcoef.ReInit(SL_vel.Dim());
        #pragma omp parallel for
        for(long i=0;i<Vcoef.Dim();i++) Vcoef[i]=SL_vel[i]+DL_vel[i];
      }else{
        if(SL_vel.Dim()) Vcoef.ReInit(SL_vel.Dim(),&SL_vel[0]);
        else if(DL_vel.Dim()) Vcoef.ReInit(DL_vel.Dim(),&DL_vel[0]);
        else Vcoef.ReInit(0);
      }
    }
    SphericalHarmonics<Real>::SHC2Grid(Vcoef, sh_order, sh_order   , S_vel);
    SphericalHarmonics<Real>::SHC2Grid(Vcoef, sh_order, sh_order_up, vel_up);
    SphericalHarmonics<Real>::SHC2Pole(Vcoef, sh_order, vel_pole);
    { // Set S_vel_up
      long Nves=vel_pole.Dim()/COORD_DIM/2;
      long Mves=2*sh_order_up*(1+sh_order_up);
      S_vel_up.ReInit(Nves*(Mves+2)*COORD_DIM);
      #pragma omp parallel for
      for(long i=0;i<Nves;i++){
        for(long k=0;k<COORD_DIM;k++){
          S_vel_up[(i*(Mves+2)+0)*COORD_DIM+k]=vel_pole[(i*COORD_DIM+k)*2+0];
          S_vel_up[(i*(Mves+2)+1)*COORD_DIM+k]=vel_pole[(i*COORD_DIM+k)*2+1];
        }
        for(long j=0;j<Mves;j++){
          for(long k=0;k<COORD_DIM;k++){
            S_vel_up[(i*(Mves+2)+(j+2))*COORD_DIM+k]=vel_up[(i*COORD_DIM+k)*Mves+j];
          }
        }
      }
    }
    near_singular0.SetSurfaceVel(&S_vel_up);
    near_singular1.SetSurfaceVel(&S_vel_up);
    pvfmm::Profile::Enable(prof_state);
    pvfmm::Profile::Toc();
  }

  NearSingular<Real>& near_singular=(trg_is_surf?near_singular0:near_singular1);
  PVFMMVec& trg_coord=(trg_is_surf?tcoord_repl:tcoord);
  if(!fmm_vel.Dim()){ // Compute far interaction
    pvfmm::Profile::Tic("FarInteraction",&comm,true);
    bool prof_state=pvfmm::Profile::Enable(false);
    fmm_vel.ReInit(trg_coord.Dim());
    PVFMMEval(&scoord_far[0],
              (qforce_single.Dim()?&qforce_single[0]:NULL),
              (qforce_double.Dim()?&qforce_double[0]:NULL),
              scoord_far.Dim()/COORD_DIM,
              &trg_coord[0], &fmm_vel[0], trg_coord.Dim()/COORD_DIM, &pvfmm_ctx, fmm_setup);
    fmm_setup=false;
    near_singular.SubtractDirect(fmm_vel);
    pvfmm::Profile::Enable(prof_state);
    pvfmm::Profile::Toc();
  }

  if(!trg_vel.Dim()){ // Compute near interaction
    pvfmm::Profile::Tic("NearInteraction",&comm,true);
    bool prof_state=pvfmm::Profile::Enable(false);
    const PVFMMVec& near_vel=near_singular();
    { // Compute trg_vel = fmm_vel + near_vel
      trg_vel.ReInit(trg_coord.Dim());
      assert(trg_vel.Dim()==fmm_vel.Dim());
      assert(trg_vel.Dim()==near_vel.Dim());
      #pragma omp parallel for
      for(size_t i=0;i<trg_vel.Dim();i++){
        trg_vel[i]=fmm_vel[i]+near_vel[i];
      }
    }
    if(trg_is_surf){ // trg_vel+=S_vel
      long Ngrid = 2*sh_order*(sh_order+1);
      long Nves = S_vel.Dim()/Ngrid/COORD_DIM;

      static PVFMMVec tmp;
      tmp.ReInit(Nves*COORD_DIM*Ngrid);
      #pragma omp parallel for
      for(long i=0;i<Nves;i++){
        for(long j=0;j<COORD_DIM;j++){
          for(long k=0;k<Ngrid;k++){
            tmp[(i*COORD_DIM+j)*Ngrid+k]=S_vel[(i*COORD_DIM+j)*Ngrid+k]+trg_vel[(i*Ngrid+k)*COORD_DIM+j];
          }
        }
      }
      trg_vel.Swap(tmp);
      { // filter
        #ifdef __SH_FILTER__
        pvfmm::Vector<Real>& V=trg_vel;
        pvfmm::Vector<Real> tmp;
        SphericalHarmonics<Real>::Grid2SHC(V,sh_order,sh_order,tmp);
        SphericalHarmonics<Real>::SHC2Grid(tmp,sh_order,sh_order,V);
        #endif
      }
    }
    pvfmm::Profile::Enable(prof_state);
    pvfmm::Profile::Toc();
  }

#ifdef __ENABLE_PVFMM_PROFILER__
  pvfmm::Profile::Toc();
  pvfmm::Profile::Enable(prof_state);
#endif

  return trg_vel;
}

template <class Real>
template <class Vec>
void StokesVelocity<Real>::operator()(Vec& vel){
  PVFMMVec trg_vel_(vel.size(),vel.begin(),false);

  const PVFMMVec& trg_vel=this->operator()();
  assert(trg_vel.Dim()==trg_vel_.Dim());
  trg_vel_=trg_vel;
}


template <class Real>
Real StokesVelocity<Real>::MonitorError(Real tol){
  static PVFMMVec force_single_orig, force_double_orig, tcoord_orig;
  pvfmm::Profile::Tic("StokesMonitor",&comm, true);
  bool trg_is_surf_orig;
  { // Save original state
    force_single_orig=force_single;
    force_double_orig=force_double;
    trg_is_surf_orig=trg_is_surf;
    if(!trg_is_surf) tcoord_orig=tcoord;
  }

  static PVFMMVec force;
  force.ReInit(scoord.Dim());
  long Ngrid = 2*sh_order*(sh_order+1);
  long N_ves = scoord.Dim()/COORD_DIM/Ngrid;
  for(size_t i=0;i<N_ves;i++){ // Set force
    for(size_t j=0;j<COORD_DIM;j++){
      for(size_t k=0;k<Ngrid;k++){
        force[(i*COORD_DIM+j)*Ngrid+k]=1.0;
      }
    }
  }

  SetDensitySL(NULL);
  SetDensityDL(&force);
  SetTrgCoord(NULL);

  const PVFMMVec& velocity=this->operator()();
  PVFMMVec&       velocity_fmm=fmm_vel;
  const PVFMMVec& velocity_near=near_singular0();
  const PVFMMVec& velocity_self=S_vel;
  bool collision=near_singular0.CheckCollision();

  if(box_size>0){ // Subtract average from velocity_fmm
    long long N_loc=velocity_fmm.Dim()/COORD_DIM, N_glb;
    Real sum_loc[COORD_DIM], sum_glb[COORD_DIM];
    for(long i=0;i<COORD_DIM;i++) sum_loc[i]=0;
    for(long i=0;i<N_loc;i++){
      for(long j=0;j<COORD_DIM;j++){
        sum_loc[j]+=velocity_fmm[i*COORD_DIM+j];
      }
    }
    MPI_Allreduce(sum_loc, sum_glb, 3, pvfmm::par::Mpi_datatype<     Real>::value(), pvfmm::par::Mpi_datatype<     Real>::sum(), comm);
    MPI_Allreduce(& N_loc, & N_glb, 1, pvfmm::par::Mpi_datatype<long long>::value(), pvfmm::par::Mpi_datatype<long long>::sum(), comm);
    for(long i=0;i<COORD_DIM;i++) sum_glb[i]/=N_glb;
    for(long i=0;i<N_loc;i++){
      for(long j=0;j<COORD_DIM;j++){
        velocity_fmm[i*COORD_DIM+j]-=sum_glb[j];
      }
    }
  }

  if(1){ // Write VTK file
    static unsigned long iter=0;
    unsigned long skip=1;
    if(iter%skip==0){
      char fname[1000];
      sprintf(fname, "vis/test1_%05d", (int)(iter/skip));
      WriteVTK(scoord, sh_order, std::max(sh_order_up,sh_order_up_self), fname, box_size, &velocity, comm);
    }
    iter++;
  }

  double norm_glb[3]={0,0,0};
  { // Compute error norm
    double norm_local[3]={0,0,0};
    for(size_t i=0;i<N_ves*Ngrid*COORD_DIM;i++){
      norm_local[0]=std::max(norm_local[0],fabs(velocity_self[i]+0.5));
      norm_local[1]=std::max(norm_local[1],fabs(velocity_near[i]    ));
      norm_local[2]=std::max(norm_local[2],fabs(velocity_fmm [i]    ));
    }
    if(collision) norm_local[1]=1e10;
    MPI_Allreduce(&norm_local, &norm_glb, 3, pvfmm::par::Mpi_datatype<Real>::value(), pvfmm::par::Mpi_datatype<Real>::max(), comm);
  }

  { // Restore original state
    if(force_single_orig.Dim()) SetDensitySL(&force_single_orig);
    else SetDensitySL(NULL);
    if(force_double_orig.Dim()) SetDensityDL(&force_double_orig);
    else SetDensityDL(NULL);
    if(!trg_is_surf_orig) SetTrgCoord(&tcoord_orig);
  }
  pvfmm::Profile::Toc();

  INFO("StokesVelocity: sh_order = "<<sh_order_up_self<<","<<sh_order_up<<"  Double-layer integration error: "<<norm_glb[0]<<' '<<norm_glb[1]<<' '<<norm_glb[2]);

  // Update sh_order_up_self, sh_order_up
  if(norm_glb[0]>tol) sh_order_up_self++;
  else if(norm_glb[0]<tol*1e-1) sh_order_up_self--;
  if(norm_glb[0]<tol){
    Real err=std::max(norm_glb[1],norm_glb[2]);
    if(err>tol) sh_order_up++;
    else if(err<tol*1e-1) sh_order_up--;
  }

  return norm_glb[0]+norm_glb[1]+norm_glb[2];
}

template <class Real>
void StokesVelocity<Real>::Test(){
  int p0=16;
  long Ngrid=(p0+1)*2*p0;
  long Nves=2;

  StokesVelocity<Real> S(p0,2*p0);
  pvfmm::Vector<Real> X (Nves*Ngrid*COORD_DIM);
  pvfmm::Vector<Real> FS(Nves*Ngrid*COORD_DIM);
  pvfmm::Vector<Real> FD(Nves*Ngrid*COORD_DIM);
  pvfmm::Vector<Real> V (Nves*Ngrid*COORD_DIM);

  { // Set coordinate values
    std::vector<Real> qx;
    { // compute legendre node points
      qx.resize(p0+1);
      std::vector<Real> qw(p0+1);
      cgqf(p0+1, 1, 0.0, 0.0, -1.0, 1.0, &qx[0], &qw[0]);
    }
    for(long k=0;k<Nves;k++){
      for(size_t i=0;i<p0+1;i++){
        Real cos_t=qx[i];
        Real sin_t=sqrt(1.0-cos_t*cos_t);
        for(size_t j=0;j<p0*2;j++){
          X [(k*COORD_DIM+0)*Ngrid+i*p0*2+j]=-(1+0.05*(k-0.5))*cos_t;
          X [(k*COORD_DIM+1)*Ngrid+i*p0*2+j]= (1+0.05*(k-0.5))*sin_t*sin(j*M_PI/p0);
          X [(k*COORD_DIM+2)*Ngrid+i*p0*2+j]= (1+0.05*(k-0.5))*sin_t*cos(j*M_PI/p0);

          FS[(k*COORD_DIM+0)*Ngrid+i*p0*2+j]=0;
          FS[(k*COORD_DIM+1)*Ngrid+i*p0*2+j]=0;
          FS[(k*COORD_DIM+2)*Ngrid+i*p0*2+j]=0;

          FD[(k*COORD_DIM+0)*Ngrid+i*p0*2+j]=1;
          FD[(k*COORD_DIM+1)*Ngrid+i*p0*2+j]=1;
          FD[(k*COORD_DIM+2)*Ngrid+i*p0*2+j]=1;
        }
      }
    }
  }

  if(1){
    S.SetTrgCoord(NULL);
    S.SetSrcCoord(X);
    S.SetDensitySL(&FS);
    S.SetDensityDL(&FD);
    pvfmm::Vector<Real> vel=S();

    WriteVTK(X, 16, 64, "test", 0.0, &vel);
  }else{
    pvfmm::Vector<Real> T;
    { // Set coordinate values
      int p0=64;
      long Ngrid=(p0+1)*2*p0;
      std::vector<Real> qx;
      { // compute legendre node points
        qx.resize(p0+1);
        std::vector<Real> qw(p0+1);
        cgqf(p0+1, 1, 0.0, 0.0, -1.0, 1.0, &qx[0], &qw[0]);
      }
      T.ReInit(Ngrid*COORD_DIM);
      for(size_t i=0;i<p0+1;i++){
        Real cos_t=qx[i];
        Real sin_t=sqrt(1.0-cos_t*cos_t);
        for(size_t j=0;j<p0*2;j++){
          T[(i*p0*2+j)*COORD_DIM+0]=-1.01*cos_t;
          T[(i*p0*2+j)*COORD_DIM+1]= 1.01*sin_t*sin(j*M_PI/p0);
          T[(i*p0*2+j)*COORD_DIM+2]= 1.01*sin_t*cos(j*M_PI/p0);
        }
      }
    }

    S.SetTrgCoord(&T);
    S.SetSrcCoord(X);
    S.SetDensitySL(&FS);
    S.SetDensityDL(&FD);
    pvfmm::Vector<Real> vel=S();

    pvfmm::Matrix<Real> M;
    M.ReInit(  T.Dim()/COORD_DIM,COORD_DIM,&  T[0],false); M=M.Transpose();
    M.ReInit(vel.Dim()/COORD_DIM,COORD_DIM,&vel[0],false); M=M.Transpose();
    WriteVTK(T, 64, 64, "test", 0.0, &vel);
  }
}


template <class Real>
void WriteVTK(const pvfmm::Vector<Real>& S, long p0, long p1, const char* fname, Real period=0, const pvfmm::Vector<Real>* v_ptr=NULL, MPI_Comm comm=MPI_COMM_WORLD){
  typedef double VTKReal;
  int data__dof=COORD_DIM;

  pvfmm::Vector<Real> X, Xp, V, Vp;
  { // Upsample X
    const pvfmm::Vector<Real>& X0=S;
    pvfmm::Vector<Real> X1;
    SphericalHarmonics<Real>::Grid2SHC(X0,p0,p0,X1);
    SphericalHarmonics<Real>::SHC2Grid(X1,p0,p1,X);
    SphericalHarmonics<Real>::SHC2Pole(X1, p0, Xp);
  }
  if(v_ptr){ // Upsample V
    const pvfmm::Vector<Real>& X0=*v_ptr;
    pvfmm::Vector<Real> X1;
    SphericalHarmonics<Real>::Grid2SHC(X0,p0,p0,X1);
    SphericalHarmonics<Real>::SHC2Grid(X1,p0,p1,V);
    SphericalHarmonics<Real>::SHC2Pole(X1, p0, Vp);
  }

  std::vector<VTKReal> point_coord;
  std::vector<VTKReal> point_value;
  std::vector< int32_t> poly_connect;
  std::vector< int32_t> poly_offset;
  { // Set point_coord, point_value, poly_connect
    size_t N_ves = X.Dim()/(2*p1*(p1+1)*COORD_DIM); // Number of vesicles
    assert(Xp.Dim() == N_ves*2*COORD_DIM);
    for(size_t k=0;k<N_ves;k++){ // Set point_coord
      Real C[COORD_DIM]={0,0,0};
      if(period>0){
        for(long l=0;l<COORD_DIM;l++) C[l]=0;
        for(size_t i=0;i<p1+1;i++){
          for(size_t j=0;j<2*p1;j++){
            for(size_t l=0;l<COORD_DIM;l++){
              C[l]+=X[j+2*p1*(i+(p1+1)*(l+k*COORD_DIM))];
            }
          }
        }
        for(size_t l=0;l<COORD_DIM;l++) C[l]+=Xp[0+2*(l+k*COORD_DIM)];
        for(size_t l=0;l<COORD_DIM;l++) C[l]+=Xp[1+2*(l+k*COORD_DIM)];
        for(long l=0;l<COORD_DIM;l++) C[l]/=2*p1*(p1+1)+2;
        for(long l=0;l<COORD_DIM;l++) C[l]=(round(C[l]/period))*period;
      }

      for(size_t i=0;i<p1+1;i++){
        for(size_t j=0;j<2*p1;j++){
          for(size_t l=0;l<COORD_DIM;l++){
            point_coord.push_back(X[j+2*p1*(i+(p1+1)*(l+k*COORD_DIM))]-C[l]);
          }
        }
      }
      for(size_t l=0;l<COORD_DIM;l++) point_coord.push_back(Xp[0+2*(l+k*COORD_DIM)]-C[l]);
      for(size_t l=0;l<COORD_DIM;l++) point_coord.push_back(Xp[1+2*(l+k*COORD_DIM)]-C[l]);
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

  std::vector<VTKReal>& coord=point_coord;
  std::vector<VTKReal>& value=point_value;
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
  vtufile<<"        <DataArray type=\"Float"<<sizeof(VTKReal)*8<<"\" NumberOfComponents=\""<<COORD_DIM<<"\" Name=\"Position\" format=\"appended\" offset=\""<<data_size<<"\" />\n";
  data_size+=sizeof(uint32_t)+coord.size()*sizeof(VTKReal);
  vtufile<<"      </Points>\n";
  //---------------------------------------------------------------------------
  if(value.size()){ // value
    vtufile<<"      <PointData>\n";
    vtufile<<"        <DataArray type=\"Float"<<sizeof(VTKReal)*8<<"\" NumberOfComponents=\""<<value.size()/pt_cnt<<"\" Name=\""<<"value"<<"\" format=\"appended\" offset=\""<<data_size<<"\" />\n";
    data_size+=sizeof(uint32_t)+value.size()*sizeof(VTKReal);
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
  block_size=coord.size()*sizeof(VTKReal); vtufile.write((char*)&block_size, sizeof(int32_t)); vtufile.write((char*)&coord  [0], coord.size()*sizeof(VTKReal));
  if(value.size()){ // value
    block_size=value.size()*sizeof(VTKReal); vtufile.write((char*)&block_size, sizeof(int32_t)); vtufile.write((char*)&value  [0], value.size()*sizeof(VTKReal));
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
  pvtufile<<"        <PDataArray type=\"Float"<<sizeof(VTKReal)*8<<"\" NumberOfComponents=\""<<COORD_DIM<<"\" Name=\"Position\"/>\n";
  pvtufile<<"      </PPoints>\n";
  if(value.size()){ // value
    pvtufile<<"      <PPointData>\n";
    pvtufile<<"        <PDataArray type=\"Float"<<sizeof(VTKReal)*8<<"\" NumberOfComponents=\""<<value.size()/pt_cnt<<"\" Name=\""<<"value"<<"\"/>\n";
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

template <class Surf>
void WriteVTK(const Surf& S, const char* fname, MPI_Comm comm=MPI_COMM_WORLD, const typename Surf::Vec_t* v_ptr=NULL, int order=-1, typename Surf::value_type period=0){
  typedef typename Surf::value_type Real;
  typedef typename Surf::Vec_t Vec;
  size_t p0=S.getShOrder();
  size_t p1=(order>0?order:p0); // upsample

  pvfmm::Vector<Real> S_, v_;
  S_.ReInit(S.getPosition().size(),(Real*)S.getPosition().begin(),false);
  if(v_ptr) v_.ReInit(v_ptr->size(),(Real*)v_ptr->begin(),false);
  WriteVTK(S_, p0, p1, fname, period, (v_ptr?&v_:NULL), comm);
}

