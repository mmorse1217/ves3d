#include <omp.h>
#include <iostream>
#include "Surface.h"
#include <profile.hpp>

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
    Real_t box_size,
    MPI_Comm c):
  near_singular(sht_up_.getShOrder(),box_size,c),
  sht_   (mats.p_   , mats.mats_p_   ),
  sht_up_(mats.p_up_, mats.mats_p_up_),
  box_size_(box_size),
  sim_par(sim_par_),
  move_pole(mats),
  S_up(NULL),
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

  sh_order   =sht_   .getShOrder();
  sh_order_up=sht_up_.getShOrder();
  assert(sim_par.sh_order==sht_.getShOrder());
  { // quadrature weights
    quad_weights_.resize(1,sh_order_up);
    quad_weights_.getDevice().Memcpy(quad_weights_.begin(),
        mats.quad_weights_p_up_,
        quad_weights_.size() * sizeof(Real_t),
        Surf_t::device_type::MemcpyDeviceToDevice);
  }

  { // Set w_sph_, w_sph_inv_, sing_quad_weights_
    w_sph_.resize(1, sh_order);
    w_sph_inv_.resize(1, sh_order);
    int np = w_sph_.getStride();

    w_sph_.getDevice().Memcpy(w_sph_.begin(), mats.w_sph_,
        np * sizeof(Real_t), device_type::MemcpyDeviceToDevice);
    xInv(w_sph_,w_sph_inv_);

    //Singular quadrature weights
    sing_quad_weights_.resize(1,sh_order);
    sing_quad_weights_.getDevice().Memcpy(sing_quad_weights_.begin(),
        mats.sing_quad_weights_, sing_quad_weights_.size() *
        sizeof(Real_t),
        device_type::MemcpyDeviceToDevice);
  }

  { // compute quadrature to find pole
    size_t p=sh_order;

    //Gauss-Legendre quadrature nodes and weights
    std::vector<Real_t> x(p+1),w(p+1);
    cgqf(p+1, 1, 0.0, 0.0, -1.0, 1.0, &x[0], &w[0]);

    std::vector<Real_t> leg((p+1)*(p+1));
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
}



template<typename Surf_t>
void StokesVelocity<Surf_t>::SetSrcCoord(const Surf_t& S_){
  self_flag=self_flag | StokesVelocity::UpdateSrcCoord;
  fmm_flag=fmm_flag | StokesVelocity::UpdateSrcCoord;
  near_singular.SetSrcCoord(src_coord_up);
  S=&S_;
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
  near_singular.SetTrgCoord(&trg_coord[0],N,trg_is_surf);
}



template<typename Surf_t>
void StokesVelocity<Surf_t>::GetPole(const Vec_t& v,  PVFMMVec_t& pvfmm_v){
  assert(v.getShOrder()==sh_order);
  size_t omp_p=omp_get_max_threads();
  size_t N_ves   =v.getNumSubs(); // Number of vesicles
  size_t M_ves   =(1+sh_order   )*(2*sh_order   );
  size_t M_ves_up=(1+sh_order_up)*(2*sh_order_up);
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
  wrk[0].resize(v.getNumSubs(), std::max(sh_order_up,sh_order));
  wrk[1].resize(v.getNumSubs(), std::max(sh_order_up,sh_order));
  v_out->resize(v.getNumSubs(), std::max(sh_order_up,sh_order));
  Resample(v, sht_, sht_up_, wrk[0], wrk[1], *v_out);

  if(pvfmm_v){
    Vec2PVFMMVec(*v_out, *pvfmm_v);
    GetPole     ( v    , *pvfmm_v);
  }
}

template<typename Surf_t>
void StokesVelocity<Surf_t>::Setup(){
  assert(S);
  size_t omp_p=omp_get_max_threads();
  if(fmm_flag & StokesVelocity::UpdateSrcCoord){ // Compute src_coord_up
    fmm_flag=fmm_flag & ~StokesVelocity::UpdateSrcCoord;

    S_up=NULL; // TODO delete old S_up? Shouldn't Surface class do that?
    if(sh_order_up==sh_order) S_up=S;
    else CHK(S->resample(sh_order_up, (Surf_t**)&S_up));
    const Vec_t& S_coord   =S   ->getPosition();
    const Vec_t& S_coord_up=S_up->getPosition();
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

      xv(S_up->getAreaElement(), qforce, qforce);
      ax<Sca_t>(quad_weights_, qforce, qforce);

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

      xv(S_up->getAreaElement(), qforce, qforce);
      ax<Sca_t>(quad_weights_, qforce, qforce);

      qforce_double_up.ReInit(N_ves*(M_ves+2)*(force_dim+COORD_DIM));
      const Vec_t* normal=&S_up->getNormal();
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
    { // Compute S_vel
      pvfmm::Profile::Tic("SelfInteraction",&comm);
      bool prof_state=pvfmm::Profile::Enable(false);
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

            ax<Sca_t>(w_sph_, force_single_out, force_single_out);
            S->getPosition().getDevice().DirectStokes(coord_out.begin(), force_single_out.begin(),
                sing_quad_weights_.begin(), np, nv, S->getPosition().begin(),
                ii * jmax + jj, ii * jmax + jj + 1, SL_vel.begin());

            ax<Sca_t>(w_sph_, force_double_out, force_double_out);
            S->getPosition().getDevice().DirectStokesDoubleLayer(coord_out.begin(), normal_out.begin(), force_double_out.begin(),
                sing_quad_weights_.begin(), np, nv, S->getPosition().begin(),
                ii * jmax + jj, ii * jmax + jj + 1, DL_vel.begin());
          }
        }
      }else{
        if(self_flag & (StokesVelocity::UpdateDensitySL | StokesVelocity::UpdateSrcCoord)){ // Single layer
          SL_vel.resize(nv, sh_order);
          if(force_single){
            Vec_t v1(nv, sh_order);
            Vec_t v2(nv, sh_order);
            Vec_t v3(nv, sh_order);

            Sca_t t1(nv, sh_order);
            ax(w_sph_inv_, S->getAreaElement(), t1);
            xv(        t1,       *force_single, v3);

            int numinputs = 2;
            const Sca_t* inputs[] = {&S->getPosition(), &v3};
            Sca_t*      outputs[] = {&              v1, &v2};
            move_pole.setOperands(inputs, numinputs, sim_par.singular_stokes);

            for(int ii=0;ii < imax; ++ii)
            for(int jj=0;jj < jmax; ++jj){
              move_pole(ii, jj, outputs);

              ax<Sca_t>(w_sph_, v2, v2);
              S->getPosition().getDevice().DirectStokes(v1.begin(), v2.begin(),
                  sing_quad_weights_.begin(), np, nv, S->getPosition().begin(),
                  ii * jmax + jj, ii * jmax + jj + 1, SL_vel.begin());
            }
          }else{
            Vec_t::getDevice().Memset(SL_vel.begin(),0,SL_vel.size()*sizeof(Real_t));
          }
        }
        if(self_flag & (StokesVelocity::UpdateDensityDL | StokesVelocity::UpdateSrcCoord)){ // Double layer
          DL_vel.resize(nv, sh_order);
          if(force_double){
            Vec_t v1(nv, sh_order);
            Vec_t v2(nv, sh_order);
            Vec_t v3(nv, sh_order);
            Vec_t v4(nv, sh_order);

            Sca_t t1(nv, sh_order);
            ax(w_sph_inv_, S->getAreaElement(), t1);
            xv(        t1,       *force_double, v3);

            int numinputs = 3;
            const Sca_t* inputs[] = {&S->getPosition(), &S->getNormal(), &v3};
            Sca_t*      outputs[] = {&              v1, &            v4, &v2};
            move_pole.setOperands(inputs, numinputs, sim_par.singular_stokes);

            for(int ii=0;ii < imax; ++ii)
            for(int jj=0;jj < jmax; ++jj){
              move_pole(ii, jj, outputs);

              ax<Sca_t>(w_sph_, v2, v2);
              S->getPosition().getDevice().DirectStokesDoubleLayer(v1.begin(), v4.begin(), v2.begin(),
                  sing_quad_weights_.begin(), np, nv, S->getPosition().begin(),
                  ii * jmax + jj, ii * jmax + jj + 1, DL_vel.begin());
            }
          }else{
            Vec_t::getDevice().Memset(DL_vel.begin(),0,DL_vel.size()*sizeof(Real_t));
          }
        }
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

  // Reading operators from file
  bool readFromFile = true;
  Mats_t mats(readFromFile, sim_par);

  //Setting the background flow
  //ShearFlow<Vec_t> vInf(sim_par.bg_flow_param);

  //StokesVelocity<Surf_t> stokes_vel(mats, sim_par, vInf, comm);
  StokesVelocity<Surf_t> stokes_vel(mats, sim_par, -5, comm);

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
            x_k[j+i*jmax]=cos_t;
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
void WriteVTK(const Surf_t& S, const char* fname, MPI_Comm comm=MPI_COMM_WORLD){
  typedef typename Surf_t::value_type Real_t;
  typedef typename Surf_t::Vec_t Vec_t;
  typedef float VTKReal_t;

  std::vector<Real_t> pole_quad;
  { // compute quadrature to find pole
    size_t p=S.getShOrder();

    //Gauss-Legendre quadrature nodes and weights
    std::vector<Real_t> x(p+1),w(p+1);
    cgqf(p+1, 1, 0.0, 0.0, -1.0, 1.0, &x[0], &w[0]);

    std::vector<Real_t> leg((p+1)*(p+1));
    LegPoly(&x[0],x.size(),p+1,&leg[0]);
    pole_quad.resize(p+1,0);
    for(size_t j=0;j<p+1;j++){
      for(size_t i=0;i<p+1;i++){
        pole_quad[i]+=leg[j*(p+1)+i]*sqrt(2.0*j+1.0);
      }
    }
    for(size_t i=0;i<p+1;i++){
      pole_quad[i]*=w[i]*0.25/p;
    }
  }

  std::vector<VTKReal_t> point_coord;
  std::vector< int32_t> poly_connect;
  std::vector< int32_t> poly_offset;
  { // Set point_coord, poly_connect
    const Vec_t& x=S.getPosition();
    size_t N_ves = x.getNumSubs(); // Number of vesicles
    size_t M_ves = x.getStride(); // Points per vesicle
    int imax(x.getGridDim().first);
    int jmax(x.getGridDim().second);
    int fLen = x.getStride();
    assert(fLen==M_ves);
    for(size_t k=0;k<N_ves;k++){
      std::vector<Real_t> pole_coord(COORD_DIM*2,0);

      const Real_t* xk=x.getSubN_begin(k)+0*fLen;
      const Real_t* yk=x.getSubN_begin(k)+1*fLen;
      const Real_t* zk=x.getSubN_begin(k)+2*fLen;
      for(size_t i=0;i<imax;i++){
        for(size_t j=0;j<jmax;j++){
          point_coord.push_back(xk[j+i*jmax]);
          point_coord.push_back(yk[j+i*jmax]);
          point_coord.push_back(zk[j+i*jmax]);

          pole_coord[0*COORD_DIM+0]+=pole_quad[imax-1-i]*xk[j+i*jmax];
          pole_coord[0*COORD_DIM+1]+=pole_quad[imax-1-i]*yk[j+i*jmax];
          pole_coord[0*COORD_DIM+2]+=pole_quad[imax-1-i]*zk[j+i*jmax];
          pole_coord[1*COORD_DIM+0]+=pole_quad[       i]*xk[j+i*jmax];
          pole_coord[1*COORD_DIM+1]+=pole_quad[       i]*yk[j+i*jmax];
          pole_coord[1*COORD_DIM+2]+=pole_quad[       i]*zk[j+i*jmax];
        }
      }
      point_coord.push_back(pole_coord[0]);
      point_coord.push_back(pole_coord[1]);
      point_coord.push_back(pole_coord[2]);
      point_coord.push_back(pole_coord[3]);
      point_coord.push_back(pole_coord[4]);
      point_coord.push_back(pole_coord[5]);
    }

    for(size_t k=0;k<N_ves;k++){
      for(size_t j=0;j<jmax;j++){
        size_t i0=     0;
        size_t i1=imax-1;
        size_t j0=((j+0)     );
        size_t j1=((j+1)%jmax);

        poly_connect.push_back((jmax*imax+2)*k + jmax*imax+0);
        poly_connect.push_back((jmax*imax+2)*k + jmax*i0+j0);
        poly_connect.push_back((jmax*imax+2)*k + jmax*i0+j1);
        poly_offset.push_back(poly_connect.size());

        poly_connect.push_back((jmax*imax+2)*k + jmax*imax+1);
        poly_connect.push_back((jmax*imax+2)*k + jmax*i1+j0);
        poly_connect.push_back((jmax*imax+2)*k + jmax*i1+j1);
        poly_offset.push_back(poly_connect.size());
      }
      for(size_t i=0;i<imax-1;i++){
        for(size_t j=0;j<jmax;j++){
          size_t i0=((i+0)     );
          size_t i1=((i+1)     );
          size_t j0=((j+0)     );
          size_t j1=((j+1)%jmax);
          poly_connect.push_back((jmax*imax+2)*k + jmax*i0+j0);
          poly_connect.push_back((jmax*imax+2)*k + jmax*i1+j0);
          poly_connect.push_back((jmax*imax+2)*k + jmax*i1+j1);
          poly_connect.push_back((jmax*imax+2)*k + jmax*i0+j1);
          poly_offset.push_back(poly_connect.size());
        }
      }
    }

  }

  int myrank, np;
  MPI_Comm_size(comm,&np);
  MPI_Comm_rank(comm,&myrank);

  std::vector<VTKReal_t>& coord=point_coord;
  std::vector<int32_t> connect=poly_connect;
  std::vector<int32_t> offset=poly_offset;

  int pt_cnt=coord.size()/COORD_DIM;
  int poly_cnt=poly_offset.size();

  //Open file for writing.
  std::stringstream vtufname;
  vtufname<<fname<<std::setfill('0')<<std::setw(6)<<myrank<<".vtp";
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
  block_size=coord   .size()*sizeof(VTKReal_t); vtufile.write((char*)&block_size, sizeof(int32_t)); vtufile.write((char*)&coord   [0], coord   .size()*sizeof(VTKReal_t));
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
  {
    // Extract filename from path.
    std::stringstream vtupath;
    vtupath<<'/'<<fname;
    std::string pathname = vtupath.str();
    unsigned found = pathname.find_last_of("/\\");
    std::string fname_ = pathname.substr(found+1);
    for(int i=0;i<np;i++) pvtufile<<"      <Piece Source=\""<<fname_<<std::setfill('0')<<std::setw(6)<<i<<".vtp\"/>\n";
  }
  pvtufile<<"  </PPolyData>\n";
  pvtufile<<"</VTKFile>\n";
  pvtufile.close();
}
