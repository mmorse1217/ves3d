#ifdef HAVE_PVFMM

#include <omp.h>
#include <iostream>
#include "Surface.h"
//#include "BgFlow.h"
#include <profile.hpp>

template<typename Surf_t>
StokesVelocity<Surf_t>::StokesVelocity(
    OperatorsMats<Arr_t> &mats,
    const Parameters<Real_t> &sim_par_,
    //const BgFlowBase<Vec_t> &bgFlow_,
    MPI_Comm c):
  move_pole(mats),
  sim_par(sim_par_),
  //bg_flow(bgFlow_),
  near_singular(c),
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

  sh_order=sim_par.sh_order;
  { // quadrature weights
    quad_weights_.resize(1,sh_order);
    quad_weights_.getDevice().Memcpy(quad_weights_.begin(),
        mats.quad_weights_,
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

  pvfmm_ctx=PVFMMCreateContext<Real_t>();
}

template<typename Surf_t>
StokesVelocity<Surf_t>::~StokesVelocity(){
  PVFMMDestroyContext<Real_t>(&pvfmm_ctx);
}



template<typename Surf_t>
void StokesVelocity<Surf_t>::SetSrcCoord(const Surf_t& S_){
  self_flag=self_flag | StokesVelocity::UpdateSrcCoord;
  fmm_flag=fmm_flag | StokesVelocity::UpdateSrcCoord;
  near_singular.SetSrcCoord(S_);
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
  near_singular.SetTrgCoord(&trg_coord[0],trg_coord.Dim()/COORD_DIM);
}

template<typename Surf_t>
void StokesVelocity<Surf_t>::SetTrgCoord(Real_t* trg_coord_, size_t N){ // TODO: change to const Real_t*
  fmm_flag=fmm_flag | StokesVelocity::UpdateTrgCoord;
  trg_coord.ReInit(N*COORD_DIM, trg_coord_);
  near_singular.SetTrgCoord(&trg_coord[0],N);
}



template<typename Surf_t>
void StokesVelocity<Surf_t>::ApplyQuadWeights(){
  Vec_t qforce;
  const int force_dim=COORD_DIM;
  size_t omp_p=omp_get_max_threads();
  if(fmm_flag & StokesVelocity::UpdateDensitySL){
    fmm_flag=fmm_flag & ~StokesVelocity::UpdateDensitySL;
    fmm_flag=fmm_flag | StokesVelocity::UpdateqDensitySL;
    if(force_single){
      size_t N_ves = force_single->getNumSubs(); // Number of vesicles
      size_t M_ves = force_single->getStride(); // Points per vesicle

      assert(S);
      qforce.resize(N_ves, sh_order);
      xv(S->getAreaElement(), *force_single, qforce);
      ax<Sca_t>(quad_weights_, qforce, qforce);

      qforce_single.ReInit(N_ves*M_ves*(force_dim));
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
            qforce_single[(i*M_ves+j)*force_dim+0]=fxk[j];
            qforce_single[(i*M_ves+j)*force_dim+1]=fyk[j];
            qforce_single[(i*M_ves+j)*force_dim+2]=fzk[j];
          }
        }
      }
      near_singular.SetDensitySL(&qforce_single);
    }else{
      qforce_single.ReInit(0);
      near_singular.SetDensitySL(NULL);
    }
  }
  if(fmm_flag & StokesVelocity::UpdateDensityDL){
    fmm_flag=fmm_flag & ~StokesVelocity::UpdateDensityDL;
    fmm_flag=fmm_flag | StokesVelocity::UpdateqDensityDL;
    if(force_double){
      size_t N_ves = force_double->getNumSubs(); // Number of vesicles
      size_t M_ves = force_double->getStride(); // Points per vesicle

      assert(S);
      qforce.resize(N_ves, sh_order);
      xv(S->getAreaElement(), *force_double, qforce);
      ax<Sca_t>(quad_weights_, qforce, qforce);

      qforce_double.ReInit(N_ves*M_ves*(force_dim+COORD_DIM));
      const Vec_t* normal=&S->getNormal();
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
            qforce_double[(i*M_ves+j)*(force_dim+COORD_DIM)+0]=fxk[j];
            qforce_double[(i*M_ves+j)*(force_dim+COORD_DIM)+1]=fyk[j];
            qforce_double[(i*M_ves+j)*(force_dim+COORD_DIM)+2]=fzk[j];
            qforce_double[(i*M_ves+j)*(force_dim+COORD_DIM)+3]=nxk[j];
            qforce_double[(i*M_ves+j)*(force_dim+COORD_DIM)+4]=nyk[j];
            qforce_double[(i*M_ves+j)*(force_dim+COORD_DIM)+5]=nzk[j];
          }
        }
      }
      near_singular.SetDensityDL(&qforce_double);
    }else{
      qforce_double.ReInit(0);
      near_singular.SetDensityDL(NULL);
    }
  }
}

template<typename Surf_t>
const StokesVelocity<Surf_t>::Vec_t& StokesVelocity<Surf_t>::SelfInteraction(bool update_self){
  if(self_flag & StokesVelocity::UpdateSurfaceVel){
    if(S_vel_ptr){ // Nothing to do; it has already been computed for us.
      near_singular.SetSurfaceVel(*S_vel_ptr);
      self_flag=StokesVelocity::UpdateNone;
    }else{ // Compute from SL and DL density
      self_flag=self_flag | StokesVelocity::UpdateDensitySL;
      self_flag=self_flag | StokesVelocity::UpdateDensityDL;
    }
  }
  if(update_self && self_flag){ // Compute self interaction
    { // Compute S_vel
      pvfmm::Profile::Tic("SelfInteraction",&comm);
      assert(S);
      int imax(S->getPosition().getGridDim().first);
      int jmax(S->getPosition().getGridDim().second);
      int np = S->getPosition().getStride();
      int nv = S->getPosition().getNumSubs();
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

          int numinputs = 4;
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
      S_vel.resize(nv, sh_order);
      axpy(1.0,SL_vel,DL_vel,S_vel);
      pvfmm::Profile::Toc();
    }
    S_vel_ptr=&S_vel;
    near_singular.SetSurfaceVel(*S_vel_ptr);
    self_flag=StokesVelocity::UpdateNone;
  }
  return *S_vel_ptr;
}

template<typename Surf_t>
const StokesVelocity<Surf_t>::PVFMMVec_t& StokesVelocity<Surf_t>::NearInteraction(bool update_near){
  if(update_near) ApplyQuadWeights();
  return near_singular(update_near);
}

template<typename Surf_t>
const StokesVelocity<Surf_t>::PVFMMVec_t& StokesVelocity<Surf_t>::FarInteraction(bool update_far){
  if(update_far && fmm_flag){ // Compute velocity with FMM
    ApplyQuadWeights();
    if(fmm_flag & StokesVelocity::UpdateSrcCoord){ // Set src_coord
      fmm_flag=fmm_flag & ~StokesVelocity::UpdateSrcCoord;

      assert(S);
      const Vec_t& x=S->getPosition();
      size_t N_ves = x.getNumSubs(); // Number of vesicles
      size_t M_ves = x.getStride(); // Points per vesicle
      assert(M_ves==x.getGridDim().first*x.getGridDim().second);
      src_coord.ReInit(N_ves*M_ves*COORD_DIM);

      size_t omp_p=omp_get_max_threads();
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
            src_coord[(i*M_ves+j)*COORD_DIM+0]=xk[j];
            src_coord[(i*M_ves+j)*COORD_DIM+1]=yk[j];
            src_coord[(i*M_ves+j)*COORD_DIM+2]=zk[j];
          }
        }
      }
    }
    fmm_vel.ReInit(trg_coord.Dim());
    PVFMMEval(&src_coord[0],
              (qforce_single.Dim()?&qforce_single[0]:NULL),
              (qforce_double.Dim()?&qforce_double[0]:NULL),
              src_coord.Dim()/COORD_DIM,
              &trg_coord[0], &fmm_vel[0], trg_coord.Dim()/COORD_DIM, &pvfmm_ctx);
    near_singular.SubtractDirect(fmm_vel);
    fmm_flag=StokesVelocity::UpdateNone;
  }
  return fmm_vel;
}



template<typename Surf_t>
StokesVelocity<Surf_t>::Real_t* StokesVelocity<Surf_t>::operator()(unsigned int flag){
  bool update_self=(flag & StokesVelocity::UpdateSelf);
  bool update_near=(flag & StokesVelocity::UpdateNear);
  bool update_far =(flag & StokesVelocity::UpdateFar );

  SelfInteraction(update_self);
  const PVFMMVec_t& near_vel=NearInteraction(update_near);
  const PVFMMVec_t& far_vel=FarInteraction(update_far);
  { // Compute far + near
    trg_vel.ReInit(trg_coord.Dim());
    assert(trg_vel.Dim()==far_vel.Dim());
    assert(trg_vel.Dim()==near_vel.Dim());
    #pragma omp parallel for
    for(size_t i=0;i<trg_vel.Dim();i++){
      trg_vel[i]=fmm_vel[i]+near_vel[i];
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

  // Background velocity
  //!@bug the time should be passed to the BgFlow handle.
  //bg_flow_(S->getPosition(), 0, x);

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



template<typename Real_t>
static void u_ref(const Real_t* coord, int n, Real_t* out){ //Analytical velocity for sphere
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

template<typename Real_t>
static void force(const Real_t* coord, int n, Real_t* out){ // Force on sphere
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
  sim_par.sh_order = 32;
  sim_par.rep_up_freq = 32;

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
  Vec_t f0(nVec, sim_par.sh_order); // force
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

      Real_t* fx_k=f0.getSubN_begin(k)+0*fLen;
      Real_t* fy_k=f0.getSubN_begin(k)+1*fLen;
      Real_t* fz_k=f0.getSubN_begin(k)+2*fLen;
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
            Real_t f[3]={0,0,0};
            force(c,1,f);
            fx_k[j+i*jmax]=f[0];
            fy_k[j+i*jmax]=f[1];
            fz_k[j+i*jmax]=f[2];
          }else{
            fx_k[j+i*jmax]=0;
            fy_k[j+i*jmax]=0;
            fz_k[j+i*jmax]=0;
          }
        }
      }
    }
  }

  // Creating surface objects
  Surf_t S(x0, mats);

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
        trg_coord[(j+i*jmax)*COORD_DIM+0]=R*cos_t;
        trg_coord[(j+i*jmax)*COORD_DIM+1]=R*sin_t*sin(j*2*M_PI/jmax);
        trg_coord[(j+i*jmax)*COORD_DIM+2]=R*sin_t*cos(j*2*M_PI/jmax);
      }
    }
  }

  stokes_vel.SetSrcCoord(S);
  //stokes_vel.SetSurfaceVel(vel_surf); // optional
  stokes_vel.SetDensitySL(&f0);
  stokes_vel.SetDensityDL(NULL);

  stokes_vel.SetTrgCoord(&trg_coord[0],trg_coord.Dim()/COORD_DIM);
  //stokes_vel.SetTrgCoord(S);

  pvfmm::Profile::Tic("NearInteractions",&comm);
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

#endif

