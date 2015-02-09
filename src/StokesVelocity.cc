#ifdef HAVE_PVFMM

#include <omp.h>
#include <iostream>
#include <PVFMMInterface.h>

template<typename Surf_t>
void StokesVelocity<Surf_t>::SetSrcCoord(const Surf_t& S_){
  near_singular.SetSrcCoord(S_);
  S=&S_;
}

template<typename Surf_t>
void StokesVelocity<Surf_t>::SetTrgCoord(const Surf_t& T_){
  size_t omp_p=omp_get_max_threads();

  const Vec_t& x=T_.getPosition();
  size_t N_ves = x.getNumSubs(); // Number of vesicles
  size_t M_ves = x.getStride(); // Points per vesicle
  assert(M_ves==x.getGridDim().first*x.getGridDim().second);

  PVFMMVec_t& pt_coord=T;
  pt_coord.ReInit(N_ves*M_ves*COORD_DIM);

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
        pt_coord[(i*M_ves+j)*COORD_DIM+0]=xk[j];
        pt_coord[(i*M_ves+j)*COORD_DIM+1]=yk[j];
        pt_coord[(i*M_ves+j)*COORD_DIM+2]=zk[j];
      }
    }
  }
  near_singular.SetTrgCoord(&T[0],T.Dim()/COORD_DIM);
}

template<typename Surf_t>
void StokesVelocity<Surf_t>::SetTrgCoord(Real_t* trg_coord, size_t N){ // TODO: change to const Real_t*
  near_singular.SetTrgCoord(trg_coord,N);
  T.ReInit(N*COORD_DIM,trg_coord);
}

template<typename Surf_t>
void StokesVelocity<Surf_t>::SetSurfaceVel(const Vec_t& S_vel_){
  near_singular.SetSurfaceVel(S_vel_);
  S_vel=&S_vel_;
}

template<typename Surf_t>
void StokesVelocity<Surf_t>::SetDensity(const Vec_t* force_single_, const Vec_t* force_double_){
  near_singular.SetDensity(force_single_,force_double_);
  force_single=force_single_;
  force_double=force_double_;
}

template<typename Surf_t>
StokesVelocity<Surf_t>::Real_t* StokesVelocity<Surf_t>::operator()(){
  trg_vel.ReInit(T.Dim());
  trg_vel.SetZero(); // TODO: Check of this is needed.

  { // Compute velocity with FMM
    size_t omp_p=omp_get_max_threads();

    PVFMMVec_t src_coord;
    PVFMMVec_t force_sl;
    PVFMMVec_t force_dl;

    { // Set src_coord
      const Vec_t& x=S->getPosition();
      size_t N_ves = x.getNumSubs(); // Number of vesicles
      size_t M_ves = x.getStride(); // Points per vesicle
      assert(M_ves==x.getGridDim().first*x.getGridDim().second);
      src_coord.ReInit(N_ves*M_ves*COORD_DIM);

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

    if(force_single){
      const Vec_t& x=*force_single;
      size_t N_ves = x.getNumSubs(); // Number of vesicles
      size_t M_ves = x.getStride(); // Points per vesicle
      assert(M_ves==x.getGridDim().first*x.getGridDim().second);
      force_sl.ReInit(N_ves*M_ves*COORD_DIM);

      #pragma omp parallel for
      for(size_t tid=0;tid<omp_p;tid++){ // Set tree pt data
        size_t a=((tid+0)*N_ves)/omp_p;
        size_t b=((tid+1)*N_ves)/omp_p;

        for(size_t i=a;i<b;i++){
          // read each component of x
          const Real_t* fxk=x.getSubN_begin(i)+0*M_ves;
          const Real_t* fyk=x.getSubN_begin(i)+1*M_ves;
          const Real_t* fzk=x.getSubN_begin(i)+2*M_ves;

          for(size_t j=0;j<M_ves;j++){
            force_sl[(i*M_ves+j)*COORD_DIM+0]=fxk[j];
            force_sl[(i*M_ves+j)*COORD_DIM+1]=fyk[j];
            force_sl[(i*M_ves+j)*COORD_DIM+2]=fzk[j];
          }
        }
      }
    }

    if(force_double){
      const Vec_t& x=*force_double;
      const Vec_t& n=S->getNormal();
      size_t N_ves = x.getNumSubs(); // Number of vesicles
      size_t M_ves = x.getStride(); // Points per vesicle
      assert(M_ves==x.getGridDim().first*x.getGridDim().second);
      assert(N_ves == n.getNumSubs() && M_ves == n.getStride());
      force_dl.ReInit(N_ves*M_ves*(COORD_DIM+COORD_DIM));

      #pragma omp parallel for
      for(size_t tid=0;tid<omp_p;tid++){ // Set tree pt data
        size_t a=((tid+0)*N_ves)/omp_p;
        size_t b=((tid+1)*N_ves)/omp_p;

        for(size_t i=a;i<b;i++){
          // read each component of x
          const Real_t* fxk=x.getSubN_begin(i)+0*M_ves;
          const Real_t* fyk=x.getSubN_begin(i)+1*M_ves;
          const Real_t* fzk=x.getSubN_begin(i)+2*M_ves;

          const Real_t* nxk=n.getSubN_begin(i)+0*M_ves;
          const Real_t* nyk=n.getSubN_begin(i)+1*M_ves;
          const Real_t* nzk=n.getSubN_begin(i)+2*M_ves;

          for(size_t j=0;j<M_ves;j++){
            force_dl[(i*M_ves+j)*(COORD_DIM+COORD_DIM)+0]=fxk[j];
            force_dl[(i*M_ves+j)*(COORD_DIM+COORD_DIM)+1]=fyk[j];
            force_dl[(i*M_ves+j)*(COORD_DIM+COORD_DIM)+2]=fzk[j];
            force_dl[(i*M_ves+j)*(COORD_DIM+COORD_DIM)+3]=nxk[j];
            force_dl[(i*M_ves+j)*(COORD_DIM+COORD_DIM)+4]=nyk[j];
            force_dl[(i*M_ves+j)*(COORD_DIM+COORD_DIM)+5]=nzk[j];
          }
        }
      }
    }

    PVFMMEval(&src_coord[0],
              (force_sl.Dim()?&force_sl[0]:NULL),
              (force_dl.Dim()?&force_dl[0]:NULL),
              src_coord.Dim()/COORD_DIM,
              &T[0], &trg_vel[0], T.Dim()/COORD_DIM, &pvfmm_ctx);
  }

  { // Add near_singular
    Real_t* near_vel=near_singular();
    for(size_t i=0;i<trg_vel.Dim();i++){
      trg_vel[i]+=near_vel[i];
    }
  }

  return &trg_vel[0];
}

template<typename Surf_t>
void StokesVelocity<Surf_t>::operator()(Vec_t& T_vel){
  size_t omp_p=omp_get_max_threads();
  (*this)();

  const Vec_t& x=T_vel;
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
      const Real_t* xk=x.getSubN_begin(i)+0*M_ves;
      const Real_t* yk=x.getSubN_begin(i)+1*M_ves;
      const Real_t* zk=x.getSubN_begin(i)+2*M_ves;

      for(size_t j=0;j<M_ves;j++){
        xk[j]=trg_vel[(i*M_ves+j)*COORD_DIM+0];
        yk[j]=trg_vel[(i*M_ves+j)*COORD_DIM+1];
        zk[j]=trg_vel[(i*M_ves+j)*COORD_DIM+2];
      }
    }
  }
}

template<typename Surf_t>
void StokesVelocity<Surf_t>::Test(){
  typedef typename Vec_t::scalars_type Sca_t;
  typedef typename Vec_t::array_type Arr_t;
  typedef OperatorsMats<Arr_t> Mats_t;
  typedef VesInteraction<Real_t> Interaction_t;
  typedef InterfacialVelocity<Surf_t, Interaction_t> IntVel_t;

  int rank, size;
  MPI_Comm comm=MPI_COMM_WORLD;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);
  size_t nVec=2/size;
  if(!nVec) nVec=1;
  srand48(rank);

  //Set parameters
  Parameters<Real_t> sim_par;
  sim_par.sh_order = 32;
  sim_par.rep_up_freq = 32;

  //Create vectors
  Vec_t x0(nVec, sim_par.sh_order); // coordinates
  Vec_t v0(nVec, sim_par.sh_order); // velocity
  Vec_t v0_ref(nVec, sim_par.sh_order); // reference velocity
  Vec_t f0(nVec, sim_par.sh_order); // force
  int fLen = x0.getStride();
  { //Set coordinate values
    int imax(x0.getGridDim().first);
    int jmax(x0.getGridDim().second);
    Real_t scal=pow((Real_t)nVec*size,1.0/3.0);

    std::vector<Real_t> qx;
    { // compute legendre node points
      size_t p=x0.getShOrder();
      qx.resize(p+1);
      std::vector<Real_t> qw(p+1);
      cgqf(p+1, 1, 0.0, 0.0, -1.0, 1.0, &qx[0], &qw[0]);
    }

    for(size_t k=0;k<nVec;k++){
      Real_t* x_k=x0.getSubN_begin(k)+0*fLen;
      Real_t* y_k=x0.getSubN_begin(k)+1*fLen;
      Real_t* z_k=x0.getSubN_begin(k)+2*fLen;

      Real_t* vx_k=v0.getSubN_begin(k)+0*fLen;
      Real_t* vy_k=v0.getSubN_begin(k)+1*fLen;
      Real_t* vz_k=v0.getSubN_begin(k)+2*fLen;
      Real_t* fx_k=f0.getSubN_begin(k)+0*fLen;
      Real_t* fy_k=f0.getSubN_begin(k)+1*fLen;
      Real_t* fz_k=f0.getSubN_begin(k)+2*fLen;

      Real_t coord[3];
      coord[0]=scal*drand48();
      coord[1]=scal*drand48();
      coord[2]=scal*drand48();
      for(size_t i=0;i<imax;i++){
        Real_t cos_t=qx[i];
        //Real_t cos_t=cos((i+1)*M_PI/(imax+1));
        Real_t sin_t=sqrt(1.0-cos_t*cos_t);
        for(size_t j=0;j<jmax;j++){

          if(0){
            x_k[j+i*jmax]=cos_t;
            y_k[j+i*jmax]=sin_t*sin(j*2*M_PI/jmax);
            z_k[j+i*jmax]=sin_t*cos(j*2*M_PI/jmax);
            if(k){
              x_k[j+i*jmax]*=0.49;
              y_k[j+i*jmax]*=0.49;
              z_k[j+i*jmax]*=0.49;
              y_k[j+i*jmax]+=1.49+1e-4;
            }
          }else{
            x_k[j+i*jmax]=coord[0]+scal*cos_t;
            y_k[j+i*jmax]=coord[1]+scal*sin_t*sin(j*2*M_PI/jmax);
            z_k[j+i*jmax]=coord[2]+scal*sin_t*cos(j*2*M_PI/jmax);
          }

          vx_k[j+i*jmax]=0;
          vy_k[j+i*jmax]=0;
          vz_k[j+i*jmax]=0;
          fx_k[j+i*jmax]=0;
          fy_k[j+i*jmax]=0;
          fz_k[j+i*jmax]=0;

          if(!rank && k==0){
            Real_t c[3]={x_k[j+i*jmax],
                       y_k[j+i*jmax],
                       z_k[j+i*jmax]};
            Real_t f[3]={0,0,0};
            force(c,1,f);
            fx_k[j+i*jmax]=f[0];
            fy_k[j+i*jmax]=f[1];
            fz_k[j+i*jmax]=f[2];
          }
        }
      }
    }
  }

  //Reading operators from file
  bool readFromFile = true;
  Mats_t mats(readFromFile, sim_par);

  //Creating objects
  COUT("Creating the surface object");
  Surf_t S(x0, mats);

  //Setting the background flow
  ShearFlow<Vec_t> vInf(sim_par.bg_flow_param);

  Interaction_t interaction(NULL);
  IntVel_t F(S, interaction, mats, sim_par, vInf);

  // Compute self interaction
  Vec_t veloc_self(nVec, sim_par.sh_order);
  F.stokes(f0,veloc_self);

  Vec_t qforce(nVec,sim_par.sh_order);
  { //Incorporating the quadrature weights and area into the force
    Sca_t quad_weights_;
    { // quadrature weights
      quad_weights_.resize(1,sim_par.sh_order);
      quad_weights_.getDevice().Memcpy(quad_weights_.begin(),
          mats.quad_weights_,
          quad_weights_.size() * sizeof(Real_t),
          Surf_t::device_type::MemcpyDeviceToDevice);
    }
    xv(S.getAreaElement(), f0, qforce);
    ax<Sca_t>(quad_weights_, qforce, qforce);
  }

  StokesVelocity<Surf_t> stokes_vel(comm);
  stokes_vel.SetSrcCoord(S);
  stokes_vel.SetTrgCoord(S);
  stokes_vel.SetSurfaceVel(veloc_self);
  stokes_vel.SetDensity(&qforce, NULL);

  pvfmm::Profile::Tic("NearInteractions",&comm);
  stokes_vel();
  pvfmm::Profile::Toc();
}

#endif
