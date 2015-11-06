#include <omp.h>
#include <iostream>
#include "Surface.h"

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
    const Parameters<Real_t> &sim_par_):
  move_pole(mats),
  sim_par(sim_par_)
{
  S=NULL;
  force_single=NULL;
  force_double=NULL;

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
}

template<typename Surf_t>
StokesVelocity<Surf_t>::~StokesVelocity(){
}


template<typename Surf_t>
void StokesVelocity<Surf_t>::SetSrcCoord(const Surf_t& S_){
  S=&S_;
}

template<typename Surf_t>
void StokesVelocity<Surf_t>::SetDensitySL(const Vec_t* force_single_){
  force_single=force_single_;
}

template<typename Surf_t>
void StokesVelocity<Surf_t>::SetDensityDL(const Vec_t* force_double_){
  force_double=force_double_;
}



template<typename Surf_t>
void StokesVelocity<Surf_t>::SetTrgCoord(const Surf_t& T_){
  T=&T_;
}




template<typename Surf_t>
StokesVelocity<Surf_t>::Real_t* StokesVelocity<Surf_t>::operator()(unsigned int flag){
  Vec_t T_vel;
  (*this)(T_vel, flag);

  Vec_t& x=T_vel;
  size_t N_ves = x.getNumSubs(); // Number of vesicles
  size_t M_ves = x.getStride(); // Points per vesicle
  assert(M_ves==x.getGridDim().first*x.getGridDim().second);

  trg_vel.resize(N_ves*M_ves*3);
  assert(N_ves*M_ves*3==trg_vel.size());

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
        trg_vel[(i*M_ves+j)*3+0]=xk[j];
        trg_vel[(i*M_ves+j)*3+1]=yk[j];
        trg_vel[(i*M_ves+j)*3+2]=zk[j];
      }
    }
  }

  return &trg_vel[0];
}

template<typename Surf_t>
void StokesVelocity<Surf_t>::operator()(Vec_t& T_vel, unsigned int flag){
  assert(S==T);
  assert(S);

  { // Compute S_vel
    int imax(S->getPosition().getGridDim().first);
    int jmax(S->getPosition().getGridDim().second);
    int np = S->getPosition().getStride();
    int nv = S->getPosition().getNumSubs();
    { // Single layer
      SL_vel.resize(nv, sh_order);
      Vec_t::getDevice().Memset(SL_vel.begin(),0,SL_vel.size()*sizeof(Real_t));
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
              sing_quad_weights_.begin(), np, np, nv, S->getPosition().begin(),
              ii * jmax + jj, ii * jmax + jj + 1, SL_vel.begin());
        }
      }
    }
    { // Double layer
      DL_vel.resize(nv, sh_order);
      Vec_t::getDevice().Memset(DL_vel.begin(),0,DL_vel.size()*sizeof(Real_t));
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
              sing_quad_weights_.begin(), np, np, nv, S->getPosition().begin(),
              ii * jmax + jj, ii * jmax + jj + 1, DL_vel.begin());
        }
      }
    }

    S_vel.resize(nv, sh_order);
    Vec_t::getDevice().Memset(S_vel.begin(),0,S_vel.size()*sizeof(Real_t));
    axpy(1.0,SL_vel,DL_vel,S_vel);
  }

  { // Compute direct interaction
    int N_ves = S->getPosition().getNumSubs();
    int M_ves = S->getPosition().getStride();

    T_vel.resize(N_ves, sh_order);
    Vec_t::getDevice().Memset(T_vel.begin(),0,T_vel.size()*sizeof(Real_t));

    if(force_single){ // Single layer
      Vec_t qforce;
      qforce.resize(N_ves, sh_order);
      xv(S->getAreaElement(), *force_single, qforce);
      ax<Sca_t>(quad_weights_, qforce, qforce);

      Vec_t tmp_vel;
      tmp_vel.resize(N_ves, sh_order);
      for(size_t s=0;s<N_ves;s++){
        Vec_t::getDevice().Memset(tmp_vel.begin(),0,tmp_vel.size()*sizeof(Real_t));

        for(size_t t=0;t<N_ves;t++) if(s!=t){
          tmp_vel.getDevice().DirectStokes(
              S->getPosition().begin()+3*M_ves*s, qforce.begin()+3*M_ves*s, (Real_t*)NULL,
              M_ves, M_ves, 1, S->getPosition().begin()+3*M_ves*t /* target */,
              0, M_ves /* number of trgs per surface */,
              tmp_vel.begin()+3*M_ves*t);
        }

        axpy(1.0,tmp_vel,T_vel,T_vel); // accumulate
      }
    }

    if(force_double){ // Double layer
      Vec_t qforce;
      qforce.resize(N_ves, sh_order);
      xv(S->getAreaElement(), *force_double, qforce);
      ax<Sca_t>(quad_weights_, qforce, qforce);

      Vec_t tmp_vel;
      tmp_vel.resize(N_ves, sh_order);
      for(size_t s=0;s<N_ves;s++){
        Vec_t::getDevice().Memset(tmp_vel.begin(),0,tmp_vel.size()*sizeof(Real_t));

        for(size_t t=0;t<N_ves;t++) if(s!=t){
          tmp_vel.getDevice().DirectStokesDoubleLayer(
              S->getPosition().begin()+3*M_ves*s, S->getNormal().begin()+3*M_ves*s, qforce.begin()+3*M_ves*s, (Real_t*)NULL,
              M_ves, M_ves, 1, S->getPosition().begin()+3*M_ves*t /* target */,
              0, M_ves /* number of trgs per surface */,
              tmp_vel.begin()+3*M_ves*t);
        }

        axpy(1.0,tmp_vel,T_vel,T_vel); // accumulate
      }
    }

  }

  axpy(1.0,S_vel,T_vel,T_vel);
}



template<typename Surf_t>
void StokesVelocity<Surf_t>::u_ref(const Real_t* coord, int n, Real_t* out){ //Analytical velocity for sphere
  Real_t R0=1.0;
  for(int i=0;i<n;i++){
    const Real_t* c=&coord[i*3];
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
    out[i*3+0]=cos_t*ur-sin_t*ut-1.0;

    Real_t r_yz=sqrt(c[1]*c[1]+c[2]*c[2]);
    out[i*3+1]=(sin_t*ur+cos_t*ut)*c[1]/r_yz;
    out[i*3+2]=(sin_t*ur+cos_t*ut)*c[2]/r_yz;

    if(fabs(r-1.0)<1e-6){
      out[i*3+0]+=0.5;
      out[i*3+1]+=0.5;
      out[i*3+2]+=0.5;
    }
  }
}

template<typename Surf_t>
void StokesVelocity<Surf_t>::force(const Real_t* coord, int n, Real_t* out){ // Force on sphere
  Real_t R0=1.0;
  for(int i=0;i<n;i++){
    const Real_t* c=&coord[i*3];
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
    out[i*3+0]=cos_t*ur-sin_t*ut;

    Real_t r_yz=sqrt(c[1]*c[1]+c[2]*c[2]);
    out[i*3+1]=(sin_t*ur+cos_t*ut)*c[1]/r_yz;
    out[i*3+2]=(sin_t*ur+cos_t*ut)*c[2]/r_yz;
  }
}

// #include <legendre_rule.hpp>
// ABT: legendre_rule seems like code from pvfmm. Not sure why we need it

template<typename Surf_t>
void StokesVelocity<Surf_t>::Test(){

  // Set parameters
  Parameters<Real_t> sim_par;
  sim_par.sh_order = 16;
  sim_par.rep_up_freq = 16;

  // Reading operators from file
  bool readFromFile = true;
  Mats_t mats(readFromFile, sim_par);

  StokesVelocity<Surf_t> stokes_vel(mats, sim_par);

  //=================================================================//

  // Create vectors
  size_t nVec=10;
  Vec_t x0(nVec, sim_par.sh_order); // coordinates
  Vec_t f0(nVec, sim_par.sh_order); // force single layer
  Vec_t f1(nVec, sim_par.sh_order); // force double layer
  std::vector<Real_t> trg_coord; Real_t R0=1.0, R1=1.6;
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
    trg_coord.resize(3*fLen*nVec);
    size_t N=trg_coord.size()/3;
    for(size_t k=0;k<nVec;k++){
      Real_t R=R0+k*(R1-R0)/nVec;

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
            x_k[j+i*jmax]=R*cos_t;
            y_k[j+i*jmax]=R*sin_t*sin(j*2*M_PI/jmax);
            z_k[j+i*jmax]=R*sin_t*cos(j*2*M_PI/jmax);

            trg_coord[(k*fLen+i*jmax+j)*3+0]=x_k[j+i*jmax];
            trg_coord[(k*fLen+i*jmax+j)*3+1]=y_k[j+i*jmax];
            trg_coord[(k*fLen+i*jmax+j)*3+2]=z_k[j+i*jmax];
          }
          if(!k){ // Set fx, fy, fz
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
  Surf_t S(x0, mats);

  stokes_vel.SetSrcCoord(S);
  stokes_vel.SetTrgCoord(S);
  stokes_vel.SetDensitySL(&f0);
  stokes_vel.SetDensityDL(&f1);

  Real_t* v=stokes_vel();
  { // Compute error
    size_t N=trg_coord.size()/3;
    std::vector<Real_t> trg_vel(v,v+N*3);
    std::vector<Real_t> ref_vel(N*3);
    u_ref(&trg_coord[0], N, &ref_vel[0]);

    size_t M=15;
    std::vector<Real_t> max_err(M,0);
    std::vector<Real_t> max_vel(M,0);
    std::vector<Real_t> err(N*3);
    for(size_t i=0;i<err.size();i++){
      err[i]=trg_vel[i]-ref_vel[i];
    }
    for(size_t j=0;j<M;j++){
      size_t a=((j+0)*err.size())/M;
      size_t b=((j+1)*err.size())/M;
      for(size_t i=a;i<b;i++){
        if(fabs(err    [i])>max_err[j]) max_err[j]=fabs(err    [i]);
        if(fabs(ref_vel[i])>max_vel[j]) max_vel[j]=fabs(ref_vel[i]);
      }
      Real_t R=R0+(a+b)*(R1-R0)/N/6.0;
      max_err[j]=log(max_err[j]/max_vel[j])/log(10.0);
      max_vel[j]=R;
    }

    std::ios::fmtflags f(std::cout.flags());
    std::cout<<std::fixed<<std::setprecision(4)<<std::setiosflags(std::ios::left);
    for(size_t i=0;i<max_err.size();i++) std::cout<<std::setw(10)<<max_err[i]<<' '; std::cout<<'\n';
    for(size_t i=0;i<max_vel.size();i++) std::cout<<std::setw(10)<<max_vel[i]<<' '; std::cout<<'\n';
    std::cout.flags(f);

  }
}
