#include <cstring>
#include <parUtils.h>
#include <vector.hpp>
#include <mortonid.hpp>

#include <fmm_pts.hpp>
#include <fmm_node.hpp>
#include <fmm_tree.hpp>

#include "Logger.h"

template <class Real_t>
static Real_t machine_eps(){
  Real_t eps=1;
  while(eps*(Real_t)0.5+(Real_t)1.0>1.0) eps*=0.5;
  return eps;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////// Kernel Function Declarations ////////////////////////

////////// Stokes Kernel //////////

template <class Real_t, class Vec_t=Real_t, Vec_t (*RSQRT_INTRIN)(Vec_t)=pvfmm::rsqrt_intrin0<Vec_t> >
void stokes_sl_m2l_uKernel(pvfmm::Matrix<Real_t>& src_coord, pvfmm::Matrix<Real_t>& src_value, pvfmm::Matrix<Real_t>& trg_coord, pvfmm::Matrix<Real_t>& trg_value){
  #define SRC_BLK 500
  size_t VecLen=sizeof(Vec_t)/sizeof(Real_t);

  //// Number of newton iterations
  size_t NWTN_ITER=0;
  if(RSQRT_INTRIN==(Vec_t (*)(Vec_t))pvfmm::rsqrt_intrin0<Vec_t,Real_t>) NWTN_ITER=0;
  if(RSQRT_INTRIN==(Vec_t (*)(Vec_t))pvfmm::rsqrt_intrin1<Vec_t,Real_t>) NWTN_ITER=1;
  if(RSQRT_INTRIN==(Vec_t (*)(Vec_t))pvfmm::rsqrt_intrin2<Vec_t,Real_t>) NWTN_ITER=2;
  if(RSQRT_INTRIN==(Vec_t (*)(Vec_t))pvfmm::rsqrt_intrin3<Vec_t,Real_t>) NWTN_ITER=3;

  Real_t nwtn_scal=1; // scaling factor for newton iterations
  for(int i=0;i<NWTN_ITER;i++){
    nwtn_scal=2*nwtn_scal*nwtn_scal*nwtn_scal;
  }
  const Real_t OOEP = 1.0/(8*nwtn_scal*pvfmm::const_pi<Real_t>());
  Vec_t inv_nwtn_scal2=pvfmm::set_intrin<Vec_t,Real_t>(1.0/(nwtn_scal*nwtn_scal));

  size_t src_cnt_=src_coord.Dim(1);
  size_t trg_cnt_=trg_coord.Dim(1);
  for(size_t sblk=0;sblk<src_cnt_;sblk+=SRC_BLK){
    size_t src_cnt=src_cnt_-sblk;
    if(src_cnt>SRC_BLK) src_cnt=SRC_BLK;
    for(size_t t=0;t<trg_cnt_;t+=VecLen){
      Vec_t tx=pvfmm::load_intrin<Vec_t>(&trg_coord[0][t]);
      Vec_t ty=pvfmm::load_intrin<Vec_t>(&trg_coord[1][t]);
      Vec_t tz=pvfmm::load_intrin<Vec_t>(&trg_coord[2][t]);

      Vec_t tvx=pvfmm::zero_intrin<Vec_t>();
      Vec_t tvy=pvfmm::zero_intrin<Vec_t>();
      Vec_t tvz=pvfmm::zero_intrin<Vec_t>();
      for(size_t s=sblk;s<sblk+src_cnt;s++){
        Vec_t dx=pvfmm::sub_intrin(tx,pvfmm::bcast_intrin<Vec_t>(&src_coord[0][s]));
        Vec_t dy=pvfmm::sub_intrin(ty,pvfmm::bcast_intrin<Vec_t>(&src_coord[1][s]));
        Vec_t dz=pvfmm::sub_intrin(tz,pvfmm::bcast_intrin<Vec_t>(&src_coord[2][s]));

        Vec_t svx       =pvfmm::bcast_intrin<Vec_t>(&src_value[0][s]);
        Vec_t svy       =pvfmm::bcast_intrin<Vec_t>(&src_value[1][s]);
        Vec_t svz       =pvfmm::bcast_intrin<Vec_t>(&src_value[2][s]);
        Vec_t inner_prod=pvfmm::bcast_intrin<Vec_t>(&src_value[3][s]);

        Vec_t r2=               pvfmm::mul_intrin(dx,dx) ;
        r2=pvfmm::add_intrin(r2,pvfmm::mul_intrin(dy,dy));
        r2=pvfmm::add_intrin(r2,pvfmm::mul_intrin(dz,dz));

        Vec_t rinv=RSQRT_INTRIN(r2);
        Vec_t rinv2=pvfmm::mul_intrin(pvfmm::mul_intrin(rinv,rinv),inv_nwtn_scal2);

        inner_prod=pvfmm::add_intrin(inner_prod,pvfmm::mul_intrin(svx,dx));
        inner_prod=pvfmm::add_intrin(inner_prod,pvfmm::mul_intrin(svy,dy));
        inner_prod=pvfmm::add_intrin(inner_prod,pvfmm::mul_intrin(svz,dz));
        inner_prod=pvfmm::mul_intrin(inner_prod,rinv2);

        tvx=pvfmm::add_intrin(tvx,pvfmm::mul_intrin(rinv,pvfmm::add_intrin(svx,pvfmm::mul_intrin(dx,inner_prod))));
        tvy=pvfmm::add_intrin(tvy,pvfmm::mul_intrin(rinv,pvfmm::add_intrin(svy,pvfmm::mul_intrin(dy,inner_prod))));
        tvz=pvfmm::add_intrin(tvz,pvfmm::mul_intrin(rinv,pvfmm::add_intrin(svz,pvfmm::mul_intrin(dz,inner_prod))));
      }
      Vec_t ooep=pvfmm::set_intrin<Vec_t,Real_t>(OOEP);

      tvx=pvfmm::add_intrin(pvfmm::mul_intrin(tvx,ooep),pvfmm::load_intrin<Vec_t>(&trg_value[0][t]));
      tvy=pvfmm::add_intrin(pvfmm::mul_intrin(tvy,ooep),pvfmm::load_intrin<Vec_t>(&trg_value[1][t]));
      tvz=pvfmm::add_intrin(pvfmm::mul_intrin(tvz,ooep),pvfmm::load_intrin<Vec_t>(&trg_value[2][t]));

      pvfmm::store_intrin(&trg_value[0][t],tvx);
      pvfmm::store_intrin(&trg_value[1][t],tvy);
      pvfmm::store_intrin(&trg_value[2][t],tvz);
    }
  }

  { // Add FLOPS
    #ifndef __MIC__
    pvfmm::Profile::Add_FLOP((long long)trg_cnt_*(long long)src_cnt_*(29+4*(NWTN_ITER)));
    #endif
  }
  #undef SRC_BLK
}

template <class T, int newton_iter=0>
void stokes_sl_m2l(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* v_trg, pvfmm::mem::MemoryManager* mem_mgr){
  #define STK_KER_NWTN(nwtn) if(newton_iter==nwtn) \
        pvfmm::generic_kernel<Real_t, 4, 3, stokes_sl_m2l_uKernel<Real_t,Vec_t, pvfmm::rsqrt_intrin##nwtn<Vec_t,Real_t> > > \
            ((Real_t*)r_src, src_cnt, (Real_t*)v_src, dof, (Real_t*)r_trg, trg_cnt, (Real_t*)v_trg, mem_mgr)
  #define STOKES_KERNEL STK_KER_NWTN(0); STK_KER_NWTN(1); STK_KER_NWTN(2); STK_KER_NWTN(3);

  if(pvfmm::mem::TypeTraits<T>::ID()==pvfmm::mem::TypeTraits<float>::ID()){
    typedef float Real_t;
    #if defined __MIC__
      #define Vec_t Real_t
    #elif defined __AVX__
      #define Vec_t __m256
    #elif defined __SSE3__
      #define Vec_t __m128
    #else
      #define Vec_t Real_t
    #endif
    STOKES_KERNEL;
    #undef Vec_t
  }else if(pvfmm::mem::TypeTraits<T>::ID()==pvfmm::mem::TypeTraits<double>::ID()){
    typedef double Real_t;
    #if defined __MIC__
      #define Vec_t Real_t
    #elif defined __AVX__
      #define Vec_t __m256d
    #elif defined __SSE3__
      #define Vec_t __m128d
    #else
      #define Vec_t Real_t
    #endif
    STOKES_KERNEL;
    #undef Vec_t
  }else{
    typedef T Real_t;
    #define Vec_t Real_t
    STOKES_KERNEL;
    #undef Vec_t
  }

  #undef STK_KER_NWTN
  #undef STOKES_KERNEL
}

template <class T>
void stokes_m2l_vol_poten(const T* coord, int n, T* out){
  for(int i=0;i<n;i++){
    const T* c=&coord[i*COORD_DIM];
    T rx_2=c[1]*c[1]+c[2]*c[2];
    T ry_2=c[0]*c[0]+c[2]*c[2];
    T rz_2=c[0]*c[0]+c[1]*c[1];
    out[n*3*0+i*3+0]=-rx_2/6; out[n*3*0+i*3+1]=      0; out[n*3*0+i*3+2]=      0;
    out[n*3*1+i*3+0]=      0; out[n*3*1+i*3+1]=-ry_2/6; out[n*3*1+i*3+2]=      0;
    out[n*3*2+i*3+0]=      0; out[n*3*2+i*3+1]=      0; out[n*3*2+i*3+2]=-rz_2/6;
    out[n*3*3+i*3+0]= c[0]/6; out[n*3*3+i*3+1]= c[1]/6; out[n*3*3+i*3+2]= c[2]/6;
  }
}


template <class Real_t, class Vec_t=Real_t, Vec_t (*RSQRT_INTRIN)(Vec_t)=pvfmm::rsqrt_intrin0<Vec_t> >
void stokes_sl_uKernel(pvfmm::Matrix<Real_t>& src_coord, pvfmm::Matrix<Real_t>& src_value, pvfmm::Matrix<Real_t>& trg_coord, pvfmm::Matrix<Real_t>& trg_value){
  #define SRC_BLK 500
  static Real_t eps=machine_eps<Real_t>();
  size_t VecLen=sizeof(Vec_t)/sizeof(Real_t);

  //// Number of newton iterations
  size_t NWTN_ITER=0;
  if(RSQRT_INTRIN==(Vec_t (*)(Vec_t))pvfmm::rsqrt_intrin0<Vec_t,Real_t>) NWTN_ITER=0;
  if(RSQRT_INTRIN==(Vec_t (*)(Vec_t))pvfmm::rsqrt_intrin1<Vec_t,Real_t>) NWTN_ITER=1;
  if(RSQRT_INTRIN==(Vec_t (*)(Vec_t))pvfmm::rsqrt_intrin2<Vec_t,Real_t>) NWTN_ITER=2;
  if(RSQRT_INTRIN==(Vec_t (*)(Vec_t))pvfmm::rsqrt_intrin3<Vec_t,Real_t>) NWTN_ITER=3;

  Real_t nwtn_scal=1; // scaling factor for newton iterations
  for(int i=0;i<NWTN_ITER;i++){
    nwtn_scal=2*nwtn_scal*nwtn_scal*nwtn_scal;
  }
  const Real_t OOEP = 1.0/(8*nwtn_scal*pvfmm::const_pi<Real_t>());
  Vec_t inv_nwtn_scal2=pvfmm::set_intrin<Vec_t,Real_t>(1.0/(nwtn_scal*nwtn_scal));

  size_t src_cnt_=src_coord.Dim(1);
  size_t trg_cnt_=trg_coord.Dim(1);
  for(size_t sblk=0;sblk<src_cnt_;sblk+=SRC_BLK){
    size_t src_cnt=src_cnt_-sblk;
    if(src_cnt>SRC_BLK) src_cnt=SRC_BLK;
    for(size_t t=0;t<trg_cnt_;t+=VecLen){
      Vec_t tx=pvfmm::load_intrin<Vec_t>(&trg_coord[0][t]);
      Vec_t ty=pvfmm::load_intrin<Vec_t>(&trg_coord[1][t]);
      Vec_t tz=pvfmm::load_intrin<Vec_t>(&trg_coord[2][t]);

      Vec_t tvx=pvfmm::zero_intrin<Vec_t>();
      Vec_t tvy=pvfmm::zero_intrin<Vec_t>();
      Vec_t tvz=pvfmm::zero_intrin<Vec_t>();
      for(size_t s=sblk;s<sblk+src_cnt;s++){
        Vec_t dx=pvfmm::sub_intrin(tx,pvfmm::bcast_intrin<Vec_t>(&src_coord[0][s]));
        Vec_t dy=pvfmm::sub_intrin(ty,pvfmm::bcast_intrin<Vec_t>(&src_coord[1][s]));
        Vec_t dz=pvfmm::sub_intrin(tz,pvfmm::bcast_intrin<Vec_t>(&src_coord[2][s]));

        Vec_t svx=             pvfmm::bcast_intrin<Vec_t>(&src_value[0][s]) ;
        Vec_t svy=             pvfmm::bcast_intrin<Vec_t>(&src_value[1][s]) ;
        Vec_t svz=             pvfmm::bcast_intrin<Vec_t>(&src_value[2][s]) ;

        Vec_t r2=               pvfmm::mul_intrin(dx,dx) ;
        r2=pvfmm::add_intrin(r2,pvfmm::mul_intrin(dy,dy));
        r2=pvfmm::add_intrin(r2,pvfmm::mul_intrin(dz,dz));
        r2=pvfmm::and_intrin(pvfmm::cmplt_intrin(pvfmm::set_intrin<Vec_t,Real_t>(eps),r2),r2);

        Vec_t rinv=RSQRT_INTRIN(r2);
        Vec_t rinv2=pvfmm::mul_intrin(pvfmm::mul_intrin(rinv,rinv),inv_nwtn_scal2);

        Vec_t inner_prod=                       pvfmm::mul_intrin(svx,dx) ;
        inner_prod=pvfmm::add_intrin(inner_prod,pvfmm::mul_intrin(svy,dy));
        inner_prod=pvfmm::add_intrin(inner_prod,pvfmm::mul_intrin(svz,dz));
        inner_prod=pvfmm::mul_intrin(inner_prod,rinv2);

        tvx=pvfmm::add_intrin(tvx,pvfmm::mul_intrin(rinv,pvfmm::add_intrin(svx,pvfmm::mul_intrin(dx,inner_prod))));
        tvy=pvfmm::add_intrin(tvy,pvfmm::mul_intrin(rinv,pvfmm::add_intrin(svy,pvfmm::mul_intrin(dy,inner_prod))));
        tvz=pvfmm::add_intrin(tvz,pvfmm::mul_intrin(rinv,pvfmm::add_intrin(svz,pvfmm::mul_intrin(dz,inner_prod))));
      }
      Vec_t ooep=pvfmm::set_intrin<Vec_t,Real_t>(OOEP);

      tvx=pvfmm::add_intrin(pvfmm::mul_intrin(tvx,ooep),pvfmm::load_intrin<Vec_t>(&trg_value[0][t]));
      tvy=pvfmm::add_intrin(pvfmm::mul_intrin(tvy,ooep),pvfmm::load_intrin<Vec_t>(&trg_value[1][t]));
      tvz=pvfmm::add_intrin(pvfmm::mul_intrin(tvz,ooep),pvfmm::load_intrin<Vec_t>(&trg_value[2][t]));

      pvfmm::store_intrin(&trg_value[0][t],tvx);
      pvfmm::store_intrin(&trg_value[1][t],tvy);
      pvfmm::store_intrin(&trg_value[2][t],tvz);
    }
  }

  { // Add FLOPS
    #ifndef __MIC__
    pvfmm::Profile::Add_FLOP((long long)trg_cnt_*(long long)src_cnt_*(29+4*(NWTN_ITER)));
    #endif
  }
  #undef SRC_BLK
}

template <class T, int newton_iter=0>
void stokes_sl(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* v_trg, pvfmm::mem::MemoryManager* mem_mgr){
  #define STK_KER_NWTN(nwtn) if(newton_iter==nwtn) \
        pvfmm::generic_kernel<Real_t, 3, 3, stokes_sl_uKernel<Real_t,Vec_t, pvfmm::rsqrt_intrin##nwtn<Vec_t,Real_t> > > \
            ((Real_t*)r_src, src_cnt, (Real_t*)v_src, dof, (Real_t*)r_trg, trg_cnt, (Real_t*)v_trg, mem_mgr)
  #define STOKES_KERNEL STK_KER_NWTN(0); STK_KER_NWTN(1); STK_KER_NWTN(2); STK_KER_NWTN(3);

  if(pvfmm::mem::TypeTraits<T>::ID()==pvfmm::mem::TypeTraits<float>::ID()){
    typedef float Real_t;
    #if defined __MIC__
      #define Vec_t Real_t
    #elif defined __AVX__
      #define Vec_t __m256
    #elif defined __SSE3__
      #define Vec_t __m128
    #else
      #define Vec_t Real_t
    #endif
    STOKES_KERNEL;
    #undef Vec_t
  }else if(pvfmm::mem::TypeTraits<T>::ID()==pvfmm::mem::TypeTraits<double>::ID()){
    typedef double Real_t;
    #if defined __MIC__
      #define Vec_t Real_t
    #elif defined __AVX__
      #define Vec_t __m256d
    #elif defined __SSE3__
      #define Vec_t __m128d
    #else
      #define Vec_t Real_t
    #endif
    STOKES_KERNEL;
    #undef Vec_t
  }else{
    typedef T Real_t;
    #define Vec_t Real_t
    STOKES_KERNEL;
    #undef Vec_t
  }

  #undef STK_KER_NWTN
  #undef STOKES_KERNEL
}

template <class Real_t, class Vec_t=Real_t, Vec_t (*RSQRT_INTRIN)(Vec_t)=pvfmm::rsqrt_intrin0<Vec_t> >
void stokes_dl_uKernel(pvfmm::Matrix<Real_t>& src_coord, pvfmm::Matrix<Real_t>& src_value, pvfmm::Matrix<Real_t>& trg_coord, pvfmm::Matrix<Real_t>& trg_value){
  #define SRC_BLK 500
  static Real_t eps=machine_eps<Real_t>();
  size_t VecLen=sizeof(Vec_t)/sizeof(Real_t);

  //// Number of newton iterations
  size_t NWTN_ITER=0;
  if(RSQRT_INTRIN==(Vec_t (*)(Vec_t))pvfmm::rsqrt_intrin0<Vec_t,Real_t>) NWTN_ITER=0;
  if(RSQRT_INTRIN==(Vec_t (*)(Vec_t))pvfmm::rsqrt_intrin1<Vec_t,Real_t>) NWTN_ITER=1;
  if(RSQRT_INTRIN==(Vec_t (*)(Vec_t))pvfmm::rsqrt_intrin2<Vec_t,Real_t>) NWTN_ITER=2;
  if(RSQRT_INTRIN==(Vec_t (*)(Vec_t))pvfmm::rsqrt_intrin3<Vec_t,Real_t>) NWTN_ITER=3;

  Real_t nwtn_scal=1; // scaling factor for newton iterations
  for(int i=0;i<NWTN_ITER;i++){
    nwtn_scal=2*nwtn_scal*nwtn_scal*nwtn_scal;
  }
  const Real_t SCAL_CONST = 3.0/(4.0*nwtn_scal*nwtn_scal*nwtn_scal*nwtn_scal*nwtn_scal*pvfmm::const_pi<Real_t>());

  size_t src_cnt_=src_coord.Dim(1);
  size_t trg_cnt_=trg_coord.Dim(1);
  for(size_t sblk=0;sblk<src_cnt_;sblk+=SRC_BLK){
    size_t src_cnt=src_cnt_-sblk;
    if(src_cnt>SRC_BLK) src_cnt=SRC_BLK;
    for(size_t t=0;t<trg_cnt_;t+=VecLen){
      Vec_t tx=pvfmm::load_intrin<Vec_t>(&trg_coord[0][t]);
      Vec_t ty=pvfmm::load_intrin<Vec_t>(&trg_coord[1][t]);
      Vec_t tz=pvfmm::load_intrin<Vec_t>(&trg_coord[2][t]);

      Vec_t tvx=pvfmm::zero_intrin<Vec_t>();
      Vec_t tvy=pvfmm::zero_intrin<Vec_t>();
      Vec_t tvz=pvfmm::zero_intrin<Vec_t>();
      for(size_t s=sblk;s<sblk+src_cnt;s++){
        Vec_t dx=pvfmm::sub_intrin(tx,pvfmm::bcast_intrin<Vec_t>(&src_coord[0][s]));
        Vec_t dy=pvfmm::sub_intrin(ty,pvfmm::bcast_intrin<Vec_t>(&src_coord[1][s]));
        Vec_t dz=pvfmm::sub_intrin(tz,pvfmm::bcast_intrin<Vec_t>(&src_coord[2][s]));

        Vec_t snx=pvfmm::bcast_intrin<Vec_t>(&src_value[0][s]) ;
        Vec_t sny=pvfmm::bcast_intrin<Vec_t>(&src_value[1][s]) ;
        Vec_t snz=pvfmm::bcast_intrin<Vec_t>(&src_value[2][s]) ;

        Vec_t svx=pvfmm::bcast_intrin<Vec_t>(&src_value[3][s]) ;
        Vec_t svy=pvfmm::bcast_intrin<Vec_t>(&src_value[4][s]) ;
        Vec_t svz=pvfmm::bcast_intrin<Vec_t>(&src_value[5][s]) ;

        Vec_t r2=               pvfmm::mul_intrin(dx,dx) ;
        r2=pvfmm::add_intrin(r2,pvfmm::mul_intrin(dy,dy));
        r2=pvfmm::add_intrin(r2,pvfmm::mul_intrin(dz,dz));
        r2=pvfmm::and_intrin(pvfmm::cmplt_intrin(pvfmm::set_intrin<Vec_t,Real_t>(eps),r2),r2);

        Vec_t rinv=RSQRT_INTRIN(r2);
        Vec_t rinv2=pvfmm::mul_intrin(rinv ,rinv );
        Vec_t rinv5=pvfmm::mul_intrin(pvfmm::mul_intrin(rinv2,rinv2),rinv);

        Vec_t r_dot_n=                    pvfmm::mul_intrin(snx,dx) ;
        r_dot_n=pvfmm::add_intrin(r_dot_n,pvfmm::mul_intrin(sny,dy));
        r_dot_n=pvfmm::add_intrin(r_dot_n,pvfmm::mul_intrin(snz,dz));

        Vec_t r_dot_f=                    pvfmm::mul_intrin(svx,dx) ;
        r_dot_f=pvfmm::add_intrin(r_dot_f,pvfmm::mul_intrin(svy,dy));
        r_dot_f=pvfmm::add_intrin(r_dot_f,pvfmm::mul_intrin(svz,dz));

        Vec_t p=pvfmm::mul_intrin(pvfmm::mul_intrin(r_dot_n,r_dot_f),rinv5);
        tvx=pvfmm::add_intrin(tvx,pvfmm::mul_intrin(dx,p));
        tvy=pvfmm::add_intrin(tvy,pvfmm::mul_intrin(dy,p));
        tvz=pvfmm::add_intrin(tvz,pvfmm::mul_intrin(dz,p));
      }
      Vec_t scal_const=pvfmm::set_intrin<Vec_t,Real_t>(SCAL_CONST);

      tvx=pvfmm::add_intrin(pvfmm::mul_intrin(tvx,scal_const),pvfmm::load_intrin<Vec_t>(&trg_value[0][t]));
      tvy=pvfmm::add_intrin(pvfmm::mul_intrin(tvy,scal_const),pvfmm::load_intrin<Vec_t>(&trg_value[1][t]));
      tvz=pvfmm::add_intrin(pvfmm::mul_intrin(tvz,scal_const),pvfmm::load_intrin<Vec_t>(&trg_value[2][t]));

      pvfmm::store_intrin(&trg_value[0][t],tvx);
      pvfmm::store_intrin(&trg_value[1][t],tvy);
      pvfmm::store_intrin(&trg_value[2][t],tvz);
    }
  }

  { // Add FLOPS
    #ifndef __MIC__
    pvfmm::Profile::Add_FLOP((long long)trg_cnt_*(long long)src_cnt_*(31+4*(NWTN_ITER)));
    #endif
  }
  #undef SRC_BLK
}

template <class T, int newton_iter=0>
void stokes_dl(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* v_trg, pvfmm::mem::MemoryManager* mem_mgr){
  #define STK_KER_NWTN(nwtn) if(newton_iter==nwtn) \
        pvfmm::generic_kernel<Real_t, 6, 3, stokes_dl_uKernel<Real_t,Vec_t, pvfmm::rsqrt_intrin##nwtn<Vec_t,Real_t> > > \
            ((Real_t*)r_src, src_cnt, (Real_t*)v_src, dof, (Real_t*)r_trg, trg_cnt, (Real_t*)v_trg, mem_mgr)
  #define STOKES_KERNEL STK_KER_NWTN(0); STK_KER_NWTN(1); STK_KER_NWTN(2); STK_KER_NWTN(3);

  if(pvfmm::mem::TypeTraits<T>::ID()==pvfmm::mem::TypeTraits<float>::ID()){
    typedef float Real_t;
    #if defined __MIC__
      #define Vec_t Real_t
    #elif defined __AVX__
      #define Vec_t __m256
    #elif defined __SSE3__
      #define Vec_t __m128
    #else
      #define Vec_t Real_t
    #endif
    STOKES_KERNEL;
    #undef Vec_t
  }else if(pvfmm::mem::TypeTraits<T>::ID()==pvfmm::mem::TypeTraits<double>::ID()){
    typedef double Real_t;
    #if defined __MIC__
      #define Vec_t Real_t
    #elif defined __AVX__
      #define Vec_t __m256d
    #elif defined __SSE3__
      #define Vec_t __m128d
    #else
      #define Vec_t Real_t
    #endif
    STOKES_KERNEL;
    #undef Vec_t
  }else{
    typedef T Real_t;
    #define Vec_t Real_t
    STOKES_KERNEL;
    #undef Vec_t
  }

  #undef STK_KER_NWTN
  #undef STOKES_KERNEL
}

template <class T>
void stokes_vol_poten(const T* coord, int n, T* out){
  for(int i=0;i<n;i++){
    const T* c=&coord[i*COORD_DIM];
    T rx_2=c[1]*c[1]+c[2]*c[2];
    T ry_2=c[0]*c[0]+c[2]*c[2];
    T rz_2=c[0]*c[0]+c[1]*c[1];
    out[n*3*0+i*3+0]=-rx_2/6; out[n*3*0+i*3+1]=      0; out[n*3*0+i*3+2]=      0;
    out[n*3*1+i*3+0]=      0; out[n*3*1+i*3+1]=-ry_2/6; out[n*3*1+i*3+2]=      0;
    out[n*3*2+i*3+0]=      0; out[n*3*2+i*3+1]=      0; out[n*3*2+i*3+2]=-rz_2/6;
  }
}


template <class Real_t>
inline const pvfmm::Kernel<Real_t>& StokesKernel<Real_t>::Kernel(){

  static const pvfmm::Kernel<Real_t> ker_m2l=pvfmm::BuildKernel<Real_t, stokes_sl_m2l<Real_t,2>                      >("stokes_m2l", 3, std::pair<int,int>(4,3),
      NULL,NULL,NULL,     NULL,    NULL,    NULL, NULL,NULL, stokes_m2l_vol_poten);

  static const pvfmm::Kernel<Real_t> ker    =pvfmm::BuildKernel<Real_t, stokes_sl    <Real_t,2>, stokes_dl<Real_t,2> >("stokes_vel", 3, std::pair<int,int>(3,3),
      NULL,NULL,NULL, &ker_m2l,&ker_m2l,&ker_m2l, NULL,NULL, stokes_vol_poten    );

  return ker;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


template<typename T>
void PVFMMBoundingBox(size_t n_src, const T* x, T* scale_xr, T* shift_xr, MPI_Comm comm){
  T& scale_x=*scale_xr;
  T* shift_x= shift_xr;

  if(n_src>0){ // Compute bounding box
    double loc_min_x[COORD_DIM];
    double loc_max_x[COORD_DIM];
    assert(n_src>0);
    for(size_t k=0;k<COORD_DIM;k++){
      loc_min_x[k]=loc_max_x[k]=x[k];
    }

    for(size_t i=0;i<n_src;i++){
      const T* x_=&x[i*COORD_DIM];
      for(size_t k=0;k<COORD_DIM;k++){
        if(loc_min_x[k]>x_[0]) loc_min_x[k]=x_[0];
        if(loc_max_x[k]<x_[0]) loc_max_x[k]=x_[0];
        ++x_;
      }
    }

    double min_x[COORD_DIM];
    double max_x[COORD_DIM];
    MPI_Allreduce(loc_min_x, min_x, COORD_DIM, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(loc_max_x, max_x, COORD_DIM, MPI_DOUBLE, MPI_MAX, comm);

    static T eps=machine_eps<T>()*64; // Points should be well within the box.
    scale_x=1.0/(max_x[0]-min_x[0]+2*eps);
    for(size_t k=0;k<COORD_DIM;k++){
      scale_x=std::min(scale_x,(T)(1.0/(max_x[k]-min_x[k]+2*eps)));
    }
    if(scale_x*0.0!=0.0) scale_x=1.0; // fix for scal_x=inf
    for(size_t k=0;k<COORD_DIM;k++){
      shift_x[k]=-min_x[k]*scale_x+eps;
    }
  }
}

template<typename T>
struct PVFMMContext{
  typedef pvfmm::FMM_Node<pvfmm::MPI_Node<T> > Node_t;
  typedef pvfmm::FMM_Pts<Node_t> Mat_t;
  typedef pvfmm::FMM_Tree<Mat_t> Tree_t;

  T box_size;
  int max_pts;
  int mult_order;
  int max_depth;
  pvfmm::BoundaryType bndry;
  const pvfmm::Kernel<T>* ker;
  MPI_Comm comm;

  typename Node_t::NodeData tree_data;
  Tree_t* tree;
  Mat_t* mat;
};

template<typename T>
void* PVFMMCreateContext(T box_size, int n, int m, int max_d,
    const pvfmm::Kernel<T>* ker,
    MPI_Comm comm){
  pvfmm::Profile::Tic("FMMContext",&comm,true);
  bool prof_state=pvfmm::Profile::Enable(false);

  // Create new context.
  PVFMMContext<T>* ctx=new PVFMMContext<T>;

  // Set member variables.
  ctx->box_size=box_size;
  ctx->max_pts=n;
  ctx->mult_order=m;
  ctx->max_depth=max_d;
  ctx->bndry=(box_size<=0?pvfmm::FreeSpace:pvfmm::Periodic);
  ctx->ker=ker;
  ctx->comm=comm;

  // Initialize FMM matrices.
  ctx->mat=new typename PVFMMContext<T>::Mat_t();
  ctx->mat->Initialize(ctx->mult_order, ctx->comm, ctx->ker);

  // Set tree_data
  ctx->tree_data.dim=COORD_DIM;
  ctx->tree_data.max_depth=ctx->max_depth;
  ctx->tree_data.max_pts=ctx->max_pts;
  { // ctx->tree_data.pt_coord=... //Set points for initial tree.
    int np, myrank;
    MPI_Comm_size(ctx->comm, &np);
    MPI_Comm_rank(ctx->comm, &myrank);

    std::vector<T> coord;
    size_t NN=ceil(pow((T)np*ctx->max_pts,1.0/3.0));
    size_t N_total=NN*NN*NN;
    size_t start= myrank   *N_total/np;
    size_t end  =(myrank+1)*N_total/np;
    for(size_t i=start;i<end;i++){
      coord.push_back(((T)((i/  1    )%NN)+0.5)/NN);
      coord.push_back(((T)((i/ NN    )%NN)+0.5)/NN);
      coord.push_back(((T)((i/(NN*NN))%NN)+0.5)/NN);
    }
    ctx->tree_data.pt_coord=coord;
  }

  // Construct tree.
  bool adap=false; // no data to do adaptive.
  ctx->tree=new typename PVFMMContext<T>::Tree_t(comm);
  ctx->tree->Initialize(&ctx->tree_data);
  ctx->tree->InitFMM_Tree(adap,ctx->bndry);

  pvfmm::Profile::Enable(prof_state);
  pvfmm::Profile::Toc();
  return ctx;
}

template<typename T>
void PVFMMDestroyContext(void** ctx){
  if(!ctx[0]) return;

  // Delete tree.
  delete ((PVFMMContext<T>*)ctx[0])->tree;

  // Delete matrices.
  delete ((PVFMMContext<T>*)ctx[0])->mat;

  // Delete context.
  delete (PVFMMContext<T>*)ctx[0];
  ctx[0]=NULL;
}

template<typename T>
void PVFMMEval(const T* src_pos, const T* sl_den, size_t n_src, T* trg_vel, void** ctx_){
  PVFMMEval<T>(src_pos, sl_den, NULL, n_src, trg_vel, ctx_);
}

template<typename T>
void PVFMMEval(const T* src_pos, const T* sl_den, const T* dl_den, size_t n_src, T* trg_vel, void** ctx_, int setup){
  PVFMMEval(src_pos, sl_den, dl_den, n_src, src_pos, trg_vel, n_src, ctx_, setup);
}

template<typename T>
void PVFMMEval(const T* src_pos, const T* sl_den, const T* dl_den, size_t n_src, const T* trg_pos, T* trg_vel, size_t n_trg, void** ctx_, int setup){
  PROFILESTART();
  long long prof_FLOPS=pvfmm::Profile::Add_FLOP(0);
  size_t omp_p=omp_get_max_threads();

  typedef pvfmm::FMM_Node<pvfmm::MPI_Node<T> > Node_t;
  typedef pvfmm::FMM_Pts<Node_t> Mat_t;
  typedef pvfmm::FMM_Tree<Mat_t> Tree_t;

  if(!ctx_[0]) ctx_[0]=PVFMMCreateContext<T>();
  PVFMMContext<T>* ctx=(PVFMMContext<T>*)ctx_[0];
  const int* ker_dim=ctx->ker->ker_dim;

  pvfmm::Profile::Tic("FMM",&ctx->comm);
  T scale_x, shift_x[COORD_DIM];
  if(ctx->box_size<=0){ // determine bounding box
    T s0, x0[COORD_DIM];
    T s1, x1[COORD_DIM];
    PVFMMBoundingBox(n_src, src_pos, &s0, x0, ctx->comm);
    PVFMMBoundingBox(n_trg, trg_pos, &s1, x1, ctx->comm);

    T c0[COORD_DIM]={(0.5-x0[0])/s0, (0.5-x0[1])/s0, (0.5-x0[2])/s0};
    T c1[COORD_DIM]={(0.5-x1[0])/s1, (0.5-x1[1])/s1, (0.5-x1[2])/s1};

    scale_x=0;
    scale_x=std::max(scale_x, fabs(c0[0]-c1[0]));
    scale_x=std::max(scale_x, fabs(c0[1]-c1[1]));
    scale_x=std::max(scale_x, fabs(c0[2]-c1[2]));
    scale_x=1.0/(scale_x+1/s0+1/s1);

    shift_x[0]=0.5-(c0[0]+c1[0])*scale_x/2.0;
    shift_x[1]=0.5-(c0[1]+c1[1])*scale_x/2.0;
    shift_x[2]=0.5-(c0[2]+c1[2])*scale_x/2.0;
  }else{
    scale_x=1.0/ctx->box_size;
    shift_x[0]=0;
    shift_x[1]=0;
    shift_x[2]=0;
  }

  pvfmm::Vector<T>  src_scal;
  pvfmm::Vector<T>  trg_scal;
  pvfmm::Vector<T> surf_scal;
  { // Set src_scal, trg_scal
    pvfmm::Vector<T>& src_scal_exp=ctx->ker->src_scal;
    pvfmm::Vector<T>& trg_scal_exp=ctx->ker->trg_scal;
    src_scal .ReInit(ctx->ker->src_scal.Dim());
    trg_scal .ReInit(ctx->ker->trg_scal.Dim());
    surf_scal.ReInit(COORD_DIM+src_scal.Dim());
    for(size_t i=0;i<src_scal.Dim();i++){
      src_scal [i]=pow(scale_x, src_scal_exp[i]);
      surf_scal[i]=scale_x*src_scal[i];
    }
    for(size_t i=0;i<trg_scal.Dim();i++){
      trg_scal[i]=pow(scale_x, trg_scal_exp[i]);
    }
    for(size_t i=src_scal.Dim();i<surf_scal.Dim();i++){
      surf_scal[i]=1;
    }
  }

  pvfmm::Vector<size_t> scatter_index;
  { // Set tree_data
    pvfmm::Vector<T>&  trg_coord=ctx->tree_data. trg_coord;
    pvfmm::Vector<T>&  src_coord=ctx->tree_data. src_coord;
    pvfmm::Vector<T>&  src_value=ctx->tree_data. src_value;
    pvfmm::Vector<T>& surf_value=ctx->tree_data.surf_value;
    pvfmm::Vector<pvfmm::MortonId> pt_mid;

    std::vector<Node_t*> nodes;
    { // Get list of leaf nodes.
      std::vector<Node_t*>& all_nodes=ctx->tree->GetNodeList();
      for(size_t i=0;i<all_nodes.size();i++){
        if(all_nodes[i]->IsLeaf() && !all_nodes[i]->IsGhost()){
          nodes.push_back(all_nodes[i]);
        }
      }
    }

    pvfmm::MortonId min_mid;
    { // Get first MortonId
      Node_t* n=ctx->tree->PreorderFirst();
      while(n!=NULL){
        if(!n->IsGhost() && n->IsLeaf()) break;
        n=ctx->tree->PreorderNxt(n);
      }
      assert(n!=NULL);
      min_mid=n->GetMortonId();
    }

    { // Set src tree_data
      { // Scatter src data
        // Compute MortonId and copy coordinates and values.
        src_coord .ReInit(       n_src            *COORD_DIM);
        src_value .ReInit(sl_den?n_src*(ker_dim[0]          ):0);
        surf_value.ReInit(dl_den?n_src*(ker_dim[0]+COORD_DIM):0);
        pt_mid    .ReInit(n_src);
        #pragma omp parallel for
        for(size_t tid=0;tid<omp_p;tid++){
          size_t a=((tid+0)*n_src)/omp_p;
          size_t b=((tid+1)*n_src)/omp_p;
          for(size_t i=a;i<b;i++){
            for(size_t j=0;j<COORD_DIM;j++){
              src_coord[i*COORD_DIM+j]=src_pos[i*COORD_DIM+j]*scale_x+shift_x[j];
              while(src_coord[i*COORD_DIM+j]< 0.0) src_coord[i*COORD_DIM+j]+=1.0;
              while(src_coord[i*COORD_DIM+j]>=1.0) src_coord[i*COORD_DIM+j]-=1.0;
            }
            pt_mid[i]=pvfmm::MortonId(&src_coord[i*COORD_DIM]);
          }
          if(src_value.Dim()) for(size_t i=a;i<b;i++){
            for(size_t j=0;j<ker_dim[0];j++){
              src_value[i*ker_dim[0]+j]=sl_den[i*ker_dim[0]+j]*src_scal[j];
            }
          }
          if(surf_value.Dim()) for(size_t i=a;i<b;i++){
            for(size_t j=0;j<ker_dim[0]+COORD_DIM;j++){
              surf_value[i*(ker_dim[0]+COORD_DIM)+j]=dl_den[i*(ker_dim[0]+COORD_DIM)+j]*surf_scal[j];
            }
          }
        }

        // Scatter src coordinates and values.
        pvfmm::par::SortScatterIndex( pt_mid  , scatter_index, ctx->comm, &min_mid);
        pvfmm::par::ScatterForward  ( pt_mid  , scatter_index, ctx->comm);
        pvfmm::par::ScatterForward  (src_coord, scatter_index, ctx->comm);
        if( src_value.Dim()) pvfmm::par::ScatterForward( src_value, scatter_index, ctx->comm);
        if(surf_value.Dim()) pvfmm::par::ScatterForward(surf_value, scatter_index, ctx->comm);
      }
      { // Set src tree_data
        std::vector<size_t> part_indx(nodes.size()+1);
        part_indx[nodes.size()]=pt_mid.Dim();
        #pragma omp parallel for
        for(size_t j=0;j<nodes.size();j++){
          part_indx[j]=std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), nodes[j]->GetMortonId())-&pt_mid[0];
        }

        if(setup){
          #pragma omp parallel for
          for(size_t j=0;j<nodes.size();j++){
            size_t n_pts=part_indx[j+1]-part_indx[j];
            if(src_value.Dim()){
              nodes[j]-> src_coord.ReInit(n_pts*( COORD_DIM),& src_coord[0]+part_indx[j]*( COORD_DIM),false);
              nodes[j]-> src_value.ReInit(n_pts*(ker_dim[0]),& src_value[0]+part_indx[j]*(ker_dim[0]),false);
            }else{
              nodes[j]-> src_coord.ReInit(0,NULL,false);
              nodes[j]-> src_value.ReInit(0,NULL,false);
            }
            if(surf_value.Dim()){
              nodes[j]->surf_coord.ReInit(n_pts*(           COORD_DIM),& src_coord[0]+part_indx[j]*(           COORD_DIM),false);
              nodes[j]->surf_value.ReInit(n_pts*(ker_dim[0]+COORD_DIM),&surf_value[0]+part_indx[j]*(ker_dim[0]+COORD_DIM),false);
            }else{
              nodes[j]->surf_coord.ReInit(0,NULL,false);
              nodes[j]->surf_value.ReInit(0,NULL,false);
            }
          }
        }else{
          #pragma omp parallel for
          for(size_t j=0;j<nodes.size();j++){
            size_t n_pts=part_indx[j+1]-part_indx[j];
            if(src_value.Dim()){
              assert(nodes[j]->src_coord.Dim()==n_pts*( COORD_DIM));
              assert(nodes[j]->src_value.Dim()==n_pts*(ker_dim[0]));
              //memcpy(&nodes[j]->src_coord[0],&src_coord[0]+part_indx[j]*( COORD_DIM),n_pts*( COORD_DIM)*sizeof(T));
              memcpy(&nodes[j]->src_value[0],&src_value[0]+part_indx[j]*(ker_dim[0]),n_pts*(ker_dim[0])*sizeof(T));
            }
            if(surf_value.Dim()){
              assert(nodes[j]->surf_coord.Dim()==n_pts*(           COORD_DIM));
              assert(nodes[j]->surf_value.Dim()==n_pts*(ker_dim[0]+COORD_DIM));
              //memcpy(&nodes[j]->surf_coord[0],& src_coord[0]+part_indx[j]*(           COORD_DIM),n_pts*(           COORD_DIM)*sizeof(T));
              memcpy(&nodes[j]->surf_value[0],&surf_value[0]+part_indx[j]*(ker_dim[0]+COORD_DIM),n_pts*(ker_dim[0]+COORD_DIM)*sizeof(T));
            }
          }
        }
      }
    }
    { // Set trg tree_data
      if(trg_pos==src_pos && n_src==n_trg){ // Scatter trg data
        trg_coord.ReInit(src_coord.Dim(),&src_coord[0],false);
      }else{
        // Compute MortonId and copy coordinates.
        trg_coord.Resize(n_trg*COORD_DIM);
        pt_mid    .ReInit(n_trg);
        #pragma omp parallel for
        for(size_t tid=0;tid<omp_p;tid++){
          size_t a=((tid+0)*n_trg)/omp_p;
          size_t b=((tid+1)*n_trg)/omp_p;
          for(size_t i=a;i<b;i++){
            for(size_t j=0;j<COORD_DIM;j++){
              trg_coord[i*COORD_DIM+j]=trg_pos[i*COORD_DIM+j]*scale_x+shift_x[j];
              while(trg_coord[i*COORD_DIM+j]< 0.0) trg_coord[i*COORD_DIM+j]+=1.0;
              while(trg_coord[i*COORD_DIM+j]>=1.0) trg_coord[i*COORD_DIM+j]-=1.0;
            }
            pt_mid[i]=pvfmm::MortonId(&trg_coord[i*COORD_DIM]);
          }
        }

        // Scatter trg coordinates.
        pvfmm::par::SortScatterIndex( pt_mid  , scatter_index, ctx->comm, &min_mid);
        pvfmm::par::ScatterForward  ( pt_mid  , scatter_index, ctx->comm);
        pvfmm::par::ScatterForward  (trg_coord, scatter_index, ctx->comm);
      }
      { // Set trg tree_data
        std::vector<size_t> part_indx(nodes.size()+1);
        part_indx[nodes.size()]=pt_mid.Dim();
        #pragma omp parallel for
        for(size_t j=0;j<nodes.size();j++){
          part_indx[j]=std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), nodes[j]->GetMortonId())-&pt_mid[0];
        }

        if(setup){
          #pragma omp parallel for
          for(size_t j=0;j<nodes.size();j++){
            size_t n_pts=part_indx[j+1]-part_indx[j];
            {
              nodes[j]-> trg_coord.ReInit(n_pts*(COORD_DIM),& trg_coord[0]+part_indx[j]*(COORD_DIM),false);
            }
          }
        }else{
          #pragma omp parallel for
          for(size_t j=0;j<nodes.size();j++){
            size_t n_pts=part_indx[j+1]-part_indx[j];
            {
              assert(nodes[j]->trg_coord.Dim()==n_pts*(COORD_DIM));
              //memcpy(&nodes[j]->trg_coord[0],&trg_coord[0]+part_indx[j]*(COORD_DIM),n_pts*(COORD_DIM)*sizeof(T));
            }
          }
        }
      }
    }
  }

  if(setup){ // Optional stuff (redistribute, adaptive refine ...)
    { //Output max tree depth.
      int np, myrank;
      MPI_Comm_size(ctx->comm, &np);
      MPI_Comm_rank(ctx->comm, &myrank);

      long nleaf=0, maxdepth=0;
      std::vector<size_t> all_nodes(MAX_DEPTH+1,0);
      std::vector<size_t> leaf_nodes(MAX_DEPTH+1,0);
      std::vector<Node_t*>& nodes=ctx->tree->GetNodeList();
      for(size_t i=0;i<nodes.size();i++){
        Node_t* n=nodes[i];
        if(!n->IsGhost()) all_nodes[n->Depth()]++;
        if(!n->IsGhost() && n->IsLeaf()){
          leaf_nodes[n->Depth()]++;
          if(maxdepth<n->Depth()) maxdepth=n->Depth();
          nleaf++;
        }
      }

      std::stringstream os1,os2;
      os1<<"All Nodes";
      for(int i=0;i<MAX_DEPTH;i++){
        int local_size=all_nodes[i];
        int global_size;
        MPI_Allreduce(&local_size, &global_size, 1, MPI_INT, MPI_SUM, ctx->comm);
        os1<<global_size<<' ';
      }
      if(!myrank) COUTDEBUG(os1.str());

      os2<<"Leaf Nodes: ";
      for(int i=0;i<MAX_DEPTH;i++){
        int local_size=leaf_nodes[i];
        int global_size;
        MPI_Allreduce(&local_size, &global_size, 1, MPI_INT, MPI_SUM, ctx->comm);
        os2<<global_size<<' ';
      }
      if(!myrank) COUTDEBUG(os2.str());

      long nleaf_glb=0, maxdepth_glb=0;
      { // MPI_Reduce
        MPI_Allreduce(&nleaf, &nleaf_glb, 1, MPI_INT, MPI_SUM, ctx->comm);
        MPI_Allreduce(&maxdepth, &maxdepth_glb, 1, MPI_INT, MPI_MAX, ctx->comm);
      }
      if(!myrank) COUTDEBUG("Number of Leaf Nodes: "<<nleaf_glb);
      if(!myrank) COUTDEBUG("Tree Depth: "<<maxdepth_glb);
    }
    bool adap=true;
    ctx->tree->InitFMM_Tree(adap,ctx->bndry);
    { //Output max tree depth.
      int np, myrank;
      MPI_Comm_size(ctx->comm, &np);
      MPI_Comm_rank(ctx->comm, &myrank);

      long nleaf=0, maxdepth=0;
      std::vector<size_t> all_nodes(MAX_DEPTH+1,0);
      std::vector<size_t> leaf_nodes(MAX_DEPTH+1,0);
      std::vector<Node_t*>& nodes=ctx->tree->GetNodeList();
      for(size_t i=0;i<nodes.size();i++){
        Node_t* n=nodes[i];
        if(!n->IsGhost()) all_nodes[n->Depth()]++;
        if(!n->IsGhost() && n->IsLeaf()){
          leaf_nodes[n->Depth()]++;
          if(maxdepth<n->Depth()) maxdepth=n->Depth();
          nleaf++;
        }
      }

      std::stringstream os1,os2;
      os1<<"All  Nodes: ";
      for(int i=0;i<MAX_DEPTH;i++){
        int local_size=all_nodes[i];
        int global_size;
        MPI_Allreduce(&local_size, &global_size, 1, MPI_INT, MPI_SUM, ctx->comm);
        os1<<global_size<<' ';
      }
      if(!myrank) COUTDEBUG(os1.str());

      os2<<"Leaf Nodes: ";
      for(int i=0;i<MAX_DEPTH;i++){
        int local_size=leaf_nodes[i];
        int global_size;
        MPI_Allreduce(&local_size, &global_size, 1, MPI_INT, MPI_SUM, ctx->comm);
        os2<<global_size<<' ';
      }
      if(!myrank) COUTDEBUG(os2.str());

      long nleaf_glb=0, maxdepth_glb=0;
      { // MPI_Reduce
        MPI_Allreduce(&nleaf, &nleaf_glb, 1, MPI_INT, MPI_SUM, ctx->comm);
        MPI_Allreduce(&maxdepth, &maxdepth_glb, 1, MPI_INT, MPI_MAX, ctx->comm);
      }
      if(!myrank) COUTDEBUG("Number of Leaf Nodes: "<<nleaf_glb);
      if(!myrank) COUTDEBUG("Tree Depth: "<<maxdepth_glb);
    }
  }

  // Setup tree for FMM.
  if(setup) ctx->tree->SetupFMM(ctx->mat);
  else ctx->tree->ClearFMMData();
  ctx->tree->RunFMM();

  //Write2File
  //ctx->tree->Write2File("output",0);

  { // Get target potential.
    pvfmm::Vector<T> trg_value;
    { // Get trg data.
      Node_t* n=NULL;
      n=ctx->tree->PreorderFirst();
      while(n!=NULL){
        if(!n->IsGhost() && n->IsLeaf()) break;
        n=ctx->tree->PreorderNxt(n);
      }
      assert(n!=NULL);

      size_t trg_size=0;
      const std::vector<Node_t*>& nodes=ctx->tree->GetNodeList();
      #pragma omp parallel for reduction(+:trg_size)
      for(size_t i=0;i<nodes.size();i++){
        if(nodes[i]->IsLeaf() && !nodes[i]->IsGhost()){
          trg_size+=nodes[i]->trg_value.Dim();
        }
      }
      trg_value.ReInit(trg_size,&n->trg_value[0]);
    }
    pvfmm::par::ScatterReverse  (trg_value, scatter_index, ctx->comm, n_trg);
    #pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){
      size_t a=((tid+0)*n_trg)/omp_p;
      size_t b=((tid+1)*n_trg)/omp_p;
      for(size_t i=a;i<b;i++){
        for(size_t k=0;k<ker_dim[1];k++){
          trg_vel[i*ker_dim[1]+k]=trg_value[i*ker_dim[1]+k]*trg_scal[k];
        }
      }
    }
  }
  pvfmm::Profile::Toc();

  prof_FLOPS=pvfmm::Profile::Add_FLOP(0)-prof_FLOPS;
  //pvfmm::Profile::print(&ctx->comm);
  //PVFMMDestroyContext<T>(ctx_);
  PROFILEEND("",prof_FLOPS);
}

template<typename T>
void PVFMM_GlobalRepart(size_t nv, size_t stride,
    const T* x, const T* tension, size_t* nvr, T** xr,
    T** tensionr, void* user_ptr){
  assert(false); // I think this doesn't work
  exit(0);

  MPI_Comm comm=MPI_COMM_WORLD;

  // Get bounding box.
  T scale_x, shift_x[COORD_DIM];
  PVFMMBoundingBox(nv*stride, x, &scale_x, shift_x, comm);

  pvfmm::Vector<pvfmm::MortonId> ves_mid(nv);
  { // Create MortonIds for vesicles.
    scale_x/=stride;
    #pragma omp parallel for
    for(size_t i=0;i<nv;i++){
      T  x_ves[3]={0,0,0};
      const T* x_=&x[i*COORD_DIM*stride];
      for(size_t j=0;j<stride;j++){
        for(size_t k=0;k<COORD_DIM;k++){
          x_ves[k]+=x_[0];
          ++x_;
        }
      }
      for(size_t k=0;k<COORD_DIM;k++){
        x_ves[k]=x_ves[k]*scale_x+shift_x[k];
        assert(x_ves[k]>0.0);
        assert(x_ves[k]<1.0);
      }
      ves_mid[i]=pvfmm::MortonId(x_ves);
    }
  }

  // Determine scatter index vector.
  pvfmm::Vector<size_t> scatter_index;
  pvfmm::par::SortScatterIndex(ves_mid, scatter_index, comm);

  // Allocate memory for output.
  nvr[0]=scatter_index.Dim();
  xr      [0]=new T[nvr[0]*stride*COORD_DIM];
  tensionr[0]=new T[nvr[0]*stride          ];

  { // Scatter x
    pvfmm::Vector<T> data(nv*stride*COORD_DIM,(T*)x);
    pvfmm::par::ScatterForward(data, scatter_index, comm);

    assert(data.Dim()==nvr[0]*stride*COORD_DIM);
    memcpy(xr[0],&data[0],data.Dim()*sizeof(T));
  }

  { // Scatter tension
    pvfmm::Vector<T> data(nv*stride,(T*)tension);
    pvfmm::par::ScatterForward(data, scatter_index, comm);

    assert(data.Dim()==nvr[0]*stride);
    memcpy(tensionr[0],&data[0],data.Dim()*sizeof(T));
  }
}

