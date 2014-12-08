#include <mpi.h>
#include <omp.h>
#include <iostream>
#include "DataIO.h"
#include "BgFlow.h"
#include "Array.h"
#include "Vectors.h"
#include "Surface.h"
#include "MovePole.h"
#include "Parameters.h"
#include "OperatorsMats.h"
#include "VesInteraction.h"
#include "PVFMMInterface.h"
#include "InterfacialVelocity.h"

#include <ompUtils.h>
#include <mpi_tree.hpp>
#include <interac_list.hpp>
#include <legendre_rule.hpp>
#include <stacktrace.h>
#include <profile.hpp>

typedef Device<CPU> DevCPU;
extern const DevCPU the_cpu_device(0);
typedef double value_type;
typedef double real;

template <class Real_t>
void LegPoly(const Real_t* x, size_t n, size_t q, Real_t* y){
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

template <class Real_t>
class QuadraticPatch{

  public:

  QuadraticPatch(){}

  QuadraticPatch(Real_t* x, int dof_){
    dof=dof_;
    coeff.resize(9*dof);
    for(size_t i=0;i<dof;i++){
      Real_t tmp[3*3];
      for(size_t j=0;j<3;j++){
        tmp[j*3+0]= x[(j*3+1)*dof+i];
        tmp[j*3+1]=(x[(j*3+2)*dof+i]-x[(j*3+0)*dof+i])*0.5;
        tmp[j*3+2]=(x[(j*3+2)*dof+i]+x[(j*3+0)*dof+i])*0.5-x[(j*3+1)*dof+i];
      }
      for(size_t j=0;j<3;j++){
        coeff[i*9+0*3+j]= tmp[1*3+j];
        coeff[i*9+1*3+j]=(tmp[2*3+j]-tmp[0*3+j])*0.5;
        coeff[i*9+2*3+j]=(tmp[2*3+j]+tmp[0*3+j])*0.5-tmp[1*3+j];
      }
    }
  }

  void eval(Real_t x, Real_t y, Real_t* val){
    Real_t x_[3]={1,x,x*x};
    Real_t y_[3]={1,y,y*y};
    for(size_t k=0;k<dof;k++) val[k]=0;
    for(size_t k=0;k<dof;k++){
      for(size_t i=0;i<3;i++){
        for(size_t j=0;j<3;j++){
          val[k]+=coeff[i+3*j+9*k]*x_[i]*y_[j];
        }
      }
    }
  }

  void grad(Real_t x, Real_t y, Real_t* val){
    Real_t x_[3]={1,x,x*x};
    Real_t y_[3]={1,y,y*y};
    Real_t x__[3]={0,1,2*x};
    Real_t y__[3]={0,1,2*y};
    for(size_t k=0;k<2*dof;k++) val[k]=0;
    for(size_t k=0;k<dof;k++){
      for(size_t i=0;i<3;i++){
        for(size_t j=0;j<3;j++){
          val[k+dof*0]+=coeff[i+3*j+9*k]*x__[i]*y_[j];
          val[k+dof*1]+=coeff[i+3*j+9*k]*x_[i]*y__[j];
        }
      }
    }
  }

  private:

  int dof;
  std::vector<Real_t> coeff;
};

template <typename Surf_t>
void StokesNearSingular(const Surf_t& S, const typename Surf_t::Vec_t &force, const typename Surf_t::Vec_t& veloc_self, real r_near,
    const size_t* trg_cnt, const size_t* trg_dsp,  real* trg_coord, real* trg_veloc){
  typedef typename Surf_t::Vec_t Vec_t;
  typedef typename Vec_t::scalars_type Sca_t;
  size_t omp_p=omp_get_max_threads();

  const int force_dim=COORD_DIM;
  const int veloc_dim=COORD_DIM;

  const Vec_t& x=S.getPosition();
  size_t N_ves = x.getNumSubs(); // Number of vesicles
  size_t M_ves = x.getStride(); // Points per vesicle
  int k0_max=x.getGridDim().first;
  int k1_max=x.getGridDim().second;
  assert(M_ves==k0_max*k1_max);

  std::vector<real> pole_quad;
  { // compute quadrature to find pole
    size_t p=x.getShOrder();

    //Gauss-Legendre quadrature nodes and weights
    std::vector<real> x(p+1),w(p+1);
    cgqf(p+1, 1, 0.0, 0.0, -1.0, 1.0, &x[0], &w[0]);

    std::vector<real> leg((p+1)*(p+1));
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

    //real sum=0;
    //for(size_t i=0;i<p+1;i++) sum+=pole_quad[i]*2*p; std::cout<<sum<<'\n';
    //for(size_t i=0;i<p+1;i++) std::cout<<pole_quad[i]*2*p<<' '; std::cout<<'\n';
  }

  size_t interp_deg=10;
  pvfmm::Matrix<real> M(interp_deg,interp_deg);
  { // matrix for computing interpolation coefficients
    for(size_t i=0;i<interp_deg;i++){
      real x=0.0;
      if(i) x=(1.0+(i-1.0)/(interp_deg-2));

      real y,yj=1.0;
      y=1.0/(1+x); // basis functions are 1/(1+x)^j
      for(size_t j=0;j<interp_deg;j++){
        M[i][j]=yj;
        yj*=y;
      }
    }
    M=M.pinv();
  }

  { // Compute trg_veloc. [[ At this point we can offload. ]]
    #pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){
      size_t a=((tid+0)*N_ves)/omp_p;
      size_t b=((tid+1)*N_ves)/omp_p;

      pvfmm::Vector<real> s_coord(M_ves*COORD_DIM);
      pvfmm::Vector<real> s_force(M_ves*force_dim);
      pvfmm::Vector<real> s_veloc(M_ves*veloc_dim);

      pvfmm::Vector<real> pole_coord(2*COORD_DIM);
      pvfmm::Vector<real> pole_veloc(2*veloc_dim);

      QuadraticPatch<real> patch_coord;
      QuadraticPatch<real> patch_veloc;
      pvfmm::Vector<real> patch_coord_(3*3*COORD_DIM);
      pvfmm::Vector<real> patch_veloc_(3*3*veloc_dim);

      for(size_t i=a;i<b;i++){ // loop over all vesicles
        { // Set s_coord, s_force, s_veloc
          // read each component of x
          const real* xk=x.getSubN_begin(i)+0*M_ves;
          const real* yk=x.getSubN_begin(i)+1*M_ves;
          const real* zk=x.getSubN_begin(i)+2*M_ves;

          // read each component of qforce
          const real* fxk=force.getSubN_begin(i)+0*M_ves;
          const real* fyk=force.getSubN_begin(i)+1*M_ves;
          const real* fzk=force.getSubN_begin(i)+2*M_ves;
          assert(force_dim==3);

          // read each component of veloc
          const real* vsxk=veloc_self.getSubN_begin(i)+0*M_ves;
          const real* vsyk=veloc_self.getSubN_begin(i)+1*M_ves;
          const real* vszk=veloc_self.getSubN_begin(i)+2*M_ves;
          assert(veloc_dim==3);

          for(size_t j=0;j<M_ves;j++){
            s_coord[j*COORD_DIM+0]=xk[j];
            s_coord[j*COORD_DIM+1]=yk[j];
            s_coord[j*COORD_DIM+2]=zk[j];

            s_force[j*force_dim+0]=fxk[j];
            s_force[j*force_dim+1]=fyk[j];
            s_force[j*force_dim+2]=fzk[j];

            s_veloc[j*veloc_dim+0]=vsxk[j];
            s_veloc[j*veloc_dim+1]=vsyk[j];
            s_veloc[j*veloc_dim+2]=vszk[j];
          }
        }
        { // Set pole values: pole_coord, pole_veloc
          pole_coord.SetZero();
          pole_veloc.SetZero();
          for(size_t k0=0;k0<k0_max;k0++){
            for(size_t k1=0;k1<k1_max;k1++){
              size_t k=k1+k0*k1_max;
              pole_coord[0*COORD_DIM+0]+=pole_quad[k0_max-1-k0]*s_coord[k*COORD_DIM+0];
              pole_coord[0*COORD_DIM+1]+=pole_quad[k0_max-1-k0]*s_coord[k*COORD_DIM+1];
              pole_coord[0*COORD_DIM+2]+=pole_quad[k0_max-1-k0]*s_coord[k*COORD_DIM+2];
              pole_coord[1*COORD_DIM+0]+=pole_quad[         k0]*s_coord[k*COORD_DIM+0];
              pole_coord[1*COORD_DIM+1]+=pole_quad[         k0]*s_coord[k*COORD_DIM+1];
              pole_coord[1*COORD_DIM+2]+=pole_quad[         k0]*s_coord[k*COORD_DIM+2];

              pole_veloc[0*veloc_dim+0]+=pole_quad[k0_max-1-k0]*s_veloc[k*veloc_dim+0];
              pole_veloc[0*veloc_dim+1]+=pole_quad[k0_max-1-k0]*s_veloc[k*veloc_dim+1];
              pole_veloc[0*veloc_dim+2]+=pole_quad[k0_max-1-k0]*s_veloc[k*veloc_dim+2];
              pole_veloc[1*veloc_dim+0]+=pole_quad[         k0]*s_veloc[k*veloc_dim+0];
              pole_veloc[1*veloc_dim+1]+=pole_quad[         k0]*s_veloc[k*veloc_dim+1];
              pole_veloc[1*veloc_dim+2]+=pole_quad[         k0]*s_veloc[k*veloc_dim+2];
            }
          }
          //std::cout<<pole_coord<<'\n';
        }

        pvfmm::Vector<real> t_coord(trg_cnt[i]*COORD_DIM, &trg_coord[trg_dsp[i]*COORD_DIM], false);
        pvfmm::Vector<real> t_veloc(trg_cnt[i]*veloc_dim, &trg_veloc[trg_dsp[i]*veloc_dim], false);
        t_veloc.SetZero();

        { // Subtract wrong near potential
          stokes_sl(&s_coord[0], M_ves, &s_force[0], 1, &t_coord[0], trg_cnt[i], &t_veloc[0], NULL);
          for(size_t j=0;j<t_veloc.Dim();j++) t_veloc[j]=-t_veloc[j];
        }

        // Add corrected near potential
        for(size_t j=0;j<trg_cnt[i];j++){ // loop over target points
          real  t_coord_j[COORD_DIM];
          real* t_veloc_j;

          { // Set t_coord, t_veloc
            t_coord_j[0]=t_coord[j*COORD_DIM+0];
            t_coord_j[1]=t_coord[j*COORD_DIM+1];
            t_coord_j[2]=t_coord[j*COORD_DIM+2];
            t_veloc_j=&t_veloc[j*veloc_dim];
          }

          { // Set t_veloc_j for surface mesh points and {continue;}
            bool surface_point=false;
            for(size_t k0=0;k0<k0_max;k0++){
              for(size_t k1=0;k1<k1_max;k1++){
                size_t k=k1+k0*k1_max;
                // Add self interaction
                if(t_coord_j[0]==s_coord[k*COORD_DIM+0])
                if(t_coord_j[1]==s_coord[k*COORD_DIM+1])
                if(t_coord_j[2]==s_coord[k*COORD_DIM+2]){
                  t_veloc_j[0]+=s_veloc[k*veloc_dim+0];
                  t_veloc_j[1]+=s_veloc[k*veloc_dim+1];
                  t_veloc_j[2]+=s_veloc[k*veloc_dim+2];
                  k0=k0_max; k1=k1_max;
                  surface_point=true;
                  break;
                }
              }
            }
            if(surface_point) continue;
          }

          { // Find nearest point on mesh and create patch
            int k0_, k1_; // mesh coordinates for nearest point
            { // Find nearest point on mesh
              real r2_min=-1.0;
              for(size_t k0=0;k0<k0_max;k0++){
                for(size_t k1=0;k1<k1_max;k1++){
                  size_t k=k1+k0*k1_max;
                  real dx=(t_coord[j*COORD_DIM+0]-s_coord[k*COORD_DIM+0]);
                  real dy=(t_coord[j*COORD_DIM+1]-s_coord[k*COORD_DIM+1]);
                  real dz=(t_coord[j*COORD_DIM+2]-s_coord[k*COORD_DIM+2]);
                  real r2=dx*dx+dy*dy+dz*dz;
                  if(r2<r2_min || r2_min<0){
                    r2_min=r2;
                    k0_=k0;
                    k1_=k1;
                  }
                }
              }
            }
            { // Create patch
              { // (0,0)
                if(k0_>0){
                  int k=((k1_+k1_max-1)%k1_max)+(k0_-1)*k1_max;
                  assert(k>=0 && k<M_ves);
                  patch_coord_[(0*3+0)*COORD_DIM+0]=s_coord[k*COORD_DIM+0];
                  patch_coord_[(0*3+0)*COORD_DIM+1]=s_coord[k*COORD_DIM+1];
                  patch_coord_[(0*3+0)*COORD_DIM+2]=s_coord[k*COORD_DIM+2];

                  patch_veloc_[(0*3+0)*veloc_dim+0]=s_veloc[k*veloc_dim+0];
                  patch_veloc_[(0*3+0)*veloc_dim+1]=s_veloc[k*veloc_dim+1];
                  patch_veloc_[(0*3+0)*veloc_dim+2]=s_veloc[k*veloc_dim+2];
                }else{
                  patch_coord_[(0*3+0)*COORD_DIM+0]=pole_coord[0*COORD_DIM+0];
                  patch_coord_[(0*3+0)*COORD_DIM+1]=pole_coord[0*COORD_DIM+1];
                  patch_coord_[(0*3+0)*COORD_DIM+2]=pole_coord[0*COORD_DIM+2];

                  patch_veloc_[(0*3+0)*veloc_dim+0]=pole_veloc[0*veloc_dim+0];
                  patch_veloc_[(0*3+0)*veloc_dim+1]=pole_veloc[0*veloc_dim+1];
                  patch_veloc_[(0*3+0)*veloc_dim+2]=pole_veloc[0*veloc_dim+2];
                }
              }
              { // (0,1)
                if(k0_>0){
                  int k=  k1_                  +(k0_-1)*k1_max;
                  assert(k>=0 && k<M_ves);
                  patch_coord_[(0*3+1)*COORD_DIM+0]=s_coord[k*COORD_DIM+0];
                  patch_coord_[(0*3+1)*COORD_DIM+1]=s_coord[k*COORD_DIM+1];
                  patch_coord_[(0*3+1)*COORD_DIM+2]=s_coord[k*COORD_DIM+2];

                  patch_veloc_[(0*3+1)*veloc_dim+0]=s_veloc[k*veloc_dim+0];
                  patch_veloc_[(0*3+1)*veloc_dim+1]=s_veloc[k*veloc_dim+1];
                  patch_veloc_[(0*3+1)*veloc_dim+2]=s_veloc[k*veloc_dim+2];
                }else{
                  patch_coord_[(0*3+1)*COORD_DIM+0]=pole_coord[0*COORD_DIM+0];
                  patch_coord_[(0*3+1)*COORD_DIM+1]=pole_coord[0*COORD_DIM+1];
                  patch_coord_[(0*3+1)*COORD_DIM+2]=pole_coord[0*COORD_DIM+2];

                  patch_veloc_[(0*3+1)*veloc_dim+0]=pole_veloc[0*veloc_dim+0];
                  patch_veloc_[(0*3+1)*veloc_dim+1]=pole_veloc[0*veloc_dim+1];
                  patch_veloc_[(0*3+1)*veloc_dim+2]=pole_veloc[0*veloc_dim+2];
                }
              }
              { // (0,2)
                if(k0_>0){
                  int k=((k1_       +1)%k1_max)+(k0_-1)*k1_max;
                  assert(k>=0 && k<M_ves);
                  patch_coord_[(0*3+2)*COORD_DIM+0]=s_coord[k*COORD_DIM+0];
                  patch_coord_[(0*3+2)*COORD_DIM+1]=s_coord[k*COORD_DIM+1];
                  patch_coord_[(0*3+2)*COORD_DIM+2]=s_coord[k*COORD_DIM+2];

                  patch_veloc_[(0*3+2)*veloc_dim+0]=s_veloc[k*veloc_dim+0];
                  patch_veloc_[(0*3+2)*veloc_dim+1]=s_veloc[k*veloc_dim+1];
                  patch_veloc_[(0*3+2)*veloc_dim+2]=s_veloc[k*veloc_dim+2];
                }else{
                  patch_coord_[(0*3+2)*COORD_DIM+0]=pole_coord[0*COORD_DIM+0];
                  patch_coord_[(0*3+2)*COORD_DIM+1]=pole_coord[0*COORD_DIM+1];
                  patch_coord_[(0*3+2)*COORD_DIM+2]=pole_coord[0*COORD_DIM+2];

                  patch_veloc_[(0*3+2)*veloc_dim+0]=pole_veloc[0*veloc_dim+0];
                  patch_veloc_[(0*3+2)*veloc_dim+1]=pole_veloc[0*veloc_dim+1];
                  patch_veloc_[(0*3+2)*veloc_dim+2]=pole_veloc[0*veloc_dim+2];
                }
              }

              { // (1,0)
                int k=((k1_+k1_max-1)%k1_max)+k0_*k1_max;
                assert(k>=0 && k<M_ves);
                patch_coord_[(1*3+0)*COORD_DIM+0]=s_coord[k*COORD_DIM+0];
                patch_coord_[(1*3+0)*COORD_DIM+1]=s_coord[k*COORD_DIM+1];
                patch_coord_[(1*3+0)*COORD_DIM+2]=s_coord[k*COORD_DIM+2];

                patch_veloc_[(1*3+0)*veloc_dim+0]=s_veloc[k*veloc_dim+0];
                patch_veloc_[(1*3+0)*veloc_dim+1]=s_veloc[k*veloc_dim+1];
                patch_veloc_[(1*3+0)*veloc_dim+2]=s_veloc[k*veloc_dim+2];
              }
              { // (1,1)
                int k=  k1_                  +k0_*k1_max;
                assert(k>=0 && k<M_ves);
                patch_coord_[(1*3+1)*COORD_DIM+0]=s_coord[k*COORD_DIM+0];
                patch_coord_[(1*3+1)*COORD_DIM+1]=s_coord[k*COORD_DIM+1];
                patch_coord_[(1*3+1)*COORD_DIM+2]=s_coord[k*COORD_DIM+2];

                patch_veloc_[(1*3+1)*veloc_dim+0]=s_veloc[k*veloc_dim+0];
                patch_veloc_[(1*3+1)*veloc_dim+1]=s_veloc[k*veloc_dim+1];
                patch_veloc_[(1*3+1)*veloc_dim+2]=s_veloc[k*veloc_dim+2];
              }
              { // (1,2)
                int k=((k1_       +1)%k1_max)+k0_*k1_max;
                assert(k>=0 && k<M_ves);
                patch_coord_[(1*3+2)*COORD_DIM+0]=s_coord[k*COORD_DIM+0];
                patch_coord_[(1*3+2)*COORD_DIM+1]=s_coord[k*COORD_DIM+1];
                patch_coord_[(1*3+2)*COORD_DIM+2]=s_coord[k*COORD_DIM+2];

                patch_veloc_[(1*3+2)*veloc_dim+0]=s_veloc[k*veloc_dim+0];
                patch_veloc_[(1*3+2)*veloc_dim+1]=s_veloc[k*veloc_dim+1];
                patch_veloc_[(1*3+2)*veloc_dim+2]=s_veloc[k*veloc_dim+2];
              }

              { // (2,0)
                if(k0_<k0_max-1){
                  int k=((k1_+k1_max-1)%k1_max)+(k0_+1)*k1_max;
                  assert(k>=0 && k<M_ves);
                  patch_coord_[(2*3+0)*COORD_DIM+0]=s_coord[k*COORD_DIM+0];
                  patch_coord_[(2*3+0)*COORD_DIM+1]=s_coord[k*COORD_DIM+1];
                  patch_coord_[(2*3+0)*COORD_DIM+2]=s_coord[k*COORD_DIM+2];

                  patch_veloc_[(2*3+0)*veloc_dim+0]=s_veloc[k*veloc_dim+0];
                  patch_veloc_[(2*3+0)*veloc_dim+1]=s_veloc[k*veloc_dim+1];
                  patch_veloc_[(2*3+0)*veloc_dim+2]=s_veloc[k*veloc_dim+2];
                }else{
                  patch_coord_[(2*3+0)*COORD_DIM+0]=pole_coord[1*COORD_DIM+0];
                  patch_coord_[(2*3+0)*COORD_DIM+1]=pole_coord[1*COORD_DIM+1];
                  patch_coord_[(2*3+0)*COORD_DIM+2]=pole_coord[1*COORD_DIM+2];

                  patch_veloc_[(2*3+0)*veloc_dim+0]=pole_veloc[1*veloc_dim+0];
                  patch_veloc_[(2*3+0)*veloc_dim+1]=pole_veloc[1*veloc_dim+1];
                  patch_veloc_[(2*3+0)*veloc_dim+2]=pole_veloc[1*veloc_dim+2];
                }
              }
              { // (2,1)
                if(k0_<k0_max-1){
                  int k=  k1_                  +(k0_+1)*k1_max;
                  assert(k>=0 && k<M_ves);
                  patch_coord_[(2*3+1)*COORD_DIM+0]=s_coord[k*COORD_DIM+0];
                  patch_coord_[(2*3+1)*COORD_DIM+1]=s_coord[k*COORD_DIM+1];
                  patch_coord_[(2*3+1)*COORD_DIM+2]=s_coord[k*COORD_DIM+2];

                  patch_veloc_[(2*3+1)*veloc_dim+0]=s_veloc[k*veloc_dim+0];
                  patch_veloc_[(2*3+1)*veloc_dim+1]=s_veloc[k*veloc_dim+1];
                  patch_veloc_[(2*3+1)*veloc_dim+2]=s_veloc[k*veloc_dim+2];
                }else{
                  patch_coord_[(2*3+1)*COORD_DIM+0]=pole_coord[1*COORD_DIM+0];
                  patch_coord_[(2*3+1)*COORD_DIM+1]=pole_coord[1*COORD_DIM+1];
                  patch_coord_[(2*3+1)*COORD_DIM+2]=pole_coord[1*COORD_DIM+2];

                  patch_veloc_[(2*3+1)*veloc_dim+0]=pole_veloc[1*veloc_dim+0];
                  patch_veloc_[(2*3+1)*veloc_dim+1]=pole_veloc[1*veloc_dim+1];
                  patch_veloc_[(2*3+1)*veloc_dim+2]=pole_veloc[1*veloc_dim+2];
                }
              }
              { // (2,2)
                if(k0_<k0_max-1){
                  int k=((k1_       +1)%k1_max)+(k0_+1)*k1_max;
                  assert(k>=0 && k<M_ves);
                  patch_coord_[(2*3+2)*COORD_DIM+0]=s_coord[k*COORD_DIM+0];
                  patch_coord_[(2*3+2)*COORD_DIM+1]=s_coord[k*COORD_DIM+1];
                  patch_coord_[(2*3+2)*COORD_DIM+2]=s_coord[k*COORD_DIM+2];

                  patch_veloc_[(2*3+2)*veloc_dim+0]=s_veloc[k*veloc_dim+0];
                  patch_veloc_[(2*3+2)*veloc_dim+1]=s_veloc[k*veloc_dim+1];
                  patch_veloc_[(2*3+2)*veloc_dim+2]=s_veloc[k*veloc_dim+2];
                }else{
                  patch_coord_[(2*3+2)*COORD_DIM+0]=pole_coord[1*COORD_DIM+0];
                  patch_coord_[(2*3+2)*COORD_DIM+1]=pole_coord[1*COORD_DIM+1];
                  patch_coord_[(2*3+2)*COORD_DIM+2]=pole_coord[1*COORD_DIM+2];

                  patch_veloc_[(2*3+2)*veloc_dim+0]=pole_veloc[1*veloc_dim+0];
                  patch_veloc_[(2*3+2)*veloc_dim+1]=pole_veloc[1*veloc_dim+1];
                  patch_veloc_[(2*3+2)*veloc_dim+2]=pole_veloc[1*veloc_dim+2];
                }
              }

              patch_coord=QuadraticPatch<real>(&patch_coord_[0],COORD_DIM);
              patch_veloc=QuadraticPatch<real>(&patch_veloc_[0],veloc_dim);
            }
          }

          std::vector<real> interp_coord(interp_deg*COORD_DIM);
          std::vector<real> interp_veloc(interp_deg*veloc_dim,0);
          { // Find nearest point on patch (first interpolation point)
            real sc[COORD_DIM];
            sc[0]=patch_coord_[(3*1+1)*COORD_DIM+0];
            sc[1]=patch_coord_[(3*1+1)*COORD_DIM+1];
            sc[2]=patch_coord_[(3*1+1)*COORD_DIM+2];

            real x=0,y=0;
            while(true){
              real dR[COORD_DIM]={t_coord_j[0]-sc[0],
                                  t_coord_j[1]-sc[1],
                                  t_coord_j[2]-sc[2]};
              real dR2=dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2];

              real sgrad[2*COORD_DIM];
              patch_coord.grad(x,y,sgrad);

              real dxdx=sgrad[0]*sgrad[0]+sgrad[1]*sgrad[1]+sgrad[2]*sgrad[2];
              real dydy=sgrad[3]*sgrad[3]+sgrad[4]*sgrad[4]+sgrad[5]*sgrad[5];
              real dxdy=sgrad[0]*sgrad[3]+sgrad[1]*sgrad[4]+sgrad[2]*sgrad[5];
              real dxdR=sgrad[0]*   dR[0]+sgrad[1]*   dR[1]+sgrad[2]*   dR[2];
              real dydR=sgrad[3]*   dR[0]+sgrad[4]*   dR[1]+sgrad[5]*   dR[2];
              if(dxdR*dxdR/dxdx+dydR*dydR/dydy < dR2*1e-6) break;

              real dx, dy;
              if(dxdy){
                dx=(dxdR/dxdy-dydR/dydy)/(dxdx/dxdy-dxdy/dydy);
                dy=(dydR/dxdy-dxdR/dxdx)/(dydy/dxdy-dxdy/dxdx);
              }else{
                dx=dxdR/dxdx;
                dy=dydR/dydy;
              }

              while(true){
                patch_coord.eval(x+dx,y+dy,sc);
                real dR_[COORD_DIM]={t_coord_j[0]-sc[0],
                                     t_coord_j[1]-sc[1],
                                     t_coord_j[2]-sc[2]};
                real dR2_=dR_[0]*dR_[0]+dR_[1]*dR_[1]+dR_[2]*dR_[2];
                if(dR2_<dR2){
                  x+=dx*2.0/3.0; y+=dy*2.0/3.0;
                  patch_coord.eval(x,y,sc);
                  break;
                }else{
                  dx*=0.5;
                  dy*=0.5;
                }
              }
            }

            // Set first interpolation point
            interp_coord[0*COORD_DIM+0]=sc[0];
            interp_coord[0*COORD_DIM+1]=sc[1];
            interp_coord[0*COORD_DIM+2]=sc[2];
            patch_veloc.eval(x,y,&interp_veloc[0*veloc_dim]);

            // For visualization of near point projections
            //t_coord[j*COORD_DIM+0]=sc[0];
            //t_coord[j*COORD_DIM+1]=sc[1];
            //t_coord[j*COORD_DIM+2]=sc[2];
          }

          { // Interpolate and find target velocity t_veloc_j
            real dR[COORD_DIM]={t_coord_j[0]-interp_coord[0],
                                t_coord_j[1]-interp_coord[1],
                                t_coord_j[2]-interp_coord[2]};
            real dR_norm=sqrt(dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2]);
            real OOdR=1.0/dR_norm;

            for(size_t l=1;l<interp_deg;l++){
              real x=(1.0+(l-1.0)/(interp_deg-2.0));
              interp_coord[l*COORD_DIM+0]=interp_coord[0*COORD_DIM+0]+dR[0]*OOdR*r_near*x;
              interp_coord[l*COORD_DIM+1]=interp_coord[0*COORD_DIM+1]+dR[1]*OOdR*r_near*x;
              interp_coord[l*COORD_DIM+2]=interp_coord[0*COORD_DIM+2]+dR[2]*OOdR*r_near*x;
            }

            // Compute velocity at interpolation points
            stokes_sl(&s_coord[0], M_ves, &s_force[0], 1, &interp_coord[1*COORD_DIM], interp_deg-1, &interp_veloc[1*veloc_dim], NULL);

            // Interpolate
            for(size_t k=0;k<veloc_dim;k++){
              pvfmm::Matrix<real> y(interp_deg,1);
              pvfmm::Matrix<real> x(interp_deg,1);
              pvfmm::Matrix<real> coeff(interp_deg,1);

              real x_=1.0/(1+dR_norm/r_near);
              for(size_t l=0;l<interp_deg;l++){
                y[l][0]=interp_veloc[l*veloc_dim+k];
                x[l][0]=(l?x[l-1][0]*x_:1.0);
              }

              coeff=M*y;
              for(size_t l=0;l<interp_deg;l++) t_veloc_j[k]+=coeff[l][0]*x[l][0];
            }
          }
        }
      }
    }
  }
}

template <typename Surf_t>
void NearInteractions(const Surf_t& S, const typename Surf_t::Vec_t &force, const typename Surf_t::Vec_t& veloc_self,
    typename Surf_t::Vec_t &velocity, const MPI_Comm& comm){
  typedef typename Surf_t::Vec_t Vec_t;

  const int force_dim=COORD_DIM;
  const int veloc_dim=COORD_DIM;

  int np, rank;
  MPI_Comm_size(comm,&np);
  MPI_Comm_rank(comm,&rank);
  size_t omp_p=omp_get_max_threads();

  const Vec_t& x=S.getPosition();
  size_t N_ves = x.getNumSubs(); // Number of vesicles
  size_t M_ves = x.getStride(); // Points per vesicle
  assert(M_ves==x.getGridDim().first*x.getGridDim().second);

  assert(N_ves==force.getNumSubs());
  assert(M_ves==force.getStride());
  assert(N_ves==velocity.getNumSubs());
  assert(M_ves==velocity.getStride());

  real r_near=0;
  pvfmm::Vector<size_t> scatter_index; // to scatter points to the interacting vesicle.
  pvfmm::Vector<size_t> pt_interac_cnt; // number of times a point has to be repeated.
  pvfmm::Vector<size_t> pt_interac_dsp;
  pvfmm::Vector<size_t> ves_interac_cnt; // number of consecutive points interacting with each ves.
  pvfmm::Vector<size_t> ves_interac_dsp;
  pvfmm::Profile::Tic("FindNear",&comm);
  { // Find near points for each vesicle
    typedef pvfmm::MPI_Node<real> Node_t;
    typedef pvfmm::MPI_Tree<Node_t> Tree_t;
    typedef Node_t::NodeData NodeData_t;

    size_t max_pts=80; // TODO: Tuning parameter
    pvfmm::BoundaryType bndry=pvfmm::FreeSpace; // Boundary condition
    size_t data_dof=(sizeof(size_t)+sizeof(real)-1)/sizeof(real); // data_dof for storing size_t type data

    size_t ves_id_offset;
    { // Get ves_id_offset
      long long disp=0;
      long long size=N_ves;
      MPI_Scan(&size, &disp, 1, MPI_LONG_LONG, MPI_SUM, comm);
      ves_id_offset=disp-size;
    }

    Tree_t pt_tree(comm);
    real scale_x, shift_x[COORD_DIM];
    pvfmm::Profile::Tic("TreeConstruct",&comm);
    { // particle tree
      bool prof_orig_state=pvfmm::Profile::Enable(false);
      NodeData_t node_data;
      node_data.dim=COORD_DIM;
      node_data.max_pts=max_pts;

      pvfmm::Profile::Tic("SetPoints",&comm);
      { // Set node_data.ves_coord, node_data.ves_value, node_data.max_depth, r_near
        pvfmm::Vector<real>& pt_coord=node_data.pt_coord;
        pvfmm::Vector<real>& pt_value=node_data.pt_value;
        pt_coord.ReInit(N_ves*M_ves*COORD_DIM);
        pt_value.ReInit(N_ves*M_ves*data_dof );

        std::vector<real> r2_ves_(omp_p);
        #pragma omp parallel for
        for(size_t tid=0;tid<omp_p;tid++){
          size_t a=((tid+0)*N_ves)/omp_p;
          size_t b=((tid+1)*N_ves)/omp_p;

          real r2_ves=0;
          real one_over_M=1.0/M_ves;
          for(size_t i=a;i<b;i++){
            // read each component of x
            const real* xk=x.getSubN_begin(i)+0*M_ves;
            const real* yk=x.getSubN_begin(i)+1*M_ves;
            const real* zk=x.getSubN_begin(i)+2*M_ves;

            real center_coord[COORD_DIM]={0,0,0};
            for(size_t j=0;j<M_ves;j++){
              center_coord[0]+=xk[j];
              center_coord[1]+=yk[j];
              center_coord[2]+=zk[j];
            }
            center_coord[0]*=one_over_M;
            center_coord[1]*=one_over_M;
            center_coord[2]*=one_over_M;

            for(size_t j=0;j<M_ves;j++){
              pt_coord[(i*M_ves+j)*COORD_DIM+0]=xk[j];
              pt_coord[(i*M_ves+j)*COORD_DIM+1]=yk[j];
              pt_coord[(i*M_ves+j)*COORD_DIM+2]=zk[j];
              ((size_t*)&pt_value[i*M_ves+j])[0]=ves_id_offset+i;

              real dx=(xk[j]-center_coord[0]);
              real dy=(yk[j]-center_coord[1]);
              real dz=(zk[j]-center_coord[2]);
              real r2=dx*dx+dy*dy+dz*dz;
              r2_ves=std::max(r2_ves,r2);
            }
          }
          r2_ves_[tid]=r2_ves;
        }

        real r_ves=0;
        { // Determine r_ves
          double r_ves_loc=0, r_ves_glb=0;
          for(size_t tid=0;tid<omp_p;tid++){
            r_ves_loc=std::max(r2_ves_[tid], r_ves_loc);
          }
          r_ves_loc=sqrt(r_ves_loc);
          MPI_Allreduce(&r_ves_loc, &r_ves_glb, 1, MPI_DOUBLE, MPI_MAX, comm);
          r_ves=r_ves_glb;
        }
        r_near=r_ves*1.0; // r_near is some function of r_ves.

        { // Determine scale_x, shift_x, node_data.max_depth
          real scale_tmp;
          PVFMMBoundingBox(N_ves*M_ves, &pt_coord[0], &scale_tmp, shift_x, comm);
          size_t pt_tree_depth=0;
          { // scale_x, pt_tree_depth
            real domain_length=1.0/scale_tmp;
            real leaf_length=r_near;
            scale_x=1.0/leaf_length;
            while(domain_length*scale_x>1.0){
              scale_x*=0.5;
              pt_tree_depth++;
            }
          }
          for(size_t j=0;j<COORD_DIM;j++){ // Update shift_x
            shift_x[j]=(shift_x[j]/scale_tmp)*scale_x;
          }
          node_data.max_depth=pt_tree_depth;
        }

        #pragma omp parallel for
        for(size_t tid=0;tid<omp_p;tid++){
          size_t a=((tid+0)*N_ves*M_ves)/omp_p;
          size_t b=((tid+1)*N_ves*M_ves)/omp_p;
          for(size_t i=a;i<b;i++){
            pt_coord[i*COORD_DIM+0]=pt_coord[i*COORD_DIM+0]*scale_x+shift_x[0];
            pt_coord[i*COORD_DIM+1]=pt_coord[i*COORD_DIM+1]*scale_x+shift_x[1];
            pt_coord[i*COORD_DIM+2]=pt_coord[i*COORD_DIM+2]*scale_x+shift_x[2];
          }
        }
      }
      pvfmm::Profile::Toc();

      pvfmm::Profile::Tic("InitTree",&comm);
      pt_tree.Initialize(&node_data);
      pvfmm::Profile::Toc();

      pvfmm::Profile::Tic("Balance21",&comm);
      pt_tree.Balance21(bndry);
      pvfmm::Profile::Toc();

      pvfmm::Profile::Tic("ConstructLET",&comm);
      pt_tree.ConstructLET(bndry);
      pvfmm::Profile::Toc();

      { // Write2File
        //pt_tree.Write2File("ptree");
      }
      pvfmm::Profile::Enable(prof_orig_state);
    }
    pvfmm::Profile::Toc();

    pvfmm::Profile::Tic("NearMapping",&comm);
    { // Compute scatter_index, ves_interac_cnt, pt_interac_cnt
      std::vector<Node_t*> node_lst; // all leaf nodes
      size_t non_ghost_range[2]; // index range of nodes which are not ghost
      { // Get node_lst
        std::vector<Node_t*>& nodes=pt_tree.GetNodeList();
        non_ghost_range[0]=nodes.size();
        non_ghost_range[1]=0;

        for(size_t i=0;i<nodes.size();i++){
          if(nodes[i]->IsLeaf() && nodes[i]->pt_coord.Dim()){
            node_lst.push_back(nodes[i]);
            if(!nodes[i]->IsGhost()){
              size_t node_id=node_lst.size()-1;
              if(non_ghost_range[0]>node_id  ) non_ghost_range[0]=node_id  ;
              if(non_ghost_range[1]<node_id+1) non_ghost_range[1]=node_id+1;
            }
          }
        }

        // Set node_id
        #pragma omp parallel for
        for(size_t tid=0;tid<omp_p;tid++){
          size_t a=(node_lst.size()*(tid+0))/omp_p;
          size_t b=(node_lst.size()*(tid+1))/omp_p;
          for(size_t i=a;i<b;i++) node_lst[i]->node_id=i;
        }
      }

      size_t n_pt; // number of points in all leaf nodes
      pvfmm::Vector<size_t> pt_cnt; // point count per leaf node
      pvfmm::Vector<size_t> pt_dsp; // point displ per leaf node
      { // Point count and displacement per leaf node.
        pt_cnt.ReInit(node_lst.size());
        pt_dsp.ReInit(node_lst.size());
        pt_dsp[0]=0;

        #pragma omp parallel for
        for(size_t i=0;i<node_lst.size();i++){
          pt_cnt[i]=node_lst[i]->pt_coord.Dim()/COORD_DIM;
        }
        pvfmm::omp_par::scan(&pt_cnt[0], &pt_dsp[0], pt_cnt.Dim());
        n_pt=pt_dsp[pt_dsp.Dim()-1]+pt_cnt[pt_cnt.Dim()-1];
      }

      pvfmm::Vector<size_t> pt_id(n_pt); // global id for each point
      pvfmm::Vector<size_t> ves_id(n_pt); // global vesicle id for each point
      pvfmm::Vector<real> pt_coord(n_pt*COORD_DIM); // coordinates for each point
      { // Create linear array of points
        #pragma omp parallel for
        for(size_t tid=0;tid<omp_p;tid++){
          size_t a,b;
          { // compute a,b
            size_t a_=(tid+0)*n_pt/omp_p;
            size_t b_=(tid+1)*n_pt/omp_p;
            a=std::lower_bound(&pt_dsp[0],&pt_dsp[0]+pt_dsp.Dim(),a_)-&pt_dsp[0];
            b=std::lower_bound(&pt_dsp[0],&pt_dsp[0]+pt_dsp.Dim(),b_)-&pt_dsp[0];
          }
          for(size_t i=a;i<b;i++){
            size_t dsp=pt_dsp[i];
            pvfmm::Vector<real  >& node_pt_coord=node_lst[i]->pt_coord;
            pvfmm::Vector<real  >&   node_ves_id=node_lst[i]->pt_value;
            pvfmm::Vector<size_t>&    node_pt_id=node_lst[i]->pt_scatter;
            for(size_t j=0;j<pt_cnt[i];j++){
              pt_coord[(dsp+j)*COORD_DIM+0]=node_pt_coord[j*COORD_DIM+0];
              pt_coord[(dsp+j)*COORD_DIM+1]=node_pt_coord[j*COORD_DIM+1];
              pt_coord[(dsp+j)*COORD_DIM+2]=node_pt_coord[j*COORD_DIM+2];
              ves_id[dsp+j]=*((size_t*)&node_ves_id[j*data_dof]);
            }
            if(node_pt_id.Dim()){
              for(size_t j=0;j<pt_cnt[i];j++){
                assert(node_pt_id[j]/M_ves==ves_id[dsp+j]);
                pt_id[dsp+j]=node_pt_id[j];
              }
            }else{
              for(size_t j=0;j<pt_cnt[i];j++){
                pt_id[dsp+j]=0;
              }
            }
          }
        }
      }

      // For each node node_lst[i], the near node indices are given by:
      // near_lst[near_dsp[i]], ..., near_lst[near_dsp[i]+near_cnt[i]-1]
      pvfmm::Vector<size_t> near_cnt;
      pvfmm::Vector<size_t> near_dsp;
      pvfmm::Vector<size_t> near_lst;
      { // Get near list for each pt_tree node.

        near_cnt.ReInit(node_lst.size()); near_cnt.SetZero();
        near_dsp.ReInit(node_lst.size()); near_dsp.SetZero();

        pvfmm::Matrix<Node_t*> near_ptr_lst;
        // Can have at most 57 neighbors (including self).
        near_ptr_lst.ReInit(node_lst.size(),57);
        near_ptr_lst.SetZero();

        pt_tree.SetColleagues(bndry);
        static pvfmm::InteracList<Node_t> interac(COORD_DIM);
        #pragma omp parallel for
        for(size_t tid=0;tid<omp_p;tid++){
          size_t a=(node_lst.size()*(tid+0))/omp_p;
          size_t b=(node_lst.size()*(tid+1))/omp_p;
          std::vector<Node_t*> interac_lst;
          for(size_t i=a;i<b;i++){
            size_t iter=0;
            interac_lst=interac.BuildList(node_lst[i],pvfmm::U0_Type);
            for(size_t j=0;j<interac_lst.size();j++){
              if(interac_lst[j] && interac_lst[j]->pt_coord.Dim()) near_ptr_lst[i][iter++]=interac_lst[j];
            }
            interac_lst=interac.BuildList(node_lst[i],pvfmm::U1_Type);
            for(size_t j=0;j<interac_lst.size();j++){
              if(interac_lst[j] && interac_lst[j]->pt_coord.Dim()) near_ptr_lst[i][iter++]=interac_lst[j];
            }
            interac_lst=interac.BuildList(node_lst[i],pvfmm::U2_Type);
            for(size_t j=0;j<interac_lst.size();j++){
              if(interac_lst[j] && interac_lst[j]->pt_coord.Dim()) near_ptr_lst[i][iter++]=interac_lst[j];
            }
            assert(iter<=57);
            near_cnt[i]=iter;
          }
        }
        pvfmm::omp_par::scan(&near_cnt[0], &near_dsp[0], near_cnt.Dim());
        near_lst.ReInit(near_cnt[near_cnt.Dim()-1]+near_dsp[near_dsp.Dim()-1]);

        #pragma omp parallel for
        for(size_t tid=0;tid<omp_p;tid++){
          size_t a=(node_lst.size()*(tid+0))/omp_p;
          size_t b=(node_lst.size()*(tid+1))/omp_p;
          for(size_t i=a;i<b;i++){
            size_t dsp=near_dsp[i];
            for(size_t j=0;j<near_cnt[i];j++){
              near_lst[dsp+j]=near_ptr_lst[i][j]->node_id;
            }
          }
        }
      }

      pvfmm::Vector<pvfmm::par::SortPair<size_t, size_t> > near_pair; // (pt_id, ves_id)
      { // construct pairs (pt_id, ves_id)
        std::vector<size_t> cnt(omp_p);
        std::vector<std::vector<pvfmm::par::SortPair<size_t, size_t> > > near_pair_(omp_p);
        real r2_near=(r_near*scale_x)*(r_near*scale_x);
        size_t flops=0;
        #pragma omp parallel for reduction(+:flops)
        for(size_t tid=0;tid<omp_p;tid++){ // TODO: optimize this
          std::vector<pvfmm::par::SortPair<size_t, size_t> >& near_pair=near_pair_[tid];

          // temporary vector to sort near points by vesicle id
          std::vector<pvfmm::par::SortPair<size_t, real*> > pair_vesid_coordptr;
          std::vector<real> near_pt_coord; // coordinates of near points
          std::vector<size_t> near_ves_id; // unique vesicle ids for near points
          std::vector<size_t> near_pt_cnt; // count of near points with vesicle id: near_ves_id[i]
          std::vector<size_t> near_pt_dsp; // displ of near points with vesicle id: near_ves_id[i]
          pvfmm::par::SortPair<size_t,real*> p_vesid_coordptr;
          pvfmm::par::SortPair<size_t,size_t> p_ptid_vesid;

          size_t a=((tid+0)*non_ghost_range[1] + (omp_p-(tid+0))*non_ghost_range[0])/omp_p;
          size_t b=((tid+1)*non_ghost_range[1] + (omp_p-(tid+1))*non_ghost_range[0])/omp_p;
          assert(non_ghost_range[1]<=node_lst.size());
          for(size_t i=a;i<b;i++){ // all non-ghost leaves

            { // set near_pt_coord, near_ves_id, near_pt_cnt, near_pt_dsp

              // Sort pair_vesid_coordptr and then use it construct near_pt_coord
              pair_vesid_coordptr.clear();
              for(size_t j=0;j<near_cnt[i];j++){ // near list
                size_t nnid=near_lst[near_dsp[i]+j]; // near node id
                size_t pt_dsp_nnid=pt_dsp[nnid];
                for(size_t k=0;k<pt_cnt[nnid];k++){
                  p_vesid_coordptr.key=ves_id[pt_dsp_nnid+k];
                  p_vesid_coordptr.data=&pt_coord[(pt_dsp_nnid+k)*COORD_DIM];
                  pair_vesid_coordptr.push_back(p_vesid_coordptr);
                }
              }
              std::sort(&pair_vesid_coordptr[0],&pair_vesid_coordptr[0]+pair_vesid_coordptr.size());

              near_ves_id.clear(); near_pt_cnt.clear(); near_pt_dsp.clear();
              near_pt_coord.resize(pair_vesid_coordptr.size()*COORD_DIM);
              for(size_t j=0;j<pair_vesid_coordptr.size();j++){
                near_pt_coord[j*COORD_DIM+0]=pair_vesid_coordptr[j].data[0];
                near_pt_coord[j*COORD_DIM+1]=pair_vesid_coordptr[j].data[1];
                near_pt_coord[j*COORD_DIM+2]=pair_vesid_coordptr[j].data[2];

                if(!near_ves_id.size() || near_ves_id.back()!=pair_vesid_coordptr[j].key){
                  near_ves_id.push_back(pair_vesid_coordptr[j].key);
                  near_pt_dsp.push_back(j);
                  near_pt_cnt.push_back(0);
                }
                near_pt_cnt.back()++;
              }
            }

            { // add pairs to near_pair
              size_t pt_cnt_i=pt_cnt[i];
              size_t pt_dsp_i=pt_dsp[i];
              for(size_t j=0;j<near_ves_id.size();j++){ // all near vesicles
                size_t near_pt_dsp_j=near_pt_dsp[j];
                size_t near_pt_cnt_j=near_pt_cnt[j];
                size_t near_ves_id_j=near_ves_id[j];
                for(size_t k=0;k<pt_cnt_i;k++){ // all points in node i
                  size_t l;
                  for(l=0;l<near_pt_cnt_j;l++){ // all near points in vesicle near_ves_id_j
                    real dx=near_pt_coord[(near_pt_dsp_j+l)*COORD_DIM+0]-pt_coord[(pt_dsp_i+k)*COORD_DIM+0];
                    real dy=near_pt_coord[(near_pt_dsp_j+l)*COORD_DIM+1]-pt_coord[(pt_dsp_i+k)*COORD_DIM+1];
                    real dz=near_pt_coord[(near_pt_dsp_j+l)*COORD_DIM+2]-pt_coord[(pt_dsp_i+k)*COORD_DIM+2];
                    real r2=dx*dx+dy*dy+dz*dz;
                    if(r2<r2_near){
                      p_ptid_vesid.key=pt_id[pt_dsp_i+k];
                      p_ptid_vesid.data=near_ves_id_j;
                      near_pair.push_back(p_ptid_vesid);
                      break;
                    }
                  }
                  flops+=(l+1)*8;
                }
              }
            }

          }
          cnt[tid]=near_pair.size();
        }
        pvfmm::Profile::Add_FLOP(flops);

        // combine results from all threads
        std::vector<size_t> dsp(omp_p); dsp[0]=0;
        pvfmm::omp_par::scan(&cnt[0], &dsp[0], cnt.size());
        near_pair.ReInit(cnt[omp_p-1]+dsp[omp_p-1]);
        #pragma omp parallel for
        for(size_t tid=0;tid<omp_p;tid++){
          memcpy(&near_pair[0]+dsp[tid],&near_pair_[tid][0],cnt[tid]*sizeof(pvfmm::par::SortPair<size_t, size_t>));
        }
      }

      { // Compute scatter_index, ves_interac_cnt, pt_interac_cnt

        // Sort near_pair
        pvfmm::Vector<pvfmm::par::SortPair<size_t, size_t> > sorted_near_pair;
        pvfmm::par::HyperQuickSort(near_pair, sorted_near_pair, comm);

        { // Set pt_interac_cnt
          assert(sorted_near_pair.Dim()>0);
          size_t n_pair=sorted_near_pair.Dim();
          size_t start_id=sorted_near_pair[0].key;
          size_t end_id=sorted_near_pair[n_pair-1].key;

          pvfmm::Vector<size_t> pt_interac_cnt_(end_id-start_id+1);
          pt_interac_cnt_.SetZero();
          #pragma omp parallel for
          for(size_t tid=0;tid<omp_p;tid++){ // set pt_interac_cnt_
            size_t a=((tid+0)*n_pair)/omp_p;
            size_t b=((tid+1)*n_pair)/omp_p;
            if(a>0) while(a<n_pair && sorted_near_pair[a-1].key==sorted_near_pair[a].key) a++;
            if(b>0) while(b<n_pair && sorted_near_pair[b-1].key==sorted_near_pair[b].key) b++;
            for(size_t i=a;i<b;){
              size_t key=sorted_near_pair[i].key;
              while(i<b && sorted_near_pair[i].key==key){
                pt_interac_cnt_[key-start_id]++;
                i++;
              }
            }
          }

          size_t num_pts=pt_interac_cnt_.Dim();
          pvfmm::Vector<size_t> mins(np);
          { // Fix last element in pt_interac_cnt_
            MPI_Allgather(&start_id, 1, pvfmm::par::Mpi_datatype<size_t>::value(),
                           &mins[0], 1, pvfmm::par::Mpi_datatype<size_t>::value(), comm);

            MPI_Status status;
            size_t send_cnt=0, recv_cnt=0;
            size_t r_partner=(rank>0?rank-1:np-1);
            size_t s_partner=(rank<np-1?rank+1:0);
            if(end_id==mins[s_partner]){
              num_pts--;
              send_cnt=pt_interac_cnt_[num_pts];
            }
            MPI_Sendrecv(&send_cnt, 1, pvfmm::par::Mpi_datatype<size_t>::value(), s_partner, 0,
                         &recv_cnt, 1, pvfmm::par::Mpi_datatype<size_t>::value(), r_partner, 0, comm, &status);
            pt_interac_cnt_[0]+=recv_cnt;
          }

          { // Repartition
            pvfmm::Vector<int> recv_size(np); recv_size.SetZero();
            pvfmm::Vector<int> send_size(np); send_size.SetZero();
            size_t p0=std::lower_bound(&mins[0],&mins[0]+np, (ves_id_offset+    0)*M_ves)-&mins[0];
            size_t p1=std::lower_bound(&mins[0],&mins[0]+np, (ves_id_offset+N_ves)*M_ves)-&mins[0];
            if(p0>0) p0--;
            for(size_t i=p0;i<p1;i++){ // set recv_size
              size_t a=(ves_id_offset+    0)*M_ves;
              size_t b=(ves_id_offset+N_ves)*M_ves;
              a           =std::max(a,mins[i  ]);
              if(i<np-1) b=std::min(b,mins[i+1]);
              recv_size[i]=b-a;
            }

            MPI_Alltoall(&recv_size[0], 1, pvfmm::par::Mpi_datatype<int>::value(),
                         &send_size[0], 1, pvfmm::par::Mpi_datatype<int>::value(), comm);
            pvfmm::Vector<int> recv_disp(np); recv_disp.SetZero();
            pvfmm::Vector<int> send_disp(np); send_disp.SetZero();
            pvfmm::omp_par::scan(&recv_size[0],&recv_disp[0],np);
            pvfmm::omp_par::scan(&send_size[0],&send_disp[0],np);
            assert(recv_size[np-1]+recv_disp[np-1]==N_ves*M_ves);
            assert(send_size[np-1]+send_disp[np-1]==num_pts);

            pt_interac_cnt.ReInit(recv_size[np-1]+recv_disp[np-1]);
            pvfmm::par::Mpi_Alltoallv_sparse<size_t>(&pt_interac_cnt_[0], &send_size[0], &send_disp[0],
                                                     &pt_interac_cnt [0], &recv_size[0], &recv_disp[0], comm);
          }
        }

        { // Compute scatter_index, ves_interac_cnt
          pvfmm::Vector<size_t> ves_id(sorted_near_pair.Dim());
          for(size_t i=0;i<ves_id.Dim();i++) ves_id[i]=sorted_near_pair[i].data;
          pvfmm::par::SortScatterIndex(ves_id, scatter_index, comm, &ves_id_offset);
          pvfmm::par::ScatterForward(ves_id, scatter_index, comm);

          { // Compute ves_interac_cnt
            ves_interac_cnt.ReInit(N_ves);
            ves_interac_cnt.SetZero();

            size_t start_id=ves_id[0];
            #pragma omp parallel for
            for(size_t tid=0;tid<omp_p;tid++){ // set pt_interac_cnt_
              size_t a=((tid+0)*ves_id.Dim())/omp_p;
              size_t b=((tid+1)*ves_id.Dim())/omp_p;
              if(a>0) while(a<ves_id.Dim() && ves_id[a-1]==ves_id[a]) a++;
              if(b>0) while(b<ves_id.Dim() && ves_id[b-1]==ves_id[b]) b++;
              for(size_t i=a;i<b;){
                size_t key=ves_id[i];
                while(i<b && ves_id[i]==key){
                  ves_interac_cnt[key-start_id]++;
                  i++;
                }
              }
            }
          }
        }
      }

      pt_interac_dsp.ReInit(pt_interac_cnt.Dim()); pt_interac_dsp[0]=0;
      pvfmm::omp_par::scan(&pt_interac_cnt[0], &pt_interac_dsp[0], pt_interac_cnt.Dim());

      ves_interac_dsp.ReInit(ves_interac_cnt.Dim()); ves_interac_dsp[0]=0;
      pvfmm::omp_par::scan(&ves_interac_cnt[0], &ves_interac_dsp[0], ves_interac_cnt.Dim());
    }
    pvfmm::Profile::Toc();
  }
  pvfmm::Profile::Toc();

  pvfmm::Profile::Tic("Interactions",&comm);
  { // Compute potential at near points
    size_t n_trg=0;
    n_trg+=pt_interac_cnt[N_ves*M_ves-1];
    n_trg+=pt_interac_dsp[N_ves*M_ves-1];

    pvfmm::Vector<real> trg_coord(n_trg*COORD_DIM);
    { // Set trg_coord
      #pragma omp parallel for
      for(size_t tid=0;tid<omp_p;tid++){
        size_t a=((tid+0)*N_ves)/omp_p;
        size_t b=((tid+1)*N_ves)/omp_p;
        for(size_t i=a;i<b;i++){
          // read each component of x
          const real* xk=x.getSubN_begin(i)+0*M_ves;
          const real* yk=x.getSubN_begin(i)+1*M_ves;
          const real* zk=x.getSubN_begin(i)+2*M_ves;
          for(size_t j=0;j<M_ves;j++){
            size_t pt_id=i*M_ves+j;
            size_t idx_offset=pt_interac_dsp[pt_id];
            for(size_t k=0;k<pt_interac_cnt[pt_id];k++){
              trg_coord[(idx_offset+k)*COORD_DIM+0]=xk[j];
              trg_coord[(idx_offset+k)*COORD_DIM+1]=yk[j];
              trg_coord[(idx_offset+k)*COORD_DIM+2]=zk[j];
            }
          }
        }
      }
    }

    // ScatterForward
    pvfmm::par::ScatterForward(trg_coord, scatter_index, comm);

    pvfmm::Vector<real> trg_veloc(scatter_index.Dim()*veloc_dim);
    { // Compute NearSingular
      pvfmm::Profile::Tic("NearSingular",&comm);
      StokesNearSingular(S, force, veloc_self, r_near,
          &ves_interac_cnt[0], &ves_interac_dsp[0], &trg_coord[0], &trg_veloc[0]);
      pvfmm::Profile::Toc();
    }

    // ScatterReverse
    pvfmm::par::ScatterReverse(trg_veloc, scatter_index, comm, n_trg);

    // Combine trg_veloc
    #pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){
      size_t a=((tid+0)*N_ves)/omp_p;
      size_t b=((tid+1)*N_ves)/omp_p;
      for(size_t i=a;i<b;i++){
        // read each component of velocity
        real* vxk=velocity.getSubN_begin(i)+0*M_ves;
        real* vyk=velocity.getSubN_begin(i)+1*M_ves;
        real* vzk=velocity.getSubN_begin(i)+2*M_ves;
        assert(veloc_dim==3);

        for(size_t j=0;j<M_ves;j++){
          size_t pt_id=i*M_ves+j;
          size_t idx_offset=pt_interac_dsp[pt_id];
          for(size_t k=0;k<pt_interac_cnt[pt_id];k++){
            vxk[j]+=trg_veloc[(idx_offset+k)*veloc_dim+0];
            vyk[j]+=trg_veloc[(idx_offset+k)*veloc_dim+1];
            vzk[j]+=trg_veloc[(idx_offset+k)*veloc_dim+2];
          }
        }
      }
    }

    if(0){ // visualization (near points)
      pvfmm::Profile::Tic("vis",&comm);
      typedef pvfmm::MPI_Node<real> Node_t;
      typedef pvfmm::MPI_Tree<Node_t> Tree_t;
      typedef Node_t::NodeData NodeData_t;

      Tree_t pt_tree(comm);
      NodeData_t node_data;
      node_data.dim=COORD_DIM;
      node_data.max_depth=MAX_DEPTH;
      node_data.max_pts=10000000000;

      // Set node_data.pt_coord, node_data.pt_value
      pvfmm::par::ScatterReverse(trg_coord, scatter_index, comm, n_trg);
      node_data.pt_coord=trg_coord;
      node_data.pt_value=trg_veloc;

      { // Map to unit cube
        real scale_x;
        real shift_x[COORD_DIM];
        pvfmm::Vector<real>& trg_coord=node_data.pt_coord;
        PVFMMBoundingBox(trg_coord.Dim()/COORD_DIM, &trg_coord[0], &scale_x, shift_x, comm);
        for(size_t i=0;i<trg_coord.Dim()/COORD_DIM;i++){
          trg_coord[i*COORD_DIM+0]=trg_coord[i*COORD_DIM+0]*scale_x+shift_x[0];
          trg_coord[i*COORD_DIM+1]=trg_coord[i*COORD_DIM+1]*scale_x+shift_x[1];
          trg_coord[i*COORD_DIM+2]=trg_coord[i*COORD_DIM+2]*scale_x+shift_x[2];
        }
      }
      pt_tree.Initialize(&node_data);
      pt_tree.Write2File("ptree");
      pvfmm::Profile::Toc();
    }
  }
  pvfmm::Profile::Toc();

  if(1){ // visualization
    pvfmm::Profile::Tic("vis",&comm);
    typedef pvfmm::MPI_Node<real> Node_t;
    typedef pvfmm::MPI_Tree<Node_t> Tree_t;
    typedef Node_t::NodeData NodeData_t;

    Tree_t pt_tree(comm);
    NodeData_t node_data;
    node_data.dim=COORD_DIM;
    node_data.max_depth=MAX_DEPTH;
    node_data.max_pts=10000000000;

    // Set node_data.pt_coord, node_data.pt_value
    node_data.pt_coord.ReInit(N_ves*M_ves*COORD_DIM);
    node_data.pt_value.ReInit(N_ves*M_ves*veloc_dim);
    #pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){
      size_t a=((tid+0)*N_ves)/omp_p;
      size_t b=((tid+1)*N_ves)/omp_p;
      for(size_t i=a;i<b;i++){
        // read each component of x
        const real* xk=x.getSubN_begin(i)+0*M_ves;
        const real* yk=x.getSubN_begin(i)+1*M_ves;
        const real* zk=x.getSubN_begin(i)+2*M_ves;

        // read each component of velocity
        const real* vxk=velocity.getSubN_begin(i)+0*M_ves;
        const real* vyk=velocity.getSubN_begin(i)+1*M_ves;
        const real* vzk=velocity.getSubN_begin(i)+2*M_ves;
        assert(veloc_dim==3);

        for(size_t j=0;j<M_ves;j++){
          size_t pt_id=i*M_ves+j;
          node_data.pt_coord[pt_id*COORD_DIM+0]=xk[j];
          node_data.pt_coord[pt_id*COORD_DIM+1]=yk[j];
          node_data.pt_coord[pt_id*COORD_DIM+2]=zk[j];

          node_data.pt_value[pt_id*veloc_dim+0]=vxk[j];
          node_data.pt_value[pt_id*veloc_dim+1]=vyk[j];
          node_data.pt_value[pt_id*veloc_dim+2]=vzk[j];
        }
      }
    }

    { // Map to unit cube
      real scale_x;
      real shift_x[COORD_DIM];
      pvfmm::Vector<real>& trg_coord=node_data.pt_coord;
      PVFMMBoundingBox(trg_coord.Dim()/COORD_DIM, &trg_coord[0], &scale_x, shift_x, comm);
      for(size_t i=0;i<trg_coord.Dim()/COORD_DIM;i++){
        trg_coord[i*COORD_DIM+0]=trg_coord[i*COORD_DIM+0]*scale_x+shift_x[0];
        trg_coord[i*COORD_DIM+1]=trg_coord[i*COORD_DIM+1]*scale_x+shift_x[1];
        trg_coord[i*COORD_DIM+2]=trg_coord[i*COORD_DIM+2]*scale_x+shift_x[2];
      }
    }
    pt_tree.Initialize(&node_data);
    pt_tree.Write2File("ptree");
    pvfmm::Profile::Toc();
  }
}

template <typename Vec_t>
void test(const MPI_Comm& comm){
  typedef typename Vec_t::scalars_type Sca_t;
  typedef typename Vec_t::array_type Arr_t;
  typedef OperatorsMats<Arr_t> Mats_t;
  typedef Surface<Sca_t, Vec_t> Surf_t;
  typedef VesInteraction<real> Interaction_t;
  typedef InterfacialVelocity<Surf_t, Interaction_t> IntVel_t;

  int rank, size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);
  size_t nVec=2/size;
  srand48(rank);

  //Set parameters
  Parameters<value_type> sim_par;
  sim_par.sh_order = 32;
  sim_par.rep_up_freq = 32;

  //Create vectors
  Vec_t x0(nVec, sim_par.sh_order); // coordinates
  Vec_t v0(nVec, sim_par.sh_order); // velocity
  Vec_t f0(nVec, sim_par.sh_order); // force
  int fLen = x0.getStride();
  { //Set coordinate values
    int imax(x0.getGridDim().first);
    int jmax(x0.getGridDim().second);
    real scal=pow((real)nVec*size,1.0/3.0);

    //std::vector<real> qx;
    //{ // compute legendre node points
    //  size_t p=x0.getShOrder();
    //  qx.resize(p+1);
    //  std::vector<real> qw(p+1);
    //  cgqf(p+1, 1, 0.0, 0.0, -1.0, 1.0, &qx[0], &qw[0]);
    //}

    for(size_t k=0;k<nVec;k++){
      real* x_k=x0.getSubN_begin(k)+0*fLen;
      real* y_k=x0.getSubN_begin(k)+1*fLen;
      real* z_k=x0.getSubN_begin(k)+2*fLen;

      real coord[3];
      coord[0]=scal*drand48();
      coord[1]=scal*drand48();
      coord[2]=scal*drand48();
      if(k==0){
        coord[0]=0.0;
        coord[1]=0.0;
        coord[2]=0.0;
      }else{
        coord[0]=0.0;
        coord[1]=1.0;
        coord[2]=0.0;
      }
      for(size_t i=0;i<imax;i++){
        //real cos_t=qx[i];
        real cos_t=cos((i+1)*M_PI/(imax+1));
        real sin_t=sqrt(1.0-cos_t*cos_t);
        for(size_t j=0;j<jmax;j++){
          x_k[j+i*jmax]=coord[0]+0.2*cos_t;
          y_k[j+i*jmax]=coord[1]+0.5*sin_t*sin(j*2*M_PI/jmax)+0.02*sin_t*cos(j*10*M_PI/jmax);
          z_k[j+i*jmax]=coord[2]+0.5*sin_t*cos(j*2*M_PI/jmax)+0.02*sin_t*sin(j*10*M_PI/jmax);
        }
      }

      real* vx_k=v0.getSubN_begin(k)+0*fLen;
      real* vy_k=v0.getSubN_begin(k)+1*fLen;
      real* vz_k=v0.getSubN_begin(k)+2*fLen;
      real* fx_k=f0.getSubN_begin(k)+0*fLen;
      real* fy_k=f0.getSubN_begin(k)+1*fLen;
      real* fz_k=f0.getSubN_begin(k)+2*fLen;
      for(size_t i=0;i<imax*jmax;i++){
        vx_k[i]=0;
        vy_k[i]=0;
        vz_k[i]=0;
        fx_k[i]=0;
        fy_k[i]=0;
        fz_k[i]=0;

        if(!rank && k<1){
          fx_k[i]=1;
          fy_k[i]=1;
          fz_k[i]=1;
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
          quad_weights_.size() * sizeof(real),
          Surf_t::device_type::MemcpyDeviceToDevice);
    }
    xv(S.getAreaElement(), f0, qforce);
    ax<Sca_t>(quad_weights_, qforce, qforce);
  }

  pvfmm::Profile::Tic("NearInteractions",&comm);
  NearInteractions(S, qforce, veloc_self, v0, comm);
  pvfmm::Profile::Toc();
}

int main(int argc, char** argv){
  typedef Vectors<real, DevCPU, the_cpu_device> VecCPU_t;
  MPI_Init(&argc, &argv);
  pvfmm::SetSigHandler();

  MPI_Comm comm=MPI_COMM_WORLD;
  pvfmm::Profile::Enable(true);
  test<VecCPU_t>(comm);
  pvfmm::Profile::print(&comm);

  MPI_Finalize();
}

