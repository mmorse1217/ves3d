#ifndef _RIGID_PARTICLE_H_
#define _RIGID_PARTICLE_H_

#include "Logger.h"
#include "Enums.h"
#include "Error.h"
#include "GMRESLinSolver.h"
#include "StokesVelocity.h"
#include "SphericalHarmonics.h"

template<typename SurfContainer>
class RigidParticle
{
  public:
      typedef typename SurfContainer::value_type value_type;
      typedef typename SurfContainer::Arr_t Arr_t;
      typedef typename SurfContainer::Sca_t Sca_t;
      typedef typename SurfContainer::Vec_t Vec_t;
      typedef StokesVelocity<value_type> Stokes_t;
      typedef typename pvfmm::Vector<value_type> PVFMMVec;
      //typedef typename Stokes_t::PVFMMVec PVFMMVec;

      RigidParticle(int sh_order, int sh_order_up, int sh_order_up_self);
      ~RigidParticle();
      void EvalPotential(int num_target_points, value_type* target_address, value_type* target_potential);
      void Solve();
      value_type* GetSamplePoints(int& num_sample_points);
      void SetBoundaryData(value_type* boundary_data_address);
    
  private:
      //SurfContainer surface_;
      Vec_t position_;  
      Vec_t far_vel_;  //rhs
      Vec_t density_;  //unknown density
      Vec_t cen_;
      Vec_t tmp1_;
      Vec_t tmp2_;
      //value_type* t_vel_; //unknown translation vel
      //value_type* r_vel_; //unknown angular vel
      //value_type* cm_;    //center of mass, pole location
      Arr_t t_vel_; //unknown translation vel
      Arr_t r_vel_; //unknown angular vel
      Arr_t cm_;    //center of mass, pole location

      GMRESLinSolver<value_type> linear_solver_gmres_;

      int sh_order_;
      int sh_order_up_;
      int sh_order_up_self_;

      Stokes_t stokes_;

      PVFMMVec DLMatrix_;
      Sca_t weights_;
      //PVFMMVec weights_;

      Error_t operator()(const Vec_t density, const Arr_t t_vel, const Arr_t r_vel, 
                         Vec_t &potential, Arr_t &force_int, Arr_t &torque_int) const;
      static Error_t JacobiImplicitApply(const GMRESLinSolver<value_type> *o, const value_type *x, value_type *y, int tmp_trash=-1);
};

#include "RigidParticle.cc"

#endif // _RIGID_PARTICLE_H_
