#ifdef HAVE_PVFMM
#ifndef _NEAR_SINGULAR_H_
#define _NEAR_SINGULAR_H_

#include <mpi.h>
#include <vector.hpp>

template<typename Surf_t>
class NearSingular{

    typedef typename Surf_t::Vec_t Vec_t;
    typedef typename Surf_t::value_type Real_t;
    typedef typename pvfmm::Vector<Real_t> PVFMMVec_t;

  public:

    NearSingular(int sh_order, Real_t box_size=-1, MPI_Comm c=MPI_COMM_WORLD);

    void SetSrcCoord(const PVFMMVec_t& src_coord);
    void SetSurfaceVel(const PVFMMVec_t* S_vel);
    void SetDensitySL(const PVFMMVec_t* qforce_single=NULL);
    void SetDensityDL(const PVFMMVec_t* qforce_double=NULL, const PVFMMVec_t* force_double=NULL);

    void SetTrgCoord(Real_t* trg_coord, size_t N);

    void SubtractDirect(PVFMMVec_t& vel_fmm);
    const PVFMMVec_t& operator()(bool update=true);

  private:

    NearSingular(const NearSingular &);
    NearSingular& operator=(const NearSingular &);

    void SetupCoordData();

    void VelocityScatter(PVFMMVec_t& trg_vel);

    struct QuadraticPatch{
      public:

      QuadraticPatch(){}

      QuadraticPatch(Real_t* x, int dof_);

      void eval(Real_t x, Real_t y, Real_t* val);

      int project(Real_t* t_coord_j, Real_t& x, Real_t&y);

      static void patch_mesh(Real_t* patch_value_, size_t sh_order, pvfmm::Matrix<Real_t>& M_interp, size_t k_, const Real_t* sx_value);

      private:

      void grad(Real_t x, Real_t y, Real_t* val);

      int dof;
      std::vector<Real_t> coeff;
    };

    struct{
      Real_t r_near;
      Real_t bbox[4]; // {s,x,y,z} : scale, shift

      PVFMMVec_t            near_trg_coord;
      pvfmm::Vector<size_t> near_trg_cnt;
      pvfmm::Vector<size_t> near_trg_dsp;

      pvfmm::Vector<size_t> near_ves_pt_id;   // vesicle point index
      pvfmm::Vector<size_t> near_trg_pt_id;   // target index at original location
      pvfmm::Vector<size_t> near_trg_scatter; // Scatter velocity to original location

      pvfmm::Vector<char>   is_surf_pt;       // If a target point is a surface point
      pvfmm::Vector<char>   is_extr_pt;       // If a target point is an exterior point

      PVFMMVec_t            proj_patch_param; // Projection patch prameter coordinates
      PVFMMVec_t            proj_coord;       // Projection coordinates (x,y,z)
    } coord_setup;

    const PVFMMVec_t* S;
    const PVFMMVec_t* qforce_single;
    const PVFMMVec_t* qforce_double;
    const PVFMMVec_t* force_double;
    const PVFMMVec_t* S_vel;
    Real_t box_size_;
    MPI_Comm comm;

    int sh_order_;
    pvfmm::Matrix<Real_t> M_patch_interp;

    PVFMMVec_t T;
    PVFMMVec_t vel_direct;
    PVFMMVec_t vel_interp;
    PVFMMVec_t vel_surfac;

    enum UpdateFlag{
      UpdateNone      = 0,
      UpdateSrcCoord  = 1,
      UpdateTrgCoord  = 2,
      UpdateDensitySL = 4,
      UpdateDensityDL = 8,
      UpdateSurfaceVel=16,
      UpdateAll       =31
    };
    unsigned int update_direct;
    unsigned int update_interp;
    unsigned int update_setup ;

    static const size_t INTERP_DEG=8;
};

#include "NearSingular.cc"

#endif // _NEAR_SINGULAR_H_
#endif
