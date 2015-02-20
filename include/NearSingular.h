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

    NearSingular(MPI_Comm c=MPI_COMM_WORLD);

    void SetSrcCoord(const Surf_t& S);
    void SetSurfaceVel(const Vec_t& S_vel);
    void SetDensitySL(const PVFMMVec_t* qforce_single=NULL);
    void SetDensityDL(const PVFMMVec_t* qforce_double=NULL);

    void SetTrgCoord(Real_t* trg_coord, size_t N);

    void SubtractDirect(PVFMMVec_t& vel_fmm);
    const PVFMMVec_t& operator()(bool update=true);

  private:

    NearSingular(const NearSingular &);
    NearSingular& operator=(const NearSingular &);

    void SetupCoordData();

    void VelocityScatter(PVFMMVec_t& trg_vel);

    struct{
      Real_t r_near;
      Real_t bbox[4]; // {s,x,y,z} : scale, shift

      PVFMMVec_t near_trg_coord;
      pvfmm::Vector<size_t> near_trg_cnt;
      pvfmm::Vector<size_t> near_trg_dsp;
      pvfmm::Vector<size_t> trg_scatter; // Scatter velocity to original location
      pvfmm::Vector<size_t> trg_pt_id;   // target index at original location
    } coord_setup;

    const Surf_t* S;
    const PVFMMVec_t* qforce_single;
    const PVFMMVec_t* qforce_double;
    const Vec_t* S_vel;

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

    MPI_Comm comm;
};

#include "NearSingular.cc"

#endif // _NEAR_SINGULAR_H_
#endif
