#ifndef _FIXEDBOUNDARY_H_
#define _FIXEDBOUNDARY_H_

#include <common/kernel3d.hpp>
#include <bdry3d/patch_surf_blended.hpp>
#include <bdry3d/patch_surf_face_map.hpp>
#include <bdry3d/patch_surf_analytic.hpp>
#include <bie3d/solver_utils.hpp>
#include <bie3d/solver_gmres_double_layer.hpp>
#include <bie3d/markgrid.hpp>
#include <bie3d/evaluator_qbkix.hpp>
#include <sampling.hpp>
#include <common/nummat.hpp>
#include <common/vtk_writer.hpp>
#include <common/stats.hpp>

using namespace Ebi;


//template<typename SurfContainer, typename Interaction>
class FixedBoundary
{
  public:
      FixedBoundary();
      ~FixedBoundary();
      void EvalPotential(int num_target_points, double* target_address, double* target_potential);
      void Solve();

      PatchSurfFaceMap* surface;
      SolverGMRESDoubleLayer* solver;
      MPI_Comm comm;

      Vec solved_density;
      Vec boundary_data;

      double* GetSamplePoints(int& num_sample_points);
      void SetBoundaryData(double* boundary_data_address);
  
      //Vec computed_potential;
      //Vec targets;
  private:
};

#include "FixedBoundary.cc"

#endif // _FIXEDBOUNDARY_H_
