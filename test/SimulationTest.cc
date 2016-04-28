/**
 * @file
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @revision $Revision$
 * @tags $Tags$
 * @date $Date$
 *
 * @brief
 */

/*
 * Copyright (c) 2014, Abtin Rahimian
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include "ves3d_simulation.h"

typedef Device<CPU> Dev;
extern const Dev cpu(0);
typedef Simulation<Dev, cpu> Sim_t;
typedef Sim_t::Param_t Param_t;

int main(int argc, char **argv)
{
#ifdef HAS_PETSC
  PetscInitialize(&argc, &argv, NULL, NULL);
#endif

  DictString_t dict;
  int nproc(1), rank(0);
#ifdef HAS_MPI
  // Adding nproc and rank to template expansion dictionary
  MPI_Comm_size(VES3D_COMM_WORLD, &nproc);
  MPI_Comm_rank(VES3D_COMM_WORLD, &rank);
#endif
  std::stringstream snp, sr;
  snp<<nproc; sr<<rank;
  dict["nprocs"] = snp.str();
  dict["rank"]   = sr.str();

  {
    Param_t params;
    real_t ts(1e-4);
    params.n_surfs              = 2;
    params.sh_order             = 6;
    params.ts                   = ts;
    params.time_horizon	    	= 2*ts;
    params.scheme               = GloballyImplicit;
    params.time_precond         = DiagonalSpectral;
    params.pseudospectral       = false;
    params.bg_flow_param        = 5e-2;
    params.checkpoint           = true;
    params.checkpoint_stride    = ts;
    params.rep_maxit            = 100;
    params.checkpoint_file_name = "SimulationTest_a_{{rank}}_{{time_idx}}.chk";
    params.init_file_name       = "precomputed/dumbbell_{{sh_order}}_{{precision}}.txt";
    params.cntrs_file_name      = "precomputed/shear_centers.txt";
    params.rep_ts               = ts;
    params.expand_templates(&dict);

    Sim_t sim1(params);
    CHK(sim1.Run());
    const Sim_t::Vec_t &xref(sim1.time_stepper()->S_->getPosition());

    INFO("Setting up the second simulation");
    params.load_checkpoint      = "test/SimulationTest_a_{{rank}}_00001.chk";
    params.expand_templates(&dict);
    Sim_t sim2(params);

    sim2.run_params()->checkpoint_file_name = "SimulationTest_b_{{rank}}_{{time_idx}}.chk";
    sim2.run_params()->time_horizon = ts;
    sim2.run_params()->expand_templates(&dict);
    CHK(sim2.Run());

    //final state of sim1 and sim2 should match
    const Sim_t::Vec_t &x(sim2.time_stepper()->S_->getPosition());
    Sim_t::Vec_t err;
    err.replicate(x);
    axpy((real_t) -1.0, xref, x, err);
    real_t maxerr = MaxAbs(err);
    ASSERT(maxerr<5e-7, "inconsistent state after start from checkpoint error="<<maxerr);
  }

  {
    //checkpoints of a and b should match
    std::string s, fa("SimulationTest_a_{{rank}}_00002.chk"), fb("SimulationTest_b_{{rank}}_00001.chk");
    expand_template(&fa, dict);
    expand_template(&fb, dict);

    std::stringstream a,b;
    Sim_t::Vec_t av,bv,e;
    DataIO::SlurpFile(fa.c_str(), a);
    while (a.good()){
      a>>s;
      if(s=="VECTORS"){
	a.seekg(-7,std::ios_base::cur);
	break;
      }
    }
    av.unpack(a, Streamable::ASCII);

    DataIO::SlurpFile(fb.c_str(), b);
    while (b.good()){
      b>>s;
      if(s=="VECTORS"){
	b.seekg(-7,std::ios_base::cur);
	break;
      }
    }
    bv.unpack(b, Streamable::ASCII);
    e.replicate(av);
    axpy((real_t) -1, av, bv, e);
    real_t maxerr = MaxAbs(e);
    ASSERT(maxerr<5e-7, "inconsistent checkpoints, error="<<maxerr);
  }

  COUT(emph<<"** SimulationTest passed **"<<emph);
#ifdef HAS_PETSC
  PetscFinalize();
#endif

  return 0;
}
