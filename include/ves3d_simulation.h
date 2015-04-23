/**
 * @file
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @revision $Revision$
 * @tags $Tags$
 * @date $Date$
 *
 * @brief Setup the simulation based on the given option
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

#ifndef _VES3D_SIMULATION_H_
#define _VES3D_SIMULATION_H_

#include "Error.h"
#include "Parameters.h"
#include "EvolveSurface.h"
#include "ves3d_common.h"
#include "Logger.h"
#include "DataIO.h"

#include <fstream>
#include <sstream>

#ifdef HAS_PETSC
#include "ParallelLinSolver_Petsc.h"
#endif

#ifdef HAVE_PVFMM
#include "PVFMMInterface.h"
#endif

template<typename DT, const DT &DEVICE>
class Simulation{
  public:
    typedef EvolveSurface<real_t, DT, DEVICE> Evolve_t;
    typedef typename Evolve_t::Params_t Param_t;
    typedef typename Evolve_t::Arr_t Arr_t;
    typedef typename Evolve_t::Vec_t Vec_t;
    typedef typename Evolve_t::Interaction_t Inter_t;
    typedef typename Evolve_t::Mats_t Mats_t;
    typedef BgFlowBase<Vec_t> Flow_t;
    typedef ParallelLinSolver<real_t> LinSol_t;

    Simulation(const Param_t &ip);
    ~Simulation();

    Error_t Run();

  private:
    const Param_t &input_params_;
    Param_t run_params_;
    std::stringstream checkpoint_data_;

    bool load_checkpoint_;
    Mats_t *Mats_;
    Flow_t *vInf_;
    LinSol_t *ksp_;
    Inter_t *interaction_;
    Evolve_t *timestepper_;

    Error_t setup_basics();
    Error_t setup_from_options();
    Error_t setup_from_checkpoint();
    Error_t cleanup_run();
    Error_t prepare_run_params();
};

#include "ves3d_simulation.cc"

#endif /* _VES3D_SIMULATION_H_ */
