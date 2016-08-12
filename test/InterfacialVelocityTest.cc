/**
 * @file
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @revision $Revision$
 * @tags $Tags$
 * @date $Date$
 *
 * @brief unit test
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

#include "InterfacialVelocity.h"
#include "ves3d_common.h"
#include "Vectors.h"
#include "VesInteraction.h"
#include "Surface.h"
#include "BgFlow.h"
#include "ParallelLinSolver_Petsc.h"
#include "VesicleProps.h"
#include "ves3d_simulation.h"

#define DT CPU
typedef Device<DT> Dev;

extern const Dev the_device(0);

int main(int argc, char** argv){
    VES3D_INITIALIZE(&argc, &argv, NULL, NULL);

    typedef Simulation<Dev, the_device> Sim_t;
    typedef Sim_t::Param_t Param_t;
    typedef Sim_t::Evolve_t Evolve_t;
    typedef Evolve_t::IntVel_t IntVel_t;
    typedef Evolve_t::value_type value_type;
    typedef Evolve_t::Sca_t Sca_t;
    typedef Evolve_t::Vec_t Vec_t;

    SET_ERR_CALLBACK(&cb_abort);

    // Parameter
    int p(16), nves(2);
    real_t tol(1e-3), ts(1e-4);
    Param_t sim_par;
    sim_par.n_surfs               = nves;
    sim_par.sh_order              = p;
    sim_par.filter_freq           = p;
    sim_par.upsample_freq         = p;
    sim_par.rep_filter_freq       = p;
    sim_par.ts                    = ts;
    sim_par.time_precond          = DiagonalSpectral;
    sim_par.excess_density        = 1e3;
    sim_par.shape_gallery_file    = "precomputed/shape_gallery_{{sh_order}}.txt";
    sim_par.vesicle_geometry_file = "precomputed/lattice_geometry_rand_spec.txt";
    sim_par.expand_templates();
    COUT(sim_par);

    {
        Vec_t v1(nves,p), v2(nves,p), v3(nves,p);
        Sca_t t1(nves,p), t2(nves,p), t3(nves,p);
        Sim_t sim(sim_par);

        sim.setup();
        Evolve_t *E(sim.time_stepper());
        E->ReinitInterfacialVelocity();
        IntVel_t *F(E->F_);

        INFO("Globally implicit update (solve for velocity)");
        F->Prepare(GloballyImplicit);

        INFO("Check linearity");
        value_type zero(0), one(1.0), two(2.0);
        Vec_t::getDevice().Memset(v1.begin(), zero, v1.mem_size());
        Sca_t::getDevice().Memset(t1.begin(), zero, t1.mem_size());
        F->ImplicitMatvecPhysical(v1,t1);
        value_type err(MaxAbs(v1)+MaxAbs(t1));
        ASSERT(err<1e-15,"matvec is not linear");

        fillRand(v1);fillRand(v1);axpy(two,v1,v2,v3);
        fillRand(t1);fillRand(t1);axpy(two,t1,t2,t3);
        F->ImplicitMatvecPhysical(v1,t1);
        F->ImplicitMatvecPhysical(v2,t2);
        F->ImplicitMatvecPhysical(v3,t3);
        axpy(two,v1,v2,v1);
        axpy(two,t1,t2,t1);

        axpy(-one,v1,v3,v1);
        axpy(-one,t1,t3,t1);
        err = MaxAbs(v1)+MaxAbs(t1);
        ASSERT(err<1e-15,"matvec is not linear");

        Vec_t dxx(nves,p), dxv(nves,p);

        sim.run_params()->solve_for_velocity=true;
        F->updateImplicit(*E->S_,ts,dxv);

        INFO("Globally implicit update (solve for position)");
        sim.run_params()->solve_for_velocity=false;
        F->updateImplicit(*E->S_,ts,dxx);

        axpy(-1.0,dxv,dxx,dxv);
        err = MaxAbs(dxv);
        COUT("Difference in solving for position or velocity: "<<err);
        ASSERT(err<tol,"large error between velocity and position solve");
    }
    VES3D_FINALIZE();
}
