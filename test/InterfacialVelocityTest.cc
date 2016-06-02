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

#define DT CPU
typedef Device<DT> Dev;

extern const Dev the_device(0);

int main(int argc, char** argv){
    PetscInitialize(&argc, &argv, NULL, NULL);

    typedef Scalars<real_t, Dev, the_device> Sca_t;
    typedef Vectors<real_t, Dev, the_device> Vec_t;
    typedef typename Sca_t::array_type Arr_t;
    typedef Surface<Sca_t,Vec_t> Sur_t;
    typedef VesInteraction<real_t> Interaction_t;
    typedef InterfacialVelocity<Sur_t, Interaction_t> IntVel_t;
    typedef OperatorsMats<Arr_t> Mats_t;
    typedef ParallelLinSolverPetsc<real_t> PSol_t;
    typedef VesicleProperties<Arr_t> VProp_t;

    SET_ERR_CALLBACK(&cb_abort);

    // Parameter
    int p(12), nvec(2);
    real_t tol(4e-5), ts(1e-2);
    Parameters<real_t> sim_par;
    sim_par.n_surfs	    = nvec;
    sim_par.sh_order	    = p;
    sim_par.filter_freq     = p;
    sim_par.upsample_freq   = p;
    sim_par.rep_filter_freq = p;
    sim_par.ts              = ts;
    sim_par.time_precond    = DiagonalSpectral;
    sim_par.init_file_name  = "precomputed/biconcave_ra85_{{sh_order}}";
    sim_par.cntrs_file_name = "precomputed/lattice_centers_rand_50.txt";
    sim_par.expand_templates();
    COUT(sim_par);

    {
        DataIO myIO;
        Vec_t x(nvec, p), dxx(nvec, p), dxv(nvec, p);
        Sca_t tension(nvec, p);
        myIO.ReadData(FullPath(sim_par.init_file_name), x, DataIO::ASCII, 0, x.getSubLength());

        Arr_t cntrs(DIM * nvec);
        myIO.ReadData( FullPath(sim_par.cntrs_file_name), cntrs );
        Populate(x, cntrs);

        //Reading Operators From
        Mats_t Mats(true /* readFromFile */, sim_par);

        //Making The Surface, And Time Stepper
        Sur_t S(p, Mats, &x);

        //Setting the background flow
        BgFlowBase<Vec_t> *vInf(NULL);
        CHK(BgFlowFactory(sim_par, &vInf));

        Interaction_t interaction(&StokesAlltoAll);
        PSol_t ksp(VES3D_COMM_WORLD);
        VProp_t vp;
        vp.setFromParams(sim_par);
        IntVel_t F(S, interaction, Mats, sim_par, vp, *vInf, &ksp);

        // testing
        // INFO("Explicit update");
        // F.updateJacobiExplicit(ts);

        // INFO("Gauss-Seidel update");
        // F.updateJacobiGaussSeidel(ts);

        INFO("Globally implicit update (solve for velocity)");
        sim_par.solve_for_velocity=true;
        F.updateImplicit(S,ts,dxv);

        INFO("Globally implicit update (solve for position)");
        sim_par.solve_for_velocity=false;
        F.updateImplicit(S,ts,dxx);

        axpy(-1.0,dxv,dxx,dxv);
        real_t err(MaxAbs(dxv));
        ASSERT(err<5e-5,"large error between velocity and position solve");

        delete vInf;
    }
    PetscFinalize();
}
