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

#include "Parameters.h"
#include <sstream>

template<typename P>
class ParametersTest
{
  public:
    bool operator()(const P &p);
    bool TestStream(const P &p);
};

template<typename P>
bool ParametersTest<P>::operator()(const P &p)
{
    bool res = TestStream(p);

    if (res)
	COUT(emph<<"Parameters test passed"<<emph);
    else
	COUT(alert<<"Parameters test failed"<<alert);

    return res;
}

template<typename P>
bool ParametersTest<P>::TestStream(const P &p)
{
    COUT(". TestStream");
    std::stringstream s1;
    p.pack(s1, P::Streamable::ASCII);

    P pc(s1, P::Streamable::ASCII);
    ASSERT(p.n_surfs == pc.n_surfs					, "incorrect n_surfs");
    ASSERT(p.sh_order == pc.sh_order					, "incorrect sh_order");
    ASSERT(p.filter_freq == pc.filter_freq				, "incorrect filter_freq");
    ASSERT(p.upsample_freq == pc.upsample_freq				, "incorrect upsample_freq");
    ASSERT(p.rep_filter_freq == pc.rep_filter_freq			, "incorrect rep_filter_freq");
    ASSERT(p.rep_upsample == pc.rep_upsample				, "incorrect rep_upsample");
    ASSERT(p.bending_modulus == pc.bending_modulus			, "incorrect bending_modulus");
    ASSERT(p.viscosity_contrast == pc.viscosity_contrast		, "incorrect viscosity_contrast");
    ASSERT(p.position_solver_iter == pc.position_solver_iter		, "incorrect position_solver_iter");
    ASSERT(p.tension_solver_iter == pc.tension_solver_iter		, "incorrect tension_solver_iter");
    ASSERT(p.position_solver_restart == pc.position_solver_restart	, "incorrect position_solver_restart");
    ASSERT(p.tension_solver_restart == pc.tension_solver_restart	, "incorrect tension_solver_restart");
    ASSERT(p.position_solver_tol == pc.position_solver_tol		, "incorrect position_solver_tol");
    ASSERT(p.tension_solver_tol == pc.tension_solver_tol		, "incorrect tension_solver_tol");
    ASSERT(p.time_horizon == pc.time_horizon				, "incorrect time_horizon");
    ASSERT(p.ts == pc.ts						, "incorrect ts");
    ASSERT(p.time_tol == pc.time_tol					, "incorrect time_tol");
    ASSERT(p.time_iter_max == pc.time_iter_max				, "incorrect time_iter_max");
    ASSERT(p.solve_for_velocity == pc.solve_for_velocity		, "incorrect solve_for_velocity");
    ASSERT(p.scheme == pc.scheme					, "incorrect scheme");
    ASSERT(p.time_precond == pc.time_precond				, "incorrect time_precond");
    ASSERT(p.bg_flow == pc. bg_flow					, "incorrect  bg_flow");
    ASSERT(p.singular_stokes == pc. singular_stokes			, "incorrect  singular_stokes");
    ASSERT(p.rep_maxit == pc.rep_maxit					, "incorrect rep_maxit");
    ASSERT(p.rep_ts == pc.rep_ts					, "incorrect rep_ts");
    ASSERT(p.rep_tol == pc.rep_tol					, "incorrect rep_tol");
    ASSERT(p.bg_flow_param == pc.bg_flow_param				, "incorrect bg_flow_param");
    ASSERT(p.interaction_upsample == pc.interaction_upsample		, "incorrect interaction_upsample");
    ASSERT(p.checkpoint == pc.checkpoint        			, "incorrect checkpoint");
    ASSERT(p.checkpoint_stride == pc.checkpoint_stride			, "incorrect checkpoint_stride");
    ASSERT(p.init_file_name == pc.init_file_name			, "incorrect init_file_name");
    ASSERT(p.cntrs_file_name == pc.cntrs_file_name			, "incorrect cntrs_file_name");
    ASSERT(p.checkpoint_file_name == pc.checkpoint_file_name		, "incorrect checkpoint_file_name");
    ASSERT(p.load_checkpoint == pc.load_checkpoint      		, "incorrect load_checkpoint");
    ASSERT(p.error_factor == pc.error_factor				, "incorrect error_factor");
    ASSERT(p.num_threads == pc.num_threads				, "incorrect num_threads");
    ASSERT(p.excess_density == pc.excess_density			, "incorrect excess_density");
    for (int ii(0); ii<DIM; ++ii)
        ASSERT(p.gravity_field[ii] == pc.gravity_field[ii]		, "incorrect gravity_field");

    return true;
}

int main(int, char**){

    int argc(15);
    char *argv[] = {"execname",
		    "--n-surfs", "5",
		    "--sh-order", "13",
		    "-o", "out.txt",
		    "-l", "a.txt",
		    "--rep-upsample",
		    "--interaction-upsample",
		    "--excess-density", "5",
		    "--gravity-field", "1.1 2e1 -3",
    };
    Parameters<double> p(argc, argv);
    ParametersTest<Parameters<double> > p_test;
    p_test(p);
}
