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

template<typename DT, const DT &DEVICE>
Simulation<DT,DEVICE>::Simulation(const Param_t &ip) :
    input_params_(ip),
    load_checkpoint_(false),
    Mats_(NULL),
    vInf_(NULL),
    ksp_(NULL),
    interaction_(NULL),
    timestepper_(NULL)
{
    CHK(prepare_run_params());
}

template<typename DT, const DT &DEVICE>
Simulation<DT,DEVICE>::~Simulation()
{
    cleanup_run();
}

template<typename DT, const DT &DEVICE>
Error_t Simulation<DT,DEVICE>::Run()
{
    setup_basics();
    if (load_checkpoint_)
        setup_from_checkpoint();
    else
        setup_from_options();

    INFO("Run options:\n"<<run_params_);
    CHK(timestepper_->Evolve());
    return ErrorEvent::Success;
}

template<typename DT, const DT &DEVICE>
Error_t Simulation<DT,DEVICE>::setup_basics(){
    if (run_params_.num_threads>0){
        INFO("Setting OMP num threads to "<<run_params_.num_threads);
        omp_set_num_threads(run_params_.num_threads);
    } else {
        INFO("OMP max threads is "<<omp_get_max_threads());
        omp_set_num_threads(omp_get_max_threads());
    }

    //Reading Operators From File
    Mats_ = new Mats_t(true /*readFromFile*/, run_params_);

    //Setting the background flow
    CHK(BgFlowFactory(run_params_, &vInf_));

#ifdef HAS_PETSC
    ksp_ = new ParallelLinSolverPetsc<real_t>(VES3D_COMM_WORLD);
#endif

#ifdef HAVE_PVFMM
    interaction_ = new Inter_t(&PVFMMEval, &PVFMMDestroyContext<real_t>);
#else
    interaction_ = new Inter_t(&StokesAlltoAll);
#endif

    return ErrorEvent::Success;
}

template<typename DT, const DT &DEVICE>
Error_t Simulation<DT,DEVICE>::setup_from_checkpoint(){

    timestepper_ = new Evolve_t(&run_params_, *Mats_, vInf_, NULL, interaction_, NULL, ksp_);
    timestepper_->unpack(checkpoint_data_, Streamable::ASCII);
    return ErrorEvent::Success;
}

template<typename DT, const DT &DEVICE>
Error_t Simulation<DT,DEVICE>::setup_from_options(){
    //Initial vesicle positions
    Vec_t x0(run_params_.n_surfs, run_params_.sh_order);

    //reading the prototype form file
    DataIO myIO;
    myIO.ReadData( FullPath(run_params_.init_file_name), x0, DataIO::ASCII, 0, x0.getSubLength());

    int nproc(1), rank(0);
#ifdef HAVE_PVFMM
    MPI_Comm_size(VES3D_COMM_WORLD, &nproc);
    MPI_Comm_rank(VES3D_COMM_WORLD, &rank);
#endif

    //reading centers file
    if (run_params_.cntrs_file_name.size()){
	INFO("Reading centers from file");
	Arr_t cntrs(DIM * run_params_.n_surfs * nproc);

	myIO.ReadData( FullPath(run_params_.cntrs_file_name), cntrs, DataIO::ASCII, 0, cntrs.size());

	INFO("Populating the initial configuration using centers");
	Arr_t my_centers(DIM * run_params_.n_surfs);
    	cntrs.getDevice().Memcpy(my_centers.begin(),
	    cntrs.begin() + rank * DIM * run_params_.n_surfs,
	    DIM * run_params_.n_surfs * sizeof(Arr_t::value_type),
	    Arr_t::device_type::MemcpyDeviceToDevice);
        Populate(x0, my_centers);
    };

    timestepper_ = new Evolve_t(&run_params_, *Mats_, vInf_, NULL, interaction_, NULL, ksp_, &x0);

    return ErrorEvent::Success;
}

template<typename DT, const DT &DEVICE>
Error_t Simulation<DT,DEVICE>::cleanup_run()
{
    INFO("Cleaning up after run");
    delete Mats_;        Mats_	  = NULL;
    delete vInf_;        vInf_	  = NULL;
    delete ksp_;         ksp_	  = NULL;
    delete interaction_; interaction_ = NULL;
    delete timestepper_; timestepper_ = NULL;

    load_checkpoint_ = false;
    checkpoint_data_.str("");
    checkpoint_data_.clear();
    return ErrorEvent::Success;
}

template<typename DT, const DT &DEVICE>
Error_t Simulation<DT,DEVICE>::prepare_run_params()
{
    if (input_params_.load_checkpoint != ""){
        std::string fname = FullPath(input_params_.load_checkpoint);
        INFO("Loading checkpoint file "<<fname);
        DataIO::SlurpFile(fname.c_str(), checkpoint_data_);
        load_checkpoint_ = true;
    } else {
        input_params_.pack(checkpoint_data_, Streamable::ASCII);
        load_checkpoint_ = false;
    }

    //consume the first part of checkpoint (see Monitor for how
    //it is dumped
    run_params_.unpack(checkpoint_data_, Streamable::ASCII);
    run_params_.overrideNonState(input_params_);

    return ErrorEvent::Success;
}
