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
    load_checkpoint_(false),
    ves_props_(NULL),
    Mats_(NULL),
    vInf_(NULL),
    ksp_(NULL),
    interaction_(NULL),
    timestepper_(NULL)
{
    CHK(prepare_run_params(ip));
}

template<typename DT, const DT &DEVICE>
Simulation<DT,DEVICE>::Simulation(int argc, char** argv, const DictString_t *dict) :
    load_checkpoint_(false),
    ves_props_(NULL),
    Mats_(NULL),
    vInf_(NULL),
    ksp_(NULL),
    interaction_(NULL),
    timestepper_(NULL)
{
    CHK(prepare_run_params(argc,argv,dict));
}

template<typename DT, const DT &DEVICE>
Simulation<DT,DEVICE>::~Simulation()
{
    cleanup_run();
}

template<typename DT, const DT &DEVICE>
Error_t Simulation<DT,DEVICE>::setup()
{
    Error_t ierr;
    ierr = setup_basics();
    if (ierr) return ierr;

    if (load_checkpoint_)
        return setup_from_checkpoint();
    else
        return setup_from_options();
}

template<typename DT, const DT &DEVICE>
Error_t Simulation<DT,DEVICE>::Run()
{
    CHK(setup());
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
    omp_set_nested(0);

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
Error_t Simulation<DT,DEVICE>::setup_from_options()
{
    DataIO io;
    int nves(run_params_.n_surfs);
    Vec_t x0(nves, run_params_.sh_order);

    // loading *all* shapes
    ASSERT(run_params_.shape_gallery_file.size()>0,"shape gallery file is required");
    std::string fname(FullPath(run_params_.shape_gallery_file));
    INFO("Reading shape gallery file "<<fname);
    std::vector<value_type> shapes;
    io.ReadDataStl(fname, shapes, DataIO::ASCII);

    int nshapes(shapes.size()/x0.getStride()/VES3D_DIM);
    INFO("Loaded "<<nshapes<<" shape(s)");

    // load centers and transformations for current mpi process
    ASSERT(run_params_.vesicle_geometry_file.size()>0,"geometry file is required");
    fname = FullPath(run_params_.vesicle_geometry_file);
    INFO("Reading geometry file "<<fname);
    std::vector<value_type> all_geo_spec;
    io.ReadDataStl(fname, all_geo_spec, DataIO::ASCII);

    int nproc(1), rank(0);
#ifdef HAVE_PVFMM
    MPI_Comm_size(VES3D_COMM_WORLD, &nproc);
    MPI_Comm_rank(VES3D_COMM_WORLD, &rank);
#endif

    /*
     * slicing the geometry data for current process. expected
     * number of fields per line is currently 8: shape_idx scale center_x
     * center_y center_z rot_z, rot_y, rot_z
     */
    int stride(8*nves);
    ASSERT(all_geo_spec.size()>=nproc*stride,"bad initialization file");
    std::vector<value_type> geo_spec(all_geo_spec.begin() + rank     * stride,
                                     all_geo_spec.begin() + (rank+1) * stride);

    //Initial vesicle position container
    INFO("Initializing the starting shapes");
    InitializeShapes(x0, shapes, geo_spec);

    //setting vesicle properties
    ves_props_ = new VProp_t();
    int nprops(VProp_t::n_props);

    if (run_params_.vesicle_props_file.size()) {
        INFO("Loading vesicle properties from file: "<<run_params_.vesicle_props_file);
        Arr_t propsf( nprops * nves * nproc);
        Arr_t props( nprops * nves * nproc);

        io.ReadData( FullPath(run_params_.vesicle_props_file), propsf,
            DataIO::ASCII, 0, propsf.size());

        //order by property (column)
        props.getDevice().Transpose(propsf.begin(),nves*nproc,nprops,props.begin());

        for (int iP(0);iP<nprops;++iP){
            typename VProp_t::container_type* prp(ves_props_->getPropIdx(iP));
            prp->resize(nves);
            prp->getDevice().Memcpy(prp->begin(),
                props.begin() + (iP*nproc + rank)*nves,
                nves * sizeof(typename VProp_t::value_type),
                DT::MemcpyDeviceToDevice);
        }
        ves_props_->update();
    } else { /* populate the properties from commandline */
        CHK(ves_props_->setFromParams(run_params_));
    }

    timestepper_ = new Evolve_t(&run_params_, *Mats_, vInf_, NULL,
        interaction_, NULL, ksp_, &x0, ves_props_);

    return ErrorEvent::Success;
}

template<typename DT, const DT &DEVICE>
Error_t Simulation<DT,DEVICE>::cleanup_run()
{
    INFO("Cleaning up after run");
    delete ves_props_;   ves_props_   = NULL;
    delete Mats_;        Mats_        = NULL;
    delete vInf_;        vInf_	      = NULL;
    delete ksp_;         ksp_	      = NULL;
    delete interaction_; interaction_ = NULL;
    delete timestepper_; timestepper_ = NULL;

    load_checkpoint_ = false;
    checkpoint_data_.str("");
    checkpoint_data_.clear();
    return ErrorEvent::Success;
}

template<typename DT, const DT &DEVICE>
Error_t Simulation<DT,DEVICE>::prepare_run_params(int argc, char **argv, const DictString_t *dict)
{
    Param_t ip;
    CHK(ip.parseInput(argc, argv, dict));
    Error_t rval(prepare_run_params(ip));

    //reparse to override the checkpoint
    if (load_checkpoint_) run_params_.parseInput(argc, argv, dict);

    return rval;
}

template<typename DT, const DT &DEVICE>
Error_t Simulation<DT,DEVICE>::prepare_run_params(const Param_t &ip)
{
    if (ip.load_checkpoint != ""){
        //std::string fname = FullPath(ip.load_checkpoint);
        std::string fname = ip.load_checkpoint;
        INFO("Loading checkpoint file "<<fname);
        DataIO::SlurpFile(fname.c_str(), checkpoint_data_);
        load_checkpoint_ = true;
    } else {
        ip.pack(checkpoint_data_, Streamable::ASCII);
        load_checkpoint_ = false;
    }

    //consume the first part of checkpoint (see Monitor for how
    //it is dumped
    run_params_.unpack(checkpoint_data_, Streamable::ASCII);

    return ErrorEvent::Success;
}
