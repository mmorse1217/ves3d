FixedBoundary::
FixedBoundary()
{
    comm = MPI_COMM_WORLD;
    PetscOptionsInsertFile(comm, NULL, "opt/morse_cases.opt", PETSC_TRUE);
    std::cout<<"fixed boundary\n";
    Options::set_value_petsc_opts("-kt","311");
    Options::set_value_petsc_opts("-dom", "1"); // exterior problem
    Options::set_value_petsc_opts("-bdtype", "2"); // blended surface
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube.wrl"); // single cube
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube.wrl");
    Options::set_value_petsc_opts("-bis3d_spacing", ".1");
    Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "2");
    Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "0");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "1");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".1");
    std::cout<<"fixed boundary set petsc options\n";
    
    surface = new PatchSurfFaceMap("BD3D_", "bd3d_");
    surface->_surface_type = PatchSurfFaceMap::BLENDED;
    surface->setFromOptions();
    surface->setup();
    surface->_coarse = true;
    surface->refine_test();
    std::cout<<"surface refine\n";
    //test_constant_boundary_data(surface.get(),true);
    //
    Options::set_value_petsc_opts("-near_interpolation_num_samples", "6");
    
    solver = new SolverGMRESDoubleLayer(surface);
    solver->set_evaluation_type(EXTRAPOLATION_AVERAGE);
    solver->_compute_refined_surface = true;
    solver->eqcoefs() = vector<double>(2,1.0);
    solver->eqcoefs()[1] = .4;
    solver->setFromOptions();
    solver->setup();
    std::cout<<"solver setup\n";

    // boundary data
    int sample_dof, pole_dof, total_num_dof;
    solver->localSize(sample_dof,pole_dof,total_num_dof);
    Petsc::create_mpi_vec(solver->mpiComm(),
            total_num_dof,
            boundary_data);
    VecSet(boundary_data, 0.);
    
    Petsc::create_mpi_vec(solver->mpiComm(),
            total_num_dof,
            solved_density);
    VecSet(solved_density, 1.);

    DblNumMat boundary_data_local = get_local_vector(1, total_num_dof, solved_density);
    for(int i=0; i<boundary_data_local._n; i++)
    {
        if(i<sample_dof)
            continue;
        boundary_data_local(0,i)=0.0;
    }
    boundary_data_local.restore_local_vector();

    /*
    Vec computed_potential;
    Vec targets;
    Vec boundary_data;
    Vec solved_density;

    // targets
    int num_samples_1d = 35;
    int num_target_points = num_samples_1d*num_samples_1d*num_samples_1d;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, num_target_points*DIM, targets);
    Petsc::create_mpi_vec(MPI_COMM_WORLD, num_target_points*DIM, computed_potential);
    VecSet(computed_potential, 100);
    DblNumMat targets_local = get_local_vector(DIM, num_target_points, targets);
    Cube target_domain = Cube(Interval(4,5),Interval(4,5),Interval(4,5));
    DblNumMat volume_samples(DIM, num_target_points, true,
            Sampling::sample_3d<Sampling::equispaced>(num_samples_1d, target_domain).data());
    for(int i=0; i<num_target_points;i++)
    {
        for(int d=0; d<DIM; d++)
        {
            targets_local(d,i) = volume_samples(d,i);
        }
    }
    targets_local.restore_local_vector();
    
    // boundary data
    int sample_dof, pole_dof, total_num_dof;
    solver->localSize(sample_dof,pole_dof,total_num_dof);
    Petsc::create_mpi_vec(solver->mpiComm(),
            total_num_dof,
            boundary_data);
    VecSet(boundary_data, 1.);
    DblNumMat boundary_data_local = get_local_vector(1, total_num_dof, boundary_data);
    for(int i=0; i<boundary_data_local._n; i++)
    {
        if(i<sample_dof)
            continue;
        boundary_data_local(0,i)=0.0;
    }
    boundary_data_local.restore_local_vector();

    Petsc::create_mpi_vec(solver->mpiComm(),
            total_num_dof,
            solved_density);
    VecSet(solved_density, 0.);

    
    //solver->fareval(targets, VAR_U, boundary_data, computed_potential);
    //solver->solve(boundary_data, solved_density);
    //std::cout<<"solver fareval\n";
    //DblNumMat computed_local = get_local_vector(1, num_target_points*DIM, computed_potential);
    //DblNumMat computed_local = get_local_vector(1, total_num_dof, solved_density);
    //for(int i=0; i<computed_local._n; i++)
    //{
    //    std::cout<<computed_local(0,i)<<std::endl;
    //}
    std::cout<<"total dof: "<<total_num_dof<<". sample_dof: "<<sample_dof<<". pole_dof: "<<pole_dof<<"\n";
    std::cout<<"num_patches: "<<solver->bdry()->num_patches()<<".\n";
    std::cout<<"num_sample_points: "<<solver->patch_samples()->local_num_sample_points()<<".\n";
    vector<int>& num_sample_points_in_patch = solver->patch_samples()->num_sample_points_in_patch();
    for(int i=0; i<num_sample_points_in_patch.size(); i++)
    {
        std::cout<<num_sample_points_in_patch[i]<<"\n";
    }
    
    std::cout<<"num_refined_patches: "<<solver->refined_surface()->num_patches()<<".\n";
    std::cout<<"num_refined_sample_points: "<<solver->refined_patch_samples()->local_num_sample_points()<<".\n";
    vector<int>& num_refined_sample_points_in_patch = solver->refined_patch_samples()->num_sample_points();
    for(int i=0; i<num_refined_sample_points_in_patch.size(); i++)
    {
        std::cout<<num_refined_sample_points_in_patch[i]<<"\n";
    }
    */
}

FixedBoundary::
~FixedBoundary()
{
}

void FixedBoundary::
EvalPotential(int num_target_points, double* target_address, double* target_potential)
{
    Vec targets;
    Vec computed_potential;
    Petsc::create_mpi_vec(MPI_COMM_WORLD, num_target_points*DIM, targets);
    Petsc::create_mpi_vec(MPI_COMM_WORLD, num_target_points*DIM, computed_potential);
    VecSet(computed_potential, 0);

    // set target
    DblNumMat targets_local = get_local_vector(1, num_target_points*DIM, targets);
    for(int i=0; i<num_target_points*DIM;i++)
    {
        targets_local(0,i) = target_address[i];
    }
    targets_local.restore_local_vector();

    solver->fareval(targets, VAR_U, solved_density, computed_potential);

    // set computed_potential to target_potential
    DblNumMat computed_local = get_local_vector(1, num_target_points*DIM, computed_potential);
    for(int i=0; i<num_target_points*DIM;i++)
    {
        target_potential[i] = computed_local(0,i);
    }
    computed_local.restore_local_vector();
}
      
void FixedBoundary::
Solve()
{
    solver->solve(boundary_data, solved_density);
}

double* FixedBoundary::
GetSamplePoints(int& num_sample_points)
{
    double* sample_points_address;
    VecGetArray(solver->patch_samples()->sample_point_3d_position(), &sample_points_address);
    
    num_sample_points = solver->patch_samples()->local_num_sample_points();

    return sample_points_address;
}

void FixedBoundary::
SetBoundaryData(double* boundary_data_address)
{
    int sample_dof, pole_dof, total_num_dof;
    solver->localSize(sample_dof,pole_dof,total_num_dof);

    VecSet(boundary_data, 0.);
    DblNumMat boundary_data_local = get_local_vector(1, total_num_dof, boundary_data);
    for(int i=0; i<sample_dof; i++)
    {
        boundary_data_local(0,i)=boundary_data_address[i];
    }
    boundary_data_local.restore_local_vector();
}
