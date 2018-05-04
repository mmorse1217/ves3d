FixedBoundary::
FixedBoundary()
{
    comm = MPI_COMM_WORLD;
    PetscOptionsInsertFile(comm, NULL, "opt/morse_cases.opt", PETSC_TRUE);
    std::cout<<"fixed boundary\n";
    Options::set_value_petsc_opts("-kt","311");
    Options::set_value_petsc_opts("-dom", "0"); // 1 exterior problem, 0 interior problem
    Options::set_value_petsc_opts("-bdtype", "2"); // blended surface
    //Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube_large.wrl"); // single cube
    //Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube_large.wrl");
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/hourglass.wrl"); // hourglass
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/hourglass.wrl");
    //Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/branch.wrl"); // branch
    //Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/branch.wrl");
    Options::set_value_petsc_opts("-bis3d_spacing", ".05");
    Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "2");
    Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "0");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "1");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".075");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".075");
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
    VecSet(solved_density, 0.);
    Petsc::create_mpi_vec(solver->mpiComm(),
            total_num_dof,
            solved_density_tmp);
    VecSet(solved_density_tmp, 1.);

    std::cout<<"total num dof: "<<total_num_dof<<"\n";
    std::cout<<"sample num dof: "<<sample_dof<<"\n";
    DblNumMat boundary_data_local = get_local_vector(1, total_num_dof, solved_density_tmp);
    for(int i=0; i<boundary_data_local._n; i++)
    {
        if(i<sample_dof)
            continue;
        boundary_data_local(0,i)=0.0;
    }
    boundary_data_local.restore_local_vector();
    
    //tri_vertices_spacing = 0.02; //cube large
    tri_vertices_spacing = 0.1; //hourglass
    //tri_vertices_spacing = 0.1; //branch
    SetTriData();

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
    if(tri_vertices)
        delete tri_vertices;
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

    NumVec<OnSurfacePoint> closest_points(num_target_points);
    solver->evaluate(targets, solved_density, computed_potential, closest_points, VAR_U);
    //solver->fareval(targets, VAR_U, solved_density, computed_potential);

    // set computed_potential to target_potential
    DblNumMat computed_local = get_local_vector(1, num_target_points*DIM, computed_potential);
    for(int i=0; i<num_target_points*DIM;i++)
    {
        target_potential[i] = computed_local(0,i);
    }
    computed_local.restore_local_vector();


    // testing
    // test evaluate with computed closest points
    VecSet(computed_potential, 0);
    solver->evaluate(targets, solved_density_tmp, computed_potential, closest_points, VAR_U, true);
    DblNumMat computed_local_1 = get_local_vector(1, num_target_points*DIM, computed_potential);
    double max_1 = -100;
    for(int i=0; i<num_target_points*DIM;i++)
    {
        //std::cout<<"evaluate with constant density: "<<computed_local_1(0,i)<<"\n";
        max_1 = std::max(max_1, std::abs(computed_local_1(0,i)-1) );
    }

    // test evaluate without computed closest points
    VecSet(computed_potential, 0);
    NumVec<OnSurfacePoint> closest_points_tmp(num_target_points);
    solver->evaluate(targets, solved_density_tmp, computed_potential, closest_points_tmp, VAR_U);
    DblNumMat computed_local_2 = get_local_vector(1, num_target_points*DIM, computed_potential);
    double max_2 = -100;
    for(int i=0; i<num_target_points*DIM;i++)
    {
        max_2 = std::max(max_2, std::abs(computed_local_2(0,i)-1) );
    }

    // test far eval
    VecSet(computed_potential, 0);
    solver->fareval(targets, VAR_U, solved_density_tmp, computed_potential);
    DblNumMat computed_local_3 = get_local_vector(1, num_target_points*DIM, computed_potential);
    double max_3 = -100;
    for(int i=0; i<num_target_points*DIM;i++)
    {
        max_3 = std::max(max_3, std::abs(computed_local_3(0,i)-1) );
    }
    std::cout<<"constant density, evaluate, with closest points, error: "<<max_1<<".\n";
    std::cout<<"constant density, evaluate, without closest points, error: "<<max_2<<".\n";
    std::cout<<"constant density, fareval, error: "<<max_3<<".\n";
    //if(max_1>0.01)
        //abort();
    //
}
      
void FixedBoundary::
Solve()
{
    VecSet(solved_density, 0.);

    static int count_i = 0;
    if(count_i == 0)
    {
        std::stringstream ss;
        ss<<std::setfill('0')<<std::setw(5)<<count_i;
        std::string fname = "solved_density_";
        fname += ss.str();
        fname += ".vtp";
        write_general_points_to_vtk(solver->patch_samples()->sample_point_3d_position(), 3, 
                fname, solved_density, "data/");
        count_i++;
    }
    
    solver->solve(boundary_data, solved_density);
    
    std::stringstream ss;
    ss<<std::setfill('0')<<std::setw(5)<<count_i;
    std::string fname = "solved_density_";
    fname += ss.str();
    fname += ".vtp";
    write_general_points_to_vtk(solver->patch_samples()->sample_point_3d_position(), 3, 
            fname, solved_density, "data/");
    count_i++;

    /*
    int sample_dof, pole_dof, total_num_dof;
    solver->localSize(sample_dof,pole_dof,total_num_dof);
    std::cout<<"total num dof: "<<total_num_dof<<"\n";
    std::cout<<"sample num dof: "<<sample_dof<<"\n";
    DblNumMat boundary_data_local = get_local_vector(1, total_num_dof, solved_density);
    double max_1 = -100;
    for(int i=0; i<total_num_dof;i++)
    {
        std::cout<<"solved density with constant rhs: "<<boundary_data_local(0,i)<<"\n";
        max_1 = std::max(max_1, std::abs(boundary_data_local(0,i)-1) );
    }
    std::cout<<"constant boundary data, solved density error: "<<max_1<<"\n";
    */
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
    
    std::cout<<"total num dof: "<<total_num_dof<<"\n";
    std::cout<<"sample num dof: "<<sample_dof<<"\n";

    VecSet(boundary_data, 0.);
    
    static int count_i = 0;
    if(count_i == 0)
    {
        std::stringstream ss;
        ss<<std::setfill('0')<<std::setw(5)<<count_i;
        std::string fname = "boundary_data_";
        fname += ss.str();
        fname += ".vtp";
        write_general_points_to_vtk(solver->patch_samples()->sample_point_3d_position(), 3, 
                fname, boundary_data, "data/");
        count_i++;
    }
    
    DblNumMat boundary_data_local = get_local_vector(1, total_num_dof, boundary_data);
    for(int i=0; i<sample_dof; i++)
    {
        boundary_data_local(0,i)=boundary_data_address[i];
    }
    boundary_data_local.restore_local_vector();

    //output boundary data
    std::stringstream ss;
    ss<<std::setfill('0')<<std::setw(5)<<count_i;
    std::string fname = "boundary_data_";
    fname += ss.str();
    fname += ".vtp";
    write_general_points_to_vtk(solver->patch_samples()->sample_point_3d_position(), 3, 
            fname, boundary_data, "data/");
    count_i++;
}

void FixedBoundary::
SetTriData()
{
    num_patches = solver->bdry()->num_patches();
    num_vertices_per_patch_1d = solver->bdry()->patch(0)->mesh_num_vertices_1d(tri_vertices_spacing);
    num_vertices_per_patch = solver->bdry()->patch(0)->mesh_num_vertices(tri_vertices_spacing);
    tri_vertices = new double[COORD_DIM*num_patches*num_vertices_per_patch];
    
    for(int pi=0; pi<num_patches; ++pi)
    {
        const auto& patch = solver->bdry()->patch(pi);
        DblNumMat vertices;
        IntNumMat triangles;
        
        patch->mesh_patch(tri_vertices_spacing, Rectangle(Interval(0.,1.), Interval(0.,1.)), vertices, triangles);

        for(int i=0; i<num_vertices_per_patch; ++i)
            for(int d=0; d<COORD_DIM; ++d)
                tri_vertices[COORD_DIM*pi*num_vertices_per_patch + COORD_DIM*i + d] = vertices(d,i);
    }
}
