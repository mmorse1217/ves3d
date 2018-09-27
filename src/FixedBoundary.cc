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
    //Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/hourglass.wrl"); // hourglass
    //Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/hourglass.wrl");
    //Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/vessel_section.wrl"); // small vessel
    //Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/vessel_section.wrl");
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/larger_vessel_section.wrl"); // large vessel
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/larger_vessel_section.wrl");
    //Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cuboid.wrl"); // hourglass
    //Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cuboid.wrl");
    //Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/branch.wrl"); // branch
    //Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/branch.wrl");
    Options::set_value_petsc_opts("-bis3d_spacing", ".1"); //hourglass
    //Options::set_value_petsc_opts("-bis3d_spacing", ".05"); //large cube
    Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "2");
    Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "0");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "0");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".075");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".075");
    std::cout<<"fixed boundary set petsc options\n";
    
    surface = new PatchSurfFaceMap("BD3D_", "bd3d_");
    surface->_surface_type = PatchSurfFaceMap::BLENDED;
    surface->setFromOptions();
    surface->setup();
    surface->_coarse = true;
    surface->refine_test();
    //surface->refine_uniform(2);
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
    tri_vertices_spacing = 0.2; // vessel branch
    //tri_vertices_spacing = 0.1; //hourglass
    //tri_vertices_spacing = 0.1; //branch
    SetTriData();

    // boundary flow
    Petsc::create_mpi_vec(solver->mpiComm(),
            total_num_dof,
            boundary_flow);
    VecSet(boundary_flow, 0.);
    Petsc::create_mpi_vec(solver->mpiComm(),
            total_num_dof,
            boundary_flow_tmp);
    VecSet(boundary_flow_tmp, 0.);

    BuildInOutLets();
    
    SetBoundaryFlow();

    FillVesicle(1);

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
        delete[] tri_vertices;
    if(surface)
        delete surface;
    if(solver)
        delete solver;
    if(inslots_min)
        delete[] inslots_min;
    if(inslots_max)
        delete[] inslots_max;
    if(inslots_count)
        delete[] inslots_count;
    Petsc::destroy_vec(boundary_data);
    Petsc::destroy_vec(solved_density);
    Petsc::destroy_vec(solved_density_tmp);
    Petsc::destroy_vec(boundary_flow);
    Petsc::destroy_vec(boundary_flow_tmp);
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
    computed_local_1.restore_local_vector();

/*
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
    computed_local_2.restore_local_vector();
*/

    // test far eval
    VecSet(computed_potential, 0);
    solver->fareval(targets, VAR_U, solved_density_tmp, computed_potential);
    DblNumMat computed_local_3 = get_local_vector(1, num_target_points*DIM, computed_potential);
    double max_3 = -100;
    for(int i=0; i<num_target_points*DIM;i++)
    {
        max_3 = std::max(max_3, std::abs(computed_local_3(0,i)-1) );
    }
    computed_local_3.restore_local_vector();
    std::cout<<"constant density, evaluate, with closest points, error: "<<max_1<<".\n";
    //std::cout<<"constant density, evaluate, without closest points, error: "<<max_2<<".\n";
    std::cout<<"constant density, fareval, error: "<<max_3<<".\n";
    Petsc::destroy_vec(targets);
    Petsc::destroy_vec(computed_potential);
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
RestoreSamplePoints(double **local_sample_points)
{
    VecRestoreArray(solver->patch_samples()->sample_point_3d_position(), local_sample_points);
}

void FixedBoundary::
SetBoundaryData(double* boundary_data_address)
{
    int sample_dof, pole_dof, total_num_dof;
    solver->localSize(sample_dof,pole_dof,total_num_dof);
    
    std::cout<<"total num dof: "<<total_num_dof<<"\n";
    std::cout<<"sample num dof: "<<sample_dof<<"\n";

    VecSet(boundary_data, 0.);
    VecAXPY(boundary_data, 1, boundary_flow);
    
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
    VecAXPY(boundary_data, 1, boundary_flow);

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

void FixedBoundary::
SetBoundaryFlow()
{
    int sample_dof, pole_dof, total_num_dof;
    solver->localSize(sample_dof,pole_dof,total_num_dof);

    double* sample_points_address;
    VecGetArray(solver->patch_samples()->sample_point_3d_position(), &sample_points_address);
    int num_sample_points = solver->patch_samples()->local_num_sample_points();
    
    DblNumMat boundary_flow_local = get_local_vector(1, total_num_dof, boundary_flow);
    DblNumMat boundary_flow_local_tmp = get_local_vector(1, total_num_dof, boundary_flow_tmp);
    //loop over all inlets
    for(int i=0; i<boundary_inlets.size(); i++)
    {
        double bboxl[3];
        double bboxc[3];
        for(int j=0; j<COORD_DIM; j++)
        {
            bboxl[j] = boundary_inlets[i].bmax[j] - boundary_inlets[i].bmin[j];
            bboxc[j] = (boundary_inlets[i].bmax[j] + boundary_inlets[i].bmin[j])/2;
        }

        double r2_max = (bboxl[0]*bboxl[0] + bboxl[1]*bboxl[1] + bboxl[2]*bboxl[2])/4;
        for(int j=0; j<num_sample_points; j++)
        {
            double sample_dis[3];
            for(int jj=0; jj<COORD_DIM; jj++)
            {
                sample_dis[jj] = sample_points_address[3*j+jj] - bboxc[jj];
            }
            double r2 = sample_dis[0]*sample_dis[0] + sample_dis[1]*sample_dis[1] + sample_dis[2]*sample_dis[2];
        
            double strength = 0;
        
            if(r2>=r2_max)
                strength = 0;
            else
                strength = std::exp(1.0/(r2/r2_max-1))/std::exp(-1.0);
        
            boundary_flow_local(0,3*j+0)+=strength * boundary_inlets[i].flow_dir[0];
            boundary_flow_local(0,3*j+1)+=strength * boundary_inlets[i].flow_dir[1];
            boundary_flow_local(0,3*j+2)+=strength * boundary_inlets[i].flow_dir[2];
        }
    }

    // do integral on boundary_flow_local with only inlets
    double* sample_points_weight;
    double* sample_points_normal;
    VecGetArray(solver->patch_samples()->sample_point_combined_weight(), &sample_points_weight);
    VecGetArray(solver->patch_samples()->sample_point_normal(), &sample_points_normal);
    double flux_inlets = 0, total_flux_inlets=0;
    double udotn = 0;
    for(int i=0; i<num_sample_points; i++)
    {
        udotn = 0;
        for(int j=0; j<COORD_DIM; j++)
            udotn += boundary_flow_local(0,3*i+j)*sample_points_normal[3*i+j];
        udotn *= sample_points_weight[i];
        flux_inlets += udotn;
    }
    MPI_Allreduce(&flux_inlets, &total_flux_inlets, 1, MPI_DOUBLE, MPI_SUM, comm);
    std::cout<<"flux_inlets: "<<flux_inlets<<std::endl;
    
    //loop over all outlets
    for(int i=0; i<boundary_outlets.size(); i++)
    {
        double bboxl[3];
        double bboxc[3];
        for(int j=0; j<COORD_DIM; j++)
        {
            bboxl[j] = boundary_outlets[i].bmax[j] - boundary_outlets[i].bmin[j];
            bboxc[j] = (boundary_outlets[i].bmax[j] + boundary_outlets[i].bmin[j])/2;
        }

        double r2_max = (bboxl[0]*bboxl[0] + bboxl[1]*bboxl[1] + bboxl[2]*bboxl[2])/4;
        for(int j=0; j<num_sample_points; j++)
        {
            double sample_dis[3];
            for(int jj=0; jj<COORD_DIM; jj++)
            {
                sample_dis[jj] = sample_points_address[3*j+jj] - bboxc[jj];
            }
            double r2 = sample_dis[0]*sample_dis[0] + sample_dis[1]*sample_dis[1] + sample_dis[2]*sample_dis[2];
        
            double strength = 0;
        
            if(r2>=r2_max)
                strength = 0;
            else
                strength = std::exp(1.0/(r2/r2_max-1))/std::exp(-1.0);
        
            boundary_flow_local_tmp(0,3*j+0)+=strength * boundary_outlets[i].flow_dir[0];
            boundary_flow_local_tmp(0,3*j+1)+=strength * boundary_outlets[i].flow_dir[1];
            boundary_flow_local_tmp(0,3*j+2)+=strength * boundary_outlets[i].flow_dir[2];
        }
    }
    // do integral on boundary_flow_local with only outlets
    double flux_outlets = 0, total_flux_outlets = 0;
    udotn = 0;
    for(int i=0; i<num_sample_points; i++)
    {
        udotn = 0;
        for(int j=0; j<COORD_DIM; j++)
            udotn += boundary_flow_local_tmp(0,3*i+j)*sample_points_normal[3*i+j];
        udotn *= sample_points_weight[i];
        flux_outlets += udotn;
    }
    MPI_Allreduce(&flux_outlets, &total_flux_outlets, 1, MPI_DOUBLE, MPI_SUM, comm);

    // restore array
    VecRestoreArray(solver->patch_samples()->sample_point_combined_weight(), &sample_points_weight);
    VecRestoreArray(solver->patch_samples()->sample_point_normal(), &sample_points_normal);
    boundary_flow_local.restore_local_vector();
    boundary_flow_local_tmp.restore_local_vector();
    VecRestoreArray(solver->patch_samples()->sample_point_3d_position(), &sample_points_address);

    // adjust strength
    VecAXPY(boundary_flow, -total_flux_inlets/total_flux_outlets, boundary_flow_tmp);

    /*
    for(int i=0; i<num_sample_points; i++)
    {
        if(fabs(sample_points_address[3*i+0]) <= 14)
        {
            boundary_flow_local(0,3*i+0)=0;
            boundary_flow_local(0,3*i+1)=0;
            boundary_flow_local(0,3*i+2)=0;
        }
        else
        {
            double r2 = sample_points_address[3*i+1]*sample_points_address[3*i+1] 
                + sample_points_address[3*i+2]*sample_points_address[3*i+2];
            double r2_max = 2.6*2.6;

            if(r2>=r2_max)
            {
                boundary_flow_local(0,3*i+0)=0;
                boundary_flow_local(0,3*i+1)=0;
                boundary_flow_local(0,3*i+2)=0;
            }
            else
            {
                boundary_flow_local(0,3*i+0)=std::exp(1.0/(r2/r2_max-1))/std::exp(-1.0);
                boundary_flow_local(0,3*i+1)=0;
                boundary_flow_local(0,3*i+2)=0;
            }
        }
        boundary_flow_local(0,3*i+0)=sample_points_address[3*i+2];
        boundary_flow_local(0,3*i+1)=0;
        boundary_flow_local(0,3*i+2)=0;
    }
    */
}

void FixedBoundary::
BuildInOutLets()
{
    // temp inlet
    BINOUT tmp_let;
    tmp_let.bmin[0] = 22 - 3.1; tmp_let.bmin[1] = -18.3 - 2; tmp_let.bmin[2] = 4.4 - 1.5;
    tmp_let.bmax[0] = 22 + 3.1; tmp_let.bmax[1] = -18.3 + 2; tmp_let.bmax[2] = 4.4 + 1.5;
    tmp_let.flow_dir[0] = -1; tmp_let.flow_dir[1] = 1; tmp_let.flow_dir[2] = 0;
    boundary_inlets.push_back(tmp_let);
    
    // temp outlet
    tmp_let.bmin[0] = -12.8 - 1.5; tmp_let.bmin[1] = 20.3 - 1.5; tmp_let.bmin[2] = -12 - 1;
    tmp_let.bmax[0] = -12.8 + 1.5; tmp_let.bmax[1] = 20.3 + 1.5; tmp_let.bmax[2] = -12 + 1;
    tmp_let.flow_dir[0] = 0; tmp_let.flow_dir[1] = 0; tmp_let.flow_dir[2] = -1;
    boundary_outlets.push_back(tmp_let);
    tmp_let.bmin[0] = -24.8 - 1; tmp_let.bmin[1] = 15.8 - 1; tmp_let.bmin[2] = -12.7 - 1;
    tmp_let.bmax[0] = -24.8 + 1; tmp_let.bmax[1] = 15.8 + 1; tmp_let.bmax[2] = -12.7 + 1;
    tmp_let.flow_dir[0] = -1; tmp_let.flow_dir[1] = 1; tmp_let.flow_dir[2] = -1;
    boundary_outlets.push_back(tmp_let);

    // compute inslots/outslots
    total_inslots = 1;
    inslots_min = new double[total_inslots*COORD_DIM];
    inslots_max = new double[total_inslots*COORD_DIM];
    inslots_count = new int[total_inslots];
    inslots_min[0] = 22 - 0.7; inslots_min[1] = -18.3 - 0.7; inslots_min[2] = 4.4 - 0.7;
    inslots_max[0] = 22 + 0.7; inslots_max[1] = -18.3 + 0.7; inslots_max[2] = 4.4 + 0.7;
    inslots_count[0] = 0;
}

void FixedBoundary::
FillVesicle(double cell_size)
{
    double* sample_points_address;
    VecGetArray(solver->patch_samples()->sample_point_3d_position(), &sample_points_address);
    int num_sample_points = solver->patch_samples()->local_num_sample_points();
    std::cout<<"num sample points: "<<num_sample_points<<std::endl;

    // get fixed boundary bounding box
    double bb_min[3], bb_max[3];
    for(int i=0; i<num_sample_points; i++)
    {
        if(i==0)
        { 
            for(int j=0; j<VES3D_DIM; j++)
            {
                bb_min[j] = sample_points_address[i*VES3D_DIM + j];
                bb_max[j] = sample_points_address[i*VES3D_DIM + j]; 
            }
        }
        else
        {
            for(int j=0; j<VES3D_DIM; j++)
            {
                bb_min[j] = std::min(sample_points_address[i*VES3D_DIM + j], bb_min[j]);
                bb_max[j] = std::max(sample_points_address[i*VES3D_DIM + j], bb_max[j]); 
            }
        }
    }

    // get cell center points
    int num_cen_pts_dim[3];
    int num_cen_pts = 1;
    for(int j=0; j<VES3D_DIM; j++)
    {
        num_cen_pts_dim[j] = floor((bb_max[j]-bb_min[j])/cell_size);
        num_cen_pts *= num_cen_pts_dim[j];
    }
    double *cen_pts = new double[VES3D_DIM*num_cen_pts];
    double *val = new double[VES3D_DIM*num_cen_pts];
    bool *is_valid = new bool[num_cen_pts];
    for(int i=0; i<num_cen_pts_dim[0]; i++)
        for(int j=0; j<num_cen_pts_dim[1]; j++)
            for(int k=0; k<num_cen_pts_dim[2]; k++)
            {
                int ind = num_cen_pts_dim[2]*num_cen_pts_dim[1]*i + num_cen_pts_dim[2]*j + k;
                cen_pts[VES3D_DIM*ind + 0] = bb_min[0]+(i+0.5)*cell_size;
                cen_pts[VES3D_DIM*ind + 1] = bb_min[1]+(j+0.5)*cell_size;
                cen_pts[VES3D_DIM*ind + 2] = bb_min[2]+(k+0.5)*cell_size;
            }

    // remove points outside of fixed boundary
    // do constant density evaluation
    int num_box = 0;
    VecSet(solved_density, 1.);
    EvalPotential(num_cen_pts, cen_pts, val);
    for(int i=0; i<num_cen_pts; i++)
    {
        if(std::abs(val[i*VES3D_DIM+0]-1)>0.2 && std::abs(val[i*VES3D_DIM+1]-1)>0.2 && std::abs(val[i*VES3D_DIM+2]-1)>0.2)
            is_valid[i] = false;
        else
        {
            is_valid[i] = true;
            num_box++;
        }
    }
    // remove cells intersecting with patches bounding box
    double *bb_min_patch = new double[VES3D_DIM*num_patches];
    double *bb_max_patch = new double[VES3D_DIM*num_patches];
    for(int i=0; i<num_patches; i++)
    {
        const double *first_point = &tri_vertices[VES3D_DIM*i*num_vertices_per_patch_1d*num_vertices_per_patch_1d];
        for(int j=0; j<VES3D_DIM; j++)
        {
            bb_min_patch[i*VES3D_DIM+j] = first_point[j];
            bb_max_patch[i*VES3D_DIM+j] = first_point[j];
        }
        for(int j=0; j<num_vertices_per_patch_1d; j++)
            for(int k=0; k<num_vertices_per_patch_1d; k++)
            {
                int point_id = num_vertices_per_patch_1d*j + k;
                for(int ii=0; ii<VES3D_DIM; i++)
                {
                    bb_min_patch[i*VES3D_DIM+ii] = std::min(bb_min_patch[i*VES3D_DIM+ii], first_point[point_id*VES3D_DIM+ii]);
                    bb_max_patch[i*VES3D_DIM+ii] = std::max(bb_max_patch[i*VES3D_DIM+ii], first_point[point_id*VES3D_DIM+ii]);
                }
            }
    }

    for(int i=0; i<num_cen_pts_dim[0]; i++)
        for(int j=0; j<num_cen_pts_dim[1]; j++)
            for(int k=0; k<num_cen_pts_dim[2]; k++)
            {
                int ind = num_cen_pts_dim[2]*num_cen_pts_dim[1]*i + num_cen_pts_dim[2]*j + k;
                if(is_valid[ind]==false)
                    continue;

                double bb_min_ves[3], bb_max_ves[3];
                for(int ii=0; ii<VES3D_DIM; ii++)
                {
                    bb_min_ves[ii] = cen_pts[VES3D_DIM*ind+ii] - 0.5*cell_size;
                    bb_max_ves[ii] = cen_pts[VES3D_DIM*ind+ii] + 0.5*cell_size; 
                }
                cen_pts[VES3D_DIM*ind + 2] = bb_min[2]+(k+0.5)*cell_size;
                                
                for(int ii=0; ii<num_patches; ii++)
                {
                    double bb_min_tmp[3], bb_max_tmp[3];
                    for(int jj=0; jj<VES3D_DIM; jj++)
                    {
                        bb_min_tmp[jj] = bb_min_patch[ii*VES3D_DIM+jj];
                        bb_max_tmp[jj] = bb_max_patch[ii*VES3D_DIM+jj];
                    }
                    if(bb_min_ves[0]<=bb_max_tmp[0] && bb_max_ves[0]>=bb_min_tmp[0] && 
                       bb_min_ves[1]<=bb_max_tmp[1] && bb_max_ves[1]>=bb_min_tmp[1] && 
                       bb_min_ves[2]<=bb_max_tmp[2] && bb_max_ves[2]>=bb_min_tmp[2]
                      )
                    {
                        is_valid[ind] = false;
                        break;
                    }
                }
            }

    
    VecRestoreArray(solver->patch_samples()->sample_point_3d_position(), &sample_points_address);
    delete[] cen_pts;
    delete[] val;
    delete[] is_valid;
    delete[] bb_min_patch;
    delete[] bb_max_patch;
}
