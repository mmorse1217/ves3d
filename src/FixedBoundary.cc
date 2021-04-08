FixedBoundary::
FixedBoundary()
{
    comm = MPI_COMM_WORLD;
    PetscOptionsInsertFile(comm, NULL, "opt/morse_cases.opt", PETSC_TRUE);
    //PetscOptionsInsertFile(comm, NULL, "weak_setup/weak_65536.opt", PETSC_TRUE);
    //PetscOptionsInsertFile(comm, NULL, "weak_setup/weak_16384.opt", PETSC_TRUE);
    //PetscOptionsInsertFile(comm, NULL, "weak_setup/weak_4096.opt", PETSC_TRUE);
    //PetscOptionsInsertFile(comm, NULL, "weak_setup/weak_1024.opt", PETSC_TRUE);
    //PetscOptionsInsertFile(comm, NULL, "strong_setup/strong_cases.opt", PETSC_TRUE);
    //std::cout<<"fixed boundary\n";
    Options::set_value_petsc_opts("-kt","311");
    Options::set_value_petsc_opts("-dom", "0"); // 1 exterior problem, 0 interior problem
    Options::set_value_petsc_opts("-bdtype", "2"); // blended surface
    //Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cube_large.wrl"); // single cube
    //Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cube_large.wrl");
    Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/hourglass.wrl"); // hourglass
    Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/hourglass.wrl");
    //Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/vessel_section_scaling_large.wrl"); // small vessel
    //Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/vessel_section_scaling_large.wrl");
    //Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cylinder5.wrl"); // small vessel
    //Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cylinder5.wrl");
    //Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/larger_vessel_section2.wrl"); // large vessel
    //Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/larger_vessel_section2.wrl");
    //Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/largest_vessel_section.wrl"); // large vessel
    //Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/largest_vessel_section.wrl");
    //Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/cuboid.wrl"); // hourglass
    //Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/cuboid.wrl");
    //Options::set_value_petsc_opts("-bd3d_filename", "wrl_files/branch.wrl"); // branch
    //Options::set_value_petsc_opts("-bd3d_meshfile", "wrl_files/branch.wrl");
    //Options::set_value_petsc_opts("-bis3d_spacing", ".1"); //hourglass, vessel
    //Options::set_value_petsc_opts("-bis3d_spacing", ".075"); //large vessel
    Options::set_value_petsc_opts("-bis3d_spacing", ".1"); // for test
    //Options::set_value_petsc_opts("-bis3d_spacing", ".05"); //large cube
    Options::set_value_petsc_opts("-bd3d_bdsurf_chttyp", "2");
    Options::set_value_petsc_opts("-bd3d_facemap_adaptive", "0");
    Options::set_value_petsc_opts("-bd3d_facemap_refinement_factor", "3");
    Options::set_value_petsc_opts("-boundary_distance_ratio", ".1");
    //Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".075");
    Options::set_value_petsc_opts("-interpolation_spacing_ratio", ".1");
    //std::cout<<"fixed boundary set petsc options\n";
    
    surface = new PatchSurfFaceMap("BD3D_", "bd3d_");
    surface->_surface_type = PatchSurfFaceMap::BLENDED;
    surface->setFromOptions();
    surface->setup();
    surface->_coarse = true;
    surface->refine_test();
    //surface->refine_uniform(2);
    //std::cout<<"surface refine\n";
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
    //std::cout<<"solver setup\n";

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

    //std::cout<<"total num dof: "<<total_num_dof<<"\n";
    //std::cout<<"sample num dof: "<<sample_dof<<"\n";
    DblNumMat boundary_data_local = get_local_vector(1, total_num_dof, solved_density_tmp);
    for(int i=0; i<boundary_data_local._n; i++)
    {
        if(i<sample_dof)
            continue;
        boundary_data_local(0,i)=0.0;
    }
    boundary_data_local.restore_local_vector();
    
    //tri_vertices_spacing = 0.02; //cube large
    //tri_vertices_spacing = 0.1; // vessel branch
    //tri_vertices_spacing = 0.125; // larger vessel branch
    //tri_vertices_spacing = 0.1; // larger vessel branch
    //tri_vertices_spacing = 0.075; // larger vessel branch ???
    //tri_vertices_spacing = 0.1; //hourglass
    //tri_vertices_spacing = 0.1; //branch
    //tri_vertices_spacing = 0.05; // weak 1024, weak 4096
    //tri_vertices_spacing = 0.1; // weak 16384
    //tri_vertices_spacing = 0.2; // weak 65536
    tri_vertices_spacing = 0.3; // strong
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
    
    //LoadDensity(false);

    /*
    std::string file_name("cell_size_tmp");
    std::ifstream data_file(file_name.c_str(), std::ios::in);
    double cell_size;
    data_file>>cell_size;
    cout<<"cell size: "<<cell_size<<"\n";
    FillVesicle(cell_size); // vessel
    data_file.close();
    exit(0);
    */
    //FillVesicle(0.7); // larger vessel
    //FillVesicle(1.13); // largest vessel
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


    /*
    // testing
    // test evaluate with computed closest points
    VecSet(computed_potential, 0);
    //solver->evaluate(targets, solved_density_tmp, computed_potential, closest_points, VAR_U, true);
    DblNumMat computed_local_1 = get_local_vector(1, num_target_points*DIM, computed_potential);
    double max_1 = -100;
    for(int i=0; i<num_target_points*DIM;i++)
    {
        //std::cout<<"evaluate with constant density: "<<computed_local_1(0,i)<<"\n";
        max_1 = std::max(max_1, std::abs(computed_local_1(0,i)-1) );
    }
    computed_local_1.restore_local_vector();
    */

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

    /*
    // test far eval
    VecSet(computed_potential, 0);
    //solver->fareval(targets, VAR_U, solved_density_tmp, computed_potential);
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
    */
    Petsc::destroy_vec(targets);
    Petsc::destroy_vec(computed_potential);
}
      
void FixedBoundary::
Solve()
{
    //VecSet(solved_density, 0.);
    solver->solve(boundary_data, solved_density);

    /*
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
    
    //solver->solve(boundary_data, solved_density);
    //solver->solve(boundary_flow, solved_density);
    
    std::stringstream ss;
    ss<<std::setfill('0')<<std::setw(5)<<count_i;
    std::string fname = "solved_density_";
    fname += ss.str();
    fname += ".vtp";
    write_general_points_to_vtk(solver->patch_samples()->sample_point_3d_position(), 3, 
            fname, solved_density, "data/");
    count_i++;
    */

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
    
    //std::cout<<"total num dof: "<<total_num_dof<<"\n";
    //std::cout<<"sample num dof: "<<sample_dof<<"\n";

    VecSet(boundary_data, 0.);
    VecAXPY(boundary_data, 1, boundary_flow);
    
    /*
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
    */
    
    DblNumMat boundary_data_local = get_local_vector(1, total_num_dof, boundary_data);
    for(int i=0; i<sample_dof; i++)
    {
        boundary_data_local(0,i)=boundary_data_address[i];
    }
    boundary_data_local.restore_local_vector();
    VecAXPY(boundary_data, 1, boundary_flow);

    /*
    //output boundary data
    std::stringstream ss;
    ss<<std::setfill('0')<<std::setw(5)<<count_i;
    std::string fname = "boundary_data_";
    fname += ss.str();
    fname += ".vtp";
    write_general_points_to_vtk(solver->patch_samples()->sample_point_3d_position(), 3, 
            fname, boundary_data, "data/");
    count_i++;
    */
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
            
        /*
        for(int d=0; d<COORD_DIM; ++d)
        {
            patch->bmin_test[d] = vertices(d, 0);
            patch->bmax_test[d] = vertices(d, 0);
        }
        */

        for(int i=0; i<num_vertices_per_patch; ++i)
            for(int d=0; d<COORD_DIM; ++d)
            {
                tri_vertices[COORD_DIM*pi*num_vertices_per_patch + COORD_DIM*i + d] = vertices(d,i);
                //patch->bmin_test[d] = std::min( vertices(d, i), patch->bmin_test[d] );
                //patch->bmax_test[d] = std::max( vertices(d, i), patch->bmax_test[d] );
            }
    }
}

void FixedBoundary::
SetBoundaryFlow()
{
    int sample_dof, pole_dof, total_num_dof;
    solver->localSize(sample_dof,pole_dof,total_num_dof);

    double* sample_points_weight;
    double* sample_points_normal;
    VecGetArray(solver->patch_samples()->sample_point_combined_weight(), &sample_points_weight);
    VecGetArray(solver->patch_samples()->sample_point_normal(), &sample_points_normal);
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
                strength = 3.0/3.0 * std::exp(1.0/(r2/r2_max-1))/std::exp(-1.0);
        
            boundary_flow_local(0,3*j+0)+=strength * boundary_inlets[i].flow_dir[0];
            boundary_flow_local(0,3*j+1)+=strength * boundary_inlets[i].flow_dir[1];
            boundary_flow_local(0,3*j+2)+=strength * boundary_inlets[i].flow_dir[2];
            //boundary_flow_local(0,3*j+0)-=strength * sample_points_normal[3*j+0];
            //boundary_flow_local(0,3*j+1)-=strength * sample_points_normal[3*j+1];
            //boundary_flow_local(0,3*j+2)-=strength * sample_points_normal[3*j+2];
        }
    }

    // do integral on boundary_flow_local with only inlets
    //double* sample_points_weight;
    //double* sample_points_normal;
    //VecGetArray(solver->patch_samples()->sample_point_combined_weight(), &sample_points_weight);
    //VecGetArray(solver->patch_samples()->sample_point_normal(), &sample_points_normal);
    double flux_inlets = 0, total_flux_inlets=0;
    double udotn = 0;
    double udotx = 0;
    double bd_vol = 0, total_bd_vol = 0;
    for(int i=0; i<num_sample_points; i++)
    {
        udotn = 0;
        udotx = 0;
        for(int j=0; j<COORD_DIM; j++){
            udotn += boundary_flow_local(0,3*i+j)*sample_points_normal[3*i+j];
            udotx += sample_points_normal[3*i+j]*sample_points_address[3*i+j]*1.0/3.0;
        }
        udotn *= sample_points_weight[i];
        udotx *= sample_points_weight[i];
        flux_inlets += udotn;
        bd_vol += udotx;
    }
    MPI_Allreduce(&flux_inlets, &total_flux_inlets, 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(&bd_vol, &total_bd_vol, 1, MPI_DOUBLE, MPI_SUM, comm);
    std::cout<<"bd vol: "<<total_bd_vol<<std::endl;
    //std::cout<<"flux_inlets: "<<flux_inlets<<std::endl;
    
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
                strength = 3.0/3.0 * std::exp(1.0/(r2/r2_max-1))/std::exp(-1.0);
        
            boundary_flow_local_tmp(0,3*j+0)+=strength * boundary_outlets[i].flow_dir[0];
            boundary_flow_local_tmp(0,3*j+1)+=strength * boundary_outlets[i].flow_dir[1];
            boundary_flow_local_tmp(0,3*j+2)+=strength * boundary_outlets[i].flow_dir[2];
            //boundary_flow_local_tmp(0,3*j+0)+=strength * sample_points_normal[3*j+0];
            //boundary_flow_local_tmp(0,3*j+1)+=strength * sample_points_normal[3*j+1];
            //boundary_flow_local_tmp(0,3*j+2)+=strength * sample_points_normal[3*j+2];
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

    /*
    for(int i=0; i<num_sample_points; i++)
    {
        boundary_flow_local(0,3*i+0)=1;
        boundary_flow_local(0,3*i+1)=0;
        boundary_flow_local(0,3*i+2)=0;
    }
    */


    // restore array
    VecRestoreArray(solver->patch_samples()->sample_point_combined_weight(), &sample_points_weight);
    VecRestoreArray(solver->patch_samples()->sample_point_normal(), &sample_points_normal);
    boundary_flow_local.restore_local_vector();
    boundary_flow_local_tmp.restore_local_vector();
    VecRestoreArray(solver->patch_samples()->sample_point_3d_position(), &sample_points_address);

    // adjust strength
    VecAXPY(boundary_flow, -total_flux_inlets/total_flux_outlets, boundary_flow_tmp);

}

void FixedBoundary::
BuildInOutLets()
{
    /*
    // vessel
    // temp inlet
    BINOUT tmp_let;
    tmp_let.bmin[0] = 46 - 1.5; tmp_let.bmin[1] = -29.3 - 1.5; tmp_let.bmin[2] = 8.5 - 1.5;
    tmp_let.bmax[0] = 46 + 1.5; tmp_let.bmax[1] = -29.3 + 1.5; tmp_let.bmax[2] = 8.5 + 1.5;
    tmp_let.flow_dir[0] = 1; tmp_let.flow_dir[1] = 0; tmp_let.flow_dir[2] = 0;
    boundary_outlets.push_back(tmp_let);
    
    // temp outlet
    tmp_let.bmin[0] = -36.5 - 1.5; tmp_let.bmin[1] = 27 - 1.5; tmp_let.bmin[2] = -19.5 - 1.5;
    tmp_let.bmax[0] = -36.5 + 1.5; tmp_let.bmax[1] = 27 + 1.5; tmp_let.bmax[2] = -19.5 + 1.5;
    tmp_let.flow_dir[0] = 1; tmp_let.flow_dir[1] = -1; tmp_let.flow_dir[2] = 1;
    boundary_inlets.push_back(tmp_let);
    tmp_let.bmin[0] = -15.3 - 2; tmp_let.bmin[1] = 34.1 - 2; tmp_let.bmin[2] = -19.3 - 1;
    tmp_let.bmax[0] = -15.3 + 2; tmp_let.bmax[1] = 34.1 + 2; tmp_let.bmax[2] = -19.3 + 1;
    tmp_let.flow_dir[0] = 0; tmp_let.flow_dir[1] = 0; tmp_let.flow_dir[2] = 1;
    boundary_inlets.push_back(tmp_let);
    // end of vessel
    */

    /*
    // vessel
    // temp inlet
    BINOUT tmp_let;
    tmp_let.bmin[0] = 46 - 1.5; tmp_let.bmin[1] = -29.3 - 1.5; tmp_let.bmin[2] = 8.5 - 1.5;
    tmp_let.bmax[0] = 46 + 1.5; tmp_let.bmax[1] = -29.3 + 1.5; tmp_let.bmax[2] = 8.5 + 1.5;
    tmp_let.flow_dir[0] = -1; tmp_let.flow_dir[1] = 0; tmp_let.flow_dir[2] = 0;
    boundary_inlets.push_back(tmp_let);
    
    // temp outlet
    tmp_let.bmin[0] = -36.5 - 1.5; tmp_let.bmin[1] = 27 - 1.5; tmp_let.bmin[2] = -19.5 - 1.5;
    tmp_let.bmax[0] = -36.5 + 1.5; tmp_let.bmax[1] = 27 + 1.5; tmp_let.bmax[2] = -19.5 + 1.5;
    tmp_let.flow_dir[0] = -1; tmp_let.flow_dir[1] = 1; tmp_let.flow_dir[2] = -1;
    boundary_outlets.push_back(tmp_let);
    tmp_let.bmin[0] = -15.3 - 2; tmp_let.bmin[1] = 34.1 - 2; tmp_let.bmin[2] = -19.3 - 1;
    tmp_let.bmax[0] = -15.3 + 2; tmp_let.bmax[1] = 34.1 + 2; tmp_let.bmax[2] = -19.3 + 1;
    tmp_let.flow_dir[0] = 0; tmp_let.flow_dir[1] = 0; tmp_let.flow_dir[2] = -1;
    boundary_outlets.push_back(tmp_let);
    // end of vessel

    // large vessel
    // temp inlet
    // inlet 1
    BINOUT tmp_let;
    tmp_let.bmin[0] = -10.6 - 1.5; tmp_let.bmin[1] = 17.3 - 1.5; tmp_let.bmin[2] = -31 - 1.5;
    tmp_let.bmax[0] = -10.6 + 1.5; tmp_let.bmax[1] = 17.3 + 1.5; tmp_let.bmax[2] = -31 + 1.5;
    tmp_let.flow_dir[0] = 0; tmp_let.flow_dir[1] = -1; tmp_let.flow_dir[2] = 0;
    for(int d_i=0; d_i<3; d_i++)
    {
        tmp_let.bmin[d_i] = 2.0*tmp_let.bmin[d_i];
        tmp_let.bmax[d_i] = 2.0*tmp_let.bmax[d_i];
    }
    boundary_inlets.push_back(tmp_let);

    // temp outlet
    // outlet 1
    tmp_let.bmin[0] = 3.7 - 0.75; tmp_let.bmin[1] = 9.6 - 1.0; tmp_let.bmin[2] = -23.4 - 1.0;
    tmp_let.bmax[0] = 3.7 + 0.75; tmp_let.bmax[1] = 9.6 + 1.0; tmp_let.bmax[2] = -23.4 + 1.0;
    tmp_let.flow_dir[0] = 1; tmp_let.flow_dir[1] = 0; tmp_let.flow_dir[2] = 0;
    for(int d_i=0; d_i<3; d_i++)
    {
        tmp_let.bmin[d_i] = 2.0*tmp_let.bmin[d_i];
        tmp_let.bmax[d_i] = 2.0*tmp_let.bmax[d_i];
    }
    boundary_outlets.push_back(tmp_let);
    // outlet 2
    tmp_let.bmin[0] = 17 - 0.75; tmp_let.bmin[1] = 2.5 - 0.75; tmp_let.bmin[2] = 0 - 0.75;
    tmp_let.bmax[0] = 17 + 0.75; tmp_let.bmax[1] = 2.5 + 0.75; tmp_let.bmax[2] = 0 + 0.75;
    tmp_let.flow_dir[0] = 0; tmp_let.flow_dir[1] = -1; tmp_let.flow_dir[2] = 1;
    for(int d_i=0; d_i<3; d_i++)
    {
        tmp_let.bmin[d_i] = 2.0*tmp_let.bmin[d_i];
        tmp_let.bmax[d_i] = 2.0*tmp_let.bmax[d_i];
    }
    boundary_outlets.push_back(tmp_let);
    // outlet 3
    tmp_let.bmin[0] = 6.1 - 0.5; tmp_let.bmin[1] = 14 - 0.5; tmp_let.bmin[2] = 1 - 0.75;
    tmp_let.bmax[0] = 6.1 + 0.5; tmp_let.bmax[1] = 14 + 0.5; tmp_let.bmax[2] = 1 + 0.75;
    tmp_let.flow_dir[0] = 1; tmp_let.flow_dir[1] = 0; tmp_let.flow_dir[2] = -1;
    for(int d_i=0; d_i<3; d_i++)
    {
        tmp_let.bmin[d_i] = 2.0*tmp_let.bmin[d_i];
        tmp_let.bmax[d_i] = 2.0*tmp_let.bmax[d_i];
    }
    boundary_outlets.push_back(tmp_let);
    // outlet 4
    tmp_let.bmin[0] = 3 - 0.5; tmp_let.bmin[1] = 28 - 0.5; tmp_let.bmin[2] = 2.7 - 0.75;
    tmp_let.bmax[0] = 3 + 0.5; tmp_let.bmax[1] = 28 + 0.5; tmp_let.bmax[2] = 2.7 + 0.75;
    tmp_let.flow_dir[0] = 0; tmp_let.flow_dir[1] = 0; tmp_let.flow_dir[2] = -1;
    for(int d_i=0; d_i<3; d_i++)
    {
        tmp_let.bmin[d_i] = 2.0*tmp_let.bmin[d_i];
        tmp_let.bmax[d_i] = 2.0*tmp_let.bmax[d_i];
    }
    boundary_outlets.push_back(tmp_let);
    // outlet 5
    tmp_let.bmin[0] = -4 - 0.75; tmp_let.bmin[1] = -21.8 - 0.75; tmp_let.bmin[2] = 2.8 - 0.75;
    tmp_let.bmax[0] = -4 + 0.75; tmp_let.bmax[1] = -21.8 + 0.75; tmp_let.bmax[2] = 2.8 + 0.75;
    tmp_let.flow_dir[0] = 0; tmp_let.flow_dir[1] = 0; tmp_let.flow_dir[2] = 1;
    for(int d_i=0; d_i<3; d_i++)
    {
        tmp_let.bmin[d_i] = 2.0*tmp_let.bmin[d_i];
        tmp_let.bmax[d_i] = 2.0*tmp_let.bmax[d_i];
    }
    boundary_outlets.push_back(tmp_let);
    // outlet 6
    tmp_let.bmin[0] = -7.4 - 0.75; tmp_let.bmin[1] = -11.1 - 0.75; tmp_let.bmin[2] = -2.5 - 0.75;
    tmp_let.bmax[0] = -7.4 + 0.75; tmp_let.bmax[1] = -11.1 + 0.75; tmp_let.bmax[2] = -2.5 + 0.75;
    tmp_let.flow_dir[0] = -1; tmp_let.flow_dir[1] = 0; tmp_let.flow_dir[2] = 1;
    for(int d_i=0; d_i<3; d_i++)
    {
        tmp_let.bmin[d_i] = 2.0*tmp_let.bmin[d_i];
        tmp_let.bmax[d_i] = 2.0*tmp_let.bmax[d_i];
    }
    boundary_outlets.push_back(tmp_let);
    // outlet 7
    tmp_let.bmin[0] = -3.3 - 0.75; tmp_let.bmin[1] = -2.8 - 0.75; tmp_let.bmin[2] = 8.7 - 0.75;
    tmp_let.bmax[0] = -3.3 + 0.75; tmp_let.bmax[1] = -2.8 + 0.75; tmp_let.bmax[2] = 8.7 + 0.75;
    tmp_let.flow_dir[0] = -1; tmp_let.flow_dir[1] = 0; tmp_let.flow_dir[2] = -1;
    for(int d_i=0; d_i<3; d_i++)
    {
        tmp_let.bmin[d_i] = 2.0*tmp_let.bmin[d_i];
        tmp_let.bmax[d_i] = 2.0*tmp_let.bmax[d_i];
    }
    boundary_outlets.push_back(tmp_let);
    // outlet 8
    tmp_let.bmin[0] = 8.4 - 0.75; tmp_let.bmin[1] = 1 - 0.75; tmp_let.bmin[2] = 9.2 - 0.75;
    tmp_let.bmax[0] = 8.4 + 0.75; tmp_let.bmax[1] = 1 + 0.75; tmp_let.bmax[2] = 9.2 + 0.75;
    tmp_let.flow_dir[0] = 0; tmp_let.flow_dir[1] = 0; tmp_let.flow_dir[2] = -1;
    for(int d_i=0; d_i<3; d_i++)
    {
        tmp_let.bmin[d_i] = 2.0*tmp_let.bmin[d_i];
        tmp_let.bmax[d_i] = 2.0*tmp_let.bmax[d_i];
    }
    boundary_outlets.push_back(tmp_let);
    // outlet 9
    tmp_let.bmin[0] = 18.5 - 0.75; tmp_let.bmin[1] = -24.7 - 0.75; tmp_let.bmin[2] = 17 - 0.75;
    tmp_let.bmax[0] = 18.5 + 0.75; tmp_let.bmax[1] = -24.7 + 0.75; tmp_let.bmax[2] = 17 + 0.75;
    tmp_let.flow_dir[0] = 0; tmp_let.flow_dir[1] = -1; tmp_let.flow_dir[2] = 0;
    for(int d_i=0; d_i<3; d_i++)
    {
        tmp_let.bmin[d_i] = 2.0*tmp_let.bmin[d_i];
        tmp_let.bmax[d_i] = 2.0*tmp_let.bmax[d_i];
    }
    boundary_outlets.push_back(tmp_let);
    // outlet 10
    tmp_let.bmin[0] = 42.3 - 0.75; tmp_let.bmin[1] = -34.2 - 0.75; tmp_let.bmin[2] = 24.6 - 0.75;
    tmp_let.bmax[0] = 42.3 + 0.75; tmp_let.bmax[1] = -34.2 + 0.75; tmp_let.bmax[2] = 24.6 + 0.75;
    tmp_let.flow_dir[0] = 1; tmp_let.flow_dir[1] = 0; tmp_let.flow_dir[2] = 0;
    for(int d_i=0; d_i<3; d_i++)
    {
        tmp_let.bmin[d_i] = 2.0*tmp_let.bmin[d_i];
        tmp_let.bmax[d_i] = 2.0*tmp_let.bmax[d_i];
    }
    boundary_outlets.push_back(tmp_let);
    // end of large vessel
    */

    // largest vessel
    // temp inlet
    // inlet 1
    BINOUT tmp_let;
    tmp_let.bmin[0] = -50 - 5; tmp_let.bmin[1] = -14 - 5; tmp_let.bmin[2] = -85 - 5;
    tmp_let.bmax[0] = -50 + 5; tmp_let.bmax[1] = -14 + 5; tmp_let.bmax[2] = -85 + 5;
    tmp_let.flow_dir[0] = 0; tmp_let.flow_dir[1] = sqrt(2)/2; tmp_let.flow_dir[2] = sqrt(2)/2;
    boundary_inlets.push_back(tmp_let);
    // inlet 2
    tmp_let.bmin[0] = -40 - 2.5; tmp_let.bmin[1] = 52 - 2.5; tmp_let.bmin[2] = -76 - 3.5;
    tmp_let.bmax[0] = -40 + 2.5; tmp_let.bmax[1] = 52 + 2.5; tmp_let.bmax[2] = -76 + 3.5;
    tmp_let.flow_dir[0] = 0; tmp_let.flow_dir[1] = -1; tmp_let.flow_dir[2] = 0;
    boundary_inlets.push_back(tmp_let);

    // temp outlet
    // outlet 1
    tmp_let.bmin[0] = -49.5 - 1.25; tmp_let.bmin[1] = -64.5 - 1.25; tmp_let.bmin[2] = -87.7 - 1.25;
    tmp_let.bmax[0] = -49.5 + 1.25; tmp_let.bmax[1] = -64.5 + 1.25; tmp_let.bmax[2] = -87.7 + 1.25;
    tmp_let.flow_dir[0] = sqrt(2)/2; tmp_let.flow_dir[1] = -sqrt(2)/2; tmp_let.flow_dir[2] = 0;
    boundary_outlets.push_back(tmp_let);
    // outlet 2
    tmp_let.bmin[0] = -76 - 3; tmp_let.bmin[1] = -47.5 - 2; tmp_let.bmin[2] = -86.5 - 2.5;
    tmp_let.bmax[0] = -76 + 3; tmp_let.bmax[1] = -47.5 + 2; tmp_let.bmax[2] = -86.5 + 2.5;
    tmp_let.flow_dir[0] = 0; tmp_let.flow_dir[1] = -1; tmp_let.flow_dir[2] = 0;
    boundary_outlets.push_back(tmp_let);
    // outlet 3
    tmp_let.bmin[0] = -95 - 3.5; tmp_let.bmin[1] = -14 - 2.5; tmp_let.bmin[2] = -86.6 - 2.5;
    tmp_let.bmax[0] = -95 + 3.5; tmp_let.bmax[1] = -14 + 2.5; tmp_let.bmax[2] = -86.6 + 2.5;
    tmp_let.flow_dir[0] = 0; tmp_let.flow_dir[1] = -1; tmp_let.flow_dir[2] = 0;
    boundary_outlets.push_back(tmp_let);
    // outlet 4
    tmp_let.bmin[0] = 6.8 - 1.5; tmp_let.bmin[1] = 110.7 - 1.5; tmp_let.bmin[2] = -26.3 - 1.5;
    tmp_let.bmax[0] = 6.8 + 1.5; tmp_let.bmax[1] = 110.7 + 1.5; tmp_let.bmax[2] = -26.3 + 1.5;
    tmp_let.flow_dir[0] = 0; tmp_let.flow_dir[1] = sqrt(2)/2; tmp_let.flow_dir[2] = sqrt(2)/2;
    boundary_outlets.push_back(tmp_let);
    // outlet 5
    tmp_let.bmin[0] = 29 - 1.5; tmp_let.bmin[1] = 69 - 1.5; tmp_let.bmin[2] = 13 - 1.5;
    tmp_let.bmax[0] = 29 + 1.5; tmp_let.bmax[1] = 69 + 1.5; tmp_let.bmax[2] = 13 + 1.5;
    tmp_let.flow_dir[0] = -1; tmp_let.flow_dir[1] = 0; tmp_let.flow_dir[2] = 0;
    boundary_outlets.push_back(tmp_let);
    // outlet 6
    tmp_let.bmin[0] = -35.5 - 1.5; tmp_let.bmin[1] = -3.3 - 2; tmp_let.bmin[2] = 122.5 - 1.5;
    tmp_let.bmax[0] = -35.5 + 1.5; tmp_let.bmax[1] = -3.3 + 2; tmp_let.bmax[2] = 122.5 + 1.5;
    tmp_let.flow_dir[0] = 0; tmp_let.flow_dir[1] = 1; tmp_let.flow_dir[2] = 0;
    boundary_outlets.push_back(tmp_let);
    // outlet 7
    tmp_let.bmin[0] = 23.5 - 1.75; tmp_let.bmin[1] = -92 - 1.75; tmp_let.bmin[2] = 79 - 1.75;
    tmp_let.bmax[0] = 23.5 + 1.75; tmp_let.bmax[1] = -92 + 1.75; tmp_let.bmax[2] = 79 + 1.75;
    tmp_let.flow_dir[0] = 0; tmp_let.flow_dir[1] = -sqrt(2)/2; tmp_let.flow_dir[2] = -sqrt(2)/2;
    boundary_outlets.push_back(tmp_let);
    // outlet 8
    tmp_let.bmin[0] = 37.8 - 1.5; tmp_let.bmin[1] = -105.2 - 1.5; tmp_let.bmin[2] = 48 - 2;
    tmp_let.bmax[0] = 37.8 + 1.5; tmp_let.bmax[1] = -105.2 + 1.5; tmp_let.bmax[2] = 48 + 2;
    tmp_let.flow_dir[0] = 0; tmp_let.flow_dir[1] = -sqrt(2)/2; tmp_let.flow_dir[2] = sqrt(2)/2;
    boundary_outlets.push_back(tmp_let);
    // outlet 9
    tmp_let.bmin[0] = 57.1 - 1; tmp_let.bmin[1] = -17.5 - 1.25; tmp_let.bmin[2] = 46.15 - 1;
    tmp_let.bmax[0] = 57.1 + 1; tmp_let.bmax[1] = -17.5 + 1.25; tmp_let.bmax[2] = 46.15 + 1;
    tmp_let.flow_dir[0] = 0; tmp_let.flow_dir[1] = 1; tmp_let.flow_dir[2] = 0;
    boundary_outlets.push_back(tmp_let);
    // end of largest vessel

    // compute inslots/outslots
    //std::string file_name("wrl_files/inboxes_cm_small"); // vessel
    //std::string file_name("wrl_files/inboxes_cm_large2"); // large vessel
    std::string file_name("wrl_files/inlet_cms_largest"); // large vessel
    std::ifstream data_file(file_name.c_str(), std::ios::in);

    vector<double> boxes_cms;
    double d;
    std::string line;

    while(getline(data_file, line))
    {
        if(line.empty()) continue;
        std::istringstream is(line);
        while(is.good())
        {
            while(is.peek()==' ') is.get();
            is >> d;
            boxes_cms.push_back(d);
        }
    }
    data_file.close();

    total_inslots = boxes_cms.size()/3;
    inslots_min = new double[total_inslots*COORD_DIM];
    inslots_max = new double[total_inslots*COORD_DIM];
    inslots_count = new int[total_inslots];
    for(int i=0; i<total_inslots; i++)
    {
        double cx = boxes_cms[3*i+0];
        double cy = boxes_cms[3*i+1];
        double cz = boxes_cms[3*i+2];
        //inslots_min[i*3+0] = cx - 0.65; inslots_min[i*3+1] = cy - 0.35; inslots_min[i*3+2] = cz - 0.35; // vessel
        //inslots_max[i*3+0] = cx + 0.65; inslots_max[i*3+1] = cy + 0.35; inslots_max[i*3+2] = cz + 0.35; // vessel
        //inslots_min[i*3+0] = cx - 0.1885*2; inslots_min[i*3+1] = cy - 0.332*2; inslots_min[i*3+2] = cz - 0.1885*2; // large vessel
        //inslots_max[i*3+0] = cx + 0.1885*2; inslots_max[i*3+1] = cy + 0.332*2; inslots_max[i*3+2] = cz + 0.1885*2; // large vessel
        inslots_min[i*3+0] = cx - 0.405; inslots_min[i*3+1] = cy - 0.725; inslots_min[i*3+2] = cz - 0.405; // largest vessel
        inslots_max[i*3+0] = cx + 0.405; inslots_max[i*3+1] = cy + 0.725; inslots_max[i*3+2] = cz + 0.405; // largest vessel
        inslots_count[i] = 0;
    }
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
    std::cout<<"bb_min_x: "<<bb_min[0]<<", bb_min_y: "<<bb_min[1]<<", bb_min_z: "<<bb_min[2]<<"\n";
    std::cout<<"bb_max_x: "<<bb_max[0]<<", bb_max_y: "<<bb_max[1]<<", bb_max_z: "<<bb_max[2]<<"\n";
    std::cout<<"num_cen_pts_x: "<<num_cen_pts_dim[0]<<", num_cen_pts_y: "<<num_cen_pts_dim[1]<<", num_cen_pts_z: "<<num_cen_pts_dim[2]<<"\n";
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
        //std::cout<<"val: "<<val[i*VES3D_DIM+0]<<", "<<val[i*VES3D_DIM+1]<<", "<<val[i*VES3D_DIM+2]<<"\n";
        if(std::abs(val[i*VES3D_DIM+0])<0.05 && std::abs(val[i*VES3D_DIM+1])<0.05 && std::abs(val[i*VES3D_DIM+2])<0.05)
            is_valid[i] = false;
        else
        {
            is_valid[i] = true;
            num_box++;
        }
    }
    // remove cells intersecting with patches bounding box
    std::cout<<"generating bounding boxes for all patches\n";
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
                for(int ii=0; ii<VES3D_DIM; ii++)
                {
                    bb_min_patch[i*VES3D_DIM+ii] = std::min(bb_min_patch[i*VES3D_DIM+ii], first_point[point_id*VES3D_DIM+ii]);
                    bb_max_patch[i*VES3D_DIM+ii] = std::max(bb_max_patch[i*VES3D_DIM+ii], first_point[point_id*VES3D_DIM+ii]);
                }
            }
    }

    std::cout<<"testing intersecting bbox\n";
    for(int i=0; i<num_cen_pts_dim[0]; i++)
        for(int j=0; j<num_cen_pts_dim[1]; j++)
            for(int k=0; k<num_cen_pts_dim[2]; k++)
            {
                int ind = num_cen_pts_dim[2]*num_cen_pts_dim[1]*i + num_cen_pts_dim[2]*j + k;
                if(is_valid[ind]==false)
                //if(true)
                    continue;

                double bb_min_ves[3], bb_max_ves[3];
                for(int ii=0; ii<VES3D_DIM; ii++)
                {
                    bb_min_ves[ii] = cen_pts[VES3D_DIM*ind+ii] - 0.5*cell_size;
                    bb_max_ves[ii] = cen_pts[VES3D_DIM*ind+ii] + 0.5*cell_size; 
                }
                if(cen_pts[VES3D_DIM*ind+0] > 40) // vessel
                //if(cen_pts[VES3D_DIM*ind+2] < -28.5) // larger vessel
                //if(cen_pts[VES3D_DIM*ind+2] < -28.5) // largest vessel
                {
                    is_valid[ind]=false;
                    continue;
                }
                                
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
                        // test distance cen point to boundary, is to close to boundary mark false
                        const double *first_point = &tri_vertices[VES3D_DIM*ii*num_vertices_per_patch_1d*num_vertices_per_patch_1d];
                        for(int i_pt=0; i_pt<num_vertices_per_patch_1d; i_pt++)
                        {
                            if(is_valid[ind]==false)
                                break;
                            for(int j_pt=0; j_pt<num_vertices_per_patch_1d; j_pt++)
                            {
                                int point_id = num_vertices_per_patch_1d*i_pt + j_pt;
                                double pt_tmp[3];
                                double dis_tmp = 0;
                                for(int i_dim=0; i_dim<VES3D_DIM; i_dim++)
                                {
                                    pt_tmp[i_dim] = cen_pts[VES3D_DIM*ind+i_dim]; 
                                    pt_tmp[i_dim] -= first_point[point_id*VES3D_DIM+i_dim];
                                    dis_tmp += pt_tmp[i_dim]*pt_tmp[i_dim];
                                }
                                dis_tmp = std::sqrt(dis_tmp);
                                if(dis_tmp<0.5*cell_size)
                                {
                                    is_valid[ind] = false;
                                    break;
                                }
                            }
                        }
                        //is_valid[ind] = false;
                        if(is_valid[ind]==false)
                            break;
                    }
                }
            }

    ofstream myfile;
    myfile.open ("example.txt");
    for(int i=0; i<num_cen_pts_dim[0]; i++)
        for(int j=0; j<num_cen_pts_dim[1]; j++)
            for(int k=0; k<num_cen_pts_dim[2]; k++)
            {
                int ind = num_cen_pts_dim[2]*num_cen_pts_dim[1]*i + num_cen_pts_dim[2]*j + k;
                if(is_valid[ind]==false)
                    continue;
                //myfile<<"0 0.17 "<<cen_pts[VES3D_DIM*ind+0]<<" "<<cen_pts[VES3D_DIM*ind+1]<<" "<<cen_pts[VES3D_DIM*ind+2]<<
                //myfile<<"0 0.10 "<<cen_pts[VES3D_DIM*ind+0]<<" "<<cen_pts[VES3D_DIM*ind+1]<<" "<<cen_pts[VES3D_DIM*ind+2]<<
                myfile<<"0 "<<cell_size*0.2<<" "<<cen_pts[VES3D_DIM*ind+0]<<" "<<cen_pts[VES3D_DIM*ind+1]<<" "<<cen_pts[VES3D_DIM*ind+2]<<
                    " "<<rand()%314/100.0<<" "<<rand()%314/100.0<<" "<<rand()%314/100.0<<"\n";
            }
    myfile.close();

    std::cout<<"clean up\n";
    
    VecRestoreArray(solver->patch_samples()->sample_point_3d_position(), &sample_points_address);
    delete[] cen_pts;
    delete[] val;
    delete[] is_valid;
    delete[] bb_min_patch;
    delete[] bb_max_patch;
}

void FixedBoundary::
LoadDensity(bool flag){
    if(flag){
        int my_rank;
        MPI_Comm_rank(comm,&my_rank);

        // TODO: hard coded
        // file name, like "boundary_density_time16_rank_21.chk"
        string filename = "boundary_density_time36_rank_"+to_string(my_rank)+".chk";
        std::ifstream data_file(filename.c_str(), std::ios::in);

        double d;
        std::string line;

        int sample_dof, pole_dof, total_num_dof;
        solver->localSize(sample_dof,pole_dof,total_num_dof);

        DblNumMat boundary_data_local = get_local_vector(1, total_num_dof, solved_density);
        int i = 0;
        while(getline(data_file, line))
        {
            if(line.empty()) continue;
            std::istringstream is(line);
            while(is.good())
            {
                while(is.peek()==' ') is.get();
                is >> d;
                boundary_data_local(0,i) = d; i++;
            }
        }
        boundary_data_local.restore_local_vector();
        data_file.close();
    }
}

void FixedBoundary::
SaveDensity(){
    static int time_count = 0;
    int my_rank;
    MPI_Comm_rank(comm,&my_rank);

    // file name
    string filename = "data/boundary_density_time"+to_string(time_count)+"_rank_"+to_string(my_rank)+".chk";
    
    std::stringstream ss;
    ss<<std::scientific<<std::setprecision(16);

    int sample_dof, pole_dof, total_num_dof;
    solver->localSize(sample_dof,pole_dof,total_num_dof);
    
    DblNumMat boundary_data_local = get_local_vector(1, total_num_dof, solved_density);
    for(int i=0; i<boundary_data_local._n; i++)
    {
        ss<<boundary_data_local(0,i)<<"\n";
    }
    boundary_data_local.restore_local_vector();

    // open file
    std::ofstream fh(filename, std::ios::out);
    if(!fh)
	CERR_LOC("Cannot open file for writing: "<<filename, "", exit(1));

    //store data
    fh<<ss.rdbuf();
    fh.close();
    time_count++;
}
