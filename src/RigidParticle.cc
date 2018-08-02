template<typename SurfContainer>
RigidParticle<SurfContainer>::
RigidParticle(SurfContainer& surface, int sh_order, int sh_order_up, int sh_order_up_self):
    sh_order_(sh_order),
    sh_order_up_(sh_order_up),
    sh_order_up_self_(sh_order_up_self),
    surface_(surface),
    stokes_(sh_order,sh_order_up,-1,1e-3)
{
    //set up all containers
    //set up stokes_
    //set up weights, dl matrix
    PVFMMVec scoord(surface_.getPosition().size(), surface_.getPosition().begin(), false);
    PVFMMVec scoord_up, scoord_shc, X_theta, X_phi, scoord_area;
    SphericalHarmonics<value_type>::Grid2SHC(scoord,     sh_order, sh_order, scoord_shc);
    SphericalHarmonics<value_type>::SHC2Grid(scoord_shc, sh_order, sh_order, scoord_up, &X_theta, &X_phi);
    COUT(0);
    { // Set scoord_norm, scoord_area
        long Mves=2*sh_order*(sh_order+1);
        long N=X_theta.Dim()/Mves/COORD_DIM;
        scoord_area.ReInit(N*Mves);
        #pragma omp parallel for
        for(long i=0;i<N;i++){
          for(long j=0;j<Mves;j++){
            value_type nx, ny, nz;
            { // Compute source normal
              value_type x_theta=X_theta[(i*COORD_DIM+0)*Mves+j];
              value_type y_theta=X_theta[(i*COORD_DIM+1)*Mves+j];
              value_type z_theta=X_theta[(i*COORD_DIM+2)*Mves+j];

              value_type x_phi=X_phi[(i*COORD_DIM+0)*Mves+j];
              value_type y_phi=X_phi[(i*COORD_DIM+1)*Mves+j];
              value_type z_phi=X_phi[(i*COORD_DIM+2)*Mves+j];

              nx=(y_theta*z_phi-z_theta*y_phi);
              ny=(z_theta*x_phi-x_theta*z_phi);
              nz=(x_theta*y_phi-y_theta*x_phi);
            }
            value_type area=sqrt(nx*nx+ny*ny+nz*nz);
            scoord_area[i*Mves+j]=area;
          }
        }
    }

    COUT(1);
    { //set weights
        PVFMMVec& qw=SphericalHarmonics<value_type>::LegendreWeights(sh_order);
        COUT(2);
        long Mves = 2*sh_order*(sh_order+1);
        long Nves = scoord_area.Dim()/Mves;
        weights_.replicate(surface_.getPosition());
        #pragma omp parallel for
        for(long i=0;i<Nves;i++){
          for(long j0=0;j0<sh_order_up+1;j0++){
            for(long j1=0;j1<sh_order_up*2;j1++){
              long j=j0*sh_order_up*2+j1;
              weights_.begin()[i*Mves+j]=scoord_area[i*Mves+j]*qw[j0];
            }
          }
        }
    }

    // set cms
    COUT("cms");
    tmp1_.replicate(surface_.getPosition());
    tmp2_.replicate(surface_.getPosition());
    tmp_scalar_.replicate(surface_.getPosition());
    surface_.getCenters(tmp1_);
    int num_particles = surface_.getPosition().getNumSubs();
    cm_.resize(num_particles*COORD_DIM);
    COUT("num of particles: "<<num_particles);
    for (int i = 0; i < num_particles; i++) {
        for (int d = 0; d < COORD_DIM; d++) 
            cm_.begin()[i*num_particles+d] = tmp1_.begin()[i*num_particles+d];
        std::cout << std::endl;
    }     


    
    // memory for linear system solver
    COUT("linsolve");
    x_host_ =   new value_type[surface_.getPosition().size()+cm_.size()*2];
    rhs_host_ = new value_type[surface_.getPosition().size()+cm_.size()*2];

    // set MKL GMRES solver
    COUT("linsolve");
    linear_solver_gmres_.SetContext(static_cast<const void*>(this));

    // set up position_, far_vel_, tmp1_, tmp2_, density_, t_vel_, r_vel_
    int num_samples = surface_.getPosition().size()/COORD_DIM;
    position_.replicate(surface_.getPosition());
    position_.getDevice().Memcpy(position_.begin(),
            surface_.getPosition().begin(), 
            num_samples*COORD_DIM*sizeof(value_type), 
            device_type::MemcpyDeviceToDevice);

    density_.replicate(surface_.getPosition());
    density_.getDevice().Memset(density_.begin(), 0, sizeof(value_type)*density_.size());
    
    far_vel_.replicate(surface_.getPosition());
    far_vel_.getDevice().Memset(far_vel_.begin(), 0, sizeof(value_type)*far_vel_.size());
    
    t_vel_.resize(num_particles*COORD_DIM);
    t_vel_.getDevice().Memset(t_vel_.begin(), 0, sizeof(value_type)*t_vel_.size());
    
    r_vel_.resize(num_particles*COORD_DIM);
    r_vel_.getDevice().Memset(r_vel_.begin(), 0, sizeof(value_type)*r_vel_.size());
}

template<typename SurfContainer>
RigidParticle<SurfContainer>::
~RigidParticle()
{
    if(x_host_)
        delete x_host_;
    if(rhs_host_)
        delete rhs_host_;
}

template<typename SurfContainer>
void RigidParticle<SurfContainer>::
EvalPotential(int num_target_points, value_type* target_address, value_type* target_potential)
{
    // set source position
    stokes_.SetSrcCoord(position_);
    // set target position
    PVFMMVec t_coord(COORD_DIM*num_target_points, target_address, false); 
    stokes_.SetTrgCoord(&t_coord);

    // set single layer density 
    stokes_.SetDensitySL(NULL);
    // set double layer density
    stokes_.SetDensityDL(&density_);

    // evaluate potential
    PVFMMVec value(COORD_DIM*num_target_points, target_potential, false);
    value = stokes_();
}
      
template<typename SurfContainer>
void RigidParticle<SurfContainer>::
Solve()
{
    size_t vsz = surface_.getPosition().size();
    size_t fsz = cm_.size();
    size_t tsz = cm_.size();

    // for locally implicit scheme call fmm for inter-particle interaction to rhs_host_
    // then add boundary data to rhs_host_
    stokes_.SetSrcCoord(position_);
    stokes_.SetTrgCoord(NULL);

    // setup far evaluation of vesicles
    stokes_.SetDensitySL(NULL);
    stokes_.SetDensityDL(&density_);

    // do far evaluation
    stokes_.FarInteraction(far_vel_);
    
    // bg flow in x direction
    for (int i = 0; i < far_vel_.size()/3; i++) {
        far_vel_.begin()[i] += 1.;
    }

    // set to zero
    std::memset(rhs_host_, 0, sizeof(value_type)*(vsz + fsz + tsz));
    // move far velocity to RHS
    std::memcpy(rhs_host_, far_vel_.begin(),
            far_vel_.size()*sizeof(value_type));
    /*
    COUT("Right hand side");
    for (int i = 0; i < (vsz + fsz + tsz); i++) {
        COUT(rhs_host_[i]);
    }
    */

    // set x_host_ to solution of last time step as initial guess
    density_.getDevice().Memcpy(x_host_,         density_.begin(), vsz*sizeof(value_type), device_type::MemcpyDeviceToHost);
    t_vel_.getDevice().Memcpy(  x_host_+vsz,     t_vel_.begin(),   fsz*sizeof(value_type), device_type::MemcpyDeviceToHost);
    r_vel_.getDevice().Memcpy(  x_host_+vsz+fsz, r_vel_.begin(),   tsz*sizeof(value_type), device_type::MemcpyDeviceToHost); 
    // TODO: parameters from setup file?
    value_type relres = 1e-7;
    int iter = 200;
    int rsrt = 200;

    int N_size = vsz+fsz+tsz;
    COUT("N_size: "<<N_size);
    COUT("vsz: "<<vsz);
    COUT("fsz: "<<fsz);
    COUT("tsz: "<<tsz);
    // solve
    int solver_ret = linear_solver_gmres_(JacobiImplicitApply, JacobiImplicitPrecond, 
            x_host_, rhs_host_, 
            relres, relres*0, N_size, iter, rsrt);
    
    // put solution in density, translation vel, rotation vel
    density_.getDevice().Memcpy(density_.begin(), x_host_,         vsz*sizeof(value_type), device_type::MemcpyHostToDevice);
    t_vel_.getDevice().Memcpy(t_vel_.begin(),     x_host_+vsz,     fsz*sizeof(value_type), device_type::MemcpyHostToDevice);
    r_vel_.getDevice().Memcpy(r_vel_.begin(),     x_host_+vsz+fsz, tsz*sizeof(value_type), device_type::MemcpyHostToDevice); 
    
    /*
    COUT("Right hand side");
    for (int i = 0; i < (vsz + fsz + tsz); i++) {
        COUT(rhs_host_[i]);
    }
    */
}

template<typename SurfContainer>
typename RigidParticle<SurfContainer>::value_type* RigidParticle<SurfContainer>::
GetSamplePoints(int& num_sample_points)
{
    value_type* sample_points_address;
    sample_points_address = surface_.getPosition().begin();
    
    num_sample_points = surface_.getPosition().size()/3;

    return sample_points_address;
    return NULL;
}

template<typename SurfContainer>
void RigidParticle<SurfContainer>::
SetBoundaryData(value_type* boundary_data_address)
{
}

template<typename SurfContainer>
Error_t RigidParticle<SurfContainer>::
operator()(const Vec_t& density, const Arr_t& t_vel, const Arr_t& r_vel,
           Vec_t &potential, Arr_t &force_int, Arr_t &torque_int) const
{

    // computes double layer integral + stokeslet + rotlet, 
    // integral for translational velocity and 
    // integral for angular velocity 
    /*
    static PVFMMVec density_coef, DL_vel;
    PVFMMVec density_dl(density.size(), (value_type*)density.begin(), false);
    PVFMMVec vel_dl(potential.size(), (value_type*)potential.begin(), false);

    long Ncoef = sh_order_*(sh_order_+2);
    long Ngrid = 2*sh_order_*(sh_order_+1);
    long nv = density_dl.Dim()/Ngrid/COORD_DIM;

    DL_vel.ReInit(0);
    DL_vel.ReInit(nv*COORD_DIM*Ncoef);

    SphericalHarmonics<value_type>::Grid2SHC(density_dl, sh_order_, sh_order_, density_coef);
    
    #pragma omp parallel
    { // mat-vec
        long tid=omp_get_thread_num();
        long omp_p=omp_get_num_threads();

        long a=(tid+0)*nv/omp_p;
        long b=(tid+1)*nv/omp_p;
        for(long i=a;i<b;i++){
            pvfmm::Matrix<value_type> Mv(1,COORD_DIM*Ncoef,&DL_vel      [i*COORD_DIM*Ncoef],false);
            pvfmm::Matrix<value_type> Mf(1,COORD_DIM*Ncoef,&density_coef[i*COORD_DIM*Ncoef],false);
            pvfmm::Matrix<value_type> M(COORD_DIM*Ncoef,COORD_DIM*Ncoef,&DLMatrix_[i*COORD_DIM*Ncoef*COORD_DIM*Ncoef],false);
            pvfmm::Matrix<value_type>::GEMM(Mv,Mf,M);
        }
    }

    SphericalHarmonics<value_type>::SHC2Grid(DL_vel, sh_order_, sh_order_, vel_dl);
    */

    // singular integration for potential
    // TODO: save double layer matrix to DLMatrix_
    stokes_.SetDensitySL(NULL);
    stokes_.SetDensityDL(&density);
    stokes_.SelfInteraction(potential);
    axpy(0.5, density, potential, potential); // y = a*x+y
    
    // negate potential, 
    axpy(-1., potential, potential); //y = a*x
    
    // add stokeslet and rotlet evaluation

    // non-singular integral of density to compute translational velocity
    density.getDevice().Reduce(density.begin(), density.getTheDim(), weights_.begin(), 
                               density.getStride(), density.getNumSubs(), force_int.begin());
    // need to divide by surface area

    // non-singular integral of density cross (x-center)to compute angular velocity
    // x - cm_ -> tmp2_
    Arr_t cm_tmp;
    cm_tmp.resize(cm_.size());
    auto garbage = static_cast<value_type*>(NULL);
    auto minus_one = static_cast<value_type>(-1.0);
    Arr_t::getDevice().axpy(minus_one, cm_.begin(), garbage, cm_.size(), cm_tmp.begin());
    Arr_t::getDevice().apx(cm_tmp.begin(),
            position_.begin(),
            position_.getStride(),
            position_.getNumSubFuncs(),
            tmp2_.begin());

    // add translational velocity component
    Arr_t::getDevice().apx(t_vel.begin(),
            potential.begin(),
            potential.getStride(),
            potential.getNumSubFuncs(),
            potential.begin());
    
    axpy(0.0, potential, tmp1_);
    Arr_t::getDevice().apx(r_vel.begin(),
            tmp1_.begin(),
            tmp1_.getStride(),
            tmp1_.getNumSubFuncs(),
            tmp1_.begin());

    GeometricCross(tmp1_, tmp2_, tmp1_);
    // add angular velocity component
    axpy(1.0, tmp1_, potential, potential);


    GeometricCross(tmp2_, density, tmp1_);
    tmp1_.getDevice().Reduce(tmp1_.begin(), tmp1_.getTheDim(), 
            weights_.begin(), tmp1_.getStride(), tmp1_.getNumSubs(), 
            torque_int.begin());
    // need to divide by surface area


    /*
    COUT("potential");
    for (int i = 0; i < potential.size(); i++) {
        COUT(potential.begin()[i]);
        
    }
    COUT("force");
    for (int i = 0; i < force_int.size(); i++) {
        COUT(force_int.begin()[i]);
        
    }
    COUT("torque");
    for (int i = 0; i < torque_int.size(); i++) {
        COUT(torque_int.begin()[i]);
        
    }
    */
    
    return ErrorEvent::Success;
}

template<typename SurfContainer>
Error_t RigidParticle<SurfContainer>::
JacobiImplicitApply(const GMRESLinSolver<value_type> *o, const value_type *x, value_type *y, int tmp_trash)
{
    const RigidParticle *F(NULL);
    o->Context((const void**) &F);
    size_t vsz = F->surface_.getPosition().size();
    size_t fsz = F->cm_.size();
    size_t tsz = F->cm_.size();

    // TODO: create a poll for memory space or use static so that not reallocate every time called
    Vec_t pot_x, pot_y;
    pot_x.replicate(F->surface_.getPosition());
    pot_y.replicate(F->surface_.getPosition());

    Arr_t f_x(fsz);
    Arr_t f_y(fsz);

    Arr_t t_x(tsz);
    Arr_t t_y(tsz);

    pot_x.getDevice().Memcpy(pot_x.begin(), x,         vsz*sizeof(value_type), device_type::MemcpyHostToDevice);
    f_x.getDevice().Memcpy(f_x.begin(),     x+vsz,     fsz*sizeof(value_type), device_type::MemcpyHostToDevice);
    t_x.getDevice().Memcpy(t_x.begin(),     x+vsz+fsz, tsz*sizeof(value_type), device_type::MemcpyHostToDevice); 

    F->operator()(pot_x, f_x, t_x, pot_y, f_y, t_y);

    pot_y.getDevice().Memcpy(y,         pot_y.begin(), vsz*sizeof(value_type), device_type::MemcpyDeviceToHost);
    f_y.getDevice().Memcpy(  y+vsz,     f_y.begin(),   fsz*sizeof(value_type), device_type::MemcpyDeviceToHost);
    t_y.getDevice().Memcpy(  y+vsz+fsz, t_y.begin(),   tsz*sizeof(value_type), device_type::MemcpyDeviceToHost); 

    /*
    COUT("y:");
    for (int i = 0; i < vsz+fsz+tsz; i++) {
        COUT("y: "<<i<<" "<<y[i]);
    }
    */

    return ErrorEvent::Success;
}

template<typename SurfContainer>
Error_t RigidParticle<SurfContainer>::
JacobiImplicitPrecond(const GMRESLinSolver<value_type> *o, const value_type *x, value_type *y)
{
    const RigidParticle *F(NULL);
    o->Context((const void**) &F);
    
    size_t vsz = F->surface_.getPosition().size();
    size_t fsz = F->cm_.size();
    size_t tsz = F->cm_.size();

    Arr_t::getDevice().Memcpy(y, x, (vsz+fsz+tsz)*sizeof(value_type), device_type::MemcpyHostToHost);
    
    return ErrorEvent::Success;
}
// add stokelet and rotlet evaluation function
