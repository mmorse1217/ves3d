template<typename SurfContainer>
RigidParticle<SurfContainer>::
RigidParticle(int sh_order, int sh_order_up, int sh_order_up_self) : 
    sh_order_(sh_order),
    sh_order_up_(sh_order_up),
    sh_order_up_self_(sh_order_up_self),
    stokes_(sh_order,sh_order_up,-1,1e-3)
{
    //set up all containers
    //set up stokes_
    //set up weights, dl matrix
}

template<typename SurfContainer>
RigidParticle<SurfContainer>::
~RigidParticle()
{
}

template<typename SurfContainer>
void RigidParticle<SurfContainer>::
EvalPotential(int num_target_points, value_type* target_address, value_type* target_potential)
{
    
}
      
template<typename SurfContainer>
void RigidParticle<SurfContainer>::
Solve()
{
    stokes_.SetSrcCoord(position_);
    stokes_.SetTrgCoord(NULL);
}

template<typename SurfContainer>
typename RigidParticle<SurfContainer>::value_type* RigidParticle<SurfContainer>::
GetSamplePoints(int& num_sample_points)
{
    //double* sample_points_address;
    //VecGetArray(solver->patch_samples()->sample_point_3d_position(), &sample_points_address);
    
    //num_sample_points = solver->patch_samples()->local_num_sample_points();

    //return sample_points_address;
    return NULL;
}

template<typename SurfContainer>
void RigidParticle<SurfContainer>::
SetBoundaryData(value_type* boundary_data_address)
{
}

template<typename SurfContainer>
Error_t RigidParticle<SurfContainer>::
operator()(const Vec_t density, const Arr_t t_vel, const Arr_t r_vel,
           Vec_t &potential, Arr_t &force_int, Arr_t &torque_int) const
{
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
    stokes_.SetDensitySL(NULL);
    stokes_.SetDensityDL(density);
    stokes_.SelfInteraction(potential);
    
    // add stokeslet and rotlet evaluation

    // integral of density
    //xv(weights_, density_, tmp1_);
    density.getDevice().Reduce(density.begin(), density.getTheDim(), weights_.begin(), 
                               density.getStride(), density.getNumSubs(), force_int.begin());

    // integral of density cross (x-center)
    // TODO: should be -cm_
    position_.getDevice().apx(cm_.begin(), position_.begin(), position_.getStride(), position_.getNumSubFuncs(), tmp2_.begin());
    GeometricCross(tmp2_, density, tmp1_);
    //xv(weights_, tmp2_, tmp1_);
    tmp1_.getDevice().Reduce(tmp1_.begin(), tmp1_.getTheDim(), weights_.begin(), 
                               tmp1_.getStride(), tmp1_.getNumSubs(), torque_int.begin());
    
}

template<typename SurfContainer>
Error_t RigidParticle<SurfContainer>::
JacobiImplicitApply(const GMRESLinSolver<value_type> *o, const value_type *x, value_type *y, int tmp_trash)
{
}

// add stokelet and rotlet evaluation function
