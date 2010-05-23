template<typename Container, typename SurfContainer>
Stokes<Container, SurfContainer>::Stokes(SurfContainer *S_in) :
    S_(S_in),
    w_sph_(S_->x_.GetShOrder(), S_->x_.GetNumVecs()),
    sing_quad_weights_(S_->x_.GetShOrder(), 1)
{
    int p = S_->x_.GetShOrder();
    int np = S_->x_.GetFunLength();
    rot_mat_.Resize(1,p + 1, make_pair(2 * p,np));//to match the rot_chunck
    all_rot_mats_.Resize(1,S_->x_.GetShOrder() + 1, make_pair(np,np));
}

template<typename Container, typename SurfContainer>
void Stokes<Container, SurfContainer>::operator()(const 
    Container &interfacial_force, const typename Container::value_type &t, 
    Container &velocity) const
{
    int p = S_->x_.GetShOrder();
    int np = S_->x_.GetFunLength();
    int nv = S_->x_.GetNumVecs();
    //int rot_chunck = 2 * p * np;
    
    typename Container::value_type alpha(1.0), beta(0.0);
    int trg_idx(0);
    
    typename SurfContainer::Sca t1(p, S_->x_.GetNumVecs());
    typename SurfContainer::Sca t2(p, S_->x_.GetNumVecs());

    xyInv(S_->w_, w_sph_, t1);
    int nvX3= S_->x_.GetNumFuns();   

    Container v1(p, S_->x_.GetNumVecs());
    Container v2(p, S_->x_.GetNumVecs());

    for(int ii=0;ii <= p; ++ii)
    {
        for(int jj=0;jj < 2 * p; ++jj)
        {
            
            CircShift(all_rot_mats_.begin() + ii * np * np, jj * np, rot_mat_);

            S_->x_.GetDevice().gemm("N", "N", &np, &nvX3, &np, 
                &alpha, rot_mat_.begin(), &np, S_->x_.begin(), &np, &beta, v1.begin(), &np);
            
            S_->x_.GetDevice().gemm("N", "N", &np, &nvX3, &np,
                &alpha, rot_mat_.begin(), &np, interfacial_force.begin(), &np, &beta, v2.begin(), &np);
            
            S_->x_.GetDevice().gemm("N", "N", &np, &nv, &np,
                &alpha, rot_mat_.begin(), &np, t1.begin(), &np, &beta, t2.begin(), &np);
            
            xy(t2, w_sph_, t2);
            xv(t2, v2, v2);
            
            S_->x_.GetDevice().DirectStokes(v1.begin(), v2.begin(), sing_quad_weights_.begin(), np, 
                nv, S_->x_.begin(), ii*2*p + jj, ii*2*p + jj + 1, velocity.begin());
        }
    }
}
