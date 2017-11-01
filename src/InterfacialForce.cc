template<typename SurfContainer>
InterfacialForce<SurfContainer>::InterfacialForce(
    const Parameters<value_type> &params,
    const VProp_t &ves_props,
    const OperatorsMats<Arr_t> &mats) :
    params_(params),
    ves_props_(ves_props),
    sht_   (mats.p_   , mats.mats_p_   ),
    sht_up_(mats.p_up_, mats.mats_p_up_),
    cen(1, 1, std::make_pair(1,1)),
    S_up(NULL)
{
    int omp_p = omp_get_max_threads();
    for(int i=0; i<omp_p; ++i)
        S_ups.push_back(NULL);
}

template<typename SurfContainer>
InterfacialForce<SurfContainer>::~InterfacialForce()
{
  if(S_up) delete S_up;
  int omp_p = omp_get_max_threads();
  for(int i=0; i<omp_p; ++i)
  {
      if(S_ups[i]) 
          delete S_ups[i];
  }
}

template<typename SurfContainer>
void InterfacialForce<SurfContainer>::bendingForce(const SurfContainer &S,
    Vec_t &Fb) const
{
    S.resample(params_.upsample_freq, &S_up); // upsample

    Vec_t Fb_up;
    Fb_up.replicate(S_up->getPosition());
    s1.replicate(S_up->getPosition());
    s2.replicate(S_up->getPosition());

    xy(S_up->getMeanCurv(), S_up->getMeanCurv(), s2);
    axpy(static_cast<typename SurfContainer::value_type>(-1.0),
         S_up->getGaussianCurv(), s2, s2);
    xy(S_up->getMeanCurv(), s2, s2);

    S_up->grad(S_up->getMeanCurv(), Fb_up);
    S_up->div(Fb_up, s1);

    axpy(static_cast<typename SurfContainer::value_type>(2), s2, s1, s1);
    xv(s1, S_up->getNormal(), Fb_up);
    av(ves_props_.bending_coeff,Fb_up, Fb_up);

    { // downsample Fb
      Vec_t wrk[2]; // TODO: Pre-allocate
      wrk[0].resize(Fb_up.getNumSubs(), params_.upsample_freq);
      wrk[1].resize(Fb_up.getNumSubs(), params_.upsample_freq);
      Fb.replicate(S.getPosition());
      Resample(Fb_up, sht_up_, sht_, wrk[0], wrk[1], Fb);
    }
}

template<typename SurfContainer>
void InterfacialForce<SurfContainer>::linearBendingForce(const SurfContainer &S,
    const Vec_t &x_new, Vec_t &Fb) const
{
    S.resample(params_.upsample_freq, &S_up); // upsample
    Vec_t x_new_up;
    { // upsample x_new
      Vec_t wrk[2]; // TODO: Pre-allocate
      wrk[0]  .resize(x_new.getNumSubs(), params_.upsample_freq);
      wrk[1]  .resize(x_new.getNumSubs(), params_.upsample_freq);
      x_new_up.resize(x_new.getNumSubs(), params_.upsample_freq);
      Resample(x_new, sht_, sht_up_, wrk[0], wrk[1], x_new_up);
    }

    Vec_t Fb_up;
    Fb_up.replicate(S_up->getPosition());
    s1.replicate(S_up->getPosition());
    s2.replicate(S_up->getPosition());

    S_up->linearizedMeanCurv(x_new_up, s1);

    xy(S_up->getMeanCurv(), S_up->getMeanCurv(), s2);
    axpy(static_cast<typename SurfContainer::value_type>(-1.0),
        S_up->getGaussianCurv(), s2, s2);
    xy(s1, s2, s2);

    S_up->grad(s1, Fb_up);
    S_up->div(Fb_up, s1);
    axpy(static_cast<typename SurfContainer::value_type>(2), s2, s1, s1);

    xv(s1, S_up->getNormal(), Fb_up);
    av(ves_props_.bending_coeff, Fb_up, Fb_up);

    { // downsample Fb
      Vec_t wrk[2]; // TODO: Pre-allocate
      wrk[0].resize(Fb_up.getNumSubs(), params_.upsample_freq);
      wrk[1].resize(Fb_up.getNumSubs(), params_.upsample_freq);
      Fb.replicate(S.getPosition());
      Resample(Fb_up, sht_up_, sht_, wrk[0], wrk[1], Fb);
    }
}

template<typename SurfContainer>
void InterfacialForce<SurfContainer>::linearBendingForcePerVesicle(const SurfContainer &S,
    const Vec_t &x_new, Vec_t &Fb, const int vesicle_i) const
{
    // TODO: add S_ups for omp threads
    int tid = omp_get_thread_num();
    S.resample(params_.upsample_freq, &S_ups[tid]); // upsample
    Vec_t x_new_up;
    { // upsample x_new
      Vec_t wrk[2]; // TODO: Pre-allocate
      wrk[0]  .resize(x_new.getNumSubs(), params_.upsample_freq);
      wrk[1]  .resize(x_new.getNumSubs(), params_.upsample_freq);
      x_new_up.resize(x_new.getNumSubs(), params_.upsample_freq);
      Resample(x_new, sht_, sht_up_, wrk[0], wrk[1], x_new_up);
    }

    Vec_t Fb_up;
    Sca_t s1, s2;
    Fb_up.replicate(S_ups[tid]->getPosition());
    s1.replicate(S_ups[tid]->getPosition());
    s2.replicate(S_ups[tid]->getPosition());

    S_ups[tid]->linearizedMeanCurv(x_new_up, s1);

    xy(S_ups[tid]->getMeanCurv(), S_ups[tid]->getMeanCurv(), s2);
    axpy(static_cast<typename SurfContainer::value_type>(-1.0),
        S_ups[tid]->getGaussianCurv(), s2, s2);
    xy(s1, s2, s2);

    S_ups[tid]->grad(s1, Fb_up);
    S_ups[tid]->div(Fb_up, s1);
    axpy(static_cast<typename SurfContainer::value_type>(2), s2, s1, s1);

    xv(s1, S_ups[tid]->getNormal(), Fb_up);
    
    // get the i_th vesicle bending coeff 
    Arr_t bending_coeff_tmp(1);
    (*bending_coeff_tmp.begin()) = (ves_props_.bending_coeff.begin()[vesicle_i]);

    av(bending_coeff_tmp, Fb_up, Fb_up);
    //av(ves_props_.bending_coeff, Fb_up, Fb_up);

    { // downsample Fb
      Vec_t wrk[2]; // TODO: Pre-allocate
      wrk[0].resize(Fb_up.getNumSubs(), params_.upsample_freq);
      wrk[1].resize(Fb_up.getNumSubs(), params_.upsample_freq);
      Fb.replicate(S.getPosition());
      Resample(Fb_up, sht_up_, sht_, wrk[0], wrk[1], Fb);
    }
}

template<typename SurfContainer>
void InterfacialForce<SurfContainer>::tensileForce(const SurfContainer &S,
    const Sca_t &tension, Vec_t &Fs) const
{
    S.resample(params_.upsample_freq, &S_up); // upsample
    Sca_t tension_up;
    { // upsample tension
      Sca_t wrk[2]; // TODO: Pre-allocate
      wrk[0]    .resize(tension.getNumSubs(), params_.upsample_freq);
      wrk[1]    .resize(tension.getNumSubs(), params_.upsample_freq);
      tension_up.resize(tension.getNumSubs(), params_.upsample_freq);
      Resample(tension, sht_, sht_up_, wrk[0], wrk[1], tension_up);
    }

    Vec_t Fs_up;
    Fs_up.replicate(S_up->getPosition());
    v1.replicate(S_up->getPosition());

    xv(S_up->getMeanCurv(), S_up->getNormal(), Fs_up);
    xv(tension_up, Fs_up, Fs_up);
    S_up->grad(tension_up, v1);
    axpy(static_cast<typename SurfContainer::value_type>(2), Fs_up, v1, Fs_up);

    { // downsample Fs
      Vec_t wrk[2]; // TODO: Pre-allocate
      wrk[0].resize(Fs_up.getNumSubs(), params_.upsample_freq);
      wrk[1].resize(Fs_up.getNumSubs(), params_.upsample_freq);
      Fs.replicate(S.getPosition());
      Resample(Fs_up, sht_up_, sht_, wrk[0], wrk[1], Fs);
    }
}

template<typename SurfContainer>
void InterfacialForce<SurfContainer>::tensileForcePerVesicle(const SurfContainer &S,
    const Sca_t &tension, Vec_t &Fs) const
{
    // TODO: add S_ups for omp
    int tid = omp_get_thread_num();
    S.resample(params_.upsample_freq, &S_ups[tid]); // upsample
    Sca_t tension_up;
    { // upsample tension
      Sca_t wrk[2]; // TODO: Pre-allocate
      wrk[0]    .resize(tension.getNumSubs(), params_.upsample_freq);
      wrk[1]    .resize(tension.getNumSubs(), params_.upsample_freq);
      tension_up.resize(tension.getNumSubs(), params_.upsample_freq);
      Resample(tension, sht_, sht_up_, wrk[0], wrk[1], tension_up);
    }

    Vec_t Fs_up, v1;
    Fs_up.replicate(S_ups[tid]->getPosition());
    v1.replicate(S_ups[tid]->getPosition());

    xv(S_ups[tid]->getMeanCurv(), S_ups[tid]->getNormal(), Fs_up);
    xv(tension_up, Fs_up, Fs_up);
    S_ups[tid]->grad(tension_up, v1);
    axpy(static_cast<typename SurfContainer::value_type>(2), Fs_up, v1, Fs_up);

    { // downsample Fs
      Vec_t wrk[2]; // TODO: Pre-allocate
      wrk[0].resize(Fs_up.getNumSubs(), params_.upsample_freq);
      wrk[1].resize(Fs_up.getNumSubs(), params_.upsample_freq);
      Fs.replicate(S.getPosition());
      Resample(Fs_up, sht_up_, sht_, wrk[0], wrk[1], Fs);
    }
}

template<typename SurfContainer>
void InterfacialForce<SurfContainer>::gravityForce(const SurfContainer &S, const Vec_t &x, Vec_t &Fg) const
{
    value_type zero(0);
    value_type one(1.0);

    if(!ves_props_.has_excess_density){
        Vec_t::getDevice().Memset(Fg.begin(), zero, Fg.mem_size());
        return;
    }

    /*
     * The traction jump due to gravity is (excess_density)*(g.x) n
     */

    // g.x
    Sca_t s;
    Vec_t w;
    s.replicate(x);
    w.replicate(x);

    ShufflePoints(x,w);
    int m(1), k(VES3D_DIM), n(x.size()/VES3D_DIM);

    Vec_t::getDevice().gemm("N","N", &m, &n, &k,
        &one, params_.gravity_field, &m,
        w.begin(), &k,
        &zero, s.begin(), &m);

    // (g.x) n
    xv(s,S.getNormal(),Fg);

    // scale by excess density
    av(ves_props_.excess_density,Fg,Fg);
}

template<typename SurfContainer>
void InterfacialForce<SurfContainer>::gravityForcePerVesicle(const SurfContainer &S, const Vec_t &x, Vec_t &Fg, 
        const int vesicle_i) const
{
    value_type zero(0);
    value_type one(1.0);

    if(!ves_props_.has_excess_density){
        Vec_t::getDevice().Memset(Fg.begin(), zero, Fg.mem_size());
        return;
    }

    /*
     * The traction jump due to gravity is (excess_density)*(g.x) n
     */

    // g.x
    Sca_t s;
    Vec_t w;
    s.replicate(x);
    w.replicate(x);

    ShufflePoints(x,w);
    int m(1), k(VES3D_DIM), n(x.size()/VES3D_DIM);

    Vec_t::getDevice().gemm("N","N", &m, &n, &k,
        &one, params_.gravity_field, &m,
        w.begin(), &k,
        &zero, s.begin(), &m);

    // (g.x) n
    xv(s,S.getNormal(),Fg);

    // scale by excess density
    // get the i_th vesicle excess_density
    Arr_t excess_density_tmp(1);
    (*excess_density_tmp.begin()) = (ves_props_.excess_density.begin()[vesicle_i]);

    av(excess_density_tmp,Fg,Fg);
    //av(ves_props_.excess_density,Fg,Fg);
}

template<typename SurfContainer>
void InterfacialForce<SurfContainer>::explicitTractionJump(const SurfContainer &S, Vec_t &F) const
{
    bendingForce(S, F);

    ftmp.replicate(S.getPosition());
    gravityForce(S, S.getPosition(), ftmp);
    axpy((value_type) 1.0, F, ftmp, F);
}

template<typename SurfContainer>
void InterfacialForce<SurfContainer>::implicitTractionJump(const SurfContainer &S, const Vec_t &x,
    const Sca_t &tension, Vec_t &F) const
{
    linearBendingForce(S, x, F);

    ftmp.replicate(x);
    gravityForce(S, x, ftmp);
    axpy((value_type) 1.0, F, ftmp, F);

    tensileForce(S, tension, ftmp);
    axpy( (value_type) 1.0, F, ftmp, F);
}

template<typename SurfContainer>
void InterfacialForce<SurfContainer>::implicitTractionJumpPerVesicle(const SurfContainer &S, const Vec_t &x,
    const Sca_t &tension, Vec_t &F, const int vesicle_i) const
{
    linearBendingForcePerVesicle(S, x, F, vesicle_i);

    Vec_t ftmp;

    ftmp.replicate(x);
    gravityForcePerVesicle(S, x, ftmp, vesicle_i);
    axpy((value_type) 1.0, F, ftmp, F);

    tensileForcePerVesicle(S, tension, ftmp);
    axpy( (value_type) 1.0, F, ftmp, F);
}
