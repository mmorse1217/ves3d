//template<typename SurfContainer>
//InterfacialForce<SurfContainer>::InterfacialForce(value_type bending_modulus) :
//    bending_modulus_(bending_modulus), S_up(NULL)
//{}

template<typename SurfContainer>
InterfacialForce<SurfContainer>::InterfacialForce(
    const Parameters<value_type> &params,
    const OperatorsMats<Arr_t> &mats) :
    bending_modulus_(params.bending_modulus),
    sht_   (mats.p_   , mats.mats_p_   ),
    sht_up_(mats.p_up_, mats.mats_p_up_),
    params_(params),
    S_up(NULL)
{}

template<typename SurfContainer>
InterfacialForce<SurfContainer>::~InterfacialForce()
{
  if(S_up) delete S_up;
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
    axpy(-bending_modulus_, Fb_up, Fb_up);

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
    axpy(-bending_modulus_, Fb_up, Fb_up);

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

