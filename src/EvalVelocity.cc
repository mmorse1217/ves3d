template<typename Scalar, typename Vector, typename StokesEvaluator>
EvalVelocity<Scalar, Vector, StokesEvaluator>::EvalVelocity(
    const StokesEvaluator &stokes, const BgFlowBase<Vector> &vInf,
    OperatorsMats<Arr_t> &mats, value_type bending_modulus) :
    stokes_(stokes),
    vInf_(vInf),
    mats_(mats),
    S_ptr_(NULL),
    Force_(bending_modulus)
{
    quad_weights_.resize(1,mats_.p_);
    quad_weights_.getDevice().Memcpy(quad_weights_.begin(),
        mats.quad_weights_, quad_weights_.size()
        * sizeof(value_type),
        quad_weights_.getDevice().MemcpyDeviceToDevice);
}

template<typename Scalar, typename Vector, typename StokesEvaluator>
EvalVelocity<Scalar, Vector, StokesEvaluator>::~EvalVelocity()
{
    delete S_ptr_;
}

template<typename Scalar, typename Vector, typename StokesEvaluator>
Error_t EvalVelocity<Scalar, Vector, StokesEvaluator>::operator()(
    const Vector &x_src, const Scalar &tension, const Vector &x_eval,
        Vector &vel)
{
    if (S_ptr_ == NULL)
        S_ptr_ = new Sur_t(x_src, mats_);
    else
        S_ptr_->setPosition(x_src);

    //Getting the interfacial Force
    Fb.replicate(x_src);
    Fs.replicate(x_src);

    Force_.bendingForce(*S_ptr_, Fb);
    Force_.tensileForce(*S_ptr_, tension, Fs);
    axpy(static_cast<value_type>(1), Fb, Fs, Fb);

    //Incorporating the quadrature weights into the density
    xv(S_ptr_->getAreaElement(), Fb, Fb);
    ax<Scalar>(quad_weights_,Fb, Fb);

    //Collecting all points to pass to FMM
    int src_size = x_src.size();
    int eval_size= x_eval.size();
    int all_size = ( src_size + eval_size ) / DIM;

    all_src.resize(1,1, std::make_pair(all_size,1));
    all_den.replicate(all_src);
    all_pot.replicate(all_src);

    ShufflePoints(x_src, Fs);
    Scalar::getDevice().Memcpy(all_src.begin(), Fs.begin(),
        src_size * sizeof(value_type),
        Scalar::getDevice().MemcpyDeviceToDevice);

    ShufflePoints(x_eval, all_den);
    Scalar::getDevice().Memcpy(all_src.begin() + src_size,
        all_den.begin(), eval_size * sizeof(value_type),
        Scalar::getDevice().MemcpyDeviceToDevice);

    ShufflePoints(Fb, Fs);
    Scalar::getDevice().Memcpy(all_den.begin(), Fs.begin(),
        src_size * sizeof(value_type),
        Scalar::getDevice().MemcpyDeviceToDevice);
    Vector::getDevice().Memset(all_den.begin() + src_size, 0,
        eval_size * sizeof(value_type));

    stokes_(all_src.begin(), all_den.begin(), all_size,
        all_pot.begin(), NULL);

    Scalar::getDevice().Memcpy(vel.begin(), all_pot.begin()
        + src_size, eval_size * sizeof(value_type),
        Scalar::getDevice().MemcpyDeviceToDevice);

    all_den.replicate(vel);
    vel.setPointOrder(PointMajor);
    ShufflePoints(vel, all_den);
    vel.setPointOrder(AxisMajor);

    //The background flow
    vInf_(x_eval, 0, vel);
    axpy(static_cast<value_type>(1), vel, all_den, vel);

    //cleaning up
    Fs.setPointOrder(AxisMajor);

    return ErrorEvent::Success;
}
