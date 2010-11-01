template<typename Scalar, typename Vector, typename StokesEvaluator>
EvalVelocity<Scalar, Vector, StokesEvaluator>::EvalVelocity(
    const StokesEvaluator &stokes, const BgFlowBase<Vector> &vInf,
    const OperatorsMats<Scalar> &mats, value_type bending_modulus) :
    stokes_(stokes),
    vInf_(vInf),
    mats_(mats),
    S_ptr_(NULL),
    Force_(bending_modulus)
{}

template<typename Scalar, typename Vector, typename StokesEvaluator>
EvalVelocity<Scalar, Vector, StokesEvaluator>::~EvalVelocity()
{
    delete S_ptr_;
}

template<typename Scalar, typename Vector, typename StokesEvaluator>
Error_t EvalVelocity<Scalar, Vector, StokesEvaluator>::operator()(
    const Vector &x_src, const Vector &x_eval, Vector &vel)
{

}

template<typename Scalar, typename Vector, typename StokesEvaluator>
Error_t EvalVelocity<Scalar, Vector, StokesEvaluator>::operator()(
    const Vector &x_src, const Scalar &tension, const Vector &x_eval, 
        Vector &vel)
{}
