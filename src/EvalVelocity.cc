template<typename Container>
Error_t EvalVelocity<Container>::EvalVelocity()
{}

template<typename Container>
Error_t EvalVelocity<Container>::~EvalVelocity()
{}

template<typename Container>
Error_t EvalVelocity<Container>::operator()(const Container &x_src, 
        const Container &x_eval, Container &vel)
{
    return Success;
}
