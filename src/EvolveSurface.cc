template<typename SurfContainer>
Monitor<SurfContainer>::Monitor() : 
    time_hor_(Parameters<value_type>::getInstance().time_horizon),
    save_flag_(Parameters<value_type>::getInstance().save_data),
    save_stride_(Parameters<value_type>::getInstance().save_stride),
    IO(SurfContainer::Vec::getDevice())
{};
    
template<typename SurfContainer>
bool Monitor<SurfContainer>::operator()(const SurfContainer &state, 
    value_type &t, value_type &dt)
{
    typename SurfContainer::Sca A, V;
    A.replicate(state.getPosition());
    V.replicate(state.getPosition());
    state.area(A);
    state.volume(V);
    
    std::cout<<t<<" "<<A[0]<<"\t"<<V[0]<<std::endl;
    return(t<time_hor_);
}

template<typename Container, typename Interaction>
EvolveSurface<Container, Interaction>::EvolveSurface(OperatorsMats<value_type> 
    &mats) : mats_(mats){}

template<typename Container, typename Interaction>
void EvolveSurface<Container, Interaction>::operator()(Container &S, 
    Interaction &Inter)
{
    value_type t(0);
    value_type dt(Parameters<value_type>::getInstance().ts);
    InterfacialVelocity<Container, Interaction> F(S, Inter, mats_);
    Monitor<Container> M;
        
    typename Container::Vec velocity;

    while ( M(S, t, dt) )
    {
        velocity.replicate(S.getPosition());
        F.updatePositionImplicit(dt);       
        F.reparam();
        
        t += dt;
    }
}
