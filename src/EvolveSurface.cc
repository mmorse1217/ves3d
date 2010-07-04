template<typename SurfContainer>
Monitor<SurfContainer>::Monitor(const Parameters<value_type> &params) : 
    time_hor_(params.time_horizon),
    save_flag_(params.save_data),
    save_stride_(params.save_stride),
    IO(SurfContainer::Vec_t::getDevice(), params.save_file_name),
    A0(-1),
    V0(-1)
{};
    
template<typename SurfContainer>
bool Monitor<SurfContainer>::operator()(const SurfContainer &state, 
    value_type &t, value_type &dt)
{
    typename SurfContainer::Sca_t area, vol;
    area.replicate(state.getPosition());
    vol.replicate(state.getPosition());
    state.area(area);
    state.volume(vol);
    
    value_type A(area.getDevice().MaxAbs(area.begin(), state.getPosition().getNumSubs()));
    value_type V( vol.getDevice().MaxAbs( vol.begin(), state.getPosition().getNumSubs()));
    
    if(A0 == -1)
    {
        A0 = A;
        V0 = V;
    }
    
    if(save_flag_ && (t/save_stride_ == static_cast<int>(t/save_stride_)))
    {
        IO.Append(state.getPosition().begin(), state.getPosition().size());
    }

    COUT(                              "\n Monitor : t            = "<<fixed<<t
        <<scientific<<setprecision(4)<<"\n           area   error = "<<abs(A/A0-1)
        <<scientific<<setprecision(4)<<"\n           volume error = "<<abs(V/V0-1)<<endl);
    return(t<time_hor_);
}

template<typename Container, typename Interaction>
EvolveSurface<Container, Interaction>::EvolveSurface(OperatorsMats<value_type, 
    device_type> &mats, const Parameters<value_type> &params) : 
    mats_(mats), params_(params){}

template<typename Container, typename Interaction>
void EvolveSurface<Container, Interaction>::operator()(Container &S, 
    Interaction &Inter)
{
    value_type t(0);
    value_type dt(params_.ts);
    IntVel_t F(S, Inter, mats_, params_);
    Monitor<Container> M(params_);
    
    Scheme_t updater;
    switch ( params_.scheme )
    {
        case Explicit:
            updater = &IntVel_t::updatePositionExplicit; 
            break;

        case SemiImplicit:
            updater = &IntVel_t::updatePositionImplicit; 
            break;
    }
    
    typename Container::Vec_t velocity;

    while ( M(S, t, dt) )
    {
        velocity.replicate(S.getPosition());
        (F.*updater)(dt);       
        F.reparam();
        t += dt;
    }
}
