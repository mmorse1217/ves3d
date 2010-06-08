template<typename Container>
map<int, Container*> GaussLegendreIntegrator<Container>::qw_map;

template<typename Container>
GaussLegendreIntegrator<Container>::GaussLegendreIntegrator() :
    IO(Container::getDevice())
{}

template<typename Container>
GaussLegendreIntegrator<Container>::~GaussLegendreIntegrator()
{
    typename map<int, Container*>::iterator it(qw_map.begin());
    for ( ; it != qw_map.end(); ++it)
        delete it->second;
}

template<typename Container>
Container* GaussLegendreIntegrator<Container>::getQuadWeight(int key) const
{
    typename map<int, Container*>::iterator it(qw_map.find(key));
    
    if(it == qw_map.end())
    {
        ///@todo this assumes to much about the container
        Container* qw = new Container(1, key);
        it = qw_map.insert(qw_map.begin(), make_pair(key, qw));
        
        char fname[500];
        sprintf(fname,"precomputed/quad_weights_%u_single.txt",key);
        IO.ReadData(fname, qw->getStride(), qw->begin());
}
    return(it->second);
}

template<typename Container>
template<typename InputContainer>
void GaussLegendreIntegrator<Container>::operator()(const InputContainer &x_in, 
    const Container &w_in, InputContainer &x_dw) const
{
    Reduce(x_in, w_in, *getQuadWeight(x_in.getShOrder()), x_dw);
}


template<typename Container>
void GaussLegendreIntegrator<Container>::operator()(const Container &w_in, 
    Container &dw) const
{
    Reduce(w_in, *getQuadWeight(w_in.getShOrder()), dw);
}
