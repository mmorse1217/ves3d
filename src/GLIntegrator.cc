template<typename Container>
std::map<int, Container*> GaussLegendreIntegrator<Container>::qw_map;

template<typename Container>
GaussLegendreIntegrator<Container>::GaussLegendreIntegrator()
{}

template<typename Container>
GaussLegendreIntegrator<Container>::~GaussLegendreIntegrator()
{
    ///@bug this defeats the purpose of static, map is cleared
    typename std::map<int, Container*>::iterator it(qw_map.begin());
    for ( ; it != qw_map.end(); ++it)
        delete it->second;
    qw_map.clear();
}

template<typename Container>
Container* GaussLegendreIntegrator<Container>::getQuadWeights(
    int key) const
{
    typename std::map<int, Container*>::iterator it(qw_map.find(key));

    if(it == qw_map.end())
    {
        Container* qw(buildQuadWeights(key));
        it = qw_map.insert(qw_map.begin(), std::make_pair(key, qw));
    }

    return(it->second);
}

template<typename Container>
Container* GaussLegendreIntegrator<Container>::buildQuadWeights(
    int shOrder) const
{
    int order = SpharmGridDim(shOrder).first;
    int d2 = SpharmGridDim(shOrder).second;

    typedef typename Container::value_type vt;

    vt* d = new vt[order];
    vt* e = new vt[order-1];
    vt* work = new vt[order * ((d2 > 2) ? d2 : 2)];
    vt* eig = new vt[order * order];
    vt tmp;

    //Calculating the quadrature points
    for (int ii(0); ii<order-1; ++ii)
    {
        tmp = 2 * (ii + 1);
        tmp *= tmp;
        tmp = static_cast<vt>(-1)/tmp;
        tmp += 1;
        tmp = sqrt(tmp);
        *(e + ii) = static_cast<vt>(.5)/tmp;
        *(d + ii) = 0;
    }
    *(d + order - 1) = 0;

     int info;
     Steqr("I", order, d, e, eig, order, work, info);

    //d holds the quadrature points

    //replicating for the grid, multiplying by the weight of the
    //uniform grid in the longitude direction and dividing by the
    //sin(latitude)
    for (int ii=0;ii<order; ++ii)
        for(int jj=0;jj<d2; ++jj)
            work[ii * d2 + jj] = (2 * eig[ii * order] * eig[ii * order]) *
                (2 * M_PI/d2 ) / sqrt(1 - d[ii] * d[ii]);

    Container* qw = new Container(1, shOrder);
    qw->getDevice().Memcpy(qw->begin(),
        work,
        order * d2 * sizeof(vt),
        Container::device_type::MemcpyHostToDevice);

    delete[] d;
    delete[] e;
    delete[] work;
    delete[] eig;

    return(qw);
}

template<typename Container>
template<typename InputContainer>
void GaussLegendreIntegrator<Container>::operator()(const InputContainer &x_in,
    const Container &w_in, InputContainer &x_dw) const
{
    Reduce(x_in, w_in, *getQuadWeights(x_in.getShOrder()), x_dw);
}


template<typename Container>
void GaussLegendreIntegrator<Container>::operator()(const Container &w_in,
    Container &dw) const
{
    Reduce(w_in, *getQuadWeights(w_in.getShOrder()), dw);
}
