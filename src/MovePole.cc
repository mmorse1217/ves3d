template<typename Container>
MovePole<Container>::MovePole(const Container &all_rot_mats, 
    Container &rot_mat, const Container &sp_mats) :
    sp_harm_mats_(sp_mats),
    all_rot_mats_(all_rot_mats),
    rot_mat_(rot_mat),
    arr_(NULL),
    num_(0),
    rot_handle_(NULL),
    alpha(1.0),
    beta(0.0),
    shc_(NULL)
{}

template<typename Container>
void MovePole<Container>::setOperands(const Container** arr, 
    int num, enum SingularStokesRot rot_scheme)
{
    arr_ = arr;
        
    if ( rot_scheme == Direct )
        rot_handle_ = &MovePole::movePoleDirectly;
    else
    {
        rot_handle_ = &MovePole::movePoleViaSpHarm;
        if( num != num_)
        {
            delete[] shc_;
            
            shc_ = new Container[num];
            num_ = num;
        }
        
        for(int ii=0; ii<num_; ++ii)
        {
            shc_[ii].replicate(*(arr_[ii]));
            wrk_.replicate(*(arr_[ii]));
        }
    }
}

template<typename Container>
void MovePole<Container>::operator()(int trg_i, int trg_j, 
    Container** results) const
{
    (this->*rot_handle_)(trg_i, trg_j, results);
}

template<typename Container>
void MovePole<Container>::movePoleDirectly(int trg_i, int trg_j, 
    Container** results) const
{
    int p  = (*arr_)->getShOrder();
    int np = (*arr_)->getStride();
    int nv = (*arr_)->getNumSubs();
    //int rot_chunck = 2 * p * np;
    
    CircShift(all_rot_mats_.begin() + trg_i * np * np, trg_j * np, rot_mat_);
    
    for(int ii=0; ii<num_; ++ii)
    {
        int nsub(arr_[ii]->getNumSubs());
        arr_[ii]->getDevice().gemm("N", "N", &np, &nsub, &np, 
            &alpha, rot_mat_.begin(), &np, arr_[ii]->begin(), 
            &np, &beta, results[ii]->begin(), &np);
    }
}
 
template<typename Container>
void MovePole<Container>::movePoleViaSpHarm(int trg_i, int trg_j, 
    Container** results) const
{
    movePoleDirectly(trg_i, trg_j, results);
}
