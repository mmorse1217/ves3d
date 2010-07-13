template<typename Container, typename Operators>
MovePole<Container, Operators>::MovePole(Operators &mats, 
    Container &sp_mats) :
    p_(mats.p_),
    np_(gridDimOf(p_).first * gridDimOf(p_).second),
    sp_harm_mats_(p_ + 1, 1, make_pair(1, p_ * (4 * p_ * p_ -  1)/3 + 4 * p_ * p_)),
    all_rot_mats_(p_ + 1, 1, make_pair(np_, np_)),
    long_rot_(2 * p_, 1, make_pair(2 * p_, 2 * p_)),
    sht_(p_, mats.mats_p_),
    rot_handle_(NULL),
    arr_(NULL),
    num_(0),
    alpha(1.0),
    beta(0.0),
    rot_mat_(p_ +1, 1, make_pair(2*p_, np_)),
    shc_(NULL)
{
    Container::getDevice().Memcpy(all_rot_mats_.begin(), mats.all_rot_mats_,
        all_rot_mats_.size() * sizeof(typename Container::value_type),
        MemcpyDeviceToDevice);

    Container::getDevice().Memcpy(sp_harm_mats_.begin(), sp_mats.begin(),
        sp_harm_mats_.size() * sizeof(typename Container::value_type),
        MemcpyDeviceToDevice);
}

template<typename Container, typename Operators>
MovePole<Container, Operators>::~MovePole()
{
    delete[] shc_;
}

template<typename Container, typename Operators>
void MovePole<Container, Operators>::setOperands(const Container** arr, 
    int num, enum SingularStokesRot rot_scheme)
{
    arr_ = arr;
        
    if ( rot_scheme == Direct )
        rot_handle_ = &MovePole::movePoleDirectly;
    else
    {
        rot_handle_ = &MovePole::movePoleViaSpHarm;
        if ( num != num_ )
        {
            delete[] shc_;
            
            shc_ = new Container[num];
        }
        
        for(int ii=0; ii<num; ++ii)
        {
            shc_[ii].replicate(*(arr_[ii]));
            wrk_.replicate(*(arr_[ii]));
            sht_.forward(*(arr_[ii]), shc_[ii], wrk_);
            this->shuffleShCoeff(wrk_, shc_[ii]);
        }
    }
    num_ = num;
}

template<typename Container, typename Operators>
void MovePole<Container, Operators>::operator()(int trg_i, int trg_j, 
    Container** results) const
{
    (this->*rot_handle_)(trg_i, trg_j, results);
}

template<typename Container, typename Operators>
void MovePole<Container, Operators>::movePoleDirectly(int trg_i, int trg_j, 
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
        Container::getDevice().gemm("N", "N", &np, &nsub, &np, 
            &alpha, rot_mat_.begin(), &np, arr_[ii]->begin(), 
            &np, &beta, results[ii]->begin(), &np);
    }
}
 
template<typename Container, typename Operators>
void MovePole<Container, Operators>::movePoleViaSpHarm(int trg_i, int trg_j, 
    Container** results) const
{       
    const value_type* rotmat(sp_harm_mats_.getSubN(trg_i));
    
    for(int ii=0; ii<num_; ++ii)
    {
        wrk_.replicate(*(results[ii]));
        shc_out.replicate(wrk_);
        
        value_type* srcPtr(shc_[ii].begin());
        value_type* resPtr(wrk_.begin());

        int nsub(arr_[ii]->getNumSubs());
        for(int jj=0; jj<=p_; ++jj)          
        {
            int matsize = 2*jj + 1 - (jj/p_);       
            Container::getDevice().gemm("N", "N", &matsize, &nsub, &matsize, 
                &alpha, rotmat, &matsize, srcPtr, &matsize, &beta, resPtr,
                &matsize);

            rotmat += matsize * matsize;
            srcPtr += matsize * nsub;
            resPtr += matsize * nsub;
        }

       unShuffleShCoeff(wrk_, shc_out);
       sht_.backward(shc_out, wrk_, *(results[ii]));
    }
}

template<typename Container, typename Operators>
void MovePole<Container, Operators>::shuffleShCoeff(const Container &in, Container &out) const
{
    const typename Container::value_type *inPtr(NULL);
    typename Container::value_type *outPtr(NULL), *head(out.begin());

    int ns = in.getNumSubs();  
    for(int ii=0; ii<= p_; ++ii)
    {
        int len = 2*ii + 1 - (ii/p_);
                
        inPtr = in.begin() + ii;       
        for(int jj=0; jj<= 2*ii-(ii/p_); ++jj)
        {
            outPtr = head + jj;
            int dist = (p_ + 1 - (jj + 1)/2);
            for(int ss=0; ss<ns; ++ss)
            {
                *outPtr = *inPtr;
                inPtr += dist;
                outPtr+= len;
            }
            inPtr--;
            inPtr += jj%2;
        }
        head += ns * len;
    }
}

template<typename Container, typename Operators>
void MovePole<Container, Operators>::unShuffleShCoeff(const Container &in, Container &out) const
{
    const typename Container::value_type *inPtr(NULL), *head(in.begin());
    typename Container::value_type *outPtr(NULL);

    int ns = in.getNumSubs();  
    for(int ii=0; ii<= p_; ++ii)
    {
        int len = 2*ii + 1 - (ii/p_);
                
        outPtr = out.begin() + ii;       
        for(int jj=0; jj<= 2*ii-(ii/p_); ++jj)
        {
            inPtr = head + jj;
            int dist = (p_ + 1 - (jj + 1)/2);
            for(int ss=0; ss<ns; ++ss)
            {
                *outPtr = *inPtr;
                outPtr += dist;
                inPtr  += len;
            }
            outPtr--;
            outPtr += jj%2;
        }
        head += ns * len;
    }
}
