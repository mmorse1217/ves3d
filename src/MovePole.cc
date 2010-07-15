template<typename Container, typename Operators>
MovePole<Container, Operators>::MovePole(Operators &mats) :
    p_(mats.p_),
    np_(gridDimOf(p_).first * gridDimOf(p_).second),
    sp_harm_mats_(gridDimOf(p_).first, 1, 
        make_pair(1, p_ * (4 * p_ * p_ -  1)/3 + 4 * p_ * p_)),
    all_rot_mats_(gridDimOf(p_).first, 1, make_pair(np_, np_)),
    longitude_rot_(gridDimOf(p_).second, 1, make_pair(1, 4*p_ - 2)),
    row_idx((int*) Container::getDevice().Malloc((4*p_ - 2) * sizeof(int))),
    col_idx((int*) Container::getDevice().Malloc((4*p_ - 2) * sizeof(int))),
    sht_(p_, mats.mats_p_),
    rot_handle_(NULL),
    arr_(NULL),
    num_(0),
    alpha(1.0),
    beta(0.0),
    rot_mat_(gridDimOf(p_).first, 1, make_pair(gridDimOf(p_).second, np_)),
    shc_(NULL)
{
    Container::getDevice().Memcpy(all_rot_mats_.begin(), mats.all_rot_mats_,
        all_rot_mats_.size() * sizeof(value_type), MemcpyDeviceToDevice);

    Container::getDevice().Memcpy(sp_harm_mats_.begin(), mats.sh_rot_mats_,
        sp_harm_mats_.size() * sizeof(value_type), MemcpyDeviceToDevice);
    
    //Generating the longitudinal rotation matrices stored in the
    //coordinate format
    size_t size(longitude_rot_.size());
    value_type* buffer = new value_type[size];
    value_type lambda;
    
    for(int ii=0;ii<longitude_rot_.getNumSubs(); ++ii)
    {
        lambda = (M_PI/p_) * ii;
        *buffer++ = 1;
        for(int jj=1;jj<p_;++jj)
        {
            *buffer++ =   cos(jj * lambda);
            *buffer++ =  -sin(jj * lambda);
            *buffer++ =   sin(jj * lambda);
            *buffer++ =   cos(jj * lambda);
        }   
        *buffer++ = cos(p_ * lambda); 
    }
    buffer -= size;
    Container::getDevice().Memcpy(longitude_rot_.begin(), buffer, 
        size * sizeof(value_type), MemcpyHostToDevice);
    delete[] buffer;
    
    //Saving the coordinates of non-zero elements in the
    //longitude_rot_. This is to be passed to the sparse matrix
    //multiplier. The indices are one-bases (FORTRAN format) for
    //compatibility reasons. Check the sparse BLAS manual for more
    //information.
    int* row = new int[4*p_ - 2];
    int* col = new int[4*p_ - 2];
    
    col[0] = row[0] = 1;
    for(int ii=1;ii<p_;++ii)
    {
        row[4*ii-3] = row[4*ii-1] = col[4*ii-3] = col[4*ii-2] = 2*ii;
        row[4*ii-2] = row[4*ii  ] = col[4*ii-1] = col[4*ii  ] = 2*ii+1;        
    }
    row[4*p_- 3] = col[4*p_- 3] = 2*p_;

    Container::getDevice().Memcpy(row_idx, row, 
        (4*p_ - 2) * sizeof(int), MemcpyHostToDevice);
    Container::getDevice().Memcpy(col_idx, col, 
        (4*p_ - 2) * sizeof(int), MemcpyHostToDevice);

    delete[] row;
    delete[] col;
}

template<typename Container, typename Operators>
MovePole<Container, Operators>::~MovePole()
{
    Container::getDevice().Free(row_idx);
    Container::getDevice().Free(col_idx);
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
            sht_.collectSameOrder(wrk_, shc_[ii]);
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
    alignMeridian(trg_j, results);
    
    const value_type* rotmat(NULL);
    value_type* srcPtr(NULL);
    value_type* resPtr(NULL);
    int nsub;
    
    for(int ii=0; ii<num_; ++ii)
    {
        wrk_.replicate(*(results[ii]));
        shc_out.replicate(wrk_);
                
        srcPtr = results[ii]->begin();
        resPtr = wrk_.begin();
        nsub = shc_[ii].getNumSubs();
        rotmat = sp_harm_mats_.getSubN(trg_i);
        
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
        
        sht_.collectSameFreq(wrk_, shc_out);
        sht_.backward(shc_out, wrk_, *(results[ii]));
    }
}

template<typename Container, typename Operators>
void MovePole<Container, Operators>::alignMeridian(int trg_j, 
    Container** results) const
{       
    ///@todo the sparse matrix matrix multiplier is the cpu
    ///version. coomm should be added to the device.
    assert( typeid(Container::getDevice()) == typeid(Device<CPU>));

    value_type* rotmat(NULL);
    value_type* srcPtr(NULL);
    value_type* resPtr(NULL);
    int nsub, matsize, nnz;
    char* matdescra = "GGFFFF";
    value_type lalpha(1), lbeta(0);

    for(int ii=0; ii<num_; ++ii)
    {
        wrk_.replicate(*(results[ii]));
        
        srcPtr = shc_[ii].begin();
        resPtr = results[ii]->begin();
        nsub = shc_[ii].getNumSubs();       
        rotmat = const_cast<value_type*>(longitude_rot_.getSubN(trg_j));
        
        for(int jj=0; jj<=p_; ++jj)
        {
            matsize = 2*jj + 1 - (jj/p_); 
            nnz = 4*jj + 1 - 3*(jj/p_);
            
            coomm("N", &matsize, &nsub, &matsize, &lalpha, matdescra, 
                rotmat, row_idx, col_idx, &nnz, srcPtr, &matsize, &lbeta, 
                resPtr, &matsize);
            
            srcPtr += matsize * nsub;
            resPtr += matsize * nsub;
        }
    }
}
