template<typename Container, typename Operators>
MovePole<Container, Operators>::MovePole(Operators &mats) :
    p_(mats.p_),
    np_(SpharmGridDim(p_).first * SpharmGridDim(p_).second),
    sp_harm_mats_(SpharmGridDim(p_).first, 1,
        std::make_pair(1, p_ * (4 * p_ * p_ -  1)/3 + 4 * p_ * p_)),
    all_rot_mats_(SpharmGridDim(p_).first, 1, std::make_pair(np_, np_)),
    longitude_rot_(SpharmGridDim(p_).second, 1, std::make_pair(1, 4*p_ - 2)),
    row_idx((int*) Container::getDevice().Malloc((4*p_ - 2) * sizeof(int))),
    col_idx((int*) Container::getDevice().Malloc((4*p_ - 2) * sizeof(int))),
    sht_(p_, mats.mats_p_),
    rot_handle_(NULL),
    last_rot_(Direct),
    arr_(NULL),
    num_(0),
    alpha(1.0),
    beta(0.0),
    eager_n_stream_(4),
    rot_mat_(SpharmGridDim(p_).first, 1, std::make_pair(SpharmGridDim(p_).second, np_)),
    shc_(NULL),
    eager_results_(NULL),
    eager_last_latitude_(-1)
{
    if ( mats.all_rot_mats_ == NULL )
        all_rot_mats_.resize(0);
    else
        Container::getDevice().Memcpy(all_rot_mats_.begin(), mats.all_rot_mats_,
            all_rot_mats_.size() * sizeof(value_type),
            Container::getDevice().MemcpyDeviceToDevice);

    if (mats.sh_rot_mats_ == NULL )
        sp_harm_mats_.resize(0);
    else
        Container::getDevice().Memcpy(sp_harm_mats_.begin(),
            mats.sh_rot_mats_,
            sp_harm_mats_.size() * sizeof(value_type),
            Container::getDevice().MemcpyDeviceToDevice);

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
        size * sizeof(value_type),
        Container::getDevice().MemcpyHostToDevice);
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
        (4*p_ - 2) * sizeof(int),
        Container::getDevice().MemcpyHostToDevice);
    Container::getDevice().Memcpy(col_idx, col,
        (4*p_ - 2) * sizeof(int),
        Container::getDevice().MemcpyHostToDevice);

    delete[] row;
    delete[] col;
}

template<typename Container, typename Operators>
MovePole<Container, Operators>::~MovePole()
{
    Container::getDevice().Free(row_idx);
    Container::getDevice().Free(col_idx);
    delete[] shc_;
    delete[] eager_results_;
}

template<typename Container, typename Operators>
void MovePole<Container, Operators>::setOperands(const Container** arr,
    int num, enum SingularStokesRot rot_scheme)
{
    arr_ = arr;

    switch (rot_scheme)
    {
        case Direct:
            rot_handle_ = &MovePole::movePoleDirectly;
            break;

        case ViaSpHarm:
            rot_handle_ = &MovePole::movePoleViaSpHarm;
            if ( num != num_ || rot_scheme != last_rot_)
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

            break;

        case DirectEagerEval:
            rot_handle_ = &MovePole::movePoleDirectly;
            eager_last_latitude_ = -1;
            if ( eager_results_ == NULL || num_ != num )
            {
                delete[] eager_results_;
                eager_results_ = new Container[num * SpharmGridDim(p_).second];
            }

            for(int ii=0; ii<SpharmGridDim(p_).second; ++ii)
                for(int jj=0; jj<num; ++jj)
                    eager_results_[ii * num + jj].replicate(*(arr_[jj]));
            eager_wrk_.resize(eager_n_stream_,1,std::make_pair(np_,np_));
            break;
    }
    last_rot_ = rot_scheme;
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
    PROFILESTART();

    if ( last_rot_ == Direct )
    {
        CircShift(all_rot_mats_.begin() + trg_i * np_ * np_, trg_j * np_, rot_mat_);

        for(int ii=0; ii<num_; ++ii)
        {
            int nsub(arr_[ii]->getNumSubs());
            Container::getDevice().gemm("N", "N", &np_, &nsub, &np_,
                &alpha, rot_mat_.begin(), &np_, arr_[ii]->begin(),
                &np_, &beta, results[ii]->begin(), &np_);
        }
    }
    else
    {
        if ( trg_i != eager_last_latitude_ )
            updateEagerResults(eager_last_latitude_ = trg_i);

        for(int ii=0; ii<num_; ++ii)
            Container::getDevice().Memcpy(results[ii]->begin(),
                eager_results_[num_ * trg_j + ii].begin(),
                arr_[ii]->size() * sizeof(value_type),
                Container::getDevice().MemcpyDeviceToDevice);
    }
    PROFILEEND("",0);
}

template<typename Container, typename Operators>
void MovePole<Container, Operators>::movePoleViaSpHarm(int trg_i, int trg_j,
    Container** results) const
{
    PROFILESTART();
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
        rotmat = sp_harm_mats_.getSubN_begin(trg_i);

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
    PROFILEEND("",0);
}

template<typename Container, typename Operators>
void MovePole<Container, Operators>::alignMeridian(int trg_j,
    Container** results) const
{
    PROFILESTART();
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
        rotmat = const_cast<value_type*>(longitude_rot_.getSubN_begin(trg_j));

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
    PROFILEEND("",0);
}

template<typename Container, typename Operators>
void MovePole<Container, Operators>::updateEagerResults(int trg_i) const
{
    int *n_sub = new int[num_];
    const value_type** src = new const value_type*[num_];

    for(int ii=0; ii<num_; ++ii)
    {
        n_sub[ii] = arr_[ii]->getNumSubs();
        src[ii] = arr_[ii]->begin();
    }

    int n_lat = SpharmGridDim(p_).first;
    int n_long = SpharmGridDim(p_).second;
    int n_res =  num_ * n_long;
    value_type **res = new value_type*[n_res];
    for(int ii=0; ii<n_res; ++ii)
        res[ii] = eager_results_[ii].begin();

    value_type** wrk = new value_type*[eager_n_stream_];
    for (int jj = 0; jj < eager_n_stream_; ++jj)
        wrk[jj] = eager_wrk_.getSubN_begin(jj);

    Container::getDevice().AggregateRotation(
        p_, n_lat, n_long, num_, n_sub,
        all_rot_mats_.begin() + trg_i * np_ * np_,
        src, wrk, res, eager_n_stream_);

    delete[] n_sub;
    delete[] src;
    delete[] res;
    delete[] wrk;
}
