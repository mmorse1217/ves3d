template<typename Container>
const typename Container::value_type SHTrans<Container>::alpha_;

template<typename Container>
const typename Container::value_type SHTrans<Container>::beta_;

template<typename Container>
SHTrans<Container>::SHTrans(int p_in) :
    device_(&Container::getDevice()),
    mats_(device_, p_in),
    p(p_in),
    dft_size(2*p),
    filter_coeff_((value_type*) device_->Malloc(p * (p + 2) * sizeof(value_type)))
{
    int ll = p * (p + 2);
    value_type *buffer = (value_type*) malloc(p * (p + 2) * sizeof(value_type));
    
    int idx = 0, len;
    for(int ii=0; ii< 2 * p; ++ii)
    {
        len = p + 1 - (ii+1)/2;
        for(int jj=0; jj < len; ++jj)
             buffer[idx++] = (len-jj)<=(p-2*p/3) ? 0 : 1;
    }
    device_->Memcpy(filter_coeff_, buffer, p *(p + 2) * sizeof(value_type), 
        MemcpyHostToDevice);
    free(buffer);
}

template<typename Container>
SHTrans<Container>::~SHTrans() 
{
    device_->Free(filter_coeff_);
}

template<typename Container>
void SHTrans<Container>::DLT(typename Container::value_type *trans, 
    const typename Container::value_type *inputs, 
    typename Container::value_type *outputs, 
    int m, int n , int k, int mf, int nf, int kf) const 
{
    for (int freq = 0; freq <= p; freq++) {
        int num_legendre_inputs = n;
        if (freq == 0 || freq == p) num_legendre_inputs = n / 2;

        device_->gemm("N", "N", &m, &num_legendre_inputs, &k, &alpha_, 
            trans, &m, inputs, &k, &beta_,outputs, &m);

        trans += m * k;
        inputs += num_legendre_inputs * k;
        outputs += m * num_legendre_inputs;
        if (mf) m--;
        if (nf) n--;
        if (kf) k--;
    }
}

template<typename Container>
void SHTrans<Container>::back(const typename Container::value_type *inputs, 
    typename Container::value_type *work_arr, int n_funs,
    typename Container::value_type *outputs, 
    typename Container::value_type *trans, 
    typename Container::value_type *dft) const 
{    
    int num_dft_inputs = n_funs * (p + 1);
    DLT(trans, inputs, outputs, p + 1, 2 * n_funs, p + 1, 0, 0, 1);
    
    device_->Transpose(outputs, dft_size, num_dft_inputs, work_arr);
    device_->gemm("T", "N", &dft_size, &num_dft_inputs, 
        &dft_size, &alpha_, dft, &dft_size,
        work_arr, &dft_size, &beta_, outputs, &dft_size);
}

template<typename Container>
void SHTrans<Container>::forward(const Container &in, Container &work,
    Container &shc) const
{
    int n_funs = in.getNumSubs();
    int num_dft_inputs = n_funs * (p + 1);
    
    device_->gemm("N", "N", &dft_size, &num_dft_inputs, &dft_size, 
        &alpha_, mats_.dft_, &dft_size,in.begin(), &dft_size, &beta_, 
        shc.begin(), &dft_size);
    device_->Transpose(shc.begin(), num_dft_inputs, dft_size, work.begin());
    DLT(mats_.dlt_, work.begin(), shc.begin(), p + 1, 2 * n_funs, p + 1, 1, 0, 0);
}


template<typename Container>
void SHTrans<Container>::backward(const Container &shc, 
    Container &work, Container &out) const
{
    back(shc.begin(), work.begin(), out.getNumSubs(), out.begin(), 
        mats_.dlt_inv_, mats_.dft_inv_);
}


template<typename Container>
void SHTrans<Container>::backward_du(const Container &shc, 
    Container &work, Container &out) const
{
    back(shc.begin(), work.begin(), out.getNumSubs(), out.begin(), 
        mats_.dlt_inv_d1_, mats_.dft_inv_);
}

template<typename Container>
void SHTrans<Container>::backward_d2u(const Container &shc, 
    Container &work, Container &out) const
{
    back(shc.begin(), work.begin(), out.getNumSubs(), out.begin(), 
        mats_.dlt_inv_d2_, mats_.dft_inv_);
}

template<typename Container>
void SHTrans<Container>::backward_dv(const Container &shc, 
    Container &work, Container &out) const
{
    back(shc.begin(), work.begin(), out.getNumSubs(), out.begin(), 
        mats_.dlt_inv_, mats_.dft_inv_d1_);
}

template<typename Container>
void SHTrans<Container>::backward_d2v(const Container &shc, 
    Container &work, Container &out) const
{
    back(shc.begin(), work.begin(), out.getNumSubs(), out.begin(), 
        mats_.dlt_inv_, mats_.dft_inv_d2_);
}

template<typename Container>
void SHTrans<Container>::backward_duv(const Container &shc, 
    Container &work, Container &out) const
{
    back(shc.begin(), work.begin(), out.getNumSubs(), out.begin(), 
        mats_.dlt_inv_d1_, mats_.dft_inv_d1_);
}

template<typename Container>
void SHTrans<Container>::Filter(const Container &in, Container &work, 
    Container &shc, Container &out) const
{
    int n_funs = in.getNumSubs();
    forward(in, work, shc);
    ScaleFreq(shc.begin(), n_funs, this->filter_coeff_, shc.begin());
    backward(shc, work, out);
}

template<typename Container>
void SHTrans<Container>::ScaleFreq(
    const typename Container::value_type *shc_in, 
    int n_funs, 
    const typename Container::value_type* scaling_coeff, 
    typename Container::value_type *shc_out) const
{
    int leg_order = p+1;
    
    device_->ax(scaling_coeff, shc_in, leg_order, n_funs, shc_out);
    scaling_coeff += leg_order;
    shc_in += n_funs * leg_order;
    shc_out += n_funs * leg_order;
    leg_order--;

    // process remaining frequencies except the last cosine
    for (; leg_order>1; leg_order--) 
    {
        // first process cosine
        device_->ax(scaling_coeff, shc_in, leg_order, n_funs, shc_out);
        scaling_coeff += leg_order;
        shc_in += n_funs * leg_order;
        shc_out += n_funs * leg_order;
        
        // then process sine
        device_->ax(scaling_coeff, shc_in, leg_order, n_funs, shc_out);
        scaling_coeff += leg_order;
        shc_in += n_funs * leg_order;
        shc_out += n_funs * leg_order;
    }
    
    // process last cosine
    device_->ax(scaling_coeff, shc_in, leg_order, n_funs, shc_out);
}

template<typename Container>
void SHTrans<Container>::FirstDerivatives(const Container &in,
    Container &work, Container &shc, Container &du, Container &dv) const
{
    forward(in, work, shc);
    backward_du(shc, work, du);
    backward_dv(shc, work, dv);
}
