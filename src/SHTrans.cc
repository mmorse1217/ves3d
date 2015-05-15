template<typename Container, typename Mats>
const typename Container::value_type SHTrans<Container, Mats>::alpha_;

template<typename Container, typename Mats>
const typename Container::value_type SHTrans<Container, Mats>::beta_;

template<typename Container, typename Mats>
SHTrans<Container, Mats>::SHTrans(int p_in, const Mats &mats, int filter_freq) :
    device_(Container::getDevice()),
    mats_(mats),
    p(p_in),
    dft_size(2*p),
    filter_coeff_((value_type*) device_.Malloc(p * (p + 2) * sizeof(value_type)))
{
    filter_freq = (filter_freq == -1) ? 2*p/3 : filter_freq;
    INFO("Initializing with p="<<p<<", filter_freq="<<filter_freq);

    value_type *buffer = (value_type*) malloc(p * (p + 2) * sizeof(value_type));
    int idx = 0, len;
    for(int ii=0; ii< 2 * p; ++ii)
    {
        len = p + 1 - (ii+1)/2;
        for(int jj=0; jj < len; ++jj)
            buffer[idx++] = (len-jj)<=(p - filter_freq) ? 0 : 1;
    }

    device_.Memcpy(filter_coeff_, buffer, p *(p + 2) * sizeof(value_type),
        device_type::MemcpyHostToDevice);

    free(buffer);
}

template<typename Container, typename Mats>
SHTrans<Container, Mats>::~SHTrans()
{
    device_.Free(filter_coeff_);
}

template<typename Container, typename Mats>
void SHTrans<Container, Mats>::DLT(value_type *trans,
    const value_type *inputs, value_type *outputs,
    int m, int n , int k, int mf, int nf, int kf) const
{
    PROFILESTART();

    for (int freq = 0; freq <= p; freq++) {
        int num_legendre_inputs = n;
        if (freq == 0 || freq == p) num_legendre_inputs = n / 2;

        device_.gemm("N", "N", &m, &num_legendre_inputs, &k, &alpha_,
            trans, &m, inputs, &k, &beta_,outputs, &m);

        trans += m * k;
        inputs += num_legendre_inputs * k;
        outputs += m * num_legendre_inputs;
        if (mf) m--;
        //if (nf) n--;
        if (kf) k--;
    }
    PROFILEEND("SHT_",0);
}

template<typename Container, typename Mats>
void SHTrans<Container, Mats>::back(const value_type *inputs,
    value_type *work_arr, int n_funs, value_type *outputs,
    value_type *trans, value_type *dft) const
{
    PROFILESTART();
    int num_dft_inputs = n_funs * (p + 1);
    DLT(trans, inputs, outputs, p + 1, 2 * n_funs, p + 1, 0, 0, 1);

    device_.Transpose(outputs, dft_size, num_dft_inputs, work_arr);

    PROFILESTART();
    device_.gemm("T", "N", &dft_size, &num_dft_inputs,
        &dft_size, &alpha_, dft, &dft_size,
        work_arr, &dft_size, &beta_, outputs, &dft_size);
    PROFILEEND("SHT_DFT_",0);

    PROFILEEND("SHT_",0);
}

template<typename Container, typename Mats>
void SHTrans<Container, Mats>::forward(const Container &in, Container &work,
    Container &shc) const
{
    PROFILESTART();
    int n_funs = in.getNumSubFuncs();
    int num_dft_inputs = n_funs * (p + 1);

    PROFILESTART();
    device_.gemm("N", "N", &dft_size, &num_dft_inputs, &dft_size,
        &alpha_, mats_.dft_, &dft_size,in.begin(), &dft_size, &beta_,
        shc.begin(), &dft_size);
    PROFILEEND("SHT_DFT_",0);

    device_.Transpose(shc.begin(), num_dft_inputs, dft_size, work.begin());
    DLT(mats_.dlt_, work.begin(), shc.begin(), p + 1, 2 * n_funs, p + 1, 1, 0, 0);

    PROFILEEND("SHT_",0);
}


template<typename Container, typename Mats>
void SHTrans<Container, Mats>::backward(const Container &shc,
    Container &work, Container &out) const
{
    back(shc.begin(), work.begin(), out.getNumSubFuncs(), out.begin(),
        mats_.dlt_inv_, mats_.dft_inv_);
}


template<typename Container, typename Mats>
void SHTrans<Container, Mats>::backward_du(const Container &shc,
    Container &work, Container &out) const
{
    back(shc.begin(), work.begin(), out.getNumSubFuncs(), out.begin(),
        mats_.dlt_inv_d1_, mats_.dft_inv_);
}

template<typename Container, typename Mats>
void SHTrans<Container, Mats>::backward_d2u(const Container &shc,
    Container &work, Container &out) const
{
    back(shc.begin(), work.begin(), out.getNumSubFuncs(), out.begin(),
        mats_.dlt_inv_d2_, mats_.dft_inv_);
}

template<typename Container, typename Mats>
void SHTrans<Container, Mats>::backward_dv(const Container &shc,
    Container &work, Container &out) const
{
    back(shc.begin(), work.begin(), out.getNumSubFuncs(), out.begin(),
        mats_.dlt_inv_, mats_.dft_inv_d1_);
}

template<typename Container, typename Mats>
void SHTrans<Container, Mats>::backward_d2v(const Container &shc,
    Container &work, Container &out) const
{
    back(shc.begin(), work.begin(), out.getNumSubFuncs(), out.begin(),
        mats_.dlt_inv_, mats_.dft_inv_d2_);
}

template<typename Container, typename Mats>
void SHTrans<Container, Mats>::backward_duv(const Container &shc,
    Container &work, Container &out) const
{
    back(shc.begin(), work.begin(), out.getNumSubFuncs(), out.begin(),
        mats_.dlt_inv_d1_, mats_.dft_inv_d1_);
}

template<typename Container, typename Mats>
void SHTrans<Container, Mats>::lowPassFilter(const Container &in, Container &work,
    Container &shc, Container &out) const
{
    int n_funs = in.getNumSubFuncs();
    forward(in, work, shc);
    ScaleFreq(shc.begin(), n_funs, this->filter_coeff_, shc.begin());
    backward(shc, work, out);
}

template<typename Container, typename Mats>
void SHTrans<Container, Mats>::ScaleFreq(
    const value_type *shc_in, int n_funs,
    const value_type* scaling_coeff, value_type *shc_out) const
{
    int leg_order = p+1;

    device_.ax(scaling_coeff, shc_in, leg_order, n_funs, shc_out);
    scaling_coeff += leg_order;
    shc_in += n_funs * leg_order;
    shc_out += n_funs * leg_order;
    leg_order--;

    // process remaining frequencies except the last cosine
    for (; leg_order>1; leg_order--)
    {
        // first process cosine
        device_.ax(scaling_coeff, shc_in, leg_order, n_funs, shc_out);
        scaling_coeff += leg_order;
        shc_in += n_funs * leg_order;
        shc_out += n_funs * leg_order;

        // then process sine
        device_.ax(scaling_coeff, shc_in, leg_order, n_funs, shc_out);
        scaling_coeff += leg_order;
        shc_in += n_funs * leg_order;
        shc_out += n_funs * leg_order;
    }

    // process last cosine
    device_.ax(scaling_coeff, shc_in, leg_order, n_funs, shc_out);
}

template<typename Container, typename Mats>
void SHTrans<Container, Mats>::FirstDerivatives(const Container &in,
    Container &work, Container &shc, Container &du, Container &dv) const
{
    forward(in, work, shc);
    backward_du(shc, work, du);
    backward_dv(shc, work, dv);
}

template<typename Container, typename Mats>
void SHTrans<Container, Mats>::collectSameOrder(const Container &in,
    Container &out) const
{
    ///@bug this code need to be moved to the device
    const typename Container::value_type *inPtr(NULL);
    typename Container::value_type *outPtr(NULL), *head(out.begin());

    int ns = in.getNumSubFuncs();

    for(int ii=0; ii<= p; ++ii)
    {
        int len = 2*ii + 1 - (ii/p);

        inPtr = in.begin() + ii;
        for(int jj=0; jj<= 2*ii-(ii/p); ++jj)
        {
            outPtr = head + jj;
            int dist = (p + 1 - (jj + 1)/2);
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

template<typename Container, typename Mats>
void SHTrans<Container, Mats>::collectSameFreq(const Container &in,
    Container &out) const
{
    ///@bug this code need to be moved to the device
    const typename Container::value_type *inPtr(NULL), *head(in.begin());
    typename Container::value_type *outPtr(NULL);

    int ns = in.getNumSubFuncs();
    for(int ii=0; ii<= p; ++ii)
    {
        int len = 2*ii + 1 - (ii/p);

        outPtr = out.begin() + ii;
        for(int jj=0; jj<= 2*ii-(ii/p); ++jj)
        {
            inPtr = head + jj;
            int dist = (p + 1 - (jj + 1)/2);
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
