/**
 * @file   SHVectors.cc
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Fri Jan 29 09:36:58 2010
 * 
 * @brief  The implementation of SHVectors.
 */

// Constructors
template<typename ScalarType> 
SHVectors<ScalarType>::SHVectors() :
    number_of_vectors_(0){}

template<typename ScalarType> 
SHVectors<ScalarType>::SHVectors(int p_in, int num_vecs_in) :
    SHScalars<ScalarType>( p_in, 3*num_vecs_in) , number_of_vectors_(num_vecs_in){}

template<typename ScalarType> 
SHVectors<ScalarType>::SHVectors(int p_in, int num_vecs_in, const ScalarType *vec_data_in) :
    SHScalars<ScalarType>( p_in, 3*num_vecs_in, vec_data_in) , number_of_vectors_(num_vecs_in){}

// Utility functions
template<typename ScalarType> 
int SHVectors<ScalarType>::GetVecLength() const
{
    return(3*SHScalars<ScalarType>::GetFunLength());
}

template<typename ScalarType> 
void SHVectors<ScalarType>::Resize(int p_in, int number_of_vectors_in)
{
    assert(SHScalars<ScalarType>::data_ == NULL);
    SHScalars<ScalarType>::Resize(p_in, 3*number_of_vectors_in);
    number_of_vectors_ = number_of_vectors_in;
}

template<typename ScalarType> 
void DotProduct(const SHVectors<ScalarType>& a_in, 
    const SHVectors<ScalarType>& b_in, SHScalars<ScalarType>& aDb_out)
{
    assert(a_in.GetDataLength() == b_in.GetDataLength());
    assert(a_in.GetDataLength() <= 3*aDb_out.GetDataLength());
    
    int funLen = a_in.GetFunLength();
    int numVec = a_in.number_of_vectors_;
    int idx;
    const ScalarType *a_data;
    const ScalarType *b_data;
    ScalarType *dot = aDb_out.data_;
    
//ScalarType *dot = new ScalarType[numVec*funLen];
    
    for(int ii=0;ii<numVec;++ii)
    {
        a_data = a_in.GetFunctionAt(3*ii);
        b_data = b_in.GetFunctionAt(3*ii);

        for(int jj=0;jj<funLen;++jj)
        {
            idx = ii*funLen +jj;
            
            dot[idx] = a_data[jj         ]*b_data[jj         ];
            dot[idx]+= a_data[jj+  funLen]*b_data[jj+  funLen];
            dot[idx]+= a_data[jj+2*funLen]*b_data[jj+2*funLen];
        }
    }
}

template<typename ScalarType> 
void CrossProduct(const SHVectors<ScalarType>& a_in, 
    const SHVectors<ScalarType>& b_in, SHVectors<ScalarType>& aCb_out)
{
    assert(a_in.GetDataLength() == b_in.GetDataLength() &&
        b_in.GetDataLength() == aCb_out.GetDataLength());
    
    int funLen = a_in.GetFunLength();
    int numVec = a_in.number_of_vectors_;
    int idx;
    const ScalarType *a_data;
    const ScalarType *b_data;
    ScalarType *cross = new ScalarType[aCb_out.GetDataLength()];//to have output the same as inputs, we need a buffer.
    
    for(int ii=0;ii<numVec;++ii)
    {
        a_data = a_in.GetFunctionAt(3*ii);
        b_data = b_in.GetFunctionAt(3*ii);

        for(int jj=0;jj<funLen;++jj)
        {
            idx = 3*ii*funLen +jj;
            
            cross[idx         ] = a_data[jj+  funLen]*b_data[jj+2*funLen] - a_data[jj+2*funLen]*b_data[jj+  funLen];
            cross[idx+  funLen] = a_data[jj+2*funLen]*b_data[jj         ] - a_data[jj         ]*b_data[jj+2*funLen];
            cross[idx+2*funLen] = a_data[jj         ]*b_data[jj+  funLen] - a_data[jj+  funLen]*b_data[jj         ];
        }
    }
    aCb_out.SetData(cross);    
    delete []cross;
}

template<typename ScalarType> 
void AxPy(const SHScalars<ScalarType>& a_in, 
    const SHVectors<ScalarType>& x_in, const SHVectors<ScalarType>& y_in, 
        SHVectors<ScalarType>& c_out)
{

    int num_vecs = x_in.number_of_vectors_;
    int fun_len = x_in.GetFunLength();
    int vec_len = 3*fun_len;
    int v_idx, s_idx;

    for(int ii=0;ii<num_vecs;++ii)
    {
        s_idx = ii*fun_len;
        v_idx = ii*vec_len;
        
        for(int jj=0;jj<fun_len;++jj)
        {
            c_out.data_[v_idx+          jj] = a_in.data_[s_idx+jj]*x_in.data_[v_idx          +jj] + y_in.data_[v_idx+          jj];
            c_out.data_[v_idx+  fun_len+jj] = a_in.data_[s_idx+jj]*x_in.data_[v_idx+  fun_len+jj] + y_in.data_[v_idx+  fun_len+jj];
            c_out.data_[v_idx+2*fun_len+jj] = a_in.data_[s_idx+jj]*x_in.data_[v_idx+2*fun_len+jj] + y_in.data_[v_idx+2*fun_len+jj];
        }
    }
}

template<typename ScalarType> 
void AxPy(const SHScalars<ScalarType>& a_in, const SHVectors<ScalarType>& x_in, 
    ScalarType y_in, SHVectors<ScalarType>& c_out)
{

    int num_vecs = x_in.number_of_vectors_;
    int fun_len = x_in.GetFunLength();
    int vec_len = 3*fun_len;
    int v_idx, s_idx;

    for(int ii=0;ii<num_vecs;++ii)
    {
        s_idx = ii*fun_len;
        v_idx = ii*vec_len;
        
        for(int jj=0;jj<fun_len;++jj)
        {
            c_out.data_[v_idx+          jj] = a_in.data_[s_idx+jj]*x_in.data_[v_idx          +jj] + y_in;
            c_out.data_[v_idx+  fun_len+jj] = a_in.data_[s_idx+jj]*x_in.data_[v_idx+  fun_len+jj] + y_in;
            c_out.data_[v_idx+2*fun_len+jj] = a_in.data_[s_idx+jj]*x_in.data_[v_idx+2*fun_len+jj] + y_in;
        }
    }
}
template<typename ScalarType> 
void xDy(const SHVectors<ScalarType>& x_in, const SHScalars<ScalarType>& y_in, 
        SHVectors<ScalarType>& c_out)
{
    int num_vecs = x_in.number_of_vectors_;
    int fun_len = x_in.GetFunLength();
    int vec_len = 3*fun_len;
    int v_idx, s_idx;

    for(int ii=0;ii<num_vecs;++ii)
    {
        s_idx = ii*fun_len;
        v_idx = ii*vec_len;
        
        for(int jj=0;jj<fun_len;++jj)
        {
            c_out.data_[v_idx+          jj] = x_in.data_[v_idx          +jj]/y_in.data_[s_idx+jj];
            c_out.data_[v_idx+  fun_len+jj] = x_in.data_[v_idx+  fun_len+jj]/y_in.data_[s_idx+jj];
            c_out.data_[v_idx+2*fun_len+jj] = x_in.data_[v_idx+2*fun_len+jj]/y_in.data_[s_idx+jj];
        }
    }
}
