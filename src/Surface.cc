/**
 * @file   Surface.cc
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @date   Tue Feb  2 13:37:21 2010
 * 
 * @brief  The implementation of the surface class.
 */

template <typename ScalarType> 
Surface<ScalarType>::Surface() :
    p_(0), number_of_surfs_(0){}

template <typename ScalarType> 
Surface<ScalarType>::Surface(int p_in, int number_of_surfs_in) :
    p_(p_in), number_of_surfs_(number_of_surfs_in), x_(p_,number_of_surfs_),
    normal_(p_,number_of_surfs_), h_(p_,number_of_surfs_), w_(p_,number_of_surfs_),
    k_(p_,number_of_surfs_), cu_(p_,number_of_surfs_), cv_(p_,number_of_surfs_){}

template <typename ScalarType> 
Surface<ScalarType>::Surface(int p_in, int number_of_surfs_in, 
    const SHVectors<ScalarType> &x_in) :
    p_(p_in), number_of_surfs_(number_of_surfs_in), x_(p_,number_of_surfs_),
    normal_(p_,number_of_surfs_), h_(p_,number_of_surfs_), w_(p_,number_of_surfs_),
    k_(p_,number_of_surfs_), cu_(p_,number_of_surfs_), cv_(p_,number_of_surfs_)
{
    SetX(x_in);
}

template <typename ScalarType> 
void Surface<ScalarType>::SetX(const SHVectors<ScalarType> &x_in)
{
    x_.SetData(x_in.data_);
    
    //update other vectors.
}
