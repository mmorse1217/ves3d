/**
 * @file   Vesicle.cc
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @date   Tue Feb  2 14:49:58 2010
 * 
 * @brief  The definition of class Vesicle.
 */


template <typename ScalarType> 
Vesicle<ScalarType>::Vesicle() :
    Surface<ScalarType>(), kappa_(NULL){}

template <typename ScalarType> 
Vesicle<ScalarType>::Vesicle(int p_in, int number_of_vesicles_in) :
    Surface<ScalarType>(p_in, number_of_vesicles_in)
{
    kappa_ = new ScalarType[Surface<ScalarType>::number_of_surfs_];
}

template <typename ScalarType> 
Vesicle<ScalarType>::Vesicle(int p_in, int number_of_vesicles_in, SHVectors<ScalarType>* x_in) :
    Surface<ScalarType>(p_in, number_of_vesicles_in, x_in)
{
    kappa_ = new ScalarType[Surface<ScalarType>::number_of_surfs_];
}

template <typename ScalarType> 
Vesicle<ScalarType>::~Vesicle()
{
    delete [] kappa_;
}







