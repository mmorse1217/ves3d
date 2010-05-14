/**
 * @file   VectorsTest.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Mar  2 20:11:43 2010
 * 
 * @brief  The tester for the Vectors class.
 */


#include "Vectors.h"
#include "HelperFuns.h"

using namespace std;

template<typename T, enum DeviceType DT>
class VectorsTest
{
  public:
    bool PerformAll();
};

template<typename T, enum DeviceType DT>
bool VectorsTest<T,DT>::PerformAll()
{
    Device<DT> dev;
    int p(12);
    int num_vecs(3);
    Vectors<T,DT> vec(&dev, p, num_vecs);
    Scalars<T,DT> sc(&dev, p, num_vecs);

    typename Scalars<T,DT>::iterator it;
    for(it = vec.begin(); it != vec.end(); ++it)
        *it = 2;

    DotProduct(vec,vec,sc);
    cout<<sc<<endl;
    
    cout<<vec<<endl;
    return true;
}

