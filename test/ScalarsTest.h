/**
 * @file   ScalarsTest.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Mar  2 09:55:22 2010
 *
 * @brief  Tester class for the Scalars class.
 */

#include "Scalars.h"
#include "HelperFuns.h"

using namespace std;

template<typename T, enum DeviceType DT>
class ScalarsTest
{
  public:
    bool PerformAll();
};

template<typename T, enum DeviceType DT>
bool ScalarsTest<T,DT>::PerformAll()
{
    Device<DT> dev;
    int p = 12;
    int num_funs(3);

//     Scalars<T,DT> sc(&dev, p, num_funs);
    
//     typename Scalars<T,DT>::iterator it;
//     for(it = sc.begin(); it != sc.end(); ++it)
//         *it = 10;

//     cout<<sc<<endl;

    Logger::Tic();
    for(int ii=0; ii<100; ++ii)
    {
        Scalar<T,DT> sc(&dev, p, num_funs);
    }
    
    cout<<"TIME : "<<Logger::Toc()<<endl;

    return true;
}
    


