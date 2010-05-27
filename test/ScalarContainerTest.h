/**
 * @file   ScalarsTest.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Mar  2 09:55:22 2010
 *
 * @brief  Tester class for the Scalars class.
 */

#include <iostream>

using namespace std;

template<typename Container>
class ScalarsTest
{
  public:
    bool PerformAll();
};

template<typename Container>
bool ScalarsTest<Container>::PerformAll()
{
    int p = 12;
    int num_funs(3);

    Container sc;
    sc.resize(num_funs, p);
    
    typename Container::iterator it;
    for(it = sc.begin(); it != sc.end(); ++it)
        *it = 10;

    cout<<sc<<endl;
//     Logger::Tic();
//     for(int ii=0; ii<100; ++ii)
//     {
//         Scalar<T,DT> sc(&dev, p, num_funs);
//     }
    
//     cout<<"TIME : "<<Logger::Toc()<<endl;

    return true;
}
    


