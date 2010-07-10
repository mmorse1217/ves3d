/**
 * @file   ScalarContaiersTest.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Mar  2 09:55:22 2010
 *
 * @brief  Tester class for the Scalars class.
 */

#include "Logger.h"
#include <iostream>
#include "Scalars.h"
#include "cstring"
#include <typeinfo>

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
    int p(13), nsubs(3);

    Container sc;
    COUT(sc);
    COUT(" Resizing to "<<p<<" and "<<nsubs<<endl);
    sc.resize(nsubs, p);
    COUT(sc);

    Container sc_cpy;
    COUT(sc_cpy);
    COUT(" Replicating :"<<endl);
    sc_cpy.replicate(sc);
    COUT(sc_cpy);
   
    typedef typename Container::value_type T;
    size_t sz = sc.size();
    T* buffer = new T[sc.size()];
    
    for(size_t ii=0; ii<sz; ++ii)
        buffer[ii] = ii;

    sc.getDevice().Memcpy(sc.begin(), buffer, sz * sizeof(T),
        MemcpyHostToDevice);

    memset(buffer, 0, sz * sizeof(T));
    
    sc.getDevice().Memcpy(buffer, sc.begin(), sz * sizeof(T),
        MemcpyDeviceToHost);

    T err = 0;
    for(size_t ii=0; ii<sz; ++ii)
        err = (err > abs(buffer[ii]-ii)) ? err : abs(buffer[ii]-ii);
    bool res(err==0);

    delete[] buffer;
    
    COUT(" Copying form to and from the container: "<<((res) ? "Pass" : "Fail")<<endl);

    if(res)
        COUT("\n\n ====================================================================\n"
            <<"  The container "<<typeid(Container).name()<<" works fine.\n"
            <<" ===================================================================="<<endl);

    sleep(.5);
    return res;
}
    


