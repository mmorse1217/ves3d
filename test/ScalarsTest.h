/**
 * @file   ScalarContaiersTest.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Mar  2 09:55:22 2010
 *
 * @brief  Tester class for the Scalars class.
 */

#include "Logger.h"
#include <cstring>
#include <iostream>
#include <typeinfo> //for typeid
#include <cstdlib>  //for abs

template<typename Container>
class ScalarsTest
{
  public:
    bool PerformAll();
    bool TestResize();
    bool TestReplicate();
};

template<typename Container>
bool ScalarsTest<Container>::PerformAll()
{
    bool test_result;
    test_result = TestResize()
        && TestReplicate();

    if (test_result){
        COUT(emph<<" *** The container "
            <<typeid(Container).name()<<" Passed ***"
            <<emph);
    } else {
        COUT(alert<<" *** The container "
            <<typeid(Container).name()<<" Failed ***"
            <<alert);
    }

    return test_result;
}

template<typename Container>
bool ScalarsTest<Container>::TestResize()
{
    COUT(". Resize");
    int p(13), nsubs(3);

    Container sc;
    sc.resize(nsubs, p);
    ASSERT(sc.getShOrder()==p,"");
    ASSERT(sc.getNumSubs()==nsubs,"");
    ASSERT(sc.size()>0,"");

    return true;
}

template<typename Container>
bool ScalarsTest<Container>::TestReplicate()
{
    COUT(". Replicate");
    int p(8), nsubs(6);
    Container sc(nsubs, p);
    Container sc_cpy;
    sc_cpy.replicate(sc);
    ASSERT(sc_cpy.getNumSubs()==nsubs,"");
    ASSERT(sc_cpy.getSubLength()==sc.getSubLength(),"");
    ASSERT(sc_cpy.getSubFuncLength()==sc.getSubFuncLength(),"");
    ASSERT(sc_cpy.size()==sc.size(),"");

    typedef typename Container::value_type T;
    size_t sz = sc.size();
    T* buffer = new T[sc.size()];

    for(size_t ii=0; ii<sz; ++ii)
        buffer[ii] = ii;

    sc.getDevice().Memcpy(sc.begin(),
                          buffer,
                          sz * sizeof(T),
                          Container::device_type::MemcpyHostToDevice);

    memset(buffer, 0, sz * sizeof(T));

    sc.getDevice().Memcpy(buffer,
                          sc.begin(),
                          sz * sizeof(T),
                          Container::device_type::MemcpyDeviceToHost);

    T err = 0;
    for(size_t ii=0; ii<sz; ++ii)
        err = (err > abs(buffer[ii]-ii)) ? err : abs(buffer[ii]-ii);
    bool res(err==0);

    delete[] buffer;

    ASSERT(res,"Copying form to and from the container");

    return res;
}
