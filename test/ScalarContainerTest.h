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

    std::string res_print = (test_result) ? " Passed" : " Failed";

    COUT("\n\n ====================================================\n"
        <<"  The container "
        <<typeid(Container).name()<<res_print<<std::endl
        <<" ===================================================="
        <<std::endl);

    return test_result;
}

template<typename Container>
bool ScalarsTest<Container>::TestResize()
{
    int p(13), nsubs(3);

    Container sc;
    COUT(sc);
    COUT(" Resizing to "<<p<<" and "<<nsubs<<std::endl);
    sc.resize(nsubs, p);
    COUT(sc);

    bool res(sc.getNumSubs() == nsubs &&
             sc.getShOrder() == p
             );

    return res;
}

template<typename Container>
bool ScalarsTest<Container>::TestReplicate()
{
    int p(8), nsubs(6);

    Container sc(nsubs, p);
    Container sc_cpy;
    COUT(sc_cpy);
    COUT(" Replicating :"<<std::endl);
    sc_cpy.replicate(sc);
    COUT(sc_cpy);

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

    COUT(" Copying form to and from the container: "
         <<((res) ? "Pass" : "Fail")
         <<std::endl);

    return res;
}
