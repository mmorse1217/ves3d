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
#include <sstream>
#include <typeinfo> //for typeid
#include <cstdlib>  //for abs

template<typename Container>
class ScalarsTest
{
  public:
    bool PerformAll();
    bool TestResize();
    bool TestReplicate();
    bool TestStream();
};

template<typename Container>
bool ScalarsTest<Container>::PerformAll()
{
    bool test_result;
    test_result = TestResize() &&
        TestReplicate() &&
	TestStream();

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

template<typename Container>
bool ScalarsTest<Container>::TestStream()
{
    COUT(". Stream");
    int p(3), nsubs(5);

    Container sc(nsubs, p, std::make_pair<int,int>(p+5, 3*p));
    sc.set_name("XYZ");
    typedef typename Container::value_type T;
    T* buffer = new T[sc.size()];

    for(size_t ii=0; ii<sc.size(); ++ii)
        buffer[ii] = ii;

    sc.getDevice().Memcpy(sc.begin(),
	buffer,
	sc.size() * sizeof(T),
	Container::device_type::MemcpyHostToDevice);
    memset(buffer, 0, sc.size() * sizeof(T));

    std::stringstream s1;
    sc.pack(s1, Container::Streamable::ASCII);
    Container cc;
    cc.unpack(s1, Container::Streamable::ASCII);

    ASSERT(cc.getShOrder()	== sc.getShOrder()	, "bad streaming");
    ASSERT(cc.getGridDim()	== sc.getGridDim()	, "bad streaming");
    ASSERT(cc.getStride()	== sc.getStride()	, "bad streaming");
    ASSERT(cc.getNumSubFuncs()	== sc.getNumSubFuncs()	, "bad streaming");
    ASSERT(cc.name()		== sc.name()		, "bad streaming");

    cc.getDevice().Memcpy(buffer,
	cc.begin(),
	cc.size() * sizeof(T),
	Container::device_type::MemcpyDeviceToHost);

    T err = 0;
    for(size_t ii=0; ii<sc.size(); ++ii)
        err = (err > abs(buffer[ii]-ii)) ? err : abs(buffer[ii]-ii);
    ASSERT(err==0, "bad streaming for the data");

    //constructor
    std::stringstream s2,s3;
    sc.pack(s2, Container::Streamable::ASCII);
    Container cx(s2,Container::Streamable::ASCII);
    cx.pack(s3, Container::Streamable::ASCII);
    ASSERT(s2.str()==s3.str(), "bad constructor from stream");

    return true;
}
