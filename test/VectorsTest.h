/**
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @revision $Revision$
 * @tags $Tags$
 * @date $Date$
 */

/*
 * Copyright (c) 2014, Abtin Rahimian
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef _VECTORSTEST_H_
#define _VECTORSTEST_H_

#include <string>
#include "Logger.h"
#include "Enums.h"

template<typename V>
class VectorsTest
{
  public:
    typedef typename V::scalars_type S;
    typedef typename S::array_type A;

    int grid_size(int p);
    bool PerformAll();
    bool TestTheDim();
    bool TestNumSubs();
    bool TestSubLength();
    bool TestResize();
    bool TestIterators();
    bool TestReplicate();
    bool TestPointOrder();
};

template<typename V>
int VectorsTest<V>::grid_size(int p)
{ return (p+1)*(2*p) ;}

template<typename V>
bool VectorsTest<V>::PerformAll()
{
    bool test_result;
    test_result = TestTheDim()
        && TestNumSubs()
        && TestSubLength()
        && TestResize()
        && TestIterators()
        && TestReplicate()
        && TestPointOrder()
        ;

    if (test_result)
        COUT(emph<<" *** Vector class passed ***"<<emph);
    else
        COUT(alert<<" *** Vector class failed ***"<<alert);
    return test_result;
}

template<typename V>
bool VectorsTest<V>::TestTheDim()
{
    COUT(". TestTheDim");
    V v(3,4);
    S* s(&v);

    bool res = true;
    ASSERT(res &= V::getTheDim() != S::getTheDim(),
        "differnt from scalar");

    ASSERT(res &= v.getTheDim() == V::getTheDim(),
        "same as static");

    ASSERT(res &= s->getTheDim() == S::getTheDim(),
        "same as static scalar");

    ASSERT(res &= v.getTheDim() != 1,
        "dim higher than 1");

    ASSERT(res &= s->getTheDim() == 1,
        "Scalar dim is 1");

    return res;
}

template<typename V>
bool VectorsTest<V>::TestNumSubs()
{
    COUT(". TestNumSubs");
    int ns(2), p(8);
    V v(ns,p);
    S* s(&v);

    bool res = true;
    ASSERT(res &= v.getNumSubs() == ns,
        "expected num subs");

    ASSERT(res &= v.getNumSubFuncs() == ns*v.getTheDim(),
        "expected num sub functions");

    ASSERT(res &= s->getNumSubs() == ns,
        "expected num subs for scalar");

    ASSERT(res &= s->getNumSubFuncs() == ns*v.getTheDim(),
        "expected num subs functions for scalar");

    return res;
}

template<typename V>
bool VectorsTest<V>::TestSubLength()
{
    COUT(". TestSubLength");
    int ns(3), p(7);
    int sl(grid_size(p));
    V v(ns,p);
    S* s(&v);

    bool res = true;

    ASSERT(res &= s->getSubLength() == v.getTheDim()*sl,
        "same as expected for scalar");

    ASSERT(res &= s->getSubFuncLength() == sl,
        "same as expected for scalar");

    ASSERT(res &= v.getSubFuncLength() ==  sl,
        "same as expected for vector");

    ASSERT(res &= v.getSubLength() == v.getTheDim()*sl,
        "same as expected for vector");

    return res;
}

template<typename V>
bool VectorsTest<V>::TestResize()
{
    COUT(". TestResize");
    int ns(3), p(7);
    int sl(grid_size(p));
    V v;
    S* s(&v);

    bool res = true;
    ASSERT(res &= v.getNumSubs() == 0,
        "same as scalar for zero");

    ASSERT(res &= v.getNumSubs() == 0,
        "expected num subs");

    ASSERT(res &= s->getNumSubs() == 0,
        "expected num subs for scalar");

    ASSERT(res &= v.getSubLength() == s->getSubLength(),
        "same as scalar for zero");

    ASSERT(res &= s->getSubLength() == 0,
        "same as expected for scalar");

    ASSERT(res &= v.getStride() == 0,
        "zeros stride vector");

    ASSERT(res &= s->getStride() == 0,
        "zeros stride scalar");

    ASSERT(res &= v.size() == 0,
        "zero size vector");

    ASSERT(res &= s->size() == 0,
        "zero size scalar");

    ASSERT(res &= v.req_arr_size() == 0,
        "zero req_arr_size vector");

    ASSERT(res &= s->req_arr_size() == 0,
        "zero req_arr_size scalar");

    //resize vector
    v.resize(ns,p);

    ASSERT(res &= v.getNumSubs() == ns,
        "expected num subs "<<ns<<" vs. "<<v.getNumSubs());

    ASSERT(res &= s->getNumSubs() == ns,
        "expected num subs for scalar");

    ASSERT(res &= s->getNumSubFuncs() == v.getTheDim()*ns,
        "expected num subs for scalar");

    ASSERT(res &= s->getSubLength() == v.getTheDim()*sl,
        "same as expected for scalar");

    ASSERT(res &= s->getSubFuncLength() == sl,
        "same as expected for scalar");

    ASSERT(res &= v.getSubLength() == v.getTheDim() * sl,
        "same as expected for vector");

    ASSERT(res &= v.getStride() == sl,
        "expected stride vector");

    ASSERT(res &= s->getStride() == sl,
        "expected stride scalar");

    ASSERT(res &= v.size() == v.getTheDim()*sl*ns,
        "expected size vector");

    ASSERT(res &= s->size() == v.getTheDim()*ns*sl,
        "expected size scalar");

    ASSERT(res &= v.req_arr_size() == v.getTheDim()*ns*sl,
        "expected req_arr_size vector");

    ASSERT(res &= s->req_arr_size() == v.getTheDim()*ns*sl,
        "expected req_arr_size scalar");

    //resize scalar
    ns=6;
    p=7;
    sl=grid_size(p);

    s->resize(ns,p);
    ASSERT(res &= v.getNumSubs() == ns,
        "expected num subs");

    ASSERT(res &= v.getNumSubFuncs() == ns*v.getTheDim(),
        "expected num subs");

    ASSERT(res &= s->getNumSubs() == ns,
        "expected num subs for scalar");

    ASSERT(res &= s->getNumSubFuncs() == ns*v.getTheDim(),
        "expected num subs for scalar");

    ASSERT(res &= s->getSubLength() == v.getTheDim()*sl,
        "same as expected for scalar");

    ASSERT(res &= v.getSubLength() == v.getTheDim() * sl,
        "same as expected for vector");

    ASSERT(res &= v.getStride() == sl,
        "expected stride vector");

    ASSERT(res &= s->getStride() == sl,
        "expected stride scalar");

    ASSERT(res &= v.size() == sl*ns*v.getTheDim(),
        "expected size vector");

    ASSERT(res &= s->size() == ns*sl*v.getTheDim(),
        "expected size scalar");

    ASSERT(res &= v.req_arr_size() == ns*sl*v.getTheDim(),
        "expected req_arr_size vector");

    ASSERT(res &= s->req_arr_size() == ns*sl*v.getTheDim(),
        "expected req_arr_size scalar");

    return res;
}

template<typename V>
bool VectorsTest<V>::TestIterators()
{
    COUT(". TestIterators");
    int ns(3), p(7);
    int sl(grid_size(p));
    V v(ns,p);
    S* s(&v);
    A* a(&v);

    bool res = true;
    ASSERT(res &= v.getSubN_begin(0) == a->begin(),
        "");

    ASSERT(res &= s->getSubN_begin(0) == a->begin(),
        "");

    ASSERT(res &= v.getSubFuncN_begin(0) == a->begin(),
        "");

    ASSERT(res &= s->getSubFuncN_begin(0) == a->begin(),
        "");

    ASSERT(res &= v.getSubN_end(0) == a->begin()+v.getTheDim()*sl,
        "");

    ASSERT(res &= s->getSubN_end(0) == a->begin()+v.getTheDim()*sl,
        "");

    ASSERT(res &= v.getSubFuncN_end(0) == a->begin()+sl,
        "");

    ASSERT(res &= s->getSubFuncN_end(0) == a->begin()+sl,
        "");

    ASSERT(res &= v.getSubN_begin(1) == v.getSubN_end(0),
        "");

    ASSERT(res &= s->getSubN_begin(1) == v.getSubN_end(0),
        "");

    ASSERT(res &= v.getSubFuncN_begin(1) == s->getSubFuncN_end(0),
        "");

    ASSERT(res &= s->getSubFuncN_begin(1) == v.getSubFuncN_end(0),
        "");

    ASSERT(res &= v.getSubN_end(1) == a->begin()+2*v.getTheDim()*sl,
        "");

    ASSERT(res &= s->getSubN_end(1) == a->begin()+2*v.getTheDim()*sl,
        "");

    ASSERT(res &= v.getSubFuncN_end(1) == a->begin()+2*sl,
        "");

    ASSERT(res &= s->getSubFuncN_end(1) == a->begin()+2*sl,
        "");

    return res;
}

template<typename V>
bool VectorsTest<V>::TestReplicate()
{
    COUT(". TestReplicate");
    int ns(3), p(7);
    int sl(grid_size(p));
    V v(ns,p),vr,vsr;
    S s(ns,p),sr,svr;

    bool res = true;
    // array->array
    ASSERT(res &= vr.size() == 0,"");
    vr.replicate(v);
    ASSERT(res &= vr.getShOrder() == v.getShOrder(),"");
    ASSERT(res &= vr.getNumSubs() == v.getNumSubs(),"");
    ASSERT(res &= vr.getGridDim() == v.getGridDim(),"");
    ASSERT(res &= vr.getStride() == v.getStride(),"");

    // scalar->scalar
    ASSERT(res &= sr.size() == 0,"");
    sr.replicate(s);
    ASSERT(res &= sr.getShOrder() == s.getShOrder(),"");
    ASSERT(res &= sr.getNumSubs() == s.getNumSubs(),"");
    ASSERT(res &= sr.getGridDim() == s.getGridDim(),"");
    ASSERT(res &= sr.getStride() == s.getStride(),"");

    // array->scalar
    ASSERT(res &= svr.size() == 0,"");
    svr.replicate(v);
    ASSERT(res &= svr.getShOrder() == v.getShOrder(),"");
    ASSERT(res &= svr.getNumSubs() == v.getNumSubs(),"");
    ASSERT(res &= svr.getGridDim() == v.getGridDim(),"");
    ASSERT(res &= svr.getStride() == v.getStride(),"");

    // scalar->array
    ASSERT(res &= vsr.size() == 0,"");
    vsr.S::replicate(s);
    ASSERT(res &= vsr.getShOrder() == s.getShOrder(),"");
    ASSERT(res &= vsr.getNumSubs() == s.getNumSubs(),
        "Check correct number of subs vc:"
        <<vsr.getNumSubs()<<" sc:"<<s.getNumSubs()
           );

    ASSERT(res &= vsr.getGridDim() == s.getGridDim(),"");
    ASSERT(res &= vsr.getStride() == s.getStride(),"");

    return res;
}

template<typename V>
bool VectorsTest<V>::TestPointOrder()
{
    COUT(". TestPointOrder");
    V v(3,4),vr;

    bool res = true;
    ASSERT(res &= v.getPointOrder() == AxisMajor,
        "default coordinate order");
    vr.replicate(v);
    ASSERT(res &= vr.getPointOrder() == AxisMajor,"");

    v.setPointOrder(PointMajor);
    ASSERT(res &= v.getPointOrder() == PointMajor,"");
    ASSERT(res &= vr.getPointOrder() == AxisMajor,"");
    vr.replicate(v);
    ASSERT(res &= vr.getPointOrder() == AxisMajor,
        "point order not replicated");

    v.setPointOrder(PointMajor);
    ASSERT(res &= v.getPointOrder() == PointMajor,"");

    return res;
}

#endif
