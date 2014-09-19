/**
 * @file
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @revision $Revision$
 * @tags $Tags$
 * @date $Date$
 *
 * @brief unit test for home-brewed tr1
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

#include "tr1.h"
#include "Logger.h"

int main(int argc, char** argv){

    //Testing shared_ptr
    COUT(" - Default constructor"<<std::endl);
    tr1::shared_ptr<int> p0,p1; // default constructor
    ASSERT(!p0,"NULL pointer");
    ASSERT(!p1,"NULL pointer");

    int *i;
    int k(3);
    {
        i = new int(5);
        int *j(new int(8));
        COUT(" - Copy constructor from bare pointer"<<std::endl);
        tr1::shared_ptr<int> p2(i);
        ASSERT(bool(p2),"non-NULL pointer");
        ASSERT(p2.get()==i,"underlying pointer");
        ASSERT(*p2==*i,"underlying pointer");

        COUT(" - Copy constructor from shared_ptr"<<std::endl);
        tr1::shared_ptr<int> p3(p2);
        ASSERT(bool(p3),"non-NULL pointer");
        ASSERT(p3.get()==i,"underlying pointer");
        ASSERT(*p3==*i,"underlying pointer");

        COUT(" - Assignment from shared_ptr"<<std::endl);
        p0 = p2;
        ASSERT(bool(p0),"non-NULL pointer");
        ASSERT(p0.get()==i,"underlying pointer");
        ASSERT(*p0==*i,"underlying pointer");

        COUT(" - Assignment from bare pointer"<<std::endl);
        p1 = j;
        ASSERT(bool(p1),"non-NULL pointer");
        ASSERT(p1.get()==j,"underlying pointer");
        ASSERT(*p1==*j,"underlying pointer");

        *i = k;
        ASSERT(*p0==k,"underlying pointer");
        ASSERT(*p1!=k,"underlying pointer");
        ASSERT(*p1==*j,"underlying pointer");
        ASSERT(*p2==k,"underlying pointer");
        ASSERT(*p3==k,"underlying pointer");

        *j = k;
        ASSERT(*p1==k,"underlying pointer");
    }

    ASSERT(*p0==k,"underlying pointer");
    ASSERT(*p1==k,"underlying pointer");

    k   = 12;
    *p0 = k;
    ASSERT(*p0==k,"underlying pointer");
    ASSERT(*i==k,"underlying pointer");

    int l(20);
    for(int j(0);j<l;++j)
        tr1::shared_ptr<int> p(p1);

    tr1::shared_ptr<int> sp[l];
    for(int j(0);j<l;++j)
        sp[j]=p1;
}
