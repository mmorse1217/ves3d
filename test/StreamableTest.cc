/**
 * @file
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @revision $Revision$
 * @tags $Tags$
 * @date $Date$
 *
 * @brief unit test
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

#include "Streamable.h"
#include "Logger.h"
#include "TestTools.h"
#include <sstream>

class Foo  : public Streamable{
  public:
    double x,y;
    int r;
    size_t n;
    float *z;

    Foo(double x, double y, int r, int n) : x(x), y(y), r(r), n(n), z(NULL)
    {
	z = new float[n];
	for (int i=0; i<n; ++i) z[i] = i;
    }

    Foo(std::istream &is, Streamable::Format format) :
	z(NULL)
    {
	unpack(is, format);
    }

    ~Foo(){delete z;}

    Error_t pack(std::ostream &os, Streamable::Format format) const
    {
	ASSERT(format==Streamable::ASCII, "BIN is not supported yet");
    	os<<"FOO\n"<<name_<<"\n"<<x<<"\n"<<y<<"\n"<<r<<"\n"<<n<<"\n";
	pack_array(os, format, z, n);
        os<<"\n/FOO\n";
	return ErrorEvent::Success;
    }

    Error_t unpack(std::istream &is, Streamable::Format format)
    {
	ASSERT(format==Streamable::ASCII, "BIN is not supported yet");
	std::string s;
	is>>s;
	ASSERT(s=="FOO", "Bad input string (missing header).");
	is>>name_>>x>>y>>r>>n;

	delete[] z;
	z = new float[n];
	unpack_array(is, format, z, n);
	is>>s;
	ASSERT(s=="/FOO", "Bad input string (missing footer).");

	return ErrorEvent::Success;
    }
};

int main(int argc, char** argv){
    VES3D_INITIALIZE(&argc,&argv,NULL,NULL);
    Foo c(1.0,10.12345678901234567e12,3,8); c.set_name("c");
    Foo d(2.0,3.0,4,3); d.set_name("d");

    std::stringstream ss;
    ss<<std::scientific<<std::setprecision(10);

    c.pack(ss, Streamable::ASCII);
    d.unpack(ss, Streamable::ASCII);

    testtools::AssertAlmostEqual(c.x, d.x, 1e-10, "bad (un)packed x");
    testtools::AssertAlmostEqual(c.y, d.y, 1e-10, "bad (un)packed y");
    testtools::AssertEqual(c.r, d.r, "bad (un)packed r");
    testtools::AssertEqual(c.n, d.n, "bad (un)packed n");

    for (int i=0; i<d.n; ++i)
    	testtools::AssertAlmostEqual<float>(c.z[i], d.z[i], 1e-7, "bad (un)packed z");

    std::stringstream s2(ss.str());
    Foo e(s2, Streamable::ASCII);
    testtools::AssertAlmostEqual(c.x, e.x, 1e-10, "bad (un)packed x");
    testtools::AssertAlmostEqual(c.y, e.y, 1e-10, "bad (un)packed y");
    testtools::AssertEqual(c.r, e.r, "bad (un)packed r");
    testtools::AssertEqual(c.n, e.n, "bad (un)packed n");

    for (int i=0; i<e.n; ++i)
    	testtools::AssertAlmostEqual<float>(c.z[i], e.z[i], 1e-7, "bad (un)packed z");

    COUT("e:\n"<<std::scientific<<std::setprecision(16)<<&e<<"-------\nc:\n"<<&c);
    VES3D_FINALIZE();
}
