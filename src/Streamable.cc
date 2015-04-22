/**
 * @file
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @revision $Revision$
 * @tags $Tags$
 * @date $Date$
 *
 * @brief
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

Error_t Streamable::set_name(const std::string n)
{
    name_=n;
    return ErrorEvent::Success;
}

std::string Streamable::name()
{
    return name_;
}

template<typename D>
Error_t Streamable::factory(std::istream &is, Format format, D **obj)
{
    ASSERT(*obj == NULL, "initial pointer should be null");
    *obj = new D(is, format);
    return ErrorEvent::Success;
}

template<typename T>
Error_t Streamable::pack_array(std::ostream &os, Format format, const T* arr, size_t n, const char *sep) const
{
    ASSERT(sep=="\n" || sep=="\t" || sep==" ", "unsupported separator");
    if(n>0) {
	for (size_t ii=0; ii+1<n; ++ii) os<<arr[ii]<<sep;
	os<<arr[n-1]; //to skip the separator for the last one
    }
    return ErrorEvent::Success;
}

template<typename T>
Error_t Streamable::unpack_array(std::istream &is, Format format, T* arr, size_t &n, const char* sep) const
{
    ASSERT(sep=="\n" || sep=="\t" || sep==" ", "unsupported separator");
    size_t N(n);
    for (n=0; n<N; ++n){
	if(!is.good()) return ErrorEvent::IOBadStream;
	is>>arr[n];
    }
    return ErrorEvent::Success;
}

std::ostream& operator<<(std::ostream& os, const Streamable *o)
{
    o->pack(os, Streamable::ASCII);
    return os;
}
