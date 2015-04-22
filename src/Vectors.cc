/**
 * @file   Vector.cc
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

template <typename T, typename DT, const DT &DEVICE>
Vectors<T, DT, DEVICE>::Vectors(size_t num_sub_vecs, int sh_order,
    std::pair<int, int> grid_dim, CoordinateOrder po) :
    Scalars<T, DT, DEVICE>(num_sub_vecs * this->vector_dim_,
	sh_order,
	grid_dim),
    point_order_(po)
{}

template <typename T, typename DT, const DT &DEVICE>
Vectors<T, DT, DEVICE>::Vectors(std::istream &is, Streamable::Format format)
{
    Vectors<T, DT, DEVICE>::unpack(is, format);
}

template <typename T, typename DT, const DT &DEVICE>
int Vectors<T, DT, DEVICE>::getTheDim()
{
    return(vector_dim_);
}

template <typename T, typename DT, const DT &DEVICE>
size_t Vectors<T, DT, DEVICE>::getNumSubs() const
{
    return(scalars_type::num_sub_funcs_ / vector_dim_);
}

template <typename T, typename DT, const DT &DEVICE>
void Vectors<T, DT, DEVICE>::setNumSubs(size_t num_subs)
{
    scalars_type::num_sub_funcs_ = num_subs*vector_dim_;
}

template <typename T, typename DT, const DT &DEVICE>
size_t Vectors<T, DT, DEVICE>::getSubLength() const
{
    return( getTheDim() * this->getSubFuncLength());
}

template <typename T, typename DT, const DT &DEVICE>
void Vectors<T, DT, DEVICE>::setPointOrder(enum CoordinateOrder new_order)
{
    ASSERT(new_order != UnknownOrder, "Can't set the order to unknown");
    point_order_ = new_order;
}

template <typename T, typename DT, const DT &DEVICE>
enum CoordinateOrder Vectors<T, DT, DEVICE>::getPointOrder() const
{
    return(point_order_);
}

template<typename T, typename DT, const DT &DEVICE>
Error_t Vectors<T, DT, DEVICE>::pack(std::ostream &os, Streamable::Format format) const
{
    ASSERT(format==Streamable::ASCII, "BIN is not supported yet");

    os<<"VECTORS\n";
    os<<"name: "<<Streamable::name_<<"\n";
    os<<"CoordinateOrder: "<<point_order_<<"\n";
    scalars_type::pack(os, format);
    os<<"/VECTORS\n";

    return ErrorEvent::Success;
}

template<typename T, typename DT, const DT &DEVICE>
Error_t Vectors<T, DT, DEVICE>::unpack(std::istream &is, Streamable::Format format)
{
    ASSERT(format==Streamable::ASCII, "BIN is not supported yet");
    std::string s, key;
    is>>s;
    ASSERT(s=="VECTORS", "Bad input string (missing header).");

    is>>key>>Streamable::name_;
    ASSERT(key=="name:", "bad key");
    is>>key>>s;
    ASSERT(key=="CoordinateOrder:", "bad key");
    setPointOrder(EnumifyCoordinateOrder(s.c_str()));
    scalars_type::unpack(is, format);
    is>>s;
    ASSERT(s=="/VECTORS", "Bad input string (missing footer).");

    return ErrorEvent::Success;
}
