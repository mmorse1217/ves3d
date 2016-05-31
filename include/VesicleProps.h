/**
 * @file
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @revision $Revision$
 * @tags $Tags$
 * @date $Date$
 *
 * @brief Container for vesicle properties
 */

/*
 * Copyright (c) 2016, Abtin Rahimian
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

#ifndef _VESICLEPROPS_H_
#define _VESICLEPROPS_H_

#include "Logger.h"
#include "Error.h"
#include "Streamable.h"

template <typename T>
struct VesicleProperties : public Streamable
{
    typedef T container_type;
    typedef typename T::value_type  value_type;
    typedef typename T::size_type   size_type;
    typedef typename T::device_type device_type;

    VesicleProperties();
    Error_t update();
    T* getPropIdx(int i); /* for loading convenience */

    Error_t pack(std::ostream &os, Streamable::Format format) const;
    Error_t unpack(std::istream &is, Streamable::Format format);

    T bending_modulus;
    T viscosity_contrast;
    T excess_density;

    // derived values
    T vel_coeff;
    T dl_coeff;
    T bending_coeff;

    bool has_contrast;
    bool has_excess_density;

    //! number of primary properties held in this struct (i.e. input)
    static int n_props;
};

#include "VesicleProps.cc"

#endif // _VESICLEPROPS_H_
