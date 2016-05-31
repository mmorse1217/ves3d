/**
 * @file
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @date   2014-09-01 20:23
 *
 * @brief Container class templated on datatype and device
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

#ifndef _ARRAY_H_
#define _ARRAY_H_

#include "Logger.h"
#include "Streamable.h"
#include <iostream> //also has size_t

/**
 * @class Array
 * @param T       Content datatype (float,double,etc.)
 * @param DT      Device type
 * @param DEVICE  Reference to the device object
 *
 * @note Device class parameter (DT) is needed because non-type
 * template parameter (DEVICE) can't have a void type and only objects
 * having the same type as in declaration can be used during
 * instantiation.
 *
 * @warning Defined iterators are not be use for direct iteration
 * because the data is on the device. They are intended for functions
 * that are aware of the device.
 */
template <typename T,       // content type
          typename DT,      // device class
          const DT &DEVICE>
class Array : public Streamable
{
  public:
    typedef T value_type;
    typedef DT device_type;
    typedef T* iterator;
    typedef const T* const_iterator;
    typedef size_t size_type;

    explicit Array(size_t size = 0);
    explicit Array(std::istream &is, Format format);

    virtual ~Array();

    static const device_type& getDevice();

    inline size_t size() const;
    inline size_t mem_size() const;

    //! Realocates new_size and copys the current content to the new
    //! location. When new_size is zero, it frees the current
    //! allocated memory
    inline void resize(size_t new_size);

    inline iterator begin();
    inline const_iterator begin() const;

    inline iterator end();
    inline const_iterator end() const;

    // From streamable class --------------------------------------------------
    // ------------------------------------------------------------------------
    virtual Error_t pack(std::ostream &os, Format format) const;
    virtual Error_t unpack(std::istream &is, Format format);

  protected:
    size_t size_;
    size_t capacity_;
    T* data_;

    //! private copy constructory to limit pass by value
    Array(Array const& rhs);
    //! private assignment operator to limit pass by value
    Array& operator=(const Array& rhs);
};

//! streaming function for array type
template<typename T, typename DT, const DT &DEVICE>
std::ostream& operator<<(std::ostream& output,
    const Array<T, DT, DEVICE> &arr);

#include "Array.cc"

#endif // _ARRAY_H_
