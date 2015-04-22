/**
 * @file
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @revision $Revision$
 * @tags $Tags$
 * @date $Date$
 *
 * @brief base class for streamable objects intended for IO or MPI
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

#ifndef _STREAMABLE_H_
#define _STREAMABLE_H_

#include "Error.h"
#include "Logger.h"
#include <iostream>

/**
 * Base class for streamable objects. Subclasses provide interface for
 * reading and writing from bindary or ascii data.
 *
 */
class Streamable{
  public:
    enum Format {BIN,ASCII};

    Streamable(const std::string name="StreamableObject") : name_(name) {}
    virtual ~Streamable(){}

    //! Name that is optionally used in the stream header.
    Error_t set_name(std::string newname);
    std::string name();

    //! converts the streamable to a stream depending on format. This
    //! appends its data to the stream. It is expected that derived
    //! classes do (1) do not change stream flags such as precision
    //! (or reset upon return) and (2) prepend/append proper heading
    //! and footing to their data.
    virtual Error_t pack(std::ostream &os, Format format) const = 0;

    //! The inverse of pack. This potentially changes the state of
    //! object.
    virtual Error_t unpack(std::istream &is, Format format) = 0;

    // ------------------------------------------------------------------------
    // unimplemented function (should be virtual ... )
    // ------------------------------------------------------------------------
    //! converts the streamable to an ascii cstring stored in the
    //! buffer. If buffer is NULL, it is allocated. The length of
    //! cstring is returned in len. If len is set when calling, it is
    //! assumed to be the length of the buffer.
    virtual Error_t pack_cstring(char** buffer, size_t &len, size_t offset) const
    { return ErrorEvent::NotImplementedError; }

    //! Similar to pack_cstring but binary
    virtual Error_t pack_bin(void** buffer, size_t &len, size_t offset) const
    { return ErrorEvent::NotImplementedError; }

    //! Unpack data from ascii buffer to self. This may change the
    //! object's state.
    virtual Error_t unpack_cstring(const char* buffer, size_t len, size_t offset)
    { return ErrorEvent::NotImplementedError; }

    //! Similar to pack_cstring but binary
    virtual Error_t unpack_bin(const void* buffer, size_t len, size_t offset)
    { return ErrorEvent::NotImplementedError; }

    // ------------------------------------------------------------------------
    // utility function
    // ------------------------------------------------------------------------

    //! pack an array of type T and size n to the stream.
    template<typename T>
    Error_t pack_array(std::ostream &os, Format format, const T* arr, size_t n, const char *sep=" ") const;

    //! uppack a stream into an array of type T. The array is expected
    //! to be of size n. This consumes the stream until the end or for
    //! n entires. If the stream terminates before extracting n
    //! entires, this returns with an error and also set n to the
    //! number of extracted elements.
    template<typename T>
    Error_t unpack_array(std::istream &is, Format format, T* arr, size_t &n, const char *sep=" ") const;

  protected:
    std::string name_;
};

//! overloading insertion operator. Note that the argument is a
//! pointer instead of the object to avoid slicing and reusing pack
//! function.
std::ostream& operator<<(std::ostream& os, const Streamable *o);

#include "Streamable.cc"

#endif //_STREAMABLE_H_
