#ifndef _ARRAY_H_
#define _ARRAY_H_

#include "Device.h"

template <typename T, enum DeviceType DT, const Device<DT> &DEVICE>
class Array
{
  public:
    typedef Device<DT> device_type;
    typedef T* iterator;
    typedef const T* const_iterator;
    typedef T value_type;
    
    explicit Array(size_t size = 0);
    virtual ~Array();

    static enum DeviceType getDeviceType();
    static const device_type& getDevice();
    
    inline size_t size() const;
    inline void resize(size_t new_size);
                
    inline iterator begin();
    inline const_iterator begin() const;

    inline iterator end();
    inline const_iterator end() const;

    bool isNan() const;
    void fillRand();
    
  protected:
    size_t size_;
    size_t capacity_;
    T* data_;

    Array(Array const& rhs);
    Array& operator=(const Array& rhs);
};

#include "Array.cc"

#endif // _ARRAY_H_
