#ifndef _SCALARS_H_
#define _SCALARS_H_

#include "Device.h"
#include <cassert>

template <typename T> 
class Scalars
{
  private:
    Device<T> *device_;
    int p_; 
    size_t stride_;
    size_t num_funs_;   
    size_t capacity_;
    T *data_;          

    Scalars(const Scalars<T> &sc_in);
    Scalars<T>& operator=(const Scalars<T> &sc_in);
  
  public:
    typedef T* iterator;
    typedef const T* const_iterator;

    inline const Device<T>* GetDevicePtr() const;
    inline int GetShOrder() const;
    inline size_t GetFunLength() const;
    inline size_t GetNumFuns() const;
    inline size_t GetDataLength() const;

    inline T& operator[](size_t idx);
    inline const T& operator[](size_t idx) const;

    inline iterator begin();
    inline const_iterator begin() const;
    
    inline iterator end();
    inline const_iterator end() const;

    ///////////////////////////////////////////////////////////////
    Scalars(Device<T> *device_in, int p_in = 0, int n_funs_in = 0);
    virtual ~Scalars();
    void Sqrt();
    T Max();
    virtual void Resize(int n_funs_in);
};

#include "Scalars.cc"

#endif //_SCALARS_H_
