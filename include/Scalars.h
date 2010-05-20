/**
 * @file   Scalars.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Jan 26 14:44:00 2010
 * 
 * @brief The header file for the Scalars class. The implementation
 * of this class is in src/Scalars.cc. 
 */

#ifndef _SCALARS_H_
#define _SCALARS_H_

#include "Device.h"
#include "BiCGStab.h"

template <typename T, enum DeviceType DT> 
class Scalars
{
  protected:
    Device<DT> *device_;
    int sh_order_;
    pair<int, int> grid_dim_;
    size_t fun_length_;
    size_t num_funs_;
    size_t capacity_;
    T *data_;

    static const int the_dim_ = 1;

    Scalars(Scalars<T,DT> const& sc_in);
    Scalars<T,DT>& operator=(const Scalars<T,DT>& rhs);

  public:
    typedef T* iterator;
    typedef const T* const_iterator;
    typedef T value_type;

    Scalars();
    Scalars(Device<DT> *device, int sh_order = 0, 
        size_t num_funs = 0 , pair<int, int> grid_dim = EMPTY_GRID);
    virtual ~Scalars();

    virtual void Resize(size_t new_num_funs, int new_sh_order = 0,
        pair<int, int> new_grid_dim = EMPTY_GRID);
    
    inline const Device<DT>* GetDevicePtr() const;
    inline int GetShOrder() const;
    inline pair<int, int> GetGridDim() const;
    inline size_t GetFunLength() const;
    inline size_t GetNumFuns() const;
    inline size_t Size() const;
    
    inline T& operator[](size_t idx);
    inline const T& operator[](size_t idx) const;
    
    inline iterator begin();
    inline const_iterator begin() const;
    
    inline iterator end();
    inline const_iterator end() const;
    
    friend enum BiCGSReturn BiCGStab<Scalars<T,DT> >(
        void (*MatVec)(Scalars<T,DT> &in, Scalars<T,DT> &out),
        Scalars<T,DT> &x, const Scalars<T,DT> &b, int &max_iter, 
        typename Scalars<T,DT>::value_type &tol,
        void (*Precond)(Scalars<T,DT> &in, Scalars<T,DT> &out));

};

#include "Scalars.cc"

#endif //_SCALARS_H_
