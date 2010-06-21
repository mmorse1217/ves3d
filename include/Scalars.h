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
#include <queue>

namespace containers
{
    template <typename T, enum DeviceType DT, const Device<DT> &DEVICE>
    class Vectors;

    template <typename T, enum DeviceType DT, const Device<DT> &DEVICE>
    class Scalars
    {
      public:
        typedef T* iterator;
        typedef const T* const_iterator;
        typedef T value_type;

        explicit Scalars(size_t num_subs = 0, int sh_order = 0, 
            pair<int, int> grid_dim = EMPTY_GRID);
        virtual ~Scalars();
        
        static const Device<DT>& getDevice();
        static inline int getTheDim();

        inline int getShOrder() const;
        inline pair<int, int> getGridDim() const;
        inline size_t getStride() const;
        inline size_t size() const;
        inline size_t getNumSubs() const;
        
        void resize(size_t new_num_subs, int new_sh_order = -1,
            pair<int, int> new_grid_dim = EMPTY_GRID);
        inline void replicate(Scalars<T, DT, DEVICE> const& sc_in);
        inline void replicate(Vectors<T, DT, DEVICE> const& vec_in);
                
        inline value_type& operator[](size_t idx);
        inline const value_type& operator[](size_t idx) const;
    
        inline iterator begin();
        inline const_iterator begin() const;

        inline iterator getSubN(size_t n);
        inline const_iterator getSubN(size_t n) const;

        inline iterator end();
        inline const_iterator end() const;

        //void* operator new(size_t size) throw(bad_alloc);
        //void operator delete(void* ptr) throw();

      protected:
        int sh_order_;
        pair<int, int> grid_dim_;
        size_t stride_;
        size_t num_subs_;
        size_t capacity_;
        T *data_;

        static const int the_dim_ = 1;
        //static queue<void* > scalars_queue_;

        Scalars(Scalars<T, DT, DEVICE> const& sc_in);
        Scalars<T, DT, DEVICE>& operator=(const Scalars<T, DT, DEVICE>& rhs);
};
    
    template<typename T, enum DeviceType DT, const Device<DT> &DEVICE>
    std::ostream& operator<<(std::ostream& output, const Scalars<T, DT, DEVICE>&sc);
    
    //template<typename T, enum DeviceType DT, const Device<DT> &DEVICE>
    //queue<void* > Scalars<T, DT, DEVICE>::scalars_queue_;

#include "Scalars.cc"
}

#endif //_SCALARS_H_
