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
#include "Array.h"

template <typename T, enum DeviceType DT, const Device<DT> &DEVICE>
class Vectors;

template <typename T, enum DeviceType DT, const Device<DT> &DEVICE>
class Scalars : public Array<T, DT, DEVICE>
{
  public:
    explicit Scalars(size_t num_subs = 0, int sh_order = -1, 
        pair<int, int> grid_dim = EMPTY_GRID);
    virtual ~Scalars();

    static inline int getTheDim();

    inline int getShOrder() const;
    inline pair<int, int> getGridDim() const;
    inline size_t getStride() const;
    inline size_t getNumSubs() const;
    inline size_t getSubLength() const;
        
    inline virtual void resize(size_t new_num_subs, int new_sh_order = -1,
        pair<int, int> new_grid_dim = EMPTY_GRID);
    inline void replicate(Scalars<T, DT, DEVICE> const& sc_in);
    inline void replicate(Vectors<T, DT, DEVICE> const& vec_in);
                
    inline T* getSubN(size_t n);
    inline const T* getSubN(size_t n) const;
    
  protected:
    int sh_order_;
    pair<int, int> grid_dim_;
    size_t stride_;
    size_t num_subs_;
    
    static const int the_dim_ = 1;

    Scalars(Scalars<T, DT, DEVICE> const& sc_in);
    Scalars<T, DT, DEVICE>& operator=(const Scalars<T, DT, DEVICE>& rhs);
};
    
template<typename T, enum DeviceType DT, const Device<DT> &DEVICE>
std::ostream& operator<<(std::ostream& output, const Scalars<T, DT, DEVICE>&sc);
    
template<typename Container>
class ShowEntries
{       
  public:
    ShowEntries(const Container &c);
    std::ostream& operator()(std::ostream &out) const;        
      
  private:
    const Container &c_;
};
    
template<typename Container>
std::ostream& operator<<(std::ostream& output, const ShowEntries<Container> &se);
 
#include "Scalars.cc"

#endif //_SCALARS_H_
