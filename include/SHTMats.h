#ifndef _SHMATS_H_
#define _SHMATS_H_

#include <cmath>
#include "Enums.h"
#include "Logger.h"

template<typename T, typename Device>
class SHTMats{
  protected:
    int sh_order_;
    std::pair<int, int> grid_dim_;
    T *data_;
    int dft_size;
    const Device &device_;

    void gen_dft_forward();
    void gen_dft_backward();
    void gen_dft_d1backward();
    void gen_dft_d2backward();

  public:
    SHTMats(const Device &dev, int sh_order,
        T *data, bool genrateMats = false,
        std::pair<int, int> gird_dim = EMPTY_GRID);
    ~SHTMats();

    inline int getShOrder() const;
    inline std::pair<int, int> getGridDim() const;

    static inline size_t getDataLength(int sh_order,
        std::pair<int, int> grid_dim = EMPTY_GRID);
    inline size_t getDataLength() const;
    inline size_t getDFTLength() const;
    inline size_t getDLTLength() const;
    inline T* getData();
    inline const Device& getDevice() const;

    T *dft_;
    T *dft_inv_;
    T *dft_inv_d1_;
    T *dft_inv_d2_;

    T *dlt_;
    T *dlt_inv_;
    T *dlt_inv_d1_;
    T *dlt_inv_d2_;
};

#include "SHTMats.cc"

#endif //_SHMATS_H_
