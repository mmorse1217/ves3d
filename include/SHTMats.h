#ifndef _SHMATS_H_
#define _SHMATS_H_

#include "Device.h"
#include "OperatorsMats.h"

template<typename T, enum DeviceType DT>
class SHTMats{
  protected:
    int sh_order_;
    pair<int, int> grid_dim_;
    const Device<DT> *device_;
    T *data_;
    
    inline size_t GetDFTLength() const;
    inline size_t GetDLTLength() const;

    void gen_dft_forward();
    void gen_dft_backward();
    void gen_dft_d1backward();
    void gen_dft_d2backward();

    int dft_size;

  public:
    SHTMats(const Device<DT> *dev, int sh_order,
        OperatorsMats<T> &mats, pair<int,int> grid_dim = EMPTY_GRID);
    ~SHTMats();
    
    inline int GetShOrder() const;
    inline pair<int, int> GetGridDim() const;
    inline const Device<DT>* GetDevicePtr() const;
    inline size_t GetDataLength() const;

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
