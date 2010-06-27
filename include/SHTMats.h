#ifndef _SHMATS_H_
#define _SHMATS_H_

template<typename T>
class SHTMats{
  protected:
    int sh_order_;
    pair<int, int> grid_dim_;
    T *data_;
    int dft_size;
        
    void gen_dft_forward();
    void gen_dft_backward();
    void gen_dft_d1backward();
    void gen_dft_d2backward();

  public:
    SHTMats(int sh_order, T *data, bool genrateMats = false, 
        pair<int, int> gird_dim = EMPTY_GRID);
    ~SHTMats();

    inline int getShOrder() const;
    inline pair<int, int> getGridDim() const;
    
    static inline size_t getDataLength(int sh_order, 
        pair<int, int> grid_dim = EMPTY_GRID);
    inline size_t getDataLength() const;
    inline size_t getDFTLength() const;
    inline size_t getDLTLength() const;
    inline T* getData();


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
