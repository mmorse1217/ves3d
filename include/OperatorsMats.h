#ifndef _OPERATORSMATS_H_
#define _OPERATORSMATS_H_

#include "DataIO.h"
#include "Parameters.h"
#include "SHTMats.h"

template <typename T, typename Device>
struct OperatorsMats
{
  private:
    const Device &device_;
    
  public:
    
    DataIO &fileIO_;
    int p_;
    int p_up_;   

    T *data_;

    SHTMats<T,Device> mats_p_;
    SHTMats<T,Device> mats_p_up_;
    T *quad_weights_;
    T *sing_quad_weights_;
    T *w_sph_;
    T *all_rot_mats_;  
    T *sh_rot_mats_;

    OperatorsMats(const Device &dev, DataIO &fileIO_in, 
        bool readFromFile, const Parameters<T> &params);
    ~OperatorsMats();
    size_t getDataLength() const;
    const Device& getDevice() const;

  private:
    OperatorsMats(const OperatorsMats& mat_in);
    OperatorsMats& operator=(const OperatorsMats& vec_in);
};

#include "OperatorsMats.cc"

#endif //_OPERATORSMATS_H_

