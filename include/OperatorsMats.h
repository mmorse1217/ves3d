#ifndef _OPERATORSMATS_H_
#define _OPERATORSMATS_H_

#include "DataIO.h"
#include "Parameters.h"
#include "SHTMats.h"

template <typename T, typename IO>
struct OperatorsMats
{
  public:
    IO &fileIO_;
    int p_;
    int p_up_;   

    T *data_;

    SHTMats<T> mats_p_;
    SHTMats<T> mats_p_up_;
    T *quad_weights_;
    T *sing_quad_weights_;
    T *w_sph_;
    T *all_rot_mats_;  

    OperatorsMats(IO &fileIO_in, bool readFromFile, const Parameters<T> &params);
    ~OperatorsMats();
    size_t getDataLength();
    
  private:
    OperatorsMats(const OperatorsMats<T, IO>& mat_in);
    OperatorsMats<T, IO>& operator=(const OperatorsMats<T, IO>& vec_in);
};

#include "OperatorsMats.cc"

#endif //_OPERATORSMATS_H_

