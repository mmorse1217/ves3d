#ifndef _OPERATORSMATS_H_
#define _OPERATORSMATS_H_

#include "SHTMats.h"
#include "DataIO.h"
#include "Parameters.h"

template <typename Container>
struct OperatorsMats
{
  public:
    typedef typename Container::value_type value_type;
    typedef typename Container::device_type device_type;

    int p_;
    int p_up_;

    Container data_;

    SHTMats<value_type,device_type> mats_p_;
    SHTMats<value_type,device_type> mats_p_up_;

    value_type *quad_weights_;
    value_type *quad_weights_p_up_;
    value_type *sing_quad_weights_;
    value_type *w_sph_;
    value_type *all_rot_mats_;
    value_type *sh_rot_mats_;

    OperatorsMats(bool readFromFile,
        const Parameters<value_type> &params);

    size_t getDataLength(const Parameters<value_type>
        &params) const;

  private:
    OperatorsMats(const OperatorsMats& mat_in);
    OperatorsMats& operator=(const OperatorsMats& vec_in);
};

#include "OperatorsMats.cc"

#endif //_OPERATORSMATS_H_
