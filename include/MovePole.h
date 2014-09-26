#ifndef _MOVEPOLE_H_
#define _MOVEPOLE_H_

#include "HelperFuns.h"
#include "SHTMats.h"

template<typename Container, typename Operators>
class MovePole
{
  private:
    typedef void (MovePole::*rotPtr_t)(int, int, Container** ) const;
    typedef typename Container::value_type value_type;
    typedef SHTMats<value_type, typename Container::device_type> Mats_t;

  public:
    MovePole(Operators &mats);
    ~MovePole();

    void setOperands(const Container** arr, int num,
        enum SingularStokesRot rot_scheme = Direct);

    inline void operator()(int trg_i, int trg_j,
        Container** results) const;

  private:
    int p_, np_;
    //operator matrices
    Container sp_harm_mats_;
    Container all_rot_mats_;
    Container longitude_rot_;
    int* row_idx;
    int* col_idx;

    SHTrans<Container, Mats_t> sht_;

    //scheme pointer
    rotPtr_t rot_handle_;

    //operands
    enum SingularStokesRot last_rot_;
    const Container** arr_;
    int num_;
    value_type alpha, beta;
    int eager_n_stream_;

    //Work space
    Container rot_mat_;
    Container* shc_;
    mutable Container wrk_, shc_out;
    Container* eager_results_;
    mutable int eager_last_latitude_;
    mutable Container eager_wrk_;

    void movePoleDirectly(int trg_i, int trg_j,
        Container** results) const;

    void movePoleViaSpHarm(int trg_i, int trg_j,
        Container** results) const;

    void alignMeridian(int trg_j, Container** results) const;
    void updateEagerResults(int trg_i) const;
};

#include "MovePole.cc"

#endif // _MOVEPOLE_H_
