#ifndef _MOVEPOLE_H_
#define _MOVEPOLE_H_

#include "HelperFuns.h"

template<typename Container>
class MovePole
{
  private:
    typedef void (MovePole::*rotPtr_t)(int, int, Container** ) const;
    typedef typename Container::value_type value_type;

  public:
    MovePole(const Container &all_rot_mats, Container &rot_mat, 
             const Container &sp_mats);
    
    void setOperands(const Container** arr, int num, 
        enum SingularStokesRot rot_scheme);

    inline void operator()(int trg_i, int trg_j, 
        Container** results) const;
    
  private:
    const Container &sp_harm_mats_;
    const Container &all_rot_mats_;
    Container &rot_mat_;
   
    const Container** arr_;
    int num_;
    rotPtr_t rot_handle_;
    value_type alpha, beta;
    
    Container* shc_;
    Container wrk_;

    void movePoleDirectly(int trg_i, int trg_j,
        Container** results) const;
    
    void movePoleViaSpHarm(int trg_i, int trg_j, 
        Container** results) const;
};

#include "MovePole.cc"

#endif // _MOVEPOLE_H_
