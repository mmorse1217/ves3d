#include <iostream>

// #include "DataIO.h"
// #include "Scalars.h"
// #include "Vectors.h"
// #include "OperatorsMats.h"
// #include "Parameters.h"
// #include "MovePole.h"

typedef float real_t;



void CircShift(const real_t x_in, int n_sub, int sub_length, 
    int shift, real_t x_out)
{
    
}

// template<typename ScalarContainer>
// inline void CircShift(const typename ScalarContainer::value_type *x_in,
//     int shift, ScalarContainer &x_out)
// {
//     int vec_length = x_out.getStride();
//     int n_vecs = x_out.getNumSubs();
//     shift = shift % vec_length;
//     shift += (shift < 0) ?  vec_length : 0;

//     int base_in, base_out;
//     for (int ii = 0; ii < n_vecs; ii++) {
//         base_out = ii * vec_length;
//         base_in = base_out + vec_length - shift;
//         x_out.getDevice().Memcpy((typename ScalarContainer::value_type*) 
//             x_out.begin() + base_out, 
//             x_in + base_in, sizeof(typename ScalarContainer::value_type) * shift, 
//             MemcpyDeviceToDevice);
//         base_in = base_out;
//         base_out += shift;
//         x_out.getDevice().Memcpy((typename ScalarContainer::value_type*) 
//             x_out.begin() + base_out, x_in + base_in, 
//             sizeof(typename ScalarContainer::value_type) * (vec_length - shift),
//             MemcpyDeviceToDevice);
//     }
// }
/////

int main(int , char** )
{
//     typedef Scalars<real, DT, the_device> Sca_t;
//     typedef Vectors<real, DT, the_device> Vec_t;
//     typedef OperatorsMats<Sca_t> Mats_t;

//     int nVec(1024), pMax(20);
//     Parameters<real> sim_par;
//     sim_par.sh_order = 24; //Dummy, anything larger than pMax
    
//     bool readFromFile = true;
//     Mats_t mats(readFromFile, sim_par);

//     for(int p=4;p<pMax; ++p)
//     {
//         Vec_t x0(nVec, p), xr(nVec, p);
//         int fLen = x0.getStride();
        
//         //DataIO myIO;
//         //char fname[] = "precomputed/dumbbell_cart12_single.txt";
//         //myIO.ReadData(fname, x0, 0, fLen * DIM);
    
//         SHTrans<Sca_t, SHTMats<real, Device<DT> > > sht(p, mats.mats_p_);
//         MovePole<Sca_t, Mats_t> move_pole(mats);
        
//         const Sca_t* inputs[] = {&x0};
//         Sca_t* outputs[] = {&xr};
    
//         //Profile
//         PROFILECLEAR();
//         PROFILESTART();
//         move_pole.setOperands(inputs, 1);
//         for(int ii=0; ii<x0.getGridDim().first; ++ii)
//             for(int jj=0; jj<x0.getGridDim().second; ++jj)
//                 move_pole(ii, jj, outputs);

//         PROFILEEND("Direct_",0);
//         PROFILEREPORT(SortTime);
//     }
}
