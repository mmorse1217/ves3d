#include<iostream>
#include<exception>
#include "Surface.h"
#include "DataIO.h"

using namespace std;
typedef float scType;

int main(int argc, char* argv[])
{
    int p = 12;
    int np = 2*p*(p+1);
    int vec_len = 3*np;
    int m = 10;
    scType ts = .01;
    
    if(argc>1)
        m = atoi(argv[1]);
    
    if(argc>2)
        ts = (scType) atof(argv[2]);
    
    // reading data
    DataIO<scType> myIO;
    scType* pos_arr = new scType[vec_len];
    myIO.ReadData("../data/ellipse_cart12",vec_len,pos_arr);
    SHVectors<scType> pos_vec(p,1,pos_arr);

    // initializing Vesicle
    Surface<scType> S(p,1,pos_vec);
    SHVectors<scType> bend_force(p,1);

    for(int tt=0;tt<m;++tt)
    {
        AxPy(S.k_, S.normal_, (scType) 0.0, bend_force);
        //Add pressure to preserve volume
        AxPy(-ts,bend_force, S.x_, pos_vec);
        S.SetX(pos_vec);
    }

    myIO.WriteData("./cart12_final",vec_len,S.x_.data_);

    delete[] pos_arr;
    return 0;
}
