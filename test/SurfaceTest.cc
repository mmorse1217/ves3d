#include<iostream>
#include<exception>
#include "Vesicle.h"
#include "DataIO.h"

using namespace std;

typedef float scType;

int main(int argc, char* argv[])
{

    int p = 12;
    int np = 2*p*(p+1);
    int vec_len = 3*np;
    DataIO<scType> myIO;
    
    // reading data
    scType* pos_arr = new scType[vec_len];
    myIO.ReadData("../data/cart12",vec_len,pos_arr);
    SHVectors<scType> pos_vec(p,1,pos_arr);

    // initializing Vesicle
    Surface<scType> S(p,1,pos_vec);

    return 0;
}
