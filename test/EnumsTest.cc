#include <iostream>
#include <sstream>

#include "Enums.h"
#include "Logger.h"
#include "TestTools.h"
#include "ves3d_common.h"

int main(int argc, char** argv)
{
    VES3D_INITIALIZE(&argc,&argv,NULL,NULL);

    COUT("\n ==============================\n"
        <<"  enums Test:"
        <<"\n ==============================\n");

    typedef std::ostringstream strm;
    typedef std::string string;

    int N(2);
    CoordinateOrder cos [] = {PointMajor, AxisMajor};
    const char* conames[] = {"PointMajor", "AxisMajor"};

    for (int i=0;i<N; ++i){
        strm stream;
        stream<<cos[i];
        string msg("testing ");
        msg += conames[i];
        testtools::AssertTrue(stream.str()==conames[i], msg ,"bad string");
    }

    SolverScheme schemes [] = {JacobiBlockExplicit, JacobiBlockGaussSeidel, JacobiBlockImplicit,
                               GloballyImplicit, UnknownScheme};

    const char* snames[] = {"JacobiBlockExplicit", "JacobiBlockGaussSeidel", "JacobiBlockImplicit",
                            "GloballyImplicit", "UnknownScheme"};

    N = 5;
    for (int i=0; i<N; ++i){
     	strm stream;
        stream<<schemes[i];
        string msg("testing ");
        msg += snames[i];
        testtools::AssertTrue(stream.str()==snames[i], msg, "bad string");
        SolverScheme sc = EnumifyScheme(snames[i]);
        testtools::AssertTrue(sc==schemes[i],msg, "bad enum");
    }

    COUT(emph<<"** EnumTest passed **"<<emph<<std::endl);
    VES3D_FINALIZE();
}
