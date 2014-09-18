#include <iostream>
#include <sstream>

#include "Enums.h"
#include "Logger.h"

using namespace std;

int main(int argc, char** argv)
{
    COUT("\n ==============================\n"
        <<"  enums Test:"
        <<"\n ==============================\n");

    typedef ostringstream stm;
    typedef string str;
    {
        stm stream;
        str str;
        stream<<PointMajor;
        COUT(stream.str()<<endl);
        ASSERT(stream.str()=="PointMajor","");
    }

    {
        stm stream;
        str str;
        stream<<AxisMajor;
        COUT(stream.str()<<endl);
        ASSERT(stream.str()=="AxisMajor","");
    }

    {
        stm stream;
        str str;
        stream<<Explicit;
        COUT(stream.str()<<endl);
        ASSERT(stream.str()=="Explicit","");
        int sc = EnumifyScheme("Explicit");
        ASSERT(sc==Explicit,"");
    }

    {
        stm stream;
        str str;
        stream<<BlockImplicit;
        COUT(stream.str()<<endl);
        ASSERT(stream.str()=="BlockImplicit","");
        int sc = EnumifyScheme("BlockImplicit");
        ASSERT(sc==BlockImplicit,"");
    }

    {
        stm stream;
        str str;
        stream<<GloballyImplicit;
        COUT(stream.str()<<endl);
        ASSERT(stream.str()=="GloballyImplicit","");
        int sc = EnumifyScheme("GloballyImplicit");
        ASSERT(sc==GloballyImplicit,"");
    }

    {
        std::pair<int, int> gd = EMPTY_GRID;
        ASSERT(gd.first == -1, "Default grid size");
        ASSERT(gd.second == -1, "Default grid size");

        int p(3);
        gd = gridDimOf(p);
        ASSERT(gd.first == p+1, "p+1 in latitude");
        ASSERT(gd.second == 2*p+2, "2*p+2 in longitude");
    }

}
