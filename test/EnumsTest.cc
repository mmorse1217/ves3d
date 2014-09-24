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
        COUT(stream.str());
        ASSERT(stream.str()=="PointMajor","");
    }

    {
        stm stream;
        str str;
        stream<<AxisMajor;
        COUT(stream.str());
        ASSERT(stream.str()=="AxisMajor","");
    }

    {
        stm stream;
        str str;
        stream<<Explicit;
        COUT(stream.str());
        ASSERT(stream.str()=="Explicit","");
        int sc = EnumifyScheme("Explicit");
        ASSERT(sc==Explicit,"");
    }

    {
        stm stream;
        str str;
        stream<<BlockImplicit;
        COUT(stream.str());
        ASSERT(stream.str()=="BlockImplicit","");
        int sc = EnumifyScheme("BlockImplicit");
        ASSERT(sc==BlockImplicit,"");
    }

    {
        stm stream;
        str str;
        stream<<GloballyImplicit;
        COUT(stream.str());
        ASSERT(stream.str()=="GloballyImplicit","");
        int sc = EnumifyScheme("GloballyImplicit");
        ASSERT(sc==GloballyImplicit,"");
    }
    COUT(emph<<"** EnumTest passed **"<<emph<<std::endl);
}
