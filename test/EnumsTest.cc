#include <iostream>
#include <sstream>
#include <cassert>

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
        stream<<CPU;
        COUT(stream.str()<<endl);
        assert(stream.str()=="CPU");
    }

    {
        stm stream;
        str str;
        stream<<GPU;
        COUT(stream.str()<<endl);
        assert(stream.str()=="GPU");
    }

    {
        stm stream;
        str str;
        stream<<MemcpyHostToHost;
        COUT(stream.str()<<endl);
        assert(stream.str()=="MemcpyHostToHost");
    }

    {
        stm stream;
        str str;
        stream<<MemcpyHostToDevice;
        COUT(stream.str()<<endl);
        assert(stream.str()=="MemcpyHostToDevice");
    }

    {
        stm stream;
        str str;
        stream<<MemcpyDeviceToHost;
        COUT(stream.str()<<endl);
        assert(stream.str()=="MemcpyDeviceToHost");
    }

    {
        stm stream;
        str str;
        stream<<MemcpyDeviceToDevice;
        COUT(stream.str()<<endl);
        assert(stream.str()=="MemcpyDeviceToDevice");
    }

    {
        stm stream;
        str str;
        stream<<PointMajor;
        COUT(stream.str()<<endl);
        assert(stream.str()=="PointMajor");
    }

    {
        stm stream;
        str str;
        stream<<AxisMajor;
        COUT(stream.str()<<endl);
        assert(stream.str()=="AxisMajor");
    }

    {
        stm stream;
        str str;
        stream<<Explicit;
        COUT(stream.str()<<endl);
        assert(stream.str()=="Explicit");
        int sc = EnumifyScheme("Explicit");
        assert(sc==Explicit);
    }

    {
        stm stream;
        str str;
        stream<<BlockImplicit;
        COUT(stream.str()<<endl);
        assert(stream.str()=="BlockImplicit");
        int sc = EnumifyScheme("BlockImplicit");
        assert(sc==BlockImplicit);
    }

    {
        stm stream;
        str str;
        stream<<GloballyImplicit;
        COUT(stream.str()<<endl);
        assert(stream.str()=="GloballyImplicit");
        int sc = EnumifyScheme("GloballyImplicit");
        assert(sc==GloballyImplicit);
    }
}
