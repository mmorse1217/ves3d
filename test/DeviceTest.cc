#include "DeviceTest.h"
#include <iostream>
#include <sstream>

int main(int argc, char** argv)
{
    VES3D_INITIALIZE(&argc,&argv,NULL,NULL);
    COUT("\n ==============================\n"
        <<"  Device Test:"
        <<"\n ==============================");

    typedef ostringstream stm;
    typedef string str;
    {
        stm stream;
        str str;
        stream<<CPU;
        COUT(stream.str());
        ASSERT(stream.str()=="CPU","CPU enum");
    }

    {
        stm stream;
        str str;
        stream<<GPU;
        COUT(stream.str());
        ASSERT(stream.str()=="GPU","GPU enum");
    }

    PROFILESTART();
    bool res;
    string sep(60,'-');
    sep = ' ' + sep;
    {
        COUT(sep<<"\n  CPU -- float\n"<<sep);
        int id = 0;
        Error_t err;
        Device<CPU> cpu(id, &err);
        COUT(" * Initialization: "<<err<<" *");
        DeviceTest<CPU,float> dvt_f(&cpu);
        res = dvt_f.PerformAll();
    }

    {
        COUT(sep<<"\n  CPU -- double\n"<<sep);
        Device<CPU> cpu;
        DeviceTest<CPU,double> dvt_d(&cpu);
        res &= dvt_d.PerformAll();
    }

#ifdef GPU_ACTIVE
    {
        {
            COUT(sep<<"\n  GPU -- float\n"<<sep);
            int id = 0;
            Error_t err;
            Device<GPU> gpu(id, &err);
            COUT(" * Initialization: "<<err<<" *");
            DeviceTest<GPU,float> dvt_gf(&gpu);
            res &= dvt_gf.PerformAll();
        }

        {
            COUT(sep<<"\n  GPU -- double\n"<<sep);
            int id = 0;
            Error_t err;
            Device<GPU> gpu(id, &err);
            COUT(" * Initialization: "<<err<<" *");
            DeviceTest<GPU,double> dvt_gf(&gpu);
            res &= dvt_gf.PerformAll();
        }
    }
#endif //GPU_ACTIVE

    string res_print = res ? "works fine." : "is broken!";
    COUT(sep<<"\n  The device class " + res_print<<endl<<sep);

    PROFILEEND("",0);
    PROFILEREPORT(SortFunName);

    VES3D_FINALIZE();
    return 0;
}
