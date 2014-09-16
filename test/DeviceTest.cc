#include <unistd.h>  //for sleep()
#include "DeviceTest.h"
#include <iostream>
#include <sstream>

int main(int argc, char** argv)
{
    COUT("\n ==============================\n"
        <<"  Device Test:"
        <<"\n ==============================\n");
    sleep(1);

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

    PROFILESTART();
    bool res;
    string sep(60,'-');
    sep = ' ' + sep;
    {
        cout<<sep<<"\n  CPU -- float\n"<<sep<<endl;
        int id = 0;
        Error_t err;
        Device<CPU> cpu(id, &err);
        cout<<" * Initialization: "<<err<<" *"<<endl;
        DeviceTest<CPU,float> dvt_f(&cpu);
        res = dvt_f.PerformAll();
        sleep(1);
    }

    {
        cout<<sep<<"\n  CPU -- double\n"<<sep<<endl;
        Device<CPU> cpu;
        DeviceTest<CPU,double> dvt_d(&cpu);
        res &= dvt_d.PerformAll();
        sleep(1);
    }

#ifdef GPU_ACTIVE
    {
        {
            cout<<sep<<"\n  GPU -- float\n"<<sep<<endl;
            int id = 0;
            Error_t err;
            Device<GPU> gpu(id, &err);
            cout<<" * Initialization: "<<err<<" *"<<endl;
            DeviceTest<GPU,float> dvt_gf(&gpu);
            res &= dvt_gf.PerformAll();
            sleep(1);
        }

        {
            cout<<sep<<"\n  GPU -- double\n"<<sep<<endl;
            int id = 0;
            Error_t err;
            Device<GPU> gpu(id, &err);
            cout<<" * Initialization: "<<err<<" *"<<endl;
            DeviceTest<GPU,double> dvt_gf(&gpu);
            res &= dvt_gf.PerformAll();
            sleep(1);
        }
    }
#endif //GPU_ACTIVE

    string res_print = res ? "works fine." : "is broken!";
    cout<<sep<<"\n  The device class " + res_print<<endl<<sep<<endl;

    PROFILEEND("",0);
    PROFILEREPORT(SortFunName);

    sleep(1);
    return 0;
}
