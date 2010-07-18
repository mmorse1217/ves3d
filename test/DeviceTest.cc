#include "DeviceTest.h"

int main(int argc, char** argv)
{
    COUT("\n ==============================\n"
        <<"  Device Test:"
        <<"\n ==============================\n");
    sleep(1);

    PROFILESTART();
    bool res;
    string sep(60,'-');
    sep = ' ' + sep;
    {
        cout<<sep<<"\n  CPU -- float\n"<<sep<<endl;
        int id = 0;
        enum DeviceError err;
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
        cout<<sep<<"\n  GPU -- float\n"<<sep<<endl;
        int id = 0;
        enum DeviceError err;
        Device<GPU> gpu(id, &err);
        cout<<" * Initialization: "<<err<<" *"<<endl;
        DeviceTest<GPU,float> dvt_gf(&gpu);
        res &= dvt_gf.PerformAll();
        sleep(.5);
    }
#endif //GPU_ACTIVE

    string res_print = res ? "works fine." : "is broken!";
    cout<<sep<<"\n  The device class " + res_print<<endl<<sep<<endl;
   
    PROFILEEND("",0);
    PROFILEREPORT(SortTime);
    
    sleep(1);
    return 0;
}
