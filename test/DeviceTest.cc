#include "DeviceTest.h"

int main(int argc, char** argv)
{
    PROFILESTART();
    bool res;
    string sep(60,'-');

    {
        cout<<sep<<"\n  CPU -- float\n"<<sep<<endl;
        int id = 0;
        enum DeviceError err;
        Device<CPU> cpu(id, &err);
        cout<<"* Initialization: "<<err<<" *"<<endl;
        DeviceTest<CPU,float> dvt_f(&cpu);
        res = dvt_f.PerformAll();
    }

    {
        cout<<sep<<"\n  CPU -- double\n"<<sep<<endl;
        Device<CPU> cpu;
        DeviceTest<CPU,double> dvt_d(&cpu);
        res &= dvt_d.PerformAll();
    }

#ifdef GPU_ACTIVE
    {
        cout<<sep<<"\n  GPU -- float\n"<<sep<<endl;
        int id = 0;
        enum DeviceError err;
        Device<GPU> gpu(id, &err);
        cout<<"* Initialization: "<<err<<" *"<<endl;
        DeviceTest<GPU,float> dvt_gf(&gpu);
        res &= dvt_gf.PerformAll();
    }
#endif //GPU_ACTIVE

    string res_print = res ? "works fine." : "is broken!";
    cout<<sep<<"\n The device class " + res_print<<endl<<sep<<endl;
   
    PROFILEEND("",0);
    PROFILEREPORT(SortTime);

    return 0;
}
