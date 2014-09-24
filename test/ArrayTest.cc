#include "Array.h"
#include <unistd.h>  //for sleep()
#include <iostream>
#include "Logger.h"
#include "Device.h"

using namespace std;

template<typename Arr>
void test_array(){

    typedef typename Arr::value_type T;
    typedef typename Arr::device_type DT;

    COUT(" -- Testing array with device "<<Arr::getDevice().type());
    {
        size_t sz(6);
        COUT(" . Test initialization");

        Arr a,b(sz);
        ASSERT(a.size()==0,"zero size");
        ASSERT(b.size()==sz,"non-zero size");

        a.resize(sz);
        ASSERT(a.size()==sz,"resize");
    }

    {
        size_t sz(6);
        COUT(" . Test iterator");

        Arr a,b(sz);
        ASSERT((a.end()-a.begin())==a.size(),"pointer arithmetic");
        ASSERT((b.end()-b.begin())==b.size(),"pointer arithmetic");

        T *c =  new T[sz];
        T *d =  new T[sz];
        for(size_t i=0;i<sz;++i)
            c[i]=i*i;

        b.getDevice().Memcpy(
            b.begin(),
            c,
            b.mem_size(),
            DT::MemcpyHostToDevice);

        b.getDevice().Memcpy(
            d,
            b.begin(),
            b.mem_size(),
            DT::MemcpyDeviceToHost);

        for(size_t i=0;i<sz;++i)
            ASSERT(d[i]==i*i,"memcpy");
    }
}

typedef Device<CPU> C;
extern const C cpu_dev(0);
typedef Array<double, C, cpu_dev> CArr;

#ifdef GPU_ACTIVE
typedef Device<CPU> G;
extern const G gpu_dev(1);
typedef Array<double, G, gpu_dev> GArr;
#endif //GPU_ACTIVE

int main(int argc, char** argv)
{
    COUT("\n ==============================\n"
        <<"  Array Test:"
        <<"\n ==============================");
    sleep(1);

    PROFILESTART();
    test_array<CArr>();
    COUT(emph<<" ** CPU passed **"<<emph);
#ifdef GPU_ACTIVE
    test_array<GArr>();
    COUT(emph<<"** GPU passed **"<<emph);
#endif //GPU_ACTIVE
    PROFILEEND("",0);
    PROFILEREPORT(SortFunName);

    sleep(1);
    return 0;
}
