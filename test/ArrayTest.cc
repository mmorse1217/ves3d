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

    cout<<" -- Testing array with device "<<Arr::getDevice().type()<<endl;
    {
        size_t sz(6);
        cout<<" . Test initialization"<<endl;

        Arr a,b(sz);
        assert(a.size()==0);
        assert(b.size()==sz);

        a.resize(sz);
        assert(a.size()==sz);
    }

    {
        size_t sz(6);
        cout<<" . Test iterator"<<endl;

        Arr a,b(sz);
        assert((a.end()-a.begin())==a.size());
        assert((b.end()-b.begin())==b.size());

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
            assert(d[i]==i*i);
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
        <<"\n ==============================\n");
    sleep(1);

    PROFILESTART();
    test_array<CArr>();
#ifdef GPU_ACTIVE
    test_array<GArr>();
#endif //GPU_ACTIVE
    PROFILEEND("",0);
    PROFILEREPORT(SortFunName);

    sleep(1);
    return 0;
}
