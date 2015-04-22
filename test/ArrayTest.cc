#include "Array.h"
#include "Logger.h"
#include "Device.h"
#include <sstream>

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

    { //streaming
        size_t sz(11);
        COUT(" . Test streaming");

        Arr a(sz),b;
        T *c =  new T[sz];
        for(size_t i=0;i<sz;++i)
            c[i]=i*i*i;

        a.getDevice().Memcpy(
            a.begin(),
            c,
            a.mem_size(),
            DT::MemcpyHostToDevice);

	a.set_name("AAA");

	std::stringstream s1,s2,s3,s4;
	std::string ref("ARRAY\nname: AAA\nsize: 11\ndata: 0 1 8 27 64 125 216 343 512 729 1000\n/ARRAY\n");
	a.pack(s1, Streamable::ASCII);

	ASSERT(s1.str()==ref, "bad stream a");
	b.unpack(s1, Streamable::ASCII);

	b.pack(s2, Streamable::ASCII);
	ASSERT(s2.str()==ref, "bad stream b");

	Arr d(s2, Streamable::ASCII);
	d.pack(s3, Streamable::ASCII);
	ASSERT(s3.str()==ref, "bad stream d");

	delete[] c;
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

    PROFILESTART();
    test_array<CArr>();
    COUT(emph<<" ** CPU passed **"<<emph);
#ifdef GPU_ACTIVE
    test_array<GArr>();
    COUT(emph<<"** GPU passed **"<<emph);
#endif //GPU_ACTIVE
    PROFILEEND("",0);
    PROFILEREPORT(SortFunName);

    return 0;
}
