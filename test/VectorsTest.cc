//#include "VectorsTest.h"

#include "Vectors.h"
#include "Logger.h"

extern const Device<CPU> cpu_dev(0);

int main(int argc, char* argv[])
{
    int p(2);
    int num_vecs(2);

    containers::Vectors<float,CPU,cpu_dev> vec(num_vecs, p);
    containers::Scalars<float,CPU,cpu_dev> sc;
    containers::Vectors<float,CPU,cpu_dev> vec2;
    
    cout<<vec<<endl;
    vec2.replicate(vec); cout<<vec2<<endl;
    sc.replicate(vec); cout<<sc<<endl;
//     vec2.resize(0);
//     vec2.replicate(sc); cout<<vec2<<endl;

    return 0;
}
