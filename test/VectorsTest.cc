//#include "VectorsTest.h"

#include "Vectors.h"

extern const Device<CPU> cpu_dev;

int main(int argc, char* argv[])
{
    int p(2);
    int num_vecs(2);

    containers::Vectors<float,CPU,cpu_dev> vec(num_vecs, p);
    containers::Scalars<float,CPU,cpu_dev> sc(num_vecs, p);
    containers::Scalars<float,CPU,cpu_dev> vec2;
    
    vec2.replicate(vec);
    
    cout<<vec<<endl;
    cout<<vec2<<endl;
        
    return 0;
}
