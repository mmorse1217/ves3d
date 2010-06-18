//#include "VectorsTest.h"

#include "Vectors.h"
#include "Logger.h"
#include "HelperFuns.h"

extern const Device<CPU> cpu_dev(0);

int main(int argc, char* argv[])
{
    int p(2);
    int num_vecs(2);

//     containers::Vectors<float,CPU,cpu_dev> vec(num_vecs, p);
//     containers::Scalars<float,CPU,cpu_dev> sc;
//     containers::Vectors<float,CPU,cpu_dev> vec2;
    
//     cout<<vec<<endl;
//     vec2.replicate(vec); cout<<vec2<<endl;
//     sc.replicate(vec); cout<<sc<<endl;
//     vec2.resize(0);
//     vec2.replicate(sc); cout<<vec2<<endl;

    int nv = 2;
    int N = 5;
    containers::Vectors<float, CPU, cpu_dev> vec(nv, 0, make_pair(N,1));
    containers::Vectors<float, CPU, cpu_dev> vecshuffle(nv, 0, make_pair(N,1));
    
    for(int ii=0; ii<vec.size(); ++ii)
        vec[ii] = ii;
    
    cout<<vec<<endl;
    ShufflePoints(vec, vecshuffle);
    cout<<vecshuffle<<endl;

    ShufflePoints(vecshuffle, vec);
    cout<<vec<<endl;

    return 0;
}
