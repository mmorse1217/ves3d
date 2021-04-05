This codebase simulates deformable vesicles in Stokes flows and resolve inter-cell collisions and scales to thousands of cores.
It is not maintained and is made available for reference purposes; we do not guarantee timely responses to feature requests or bug fixes. 
It contributed to the following publications:

 - Dhairya Malhotra, Abtin Rahimian, Denis Zorin, and George Biros. "A parallel algorithm for long-timescale simulation of concentrated vesicle suspensions in three dimensions." (2017) [Preprint](https://cims.nyu.edu/~malhotra/files/pubs/ves3d.pdf) (2017).
 - Libin Lu, Abtin Rahimian, and Denis Zorin. "Parallel contact-aware simulations of deformable particles in 3D Stokes flow." (2018) [arXiv preprint arXiv:1812.04719](https://arxiv.org/pdf/1812.04719)
 - Libin Lu, Matthew J. Morse, Abtin Rahimian, Georg Stadler, and Denis Zorin.  "Scalable simulation of realistic volume fraction red blood cell flows through vascular networks." [Supercomputing 2019](https://dl.acm.org/doi/pdf/10.1145/3295500.3356203?casa_token=MSlWLB6jRGQAAAAA:YvpdvZpeUcNF4kj0_BjyVZMg0GyUSeqt9KE-Xb4EbhLTOu6EtUJJrZ7ixNrpqm3tL2S_C9XjwYvA)


# 3D particulate flow simulation code. #
---

## Before build: ##
---
1. Copy the ves3d project files to the install location and make an
    environment variable pointing to that location:
```
#!shell
    $ cd ves3d-cxx
    $ export VES3D_DIR=`pwd`
```

1. In makefile.in.files directory, make a file for the machine you're
installing VES3D. You can copy one of the existing machine files, for
example
```
#!shell
    $ cd makefile.in.files
    $ cp makefile.octane makefile.`hostname -s`
```
You can alternatively use one of the existing host files by passing
the VES3D_PLATFORM variable to make.

1. In the platform specific makefile you created in the previous step,
    set the correct include paths for libraries

1. If the makefile corresponding to your choice of compiler does
    not exist, you also need to write that makefile (with the correct
    extension)

 1. Optionally add makefile.<hostname> to hg for future reference

## Building: ##
---
 1. cd to VES3D_DIR and
```
#!shell

    $ make
    $ make test
    $ make check
```

 1. If the code compiles fine and tests pass, run
```
#!shell

    $ make install
```
## Misc: ##
---
1. Code-specific compiler options are: VERBOSE, PROFILING, VES3D_TESTING, VES3D_USE_GPU, VES3D_USE_PVFMM
1. If you want the GPU code compiled set "VES3D_USE_GPU=yes" in makefile.<hostname>
1. Fully functional revisions of code are tagged by vYY.ID (YY is year and ID is an integer)
