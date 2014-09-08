#ifdef HAVE_PVFMM
#ifndef _PVFMM_INTERFACE_H_
#define _PVFMM_INTERFACE_H_

#include <mpi.h>
#include <kernel.hpp>
#include <mpi_tree.hpp>

template<typename T>
void* PVFMMCreateContext(int n=400, int m=10, int max_d=20,
    pvfmm::BoundaryType bndry=pvfmm::FreeSpace,
    const pvfmm::Kernel<T>* ker=&pvfmm::ker_stokes_vel,
    const pvfmm::Kernel<T>* aux_ker=NULL,
    MPI_Comm comm=MPI_COMM_WORLD);

template<typename T>
void PVFMMDestroyContext(void** ctx);

/**
 * Stokes single-layer interactions between np  particles at locations all_pos
 * and source density all_den. The output velocity is written to all_pot. If
 * ctx_[0]==NULL, a new context is created and stored at ctx_[0].
 */
template<typename T>
void PVFMMEval(const T* all_pos, const T* all_den, size_t np, T* all_pot, void** ctx_);

/**
 * Repartition nv vesicles by the MortonId of their center-of-mass between
 * processors in MPI_COMM_WORLD.
 */
template<typename T>
void PVFMM_GlobalRepart(size_t nv, size_t stride,
    const T* x, const T* tension, size_t* nvr, T** xr,
    T** tensionr, void* user_ptr);

#include "PVFMMInterface.cc"

#endif // _PVFMM_INTERFACE_H_
#endif

