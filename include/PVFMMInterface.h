#ifdef HAVE_PVFMM
#ifndef _PVFMM_INTERFACE_H_
#define _PVFMM_INTERFACE_H_

#include <mpi.h>
#include <kernel.hpp>
#include <mpi_tree.hpp>

///////////////////////// Kernel Function Declarations ////////////////////////

template <class T>
void stokes_sl_m2l(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr);

const pvfmm::Kernel<double> ker_stokes_m2l=pvfmm::BuildKernel<double, stokes_sl_m2l>("stokes_m2l", 3, std::pair<int,int>(4,3));

template <class T>
void stokes_sl(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr);

template <class T>
void stokes_dl(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr);

const pvfmm::Kernel<double> ker_stokes=pvfmm::BuildKernel<double, stokes_sl, stokes_dl>("stokes_vel", 3, std::pair<int,int>(3,3),
    NULL, NULL, NULL, &ker_stokes_m2l, &ker_stokes_m2l, &ker_stokes_m2l, NULL, NULL);

///////////////////////////////////////////////////////////////////////////////

template<typename T>
void* PVFMMCreateContext(int n=800, int m=10, int max_d=20,
    pvfmm::BoundaryType bndry=pvfmm::FreeSpace,
    const pvfmm::Kernel<T>* ker=&ker_stokes,
    MPI_Comm comm=MPI_COMM_WORLD);

template<typename T>
void PVFMMDestroyContext(void** ctx);

/**
 * Stokes single-layer interactions between np  particles at locations all_pos
 * and source density all_den. The output velocity is written to all_pot. If
 * ctx_[0]==NULL, a new context is created and stored at ctx_[0].
 */
template<typename T>
void PVFMMEval(const T* all_pos, const T* sl_den, size_t np, T* all_pot, void** ctx_);

template<typename T>
void PVFMMEval(const T* all_pos, const T* sl_den, const T* dl_den, size_t np, T* all_pot, void** ctx_);

template<typename T>
void PVFMMEval(const T* src_pos, const T* sl_den, const T* dl_den, size_t n_src, const T* trg_pos, T* trg_vel, size_t n_trg, void** ctx_);

/**
 * Repartition nv vesicles by the MortonId of their center-of-mass between
 * processors in MPI_COMM_WORLD.
 */
template<typename T>
void PVFMM_GlobalRepart(size_t nv, size_t stride,
    const T* x, const T* tension, size_t* nvr, T** xr,
    T** tensionr, void* user_ptr);

/**
 * Determine bounding box.
 */
template<typename T>
void PVFMMBoundingBox(size_t np, const T* x, T* scale_xr, T* shift_xr, MPI_Comm comm=MPI_COMM_WORLD);

#include "PVFMMInterface.cc"

#endif // _PVFMM_INTERFACE_H_
#endif

