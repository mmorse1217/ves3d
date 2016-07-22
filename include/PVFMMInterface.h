#ifndef _PVFMM_INTERFACE_H_
#define _PVFMM_INTERFACE_H_

#include <mpi.h>
#include <kernel.hpp>
#include <mpi_tree.hpp>

///////////////////////// Kernel Function Declarations ////////////////////////

template <class T>
void stokes_sl_m2l(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr);

template <class T>
void stokes_m2l_vol_poten(const T* coord, int n, T* out);


template <class T>
void stokes_sl(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr);

template <class T>
void stokes_dl(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr);

template <class T>
void stokes_vol_poten(const T* coord, int n, T* out);

template <class Real_t>
struct StokesKernel{
  inline static const pvfmm::Kernel<Real_t>& Kernel();
};

///////////////////////////////////////////////////////////////////////////////

template<typename T>
void* PVFMMCreateContext(T box_size=-1, int n=1000, int m=10, int max_d=MAX_DEPTH,
    const pvfmm::Kernel<T>* ker=&StokesKernel<T>::Kernel(),
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
void PVFMMEval(const T* all_pos, const T* sl_den, const T* dl_den, size_t np, T* all_pot, void** ctx_, int setup=1);

template<typename T>
void PVFMMEval(const T* src_pos, const T* sl_den, const T* dl_den, size_t n_src, const T* trg_pos, T* trg_vel, size_t n_trg, void** ctx_, int setup=1);

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

