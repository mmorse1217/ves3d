#ifndef _VESINTERACTION_H_
#define _VESINTERACTION_H_

#include <typeinfo>
#include "enums.h"
#include <omp.h>
/**
 * The gateway function between the local code and FMM code. This
 * class takes care of multi-threads, organizing date, copying data to
 * the host, etc. It is a template to avoid the requirement that the
 * <tt>InteractionFun_t</tt> match the data type of
 * <tt>VecContainer</tt> (since template typedefs are not legal in
 * C++).
 */
template<typename T>
class VesInteraction
{
  public:
    ///The pointer-type to the external FMM function handle.
    typedef void(*InteractionFun_t)(const T*, const T*, size_t, T*, void*);    
 
    /** 
     * @param interaction_handle The function pointer to the FMM
     * code. When set to <tt>NULL</tt>, no interaction is performed.
     * @param num_threads The expected number of threads this class
     * need to take care of.
     */
    explicit VesInteraction(InteractionFun_t interaction_handle = NULL,
        int num_threads = omp_get_max_threads());
    ~VesInteraction();

    /** 
     * The function called inside the local code to perform
     * interaction. Each thread called this method to perform
     * interaction with other threads and/or other MPI processes.
     *
     * @param position The Cartesian coordinate of the points.
     * @param density The density at each point.
     * @param potential (return value) The calculated potential on the points. 
     * @param usr_ptr The used defined pointer, that may be required by the FMM code.
     *
     * @return Enum type defined in enums.h indicating the outcome of interaction.
     */
    template<typename VecContainer>
    InteractionReturn operator()(const VecContainer &position, 
        VecContainer &density, VecContainer &potential,
        void* usr_ptr = NULL) const;

  private:
    InteractionFun_t interaction_handle_;
    int num_threads_;
    size_t* each_thread_np_;
    size_t* each_thread_idx_;

    mutable size_t np_;
    mutable size_t containers_capacity_;

    mutable T* all_pos_;
    mutable T* all_den_;
    mutable T* all_pot_;
    
    size_t getCpyDestIdx(size_t this_thread_np) const;
    void checkContainersSize() const;
};

#include "VesInteraction.cc"

#endif // _VESINTERACTION_H_
