#ifndef _MONITOR_H_
#define _MONITOR_H_

#include "Surface.h"

template<typename SurfContainer>
class Monitor
{
  private:
    typedef typename SurfContainer::value_type value_type;
    typedef typename SurfContainer::Sca_t::device_type device_type;
    value_type time_hor_;
    bool save_flag_;
    value_type save_stride_;
    DataIO IO;
    value_type A0, V0;
    
  public:
    Monitor(const Parameters<value_type> &params);
    bool operator()(const SurfContainer &state, 
        value_type &t, value_type &dt, void* user = NULL);
};

// enum MntrRetFlg {Continue  = 101, 
//                  Terminate = 102};
                    
// struct MntrOpts
// {
//     bool save_centers_;
//     bool save_shapes_;
//     int save_freq_;
//     int area_inc_fac_;
//     bool verbose;
// };
    
// template <typename T>
// class Monitor
// {
//   public:
//     Monitor(MntrOpts opts_in, DataIO<T> &fileIO_in);
    
//     enum MntrRetFlg Quiz(Surface<T> &surf_in, int step_idx);

//   private:
//     MntrOpts opts;
//     DataIO<T> &fileIO;
//     T max_init_area_;
//     Monitor(Monitor<T> &mnt);
// };

// template <typename T>
// Monitor<T>::Monitor(MntrOpts opts_in, DataIO<T> &fileIO_in) :
//     opts(opts_in),
//     fileIO(fileIO_in),
//     max_init_area_(-1)
// {}

// template <typename T>
// enum MntrRetFlg Monitor<T>::Quiz(Surface<T> &surf_in, int step_idx)
// {

//     if(opts.verbose) 
//     {
//         cout<<" Time step   : "<<step_idx<<endl;
//         cout<<"   . n_surfs_ : "<<surf_in.params_.n_surfs_<<endl;
//     }

//     if(opts.save_shapes_ && (step_idx + 1)%opts.save_freq_ == 0)
//     {
//         fileIO.Append(surf_in.x_.data_, surf_in.x_.GetDataLength());
//     }
//     else if(opts.save_centers_)
//     {
//         int cen_len = 3 * surf_in.params_.n_surfs_;
//         T* cnts = (T*) malloc(3 * surf_in.params_.n_surfs_ * sizeof(T));
//         fileIO.Append(surf_in.GetCenters(cnts), cen_len);
//     }
    
//     T m_area = surf_in.Area();
//     if(opts.verbose) cout<<"  . Max area : "<<m_area<<endl;
//     if(max_init_area_ < 0)
//         max_init_area_ = m_area;
//     else
//         if(isnan(m_area) || m_area > opts.area_inc_fac_ * max_init_area_)
//         {
//             cerr<<" ============================== "<<endl;
//             cerr<<" The time stepper has diverged. "<<endl;
//             cerr<<" ============================== "<<endl;
//             return(Terminate);
//         }
    

//     if(step_idx == 2) surf_in.Resize(1);
//     if(step_idx == 100) surf_in.Resize(2);

//     return(Continue);
// }


#endif //_MONITOR_H_
