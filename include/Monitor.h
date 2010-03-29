#ifndef _MONITOR_H_
#define _MONITOR_H_

#include "Surface.h"

enum MntrRetFlg {Continue  = 101, 
                 Terminate = 102};
                    
struct MntrOpts
{
    bool save_centers_;
    bool save_shapes_;
    int save_freq_;
    int area_inc_fac_;
};
    
template <typename T>
class Monitor
{
  public:
    Monitor(MntrOpts opts_in, DataIO<T> &fileIO_in);
    
    enum MntrRetFlg Quiz(Surface<T> &surf_in, int step_idx);

  private:
    MntrOpts opts;
    DataIO<T> &fileIO;
    T max_init_area_;
    Monitor(Monitor<T> &mnt);
};

template <typename T>
Monitor<T>::Monitor(MntrOpts opts_in, DataIO<T> &fileIO_in) :
    opts(opts_in),
    fileIO(fileIO_in),
    max_init_area_(-1)
{}

template <typename T>
enum MntrRetFlg Monitor<T>::Quiz(Surface<T> &surf_in, int step_idx)
{

    cout<<step_idx<<endl;
    if(opts.save_shapes_ && (step_idx + 1)%opts.save_freq_ == 0)
    {
        fileIO.Append(surf_in.x_.data_, surf_in.x_.GetDataLength());
    }
    else if(opts.save_centers_)
    {
        int cen_len = 3 * surf_in.params_.n_surfs_;
        T* cnts = (T*) malloc(3 * surf_in.params_.n_surfs_ * sizeof(T));
        fileIO.Append(surf_in.GetCenters(cnts), cen_len);
    }
    
    T m_area = surf_in.Area();

    if(max_init_area_ < 0)
        max_init_area_ = m_area;
    else
        if(isnan(m_area) || m_area > opts.area_inc_fac_ * max_init_area_)
        {
            cerr<<"The time stepper has diverged"<<endl;
            return(Terminate);
        }
    
    return(Continue);
}


#endif //_MONITOR_H_
