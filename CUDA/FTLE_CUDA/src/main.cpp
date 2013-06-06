#include <math.h>
#include <iostream>
#include <vector>
#include "global_var.hpp"
//#include "GPU_info.cu"
typedef double value_type;
double t0=0.0;
int RUN (double t0);
int configuration_load ();
int tracers (double Ki, std::vector<double> obs_times);
//GPU_info();

int main ()
{
    //GPU_info();
    int error=configuration_load();
    double dt=2*M_PI/(n_frames);
    std::vector<double> obs_times (n_frames,0);
    if (!error)
        {
            for (int i=0; i<n_frames; i++) //Frame t0=0 = Frame t0=2*pi -> n_frame-1
            {
                std::cout<<"\nframe "<<i+1<<"\t t0="<<t0<<'\n';
                RUN(t0);
                std::cout<<"frame created\n";
                obs_times[i]=t0;
                t0=t0+dt;
            }
            if (tracers_flag)
            {
                std::cout<<"\nTracers computation\n";
                double Ki=1+ecc*cos(t0);
                tracers(Ki, obs_times);
            }
    return 0;
        }
    else return 1;

}
