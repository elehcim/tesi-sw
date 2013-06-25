#include <math.h>
#include <iostream>
#include <vector>
#include "global_var.hpp"
#include <stdio.h>

typedef double value_type;
int RUN (double t0);
int configuration_load (std::string config_file);
std::string command_line_load (int ac, char* av[]);
int main (int ac, char* av[])
{
    int error;
    std::string config_file;
    config_file=command_line_load(ac, av);
    error=configuration_load(config_file);
    double dt=(t0+tf)/(n_frames);
    if (!error)
    {
        for (int i=0; i<n_frames; i++) //Frame t0=0 = Frame t0=2*pi -> n_frame-1
        {
            std::cout<<"\nframe "<<i+1<<"\t t0="<<t0<<'\n';
            RUN(t0);
            std::cout<<"frame created\n";
            t0=t0+dt;
        }
        #ifdef _WIN32
        printf("Press Enter key to end");
        getchar();
        #endif // Only for Windows
        return 0;
    }
    else
    {
        #ifdef _WIN32
        printf("Press Enter key to end");
        getchar();
        #endif // Only for Windows
        return 1;
    }

}
