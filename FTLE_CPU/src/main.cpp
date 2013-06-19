#include <math.h>
#include <iostream>
#include <vector>
#include "global_var.hpp"
#include <stdio.h>



typedef double value_type;
int RUN (double t0);
int configuration_load ();

int main ()
{
    int error=configuration_load();
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
