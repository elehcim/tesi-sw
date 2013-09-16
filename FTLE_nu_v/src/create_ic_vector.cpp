#include "type_definitions.hpp"

int create_ic_vector (double c_max, double c_min, int nc, double1d& c0)
{
    double dc=0;
    c0[0]=c_min;
    if (nc!=1){dc=(c_max-c_min)/(nc-1);
    for (int i=1;  i<nc; i++)
    {
        c0[i]=c0[i-1]+dc;
    }
    }
    return 0;
}
