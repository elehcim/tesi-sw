#include <vector>
#include <fstream>

typedef std::vector<double> double1d;

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

int create_ic_vector (double C0, double dc, double1d& c0)
{
    c0[0]=C0-dc;
    c0[1]=C0;
    c0[2]=C0+dc;
    return 0;
}
