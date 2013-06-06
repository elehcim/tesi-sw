#include <math.h>
double Omega (double x, double y, double mu)
{
    double r1=sqrt((x+mu)*(x+mu)+y*y);
    double r2=sqrt((x-1+mu)*(x-1+mu)+y*y);
    double om=0.5*(x*x+y*y)+(1-mu)/r1 + mu/r2 + 0.5*(mu*(1-mu));
    return om;
}
