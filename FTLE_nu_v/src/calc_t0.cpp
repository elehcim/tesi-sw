#include <math.h>
#include <fstream>
double calc_t0(double nu)
{
    double delta_t, eps, n, f, M, E, E_new, omega;
    extern double ecc, a_jup, GM, pi;
    omega=2*pi/(365*24*3600);
    delta_t=nu/omega;
    eps=0.01;
    n=sqrt(GM/(a_jup*a_jup*a_jup));
    M=n*(delta_t);
    E=M;
    printf("ciclo E=%f\n",E);
    int i=0;
    do
    {
        E_new=E+(M-E+ecc*sin(E))/(1-ecc*cos(E));
        i=i+1;
    }
    while(E_new-E>eps && i<500);
    //f=2*atan(sqrt((1+ecc)/(1-ecc))*tan(E_new/2));
    f=2*atan2(sqrt(1-ecc)*cos(E_new/2),sqrt(1+ecc)*sin(E_new/2));
    if (f<0)
    {
        f=f+2*pi;
    }
    printf("f=%f\n",f);
    return f;
}
