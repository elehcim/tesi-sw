#include <math.h>
#include <fstream>
double calc_t0(double nu)
{
    double M, E, E_new, delta_t, eps, n, f, omega;
    //double delta_t, eps, n, f, omega;
    extern double ecc, a_jup, GM, pi;
    omega=2*pi/(365*24*3600);
    delta_t=nu/omega;
    eps=0.001;
    n=sqrt(GM/(a_jup*a_jup*a_jup));
    M=n*(delta_t);
    E=M;
    E_new=E;
    //printf("ciclo E=%f\n",E);
    int i=0;
    do
    {
        E=E_new;
        E_new=E+(M-E+ecc*sin(E))/(1-ecc*cos(E));
        i=i+1;
        //printf("E    =   %.10e \n",E);
        //printf("E_new=   %.10e \n",E_new);
    }
    while(E_new-E>eps && i<10);
    //f=2*atan(sqrt((1+ecc)/(1-ecc))*tan(E_new/2));
    //printf("2st arg = %.10e \n",sqrt(1-ecc)*cos(E_new/2));
    //printf("1nd arg = %.10e \n",sqrt(1+ecc)*sin(E_new/2));
    f=2*atan2(sqrt(1+ecc)*sin(E_new/2),sqrt(1-ecc)*cos(E_new/2));
    if (f<0)
    {
        f=f+2*pi;
        printf("ho cambiato segno a f\n");
    }
    //printf("f=%f\n",f);
    return f;
}
