#include <iostream>
#include <stdio.h>
#include <math.h>
#include "Lagrangian_points.h"
double *P5_roots(double *a);
struct Lagrangian_points* Lagr(double mu)
{
Lagrangian_points Lagrp;
double c[6]={-mu, 2*mu, -mu, 3-2*mu, mu-3, 1}; // ATTENTION!: Coefficients are listed from the lowest to the highest order!!
double *z;
z=P5_roots(c);
for (int i=0; i<5; i++)
    {
    if (z[2*i+1]==0)
        {
        Lagrp.xl1=1-mu-z[2*i];
        }
    }
c[0]=-mu; c[1]=-2*mu; c[2]=-mu; c[3]=3-2*mu; c[4]=3-mu; c[5]=1; //ATT!: order of coefficients (see above)
z=P5_roots(c);
for (int i=0; i<5; i++)
    {
        if (z[2*i+1]==0)
        {
        Lagrp.xl2=1-mu+z[2*i];
        }
    }
c[0]=mu-1; c[1]=2*mu-2; c[2]=mu-1; c[3]=1+2*mu; c[4]=2+mu; c[5]=1; //ATT!: order of coefficients (see above)
z=P5_roots(c);
for (int i=0; i<5; i++)
    {
        if (z[2*i+1]==0)
        {
        Lagrp.xl3=-mu-z[2*i];
        }
    }
Lagrp.xl4=0.5-mu;
Lagrp.xl5=Lagrp.xl4;
Lagrp.yl1=0;
Lagrp.yl2=0;
Lagrp.yl3=0;
Lagrp.yl4=sqrt(3.0)/2;
Lagrp.yl5=-sqrt(3.0)/2;
struct Lagrangian_points *p;
p=&Lagrp;
return p;
}
