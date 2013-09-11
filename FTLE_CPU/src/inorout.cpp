#include "Lagrangian_points.h"
#include <math.h>
#include <stdio.h>
#include <iostream>
Lagrangian_points Lagr (double mu);
bool inorout(double x_0, double y_0, double mu)
{
    double d1, d2, d1l3, d2l2;
    Lagrangian_points Lagrp;
    Lagrp=Lagr(mu);
    //printf("xl2=%f\t xl3=%f\n",Lagrp.xl2, Lagrp.xl3);
    d1=sqrt((x_0+mu)*(x_0+mu)+(y_0*y_0));
    d2=sqrt((x_0-1+mu)*(x_0-1+mu)+(y_0*y_0));
    d1l3=-mu-Lagrp.xl3;
    d2l2=Lagrp.xl2-1+mu;

    if ((d1>d1l3) && (d2>d2l2))
    {
        //std::cout<<"fuori";
        return 0;
    }
    else
    {
        return 1;
    }
}
