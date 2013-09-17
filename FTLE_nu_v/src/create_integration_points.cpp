#include "type_definitions.hpp"

int create_integration_points(double x, double y, double vx, double vy, double dx, double dy, double dvx, double dvy,
                              double1d& X1, double1d& X2, double1d& X3, double1d& X4, double1d& X5, double1d& X6, double1d& X7, double1d& X8)
{
    X1[0]=x+dx;
    X1[1]=y;
    X1[2]=vx;
    X1[3]=vy;

    X2[0]=x-dx;
    X2[1]=y;
    X2[2]=vx;
    X2[3]=vy;

    X3[0]=x;
    X3[1]=y+dy;
    X3[2]=vx;
    X3[3]=vy;

    X4[0]=x;
    X4[1]=y-dy;
    X4[2]=vx;
    X4[3]=vy;

    X5[0]=x;
    X5[1]=y;
    X5[2]=vx+dvx;
    X5[3]=vy;

    X6[0]=x;
    X6[1]=y;
    X6[2]=vx-dvx;
    X6[3]=vy;

    X7[0]=x;
    X7[1]=y;
    X7[2]=vx;
    X7[3]=vy+dvy;

    X8[0]=x;
    X8[1]=y;
    X8[2]=vx;
    X8[3]=vy-dvy;

    return 0;
}
