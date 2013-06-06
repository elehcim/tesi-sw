#include "global_var.hpp"
#include <math.h>
#include <fstream>
#include <vector>
#include <thrust/device_vector.h>

double Omega (double x, double y, double mu);

typedef double value_type;
typedef thrust::device_vector< value_type > state_type;

typedef std::vector<double> double1d;
typedef std::vector<double1d> double2d;
typedef std::vector<double2d> double3d;
typedef std::vector<double3d> double4d;

int print_tracers(FILE* stream, state_type X, int dim, double Tf, double4d filter)
{
    int c=0;
    fprintf(stream, "### t=%.2f ###\n",Tf);
    for (int i=0; i<nxt; i++)
    {
        for (int j=0; j<nyt; j++)
        {
            for (int k=0; k<n1t; k++)
            {
                for (int l=0; l<n2t; l++)
                {
                    if (filter[i][j][k][l])
                    {
                        value_type x = X[c];
                        value_type y = X[c+dim];
                        value_type vx = X[c+2*dim];
                        value_type vy = X[c+3*dim];
                        value_type e = 0.5*((X[c+2*dim]*X[c+2*dim])+
                                                  (X[c+3*dim]*X[c+3*dim]))-
                                                  Omega(X[c],X[c+dim],mu)/(1+ecc*cos(Tf));
                        c++;
                        fprintf(stream, "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", x, y, vx, vy, e);
                    }
                }
            }
        }
    }
    return 0;
}
