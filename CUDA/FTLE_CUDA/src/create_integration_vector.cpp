#include "global_var.hpp"
#include <vector>
#include <thrust/device_vector.h>

typedef std::vector<double> double1d;
typedef std::vector<double1d> double2d;
typedef std::vector<double2d> double3d;
typedef std::vector<double3d> double4d;


int create_integration_vector (thrust::device_vector< double >& X_state, int dim, double1d x_0, double1d y_0, double1d vx_0, double1d vy_0, double1d e_0, double4d& filter_integration, double4d& filter_ftle)
{
    int c=0;
    for (int i=0; i<nx; i++)
        {
            for (int j=0; j<ny; j++)
            {
                for (int k=0; k<n1; k++)
                {
                    for(int l=0; l<n2; l++)
                    {
                        if (filter_integration[i][j][k][l])
                        {
                            X_state[c]=x_0[i];
                            X_state[c+dim]=y_0[j];
                            if (flags[2]==0){X_state[c+2*dim]=vx_0[c];
                                            X_state[c+3*dim]=vy_0[k];
                                            X_state[c+4*dim]=e_0[l];}
                            if (flags[3]==0){X_state[c+2*dim]=vx_0[k];
                                            X_state[c+3*dim]=vy_0[c];
                                            X_state[c+4*dim]=e_0[l];}
                            if (flags[4]==0){X_state[c+2*dim]=vx_0[k];
                                            X_state[c+3*dim]=vy_0[l];
                                            X_state[c+4*dim]=e_0[c];}
                            filter_ftle[i][j][k][l]=1;
                            c=c+1;
                        }
                    }
                }
            }
        }
    return 0;
}

int create_integration_vector_t (thrust::device_vector< double >& X_state, int dim, double1d x_0, double1d y_0, double1d vx_0, double1d vy_0, double1d e_0, double4d& filter_integration)
{
    int c=0;
    for (int i=0; i<nxt; i++)
        {
            for (int j=0; j<nyt; j++)
            {
                for (int k=0; k<n1t; k++)
                {
                    for(int l=0; l<n2t; l++)
                    {
                        if (filter_integration[i][j][k][l])
                        {
                            X_state[c]=x_0[i];
                            X_state[c+dim]=y_0[j];
                            if (flags_t[2]==0){X_state[c+2*dim]=vx_0[c];
                                            X_state[c+3*dim]=vy_0[k];
                                            X_state[c+4*dim]=e_0[l];}
                            if (flags_t[3]==0){X_state[c+2*dim]=vx_0[k];
                                            X_state[c+3*dim]=vy_0[c];
                                            X_state[c+4*dim]=e_0[l];}
                            if (flags_t[4]==0){X_state[c+2*dim]=vx_0[k];
                                            X_state[c+3*dim]=vy_0[l];
                                            X_state[c+4*dim]=e_0[c];}
                            c=c+1;
                        }
                    }
                }
            }
        }
    return 0;
}
