#include "global_var.hpp"
#include <vector>
#include <math.h>
#include <fstream>

#include "type_definitions.hpp"


#ifdef _WIN32
#define copysign  _copysign
#define isnan _isnan
#endif // WIN32

double Omega (double x, double y, double mu);

int create_missing_vector(int id, double1d& x_0, double1d& y_0, double1d& vx_0, double1d& vy_0, double1d& e_0, bool4d& filter_integration, double Ki )
{
    int c=0;
    double missing_coordinate;
    double1d missing_coordinate_vec(n_tot);
    for (int i=0; i<nx; i++)
        {
            for (int j=0; j<ny; j++)
            {
                for (int k=0; k<n1; k++)
                {
                    for(int l=0; l<n2; l++)
                    {
                        if (flags[2]==0) /* missing coordinate is vx */
                        {
                            missing_coordinate=-copysign(1.00, y_0[j])*sqrt(2*Omega(x_0[i],y_0[j],mu)/Ki+2*e_0[l]-vy_0[k]*vy_0[k]);
                        }
                        if (flags[3]==0) /* missing coordinate is vy */
                        {
                            missing_coordinate=copysign(1.00, x_0[i])*sqrt(2*Omega(x_0[i],y_0[j],mu)/Ki+2*e_0[l]-vx_0[k]*vx_0[k]);
                        }
                        if (flags[4]==0) /* missing coordinate is e */
                        {
                            missing_coordinate=0.5*(vx_0[k]*vx_0[k]+vy_0[l]*vy_0[l])-Omega(x_0[i],y_0[j],mu)/Ki;
                        }

                        if (!isnan(missing_coordinate))
                        {
                            if (( (id==22011 || id==22101 || id==22110) && !(i!=1 && j!=1)) ||
                                ( (id==21021 || id==21201 || id==21210) && !(i!=1 && k!=1)) ||
                                ( (id==21012 || id==21102 || id==21120) && !(i!=1 && l!=1)) ||
                                ( (id==12021 || id==12201 || id==12210) && !(j!=1 && k!=1)) ||
                                ( (id==12012 || id==12102 || id==12120) && !(j!=1 && l!=1)) ||
                                ( (id==11022 || id==11202 || id==11220) && !(k!=1 && l!=1)) )
                            {
                                missing_coordinate_vec[c]=missing_coordinate;
                                filter_integration[i][j][k][l]=1;
                                c=c+1;
                            }
                            if (vis_count>2)
                            {
                                missing_coordinate_vec[c]=missing_coordinate;
                                filter_integration[i][j][k][l]=1;
                                c=c+1;
                            }
                        }
                    }
                }
            }
        }
        if (flags[2]==0){vx_0=missing_coordinate_vec;}
        if (flags[3]==0){vy_0=missing_coordinate_vec;}
        if (flags[4]==0){e_0=missing_coordinate_vec;}
    return c;
}

// Create missing vector - OpenMP approach
int create_missing_vector(int id, double1d& x_0, double1d& y_0, double1d& vx_0, double1d& vy_0, double1d& e_0, bool4d& filter_integration, int4d& counter,  double Ki )
{
    int c=0;
    double missing_coordinate;
    double1d missing_coordinate_vec(n_tot);
    for (int i=0; i<nx; i++)
        {
            for (int j=0; j<ny; j++)
            {
                for (int k=0; k<n1; k++)
                {
                    for(int l=0; l<n2; l++)
                    {
                        if (flags[2]==0) /* missing coordinate is vx */
                        {
                            missing_coordinate=-copysign(1.00, y_0[j])*sqrt(2*Omega(x_0[i],y_0[j],mu)/Ki+2*e_0[l]-vy_0[k]*vy_0[k]);
                        }
                        if (flags[3]==0) /* missing coordinate is vy */
                        {
                            missing_coordinate=copysign(1.00, x_0[i])*sqrt(2*Omega(x_0[i],y_0[j],mu)/Ki+2*e_0[l]-vx_0[k]*vx_0[k]);
                        }
                        if (flags[4]==0) /* missing coordinate is e */
                        {
                            missing_coordinate=0.5*(vx_0[k]*vx_0[k]+vy_0[l]*vy_0[l])-Omega(x_0[i],y_0[j],mu)/Ki;
                        }

                        if (!isnan(missing_coordinate))
                        {
                            if (( (id==22011 || id==22101 || id==22110) && !(i!=1 && j!=1)) ||
                                ( (id==21021 || id==21201 || id==21210) && !(i!=1 && k!=1)) ||
                                ( (id==21012 || id==21102 || id==21120) && !(i!=1 && l!=1)) ||
                                ( (id==12021 || id==12201 || id==12210) && !(j!=1 && k!=1)) ||
                                ( (id==12012 || id==12102 || id==12120) && !(j!=1 && l!=1)) ||
                                ( (id==11022 || id==11202 || id==11220) && !(k!=1 && l!=1)) )
                            {
                                missing_coordinate_vec[c]=missing_coordinate;
                                filter_integration[i][j][k][l]=1;
                                counter[i][j][k][l]=c;
                                c++;
                            }
                            if (vis_count>2)
                            {
                                missing_coordinate_vec[c]=missing_coordinate;
                                filter_integration[i][j][k][l]=1;
                                counter[i][j][k][l]=c;
                                c++;
                            }
                        }
                    }
                }
            }
        }
        if (flags[2]==0){vx_0=missing_coordinate_vec;}
        if (flags[3]==0){vy_0=missing_coordinate_vec;}
        if (flags[4]==0){e_0=missing_coordinate_vec;}
    return c;
}
