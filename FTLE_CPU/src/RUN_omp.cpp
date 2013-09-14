//General libraries
#include <time.h>
#include <iostream> // cout
#include <math.h>
#include <fstream> // files and printf
#include <stdio.h>
#include <string.h>
#include <sstream> // write on string

//Custom headers
#include "Lagrangian_points.h"
#include "global_var.hpp"

//GSL libraries
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>

//Boost libraries
#include <boost/numeric/odeint.hpp>
#include <boost/program_options.hpp>

//Custom headers
#include "type_definitions.hpp"

//OMP libraries
#include "omp.h"

using namespace boost::numeric::odeint;

typedef double value_type;
typedef std::vector< double > state_type;
//FIXME l'aggiunta di custom_error_checker sballa tutto quando non voglio usare il thrust
typedef runge_kutta_fehlberg78< state_type > error_stepper_type;

/* Functions declarations */
struct Lagrangian_points* Lagr(double mu);
double Omega (double x, double y, double mu);
int id_calc(int *arr);
int create_ic_vector (double c_max, double c_min, int nc, std::vector<double>& c0);
int create_ic_vector (double C0, double dc, std::vector<double>& c0);
int create_missing_vector(int id, double1d& x_0, double1d& y_0, double1d& vx_0, double1d& vy_0, double1d& e_0, bool4d& filter_integration, int4d& counter,  double Ki );
int launch_matlab( char *file_name);
int write_file(char* file_header, char* file_name, int id, double1d x_0, double1d y_0, double1d vx_0, double1d vy_0, double1d e_0, double4d ftle);

void er3bp( const state_type &x , state_type &dxdt , const double  t  );
void cr3bp_tt( const state_type &x , state_type &dxdt , const double  t  );

struct collision_check
{
    bool &m_collision_flag;
    double m_Tf;

    collision_check( bool &collision_flag , double Tf) : m_collision_flag(collision_flag), m_Tf(Tf) {}

    void operator () (const state_type &x, double &t)
    {
        if ( ((sqrt( (x[0]-mu)*(x[0]-mu)+x[1]*x[1])<d1) || (sqrt( (x[0]-1+mu)*(x[0]-1+mu)+x[1]*x[1])<d2)) && m_collision_flag==0)
        {
            m_collision_flag=1;
            t=m_Tf;
        }
    }
};

struct intersection_check
{
    int &m_count;
    double m_Tf;
    double m_intersection;
    double &m_previous;

    intersection_check( int &intersection_count, double Tf, double &previous_value, double intersection_value) :
        m_count(intersection_count), m_Tf(Tf), m_previous(previous_value), m_intersection(intersection_value) {}

    void operator () (const state_type &x, double &t)
    /* Computing intersections in both directions */
    {
        if (flags[0]==2)
        {
            if ((x[0]>m_intersection && m_previous<m_intersection) || (x[0]<m_intersection && m_previous>m_intersection))
            {
                m_count++;
            }
            m_previous=x[0];
        }
        else
        {
            if ((x[1]>m_intersection && m_previous<m_intersection) || (x[1]<m_intersection && m_previous>m_intersection))
            {
                m_count++;
            }
            m_previous=x[1];
        }
        if ( m_count==n_iterations )
        {
            t=m_Tf;
        }
    }
};

struct collision_intersection_check
{
    bool &m_collision_flag;
    int &m_count;
    double m_Tf;
    double m_intersection;
    double &m_previous;

    collision_intersection_check( bool &collision_flag, int &intersection_count , double Tf, double &previous_value, double intersection_value) :
        m_collision_flag(collision_flag), m_count(intersection_count), m_Tf(Tf), m_previous(previous_value), m_intersection(intersection_value) {}

    void operator () (const state_type &x, double &t)
    {
        if ( ((sqrt( (x[0]-mu)*(x[0]-mu)+x[1]*x[1])<d1) || (sqrt( (x[0]-1+mu)*(x[0]-1+mu)+x[1]*x[1])<d2)) && m_collision_flag==0)
        {
            m_collision_flag=1;
            t=m_Tf;
        }
        if (flags[0]==2)
        {
            if ((x[0]>m_intersection && m_previous<m_intersection) || (x[0]<m_intersection && m_previous>m_intersection))
            {
                m_count++;
            }
            m_previous=x[0];
        }
        else
        {
            if ((x[1]>m_intersection && m_previous<m_intersection) || (x[1]<m_intersection && m_previous>m_intersection))
            {
                m_count++;
            }
            m_previous=x[1];
        }
        if ( m_count==n_iterations )
        {
            t=m_Tf;
        }
    }
};

int RUN (double t0)
{
    double tempo;

    // Computation parameters
    double Tf=t0+DT;
    double Ki;
    Ki=1+ecc*cos(t0);

    int id=id_calc(flags);

    double4d ftle(nx, double3d(ny, double2d(n1, double1d(n2,0) ) ));

    bool4d filter_integration(nx, bool3d(ny, bool2d(n1, bool1d(n2,0) ) ));
    bool4d filter_ftle(nx, bool3d(ny, bool2d(n1, bool1d(n2,0) ) ));
    bool4d filter_collisions(nx, bool3d(ny, bool2d(n1, bool1d(n2,0) ) ));

    int4d counter(nx, int3d(ny, int2d(n1, int1d(n2,0) ) ));

    pp4d pf(nx, pp3d(ny, pp2d(n1, pp1d(n2) ) ));

    printf("sys dim=%i, %i, %i, %i\n", (int)ftle.size(), (int)ftle[1].size(), (int)ftle[1][1].size(), (int)ftle[1][1][1].size());

    double1d x_0(nx);
    double1d y_0(ny);
    double1d vx_0(nvx);
    double1d vy_0(nvy);
    double1d e_0(ne);
    char file_name[150]=""; //String that will become the file name
    char s_temp[50]="";
    sprintf(file_name,"ftle_ell_mu=%.4f_ecc=%.2f",mu,ecc);
    char file_header[300]=""; //String that will be printed in the file
    sprintf(file_header,"mu=%.12f\necc=%.12f\nDT=%.2f\nt0=%.2f\nn_frames=%i\nid=%i\nL=%.12f\nd1=%.12f\nd2=%.12f\nTh=%.12f\n",mu,ecc,DT,t0,n_frames,id,L,d1,d2,Th);

//Create vectors of initial conditions
    //Create x_0
    /*
    if x is a visualization variable the initial condition vector can be computed from x_max, x_min and nx.
    If it is a fixed variable, determine dx by comparison and compute initial condition vector using dx and the middle value X_0.
    */
    if (flags[0]==1)
    {
        create_ic_vector(x_max, x_min, nx, x_0);
        sprintf(s_temp,"nx=%i",nx);
        strcat(file_name,s_temp);
        strcat(file_header,s_temp);
        strcat(file_header,"\n");
    }
    else
    {
        if (flags[0]==2)
        {
            sprintf(s_temp,"x0=%.2f",X_0);
            strcat(file_name,s_temp);
            sprintf(s_temp,"x_0=%.6f",X_0);
            strcat(file_header,s_temp);
            strcat(file_header,"\n");
            if (flags[1]==1)
            {
                dx=dy;
            }
            else
            {
                if (flags[2]==1)
                {
                    dx=dvx;
                }
                else
                {
                    if (flags[3]==1)
                    {
                        dx=dvy;
                    }
                    else
                    {
                        dx=de;
                    }
                }
            }
            create_ic_vector(X_0, dx, x_0);
        }
    }

//Create y_0
    if (flags[1]==1)
    {
        create_ic_vector(y_max, y_min, ny, y_0);
        sprintf(s_temp,"ny=%i",ny);
        strcat(file_name,s_temp);
        strcat(file_header,s_temp);
        strcat(file_header,"\n");
    }
    else
    {
        if (flags[1]==2)
        {
            sprintf(s_temp,"y0=%.2f",Y_0);
            strcat(file_name,s_temp);
            sprintf(s_temp,"y_0=%.6f",Y_0);
            strcat(file_header,s_temp);
            strcat(file_header,"\n");
            if (flags[0]==1)
            {
                dy=dx;
            }
            else
            {
                if (flags[2]==1)
                {
                    dy=dvx;
                }
                else
                {
                    if (flags[3]==1)
                    {
                        dy=dvy;
                    }
                    else
                    {
                        dy=de;
                    }
                }
            }
            create_ic_vector(Y_0, dy, y_0);
        }
    }

//Create vx_0
    if (flags[2]==1)
    {
        create_ic_vector(vx_max, vx_min, nvx, vx_0);
        sprintf(s_temp,"nvx=%i",nvx);
        strcat(file_name,s_temp);
        strcat(file_header,s_temp);
        strcat(file_header,"\n");
    }
    else
    {
        if (flags[2]==2)
        {
            sprintf(s_temp,"vx0=%.2f",VX_0);
            strcat(file_name,s_temp);
            sprintf(s_temp,"vx_0=%.6f",VX_0);
            strcat(file_header,s_temp);
            strcat(file_header,"\n");
            if (flags[0]==1)
            {
                dvx=dx;
            }
            else
            {
                if (flags[1]==1)
                {
                    dvx=dy;
                }
                else
                {
                    if (flags[3]==1)
                    {
                        dvx=dvy;
                    }
                    else
                    {
                        dvx=de;
                    }
                }
            }
            create_ic_vector(VX_0, dvx, vx_0);
        }
    }

//Create vy_0
    if (flags[3]==1)
    {
        create_ic_vector(vy_max, vy_min, nvy, vy_0);
        sprintf(s_temp,"nvy=%i",nvy);
        strcat(file_name,s_temp);
        strcat(file_header,s_temp);
        strcat(file_header,"\n");
    }
    else
    {
        if (flags[3]==2)
        {
            sprintf(s_temp,"vy0=%.2f",VY_0);
            strcat(file_name,s_temp);
            sprintf(s_temp,"vy_0=%.6f",VY_0);
            strcat(file_header,s_temp);
            strcat(file_header,"\n");
            if (flags[0]==1)
            {
                dvy=dx;
            }
            else
            {
                if (flags[1]==1)
                {
                    dvy=dy;
                }
                else
                {
                    if (flags[2]==1)
                    {
                        dvy=dvx;
                    }
                    else
                    {
                        dvy=de;
                    }
                }
            }
            create_ic_vector(VY_0, dvy, vy_0);
        }
    }

//Create e_0
    if (flags[4]==1)
    {
        create_ic_vector(e_max, e_min, ne, e_0);
        sprintf(s_temp,"ne=%i",ne);
        strcat(file_name,s_temp);
        strcat(file_header,s_temp);
        strcat(file_header,"\n");
    }
    else
    {
        if (flags[4]==2)
        {
            sprintf(s_temp,"e0=%.2f",E_0);
            strcat(file_name,s_temp);
            sprintf(s_temp,"e_0=%.6f",E_0);
            strcat(file_header,s_temp);
            strcat(file_header,"\n");
            if (flags[0]==1)
            {
                de=dx;
            }
            else
            {
                if (flags[1]==1)
                {
                    de=dy;
                }
                else
                {
                    if (flags[2]==1)
                    {
                        de=dvx;
                    }
                    else
                    {
                        de=dvy;
                    }
                }
            }
            create_ic_vector(E_0, de, e_0);
        }
    }

// Create missing vector
    if (flags[0]==0)
    {
        std::cout<<"Please provide a value for both x and y";    /* x and y must always be provided as gridded or fixed variables*/
        return 1;
    }
    if (flags[1]==0)
    {
        std::cout<<"Please provide a value for both x and y";    /* x and y must always be provided as gridded or fixed variables*/
        return 1;
    }

    int dim;
    dim=create_missing_vector(id, x_0, y_0, vx_0, vy_0, e_0, filter_integration, counter, Ki);

    double perc=((double)dim/n_tot)*100;
    printf("Points to be integrated=%i\t %.1f %% n_tot\n",dim,perc);

    double steps;
    filter_ftle=filter_integration;
    //Integration
    clock_t start=clock();
    double in_step=1e-3;
    if (DT<0) {in_step=-1e-3;}
    int n_threads;
    time_t t_in=time(0);
    printf("number of threads=%i\n",n_cores);
    int failed_count=0;
    int done_count=0;
    double done_perc;
    #ifndef _WIN32
    system("setterm -cursor off");
    #endif // _WIN32
    omp_set_nested(1);
    #pragma omp parallel for num_threads (n_cores)
        for (int i=0; i<nx; i++)
        {
            #pragma omp parallel for num_threads (n_cores)
            for (int j=0; j<ny; j++)
            {
                for (int k=0; k<n1; k++)
                {
                    for (int l=0; l<n2; l++)
                    {
                        pf[i][j][k][l].x=0;
                        pf[i][j][k][l].y=0;
                        pf[i][j][k][l].vx=0;
                        pf[i][j][k][l].vy=0;
                        pf[i][j][k][l].e=0;
                        if (filter_integration[i][j][k][l])
                        {
                            std::vector<value_type> X(4,0);
                            int c=counter[i][j][k][l];
                            filter_collisions[i][j][k][l]=0;
                            bool collision_flag=filter_collisions[i][j][k][l];
                            X[0]=x_0[i];
                            X[1]=y_0[j];
                            if (flags[2]==0)
                            {
                                X[2]=vx_0[c];
                                X[3]=vy_0[k];
                            }
                            if (flags[3]==0)
                            {
                                X[2]=vx_0[k];
                                X[3]=vy_0[c];
                            }
                            if (flags[4]==0)
                            {
                                X[2]=vx_0[k];
                                X[3]=vy_0[l];
                            }
                            if (field_type=="FTLE")
                            {
                                if (distance_flag)
                                {
                                    steps=integrate_adaptive(make_controlled(abs_tol, rel_tol, error_stepper_type()), cr3bp_tt, X, t0, Tf, in_step, collision_check(collision_flag, Tf));
                                }
                                else    //FIXME check for actual number of intersection for each point
                                {
                                    steps=integrate_adaptive(make_controlled(abs_tol, rel_tol, error_stepper_type()), cr3bp_tt, X, t0, Tf, in_step);
                                }
                            }
                            else
                            {
                                int intersection_count=0;
                                double previous_value;
                                double intersection_value;
                                if (flags[0]==2){intersection_value=X_0;
                                previous_value=X_0;}
                                if (flags[1]==2){intersection_value=Y_0;
                                previous_value=Y_0;}
                                if (distance_flag)
                                {
                                    steps=integrate_adaptive(make_controlled(abs_tol, rel_tol, error_stepper_type()), er3bp, X, t0, Tf, in_step, collision_intersection_check(collision_flag,intersection_count, Tf, previous_value, intersection_value));
                                    if(intersection_count<n_iterations){failed_count++;}
                                }
                                else
                                {
                                    steps=integrate_adaptive(make_controlled(abs_tol, rel_tol, error_stepper_type()), er3bp, X, t0, Tf, in_step, intersection_check(intersection_count, Tf, previous_value, intersection_value));
                                    if(intersection_count<n_iterations){failed_count++;}
                                }
                            }

                            pf[i][j][k][l].x=X[0];
                            pf[i][j][k][l].y=X[1];
                            pf[i][j][k][l].vx=X[2];
                            pf[i][j][k][l].vy=X[3];
                            pf[i][j][k][l].e=0.5*((X[2]*X[2])+(X[3]*X[3]))-Omega(X[0],X[1],mu)/(1+ecc*cos(Tf));
                            filter_collisions[i][j][k][l]=collision_flag;

                            #pragma omp critical
                            {
                                done_count++;
                                done_perc=((double)done_count/dim)*100;
                                if ( floor(done_perc) != floor(done_perc+0.1) )
                                {
                                printf("progress: %.1f %% \r",done_perc);
                                }
                            }

                        }
                        /*
                        Do not compute ftle on the border of the phase space
                        */
                        if (i==0 || i==nx-1 || j==0 || j==ny-1 || k==0 || k==n1-1 || l==0 || l==n2-1)
                        {
                            filter_ftle[i][j][k][l]=0;
                        }
                        /*
                        Do not compute ftle for points in the neighbourhood of a not integrated point
                        */
                        else
                        {
                            if (filter_integration[i][j][k][l]==0 || filter_collisions[i][j][k][l])
                            {
                                filter_ftle[i][j][k][l]=0;
                                filter_ftle[i+1][j][k][l]=0;
                                filter_ftle[i-1][j][k][l]=0;
                                filter_ftle[i][j+1][k][l]=0;
                                filter_ftle[i][j-1][k][l]=0;
                                filter_ftle[i][j][k+1][l]=0;
                                filter_ftle[i][j][k-1][l]=0;
                                filter_ftle[i][j][k][l+1]=0;
                                filter_ftle[i][j][k][l-1]=0;
                            }
                        }
                    }
                }
            }
        }
    #ifndef _WIN32
    system("setterm -cursor on");
    #endif // _WIN32
    clock_t end_t=clock();
    tempo=((double)(end_t-start))/CLOCKS_PER_SEC;
    time_t t_fin=time(0);
    time_t t_tot=t_fin-t_in;
    std::cout<<"\nintegration time "<<t_tot<<" s\n";
    if(failed_count){printf("%i points haven\'t reached %i intersections\n",failed_count,n_iterations);}
    sprintf(s_temp, "integration_time=%i\nn_cores=%i\n", t_tot, n_cores);
    strcat(file_header, s_temp);

    //FTLE computation

    gsl_matrix *dphi=gsl_matrix_alloc (4,4);
    gsl_matrix *C=gsl_matrix_alloc (4,4);
    double l_max;
    gsl_vector *eval = gsl_vector_alloc (4);
    gsl_eigen_symm_workspace * w = gsl_eigen_symm_alloc (4);

    for (int i=0; i<nx; i++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int k=0; k<n1; k++)
            {
                for (int l=0; l<n2; l++)
                {
                    if (filter_ftle[i][j][k][l])
                    {
                        gsl_matrix_set(dphi,0,0,(pf[i+1][j][k][l].x-pf[i-1][j][k][l].x)/(2*dx));
                        gsl_matrix_set(dphi,0,1,(pf[i][j+1][k][l].x-pf[i][j-1][k][l].x)/(2*dy));
                        gsl_matrix_set(dphi,1,0,(pf[i+1][j][k][l].y-pf[i-1][j][k][l].y)/(2*dx));
                        gsl_matrix_set(dphi,1,1,(pf[i][j+1][k][l].y-pf[i][j-1][k][l].y)/(2*dy));

                        if (flags[2]==0) //k-->vy   l-->e
                        {
                            gsl_matrix_set(dphi,0,2,(pf[i][j][k+1][l].x-pf[i][j][k-1][l].x)/(2*dvy));
                            gsl_matrix_set(dphi,0,3,(pf[i][j][k][l+1].x-pf[i][j][k][l-1].x)/(2*de));

                            gsl_matrix_set(dphi,1,2,(pf[i][j][k+1][l].y-pf[i][j][k-1][l].y)/(2*dvy));
                            gsl_matrix_set(dphi,1,3,(pf[i][j][k][l+1].y-pf[i][j][k][l-1].y)/(2*de));

                            gsl_matrix_set(dphi,2,0,(pf[i+1][j][k][l].vy-pf[i-1][j][k][l].vy)/(2*dx));
                            gsl_matrix_set(dphi,2,1,(pf[i][j+1][k][l].vy-pf[i][j-1][k][l].vy)/(2*dy));
                            gsl_matrix_set(dphi,2,2,(pf[i][j][k+1][l].vy-pf[i][j][k-1][l].vy)/(2*dvy));
                            gsl_matrix_set(dphi,2,3,(pf[i][j][k][l+1].vy-pf[i][j][k][l-1].vy)/(2*de));

                            gsl_matrix_set(dphi,3,0,(pf[i+1][j][k][l].e-pf[i-1][j][k][l].e)/(2*dx));
                            gsl_matrix_set(dphi,3,1,(pf[i][j+1][k][l].e-pf[i][j-1][k][l].e)/(2*dy));
                            gsl_matrix_set(dphi,3,2,(pf[i][j][k+1][l].e-pf[i][j][k-1][l].e)/(2*dvy));
                            gsl_matrix_set(dphi,3,3,(pf[i][j][k][l+1].e-pf[i][j][k][l-1].e)/(2*de));
                        }

                        if (flags[3]==0)
                        {
                            gsl_matrix_set(dphi,0,2,(pf[i][j][k+1][l].x-pf[i][j][k-1][l].x)/(2*dvx));
                            gsl_matrix_set(dphi,0,3,(pf[i][j][k][l+1].x-pf[i][j][k][l-1].x)/(2*de));

                            gsl_matrix_set(dphi,1,2,(pf[i][j][k+1][l].y-pf[i][j][k-1][l].y)/(2*dvx));
                            gsl_matrix_set(dphi,1,3,(pf[i][j][k][l+1].y-pf[i][j][k][l-1].y)/(2*de));

                            gsl_matrix_set(dphi,2,0,(pf[i+1][j][k][l].vx-pf[i-1][j][k][l].vx)/(2*dx));
                            gsl_matrix_set(dphi,2,1,(pf[i][j+1][k][l].vx-pf[i][j-1][k][l].vx)/(2*dy));
                            gsl_matrix_set(dphi,2,2,(pf[i][j][k+1][l].vx-pf[i][j][k-1][l].vx)/(2*dvx));
                            gsl_matrix_set(dphi,2,3,(pf[i][j][k][l+1].vx-pf[i][j][k][l-1].vx)/(2*de));

                            gsl_matrix_set(dphi,3,0,(pf[i+1][j][k][l].e-pf[i-1][j][k][l].e)/(2*dx));
                            gsl_matrix_set(dphi,3,1,(pf[i][j+1][k][l].e-pf[i][j-1][k][l].e)/(2*dy));
                            gsl_matrix_set(dphi,3,2,(pf[i][j][k+1][l].e-pf[i][j][k-1][l].e)/(2*dvx));
                            gsl_matrix_set(dphi,3,3,(pf[i][j][k][l+1].e-pf[i][j][k][l-1].e)/(2*de));
                        }

                        if (flags[4]==0)
                        {
                            gsl_matrix_set(dphi,0,2,(pf[i][j][k+1][l].x-pf[i][j][k-1][l].x)/(2*dvx));
                            gsl_matrix_set(dphi,0,3,(pf[i][j][k][l+1].x-pf[i][j][k][l-1].x)/(2*dvy));

                            gsl_matrix_set(dphi,1,2,(pf[i][j][k+1][l].y-pf[i][j][k-1][l].y)/(2*dvx));
                            gsl_matrix_set(dphi,1,3,(pf[i][j][k][l+1].y-pf[i][j][k][l-1].y)/(2*dvy));

                            gsl_matrix_set(dphi,2,0,(pf[i+1][j][k][l].vx-pf[i-1][j][k][l].vx)/(2*dx));
                            gsl_matrix_set(dphi,2,1,(pf[i][j+1][k][l].vx-pf[i][j-1][k][l].vx)/(2*dy));
                            gsl_matrix_set(dphi,2,2,(pf[i][j][k+1][l].vx-pf[i][j][k-1][l].vx)/(2*dvx));
                            gsl_matrix_set(dphi,2,3,(pf[i][j][k][l+1].vx-pf[i][j][k][l-1].vx)/(2*dvy));

                            gsl_matrix_set(dphi,3,0,(pf[i+1][j][k][l].vy-pf[i-1][j][k][l].vy)/(2*dx));
                            gsl_matrix_set(dphi,3,1,(pf[i][j+1][k][l].vy-pf[i][j-1][k][l].vy)/(2*dy));
                            gsl_matrix_set(dphi,3,2,(pf[i][j][k+1][l].vy-pf[i][j][k-1][l].vy)/(2*dvx));
                            gsl_matrix_set(dphi,3,3,(pf[i][j][k][l+1].vy-pf[i][j][k][l-1].vy)/(2*dvy));
                        }

                        /* gsl_blas_dgemm(CblasTrans/CblasNoTrans,CblasTrans/CblasNoTrans,alpha,A,B,beta,C)*/ // alpha*A*B+beta*C
                        gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,dphi,dphi,0,C);

                        gsl_eigen_symm (C, eval, w);

                        gsl_vector_set(eval,0,fabs(gsl_vector_get(eval,0)));
                        gsl_vector_set(eval,1,fabs(gsl_vector_get(eval,1)));
                        gsl_vector_set(eval,2,fabs(gsl_vector_get(eval,2)));
                        gsl_vector_set(eval,3,fabs(gsl_vector_get(eval,3)));
                        l_max=gsl_vector_max(eval);
                        ftle[i][j][k][l]=1/fabs(DT)*log(sqrt(l_max));
                    }
                }
            }
        }
    }

    /*
    Output file structure:

    coordinate1   coordinate2   (coordinate3)  (coordinate4)    ftle

    Number of coordinates varies with the number of gridded variables. Coordinates always appear in this order: x, y, vx, vy, e.
    */

    sprintf(s_temp,"t0=%.2f",t0);
    strcat(file_name,s_temp);
    time_t now = time(0);
    strftime(s_temp,50,"_%H%M%S",localtime(&now));
    strcat(file_name,s_temp);
    strftime(s_temp,50,"computation_date=\'%Y-%m-%d %H:%M:%S\'\n",localtime(&now));
    strcat(file_header, s_temp);
    sprintf(s_temp,"method=\'%s\'\n",field_type.c_str());
    strcat(file_header,s_temp);
    if (field_type=="FILE")
    {
        sprintf(s_temp,"n_iterations=%i\n",n_iterations);
        strcat(file_header,s_temp);
    }
    if (name_flag)
    {
        sprintf(file_name,"%s",custom_file_name.c_str());
        if (n_frames>1)
        {
            sprintf(s_temp,"_t0=%.2f",t0);
            strcat(file_name,s_temp);
        }
    }
    strcat(file_name,".txt");
    sprintf(s_temp,"### Beginning of data ###\n");
    strcat(file_header,s_temp);

    write_file(file_header, file_name, id, x_0, y_0, vx_0, vy_0, e_0, ftle);

    if (matlab_flag)
    {
        launch_matlab(file_name);
    }
    return 0;
}
