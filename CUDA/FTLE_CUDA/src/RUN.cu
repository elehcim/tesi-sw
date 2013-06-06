
//General libraries
#include <time.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <sstream>

//Custom headers
#include "Phase_point.h"
#include "Lagrangian_points.h"
#include "global_var.hpp"
#include "er3bp_system.hpp"

//GSL libraries
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>


//Boost libraries
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/thrust/thrust_algebra.hpp>
#include <boost/numeric/odeint/external/thrust/thrust_operations.hpp>
#include <boost/numeric/odeint/external/thrust/thrust_resize.hpp>
#include <boost/program_options.hpp>

//Thrust libraries
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/remove.h>


using namespace boost::numeric::odeint;

typedef double value_type;
typedef thrust::device_vector< value_type > state_type;
typedef thrust::device_vector< size_t > index_vector_type;
//typedef runge_kutta_fehlberg78< state_type , value_type , state_type , value_type , thrust_algebra , thrust_operations > stepper_type;
typedef runge_kutta_dopri5< state_type , value_type , state_type , value_type , thrust_algebra , thrust_operations > stepper_type;
//typedef bulirsch_stoer< state_type , value_type , state_type , value_type , thrust_algebra , thrust_operations > stepper_type; //FIXME
//typedef runge_kutta_cash_karp54< state_type , value_type , state_type , value_type , thrust_algebra , thrust_operations > stepper_type;


//typedef runge_kutta_dopri5< value_type , state_type , value_type , state_type , thrust_algebra , thrust_operations > stepper;
typedef controlled_runge_kutta< stepper_type , custom_error_checker > controlled_stepper ;
//typedef dense_output_runge_kutta< controlled_stepper > dense_output_stepper;


typedef thrust::device_vector< int> int_vec_type;
    // Integrator precision
double abs_tol=1e-6;
double rel_tol=1e-6;

typedef std::vector<double> double1d;
typedef std::vector<double1d> double2d;
typedef std::vector<double2d> double3d;
typedef std::vector<double3d> double4d;

typedef std::vector<phase_point> pp1d;
typedef std::vector<pp1d> pp2d;
typedef std::vector<pp2d> pp3d;
typedef std::vector<pp3d> pp4d;

/* Functions declarations */
struct Lagrangian_points* Lagr(double mu);
double Omega (double x, double y, double mu);
int id_calc(int *arr);
int create_ic_vector (double c_max, double c_min, int nc, std::vector<double>& c0);
int create_ic_vector (double C0, double dc, std::vector<double>& c0);
int create_missing_vector(int id, double1d& x_0, double1d& y_0, double1d& vx_0, double1d& vy_0, double1d& e_0, double4d& filter_integration, double Ki );
int create_integration_vector (thrust::device_vector< double >& X_state, int dim, double1d x_0, double1d y_0, double1d vx_0, double1d vy_0, double1d e_0, double4d& filter_integration, double4d& filter_ftle);

int RUN (double t0)
{
double tempo;

  // Computation parameters
double Tf=t0+DT;
double Ki;
Ki=1+ecc*cos(t0);

int id=id_calc(flags);

double4d ftle(nx, double3d(ny, double2d(n1, double1d(n2,0) ) ));
double4d filter_integration(nx, double3d(ny, double2d(n1, double1d(n2,0) ) ));
double4d filter_ftle(nx, double3d(ny, double2d(n1, double1d(n2,0) ) ));
double4d filter_collisions(nx, double3d(ny, double2d(n1, double1d(n2,0) ) ));

pp4d pf(nx, pp3d(ny, pp2d(n1, pp1d(n2) ) ));

printf("ftle dim=%i, %i, %i, %i\n",ftle.size(),ftle[1].size(), ftle[1][1].size(), ftle[1][1][1].size());
//printf("id=%i\n",id);

double1d x_0(nx);
double1d y_0(ny);
double1d vx_0(nvx);
double1d vy_0(nvy);
double1d e_0(ne);
char file_name[100]=""; //String that will become the file name
char s_temp[50]="";
sprintf(file_name,"ftle_ell_mu=%.4f_ecc=%.2f",mu,ecc);
char file_header[200]=""; //String that will be printed in the file
sprintf(file_header,"mu=%.4f\necc=%.4f\nDT=%.2f\nt0=%.2f\nn_frames=%i\nid=%i\n",mu,ecc,DT,t0,n_frames,id);

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
    //printf("dx=%.4f\n",dx);
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
        if (flags[1]==1){dx=dy;}
        else{   if (flags[2]==1){dx=dvx;} //FIXME dimensioni diverse!
                else{   if (flags[3]==1){dx=dvy;} //FIXME dimensioni diverse!!
                        else {dx=de;}   //FIXME dimensioni diverse!!
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
        if (flags[0]==1){dy=dx;}
        else{   if (flags[2]==1){dy=dvx;} //FIXME dimensioni diverse!
                else{   if (flags[3]==1){dy=dvy;} //FIXME dimensioni diverse!!
                        else {dy=de;}   //FIXME dimensioni diverse!!
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
        if (flags[0]==1){dvx=dx;}   //FIXME dimensioni diverse!
        else{   if (flags[1]==1){dvx=dy;} //FIXME dimensioni diverse!
                else{   if (flags[3]==1){dvx=dvy;}
                        else {dvx=de;}   //FIXME dimensioni diverse!!
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
        if (flags[0]==1){dvy=dx;}   //FIXME dimensioni diverse!!
        else{   if (flags[1]==1){dvy=dy;} //FIXME dimensioni diverse!
                else{   if (flags[2]==1){dvy=dvx;}
                        else {dvy=de;}   //FIXME dimensioni diverse!!
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
        if (flags[0]==1){de=dx;} //FIXME dimensioni diverse!!
        else{   if (flags[1]==1){de=dy;} //FIXME dimensioni diverse!
                else{   if (flags[2]==1){de=dvx;} //FIXME dimensioni diverse!!
                        else {de=dvy;}   //FIXME dimensioni diverse!!
                    }
            }
        create_ic_vector(E_0, de, e_0);
    }
}

// Create missing vector
int c=0;
if (flags[0]==0){std::cout<<"Please provide a value for both x and y"; return 1;} /* x and y must always be provided as gridded or fixed variables*/
if (flags[1]==0){std::cout<<"Please provide a value for both x and y"; return 1;} /* x and y must always be provided as gridded or fixed variables*/

int dim;
dim=create_missing_vector(id, x_0, y_0, vx_0, vy_0, e_0, filter_integration, Ki);

    //Construct vector for integration
double perc=((double)dim/n_tot)*100;
printf("Points to be integrated=%i\t %.1f %% n_tot\n",dim,perc);

thrust::device_vector< double > X_state(5*dim,0); //5: x,y,vx,vy,e
create_integration_vector (X_state, dim, x_0, y_0, vx_0, vy_0, e_0, filter_integration, filter_ftle);

int_vec_type Collisions(dim,0);
er3bp_system sistema (dim, Collisions);
clock_t start=clock();
double steps;

//steps=integrate_adaptive( make_controlled(abs_tol , rel_tol, stepper_type() ), sistema , X_state , t0 , Tf , 1e-6);
//dense_output_stepper s(controlled_stepper(custom_controller( Collisions )));
//steps=integrate_adaptive( dense_output_stepper, sistema , X_state , t0 , Tf , 1e-6);
controlled_stepper s(Collisions);
steps=integrate_adaptive( s, sistema , X_state , t0 , Tf , 1e-6);

//std::cout<<steps<<'\n';
clock_t end_t=clock();
tempo=((double)(end_t-start))/CLOCKS_PER_SEC;
std::cout<<"integration time "<<tempo<<'\n';

    //Deconstruct integration results
start=clock();
c=0;
for (int i=0; i<nx; i++)
    {
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
                        if (Collisions[c]==0)
                        {
                            pf[i][j][k][l].x=X_state[c];
                            pf[i][j][k][l].y=X_state[c+dim];
                            pf[i][j][k][l].vx=X_state[c+2*dim];
                            pf[i][j][k][l].vy=X_state[c+3*dim];
                            pf[i][j][k][l].e=0.5*((X_state[c+2*dim]*X_state[c+2*dim])+
                                                  (X_state[c+3*dim]*X_state[c+3*dim]))-
                                                  Omega(X_state[c],X_state[c+dim],mu)/(1+ecc*cos(Tf));
                        }
                        filter_collisions[i][j][k][l]=Collisions[c];
                        c=c+1;
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
                        if (filter_integration[i][j][k][l]==0 || filter_collisions[i][j][k][l]==1)
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

end_t=clock();
tempo=((double)(end_t-start))/CLOCKS_PER_SEC;
std::cout<<"Deconstruction time "<<tempo<<'\n';

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
                    //ftle[i][j][k][l]=0;
                    //std::cout<<i<<'\t'<<j<<'\t'<<k<<'\t'<<l<<'\n';
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
FILE *ftle_stream;
sprintf(s_temp,"t0=%.2f",t0);
strcat(file_name,s_temp);
strcat(file_name,".txt");
sprintf(s_temp,"### Beginning of FTLE data ###\n");
strcat(file_header,s_temp);
ftle_stream=fopen(file_name,"w");
fprintf(ftle_stream,file_header);

if (id==11220 || id==11202 || id==11022)
{
    for (int i=0; i<nx; i++)
    {
        for (int j=0; j<ny; j++)
        {
            fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", x_0[i], y_0[j], ftle[i][j][1][1]);
        }
    }
}

if(id==12120 || id==12102)
{
    for (int i=0; i<nx; i++)
    {
        for (int j=0; j<nvx; j++)
        {
            fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", x_0[i], vx_0[j], ftle[i][1][j][1]);
        }
    }
}

if(id==12210)
{
    for (int i=0; i<nx; i++)
        {
            for (int j=0; j<nvy; j++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", x_0[i], vy_0[j], ftle[i][1][1][j]);
            }
        }
}

if(id==12012)
{
    for (int i=0; i<nx; i++)
        {
            for (int j=0; j<nvy; j++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", x_0[i], vy_0[j], ftle[i][1][j][1]);
            }
        }
}

if(id==12201 || id==12021)
{
    for (int i=0; i<nx; i++)
        {
            for (int j=0; j<ne; j++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", x_0[i], e_0[j], ftle[i][1][1][j]);
            }
        }
}

if (id==21120 || id==21102)
{
    for (int i=0; i<ny; i++)
    {
        for (int j=0; j<nvx; j++)
        {
            fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", y_0[i], vx_0[j], ftle[1][i][j][1]);
        }
    }
}

if (id==21012)
{
    for (int i=0; i<ny; i++)
    {
        for (int j=0; j<nvy; j++)
        {
            fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", y_0[i], vy_0[j], ftle[1][i][j][1]);
        }
    }
}

if (id==21210)
{
    for (int i=0; i<ny; i++)
    {
        for (int j=0; j<nvy; j++)
        {
            fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", y_0[i], vy_0[j], ftle[1][i][1][j]);
        }
    }
}

if (id==21201 || id==21021)
{
    for (int i=0; i<ny; i++)
    {
        for (int j=0; j<ne; j++)
        {
            fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", y_0[i], e_0[j], ftle[1][i][1][j]);
        }
    }
}

if (id==22110)
{
    for (int i=0; i<nvx; i++)
    {
        for (int j=0; j<nvy; j++)
        {
            fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", vx_0[i], vy_0[j], ftle[1][1][i][j]);
        }
    }
}

if (id==22101)
{
    for (int i=0; i<nvx; i++)
    {
        for (int j=0; j<ne; j++)
        {
            fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", vx_0[i], e_0[j], ftle[1][1][i][j]);
        }
    }
}

if (id==22011)
{
    for (int i=0; i<nvy; i++)
    {
        for (int j=0; j<ne; j++)
        {
            fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", vx_0[i], e_0[j], ftle[1][1][i][j]);
        }
    }
}
    //3D Cases
if (id==11120 || id==11102)
{
    for (int i=0; i<nx; i++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int k=0; k<nvx; k++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\n", x_0[i], y_0[j], vx_0[k], ftle[i][j][k][1]);
            }
        }
    }
}

if (id==11210)
{
    for (int i=0; i<nx; i++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int k=0; k<nvy; k++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\n", x_0[i], y_0[j], vy_0[k], ftle[i][j][1][k]);
            }
        }
    }
}

if (id==11012)
{
    for (int i=0; i<nx; i++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int k=0; k<nvy; k++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\n", x_0[i], y_0[j], vy_0[k], ftle[i][j][k][1]);
            }
        }
    }
}

if (id==11021 || id==11201)
{
    for (int i=0; i<nx; i++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int k=0; k<ne; k++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\n", x_0[i], y_0[j], e_0[k], ftle[i][j][1][k]);
            }
        }
    }
}

if (id==12110)
{
    for (int i=0; i<nx; i++)
    {
        for (int j=0; j<nvx; j++)
        {
            for (int k=0; k<nvy; k++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\n", x_0[i], vx_0[j], vy_0[k], ftle[i][1][j][k]);
            }
        }
    }
}

if (id==12101)
{
    for (int i=0; i<nx; i++)
    {
        for (int j=0; j<nvx; j++)
        {
            for (int k=0; k<ne; k++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\n", x_0[i], vx_0[j], e_0[k], ftle[i][1][j][k]);
            }
        }
    }
}

if (id==12011)
{
    for (int i=0; i<nx; i++)
    {
        for (int j=0; j<nvy; j++)
        {
            for (int k=0; k<ne; k++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\n", x_0[i], vy_0[j], e_0[k], ftle[i][1][j][k]);
            }
        }
    }
}

if (id==12011)
{
    for (int i=0; i<nx; i++)
    {
        for (int j=0; j<nvy; j++)
        {
            for (int k=0; k<ne; k++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\n", x_0[i], vy_0[j], e_0[k], ftle[i][1][j][k]);
            }
        }
    }
}

if (id==21110)
{
    for (int i=0; i<ny; i++)
    {
        for (int j=0; j<nvx; j++)
        {
            for (int k=0; k<nvy; k++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\n", y_0[i], vx_0[j], vy_0[k], ftle[1][i][j][k]);
            }
        }
    }
}

if (id==21101)
{
    for (int i=0; i<ny; i++)
    {
        for (int j=0; j<nvx; j++)
        {
            for (int k=0; k<ne; k++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\n", y_0[i], vx_0[j], e_0[k], ftle[1][i][j][k]);
            }
        }
    }
}

if (id==21011)
{
    for (int i=0; i<ny; i++)
    {
        for (int j=0; j<nvy; j++)
        {
            for (int k=0; k<ne; k++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\n", y_0[i], vy_0[j], e_0[k], ftle[1][i][j][k]);
            }
        }
    }
}

 //4D case
if (id==11110)
{
    for (int i=0; i<nx; i++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int k=0; k<nvx; k++)
            {
                for (int l=0; l<nvy; l++)
                {
                    fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\t%.12f\n", x_0[i], y_0[j], vx_0[k], vy_0[l], ftle[i][j][k][l]);
                }
            }
        }
    }
}

if (id==11101)
{
    for (int i=0; i<nx; i++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int k=0; k<nvx; k++)
            {
                for (int l=0; l<ne; l++)
                {
                    fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\t%.12f\n", x_0[i], y_0[j], vx_0[k], e_0[l], ftle[i][j][k][l]);
                }
            }
        }
    }
}

if (id==11011)
{
    for (int i=0; i<nx; i++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int k=0; k<nvy; k++)
            {
                for (int l=0; l<ne; l++)
                {
                    fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\t%.12f\n", x_0[i], y_0[j], vy_0[k], e_0[l], ftle[i][j][k][l]);
                }
            }
        }
    }
}

fclose(ftle_stream);

X_state.clear();
X_state.shrink_to_fit();
return 0;
}
