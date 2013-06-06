//General libraries

#include <time.h>
#include <iostream>
#include <math.h>
#include <stdio.h>

//Custom headers
#include "Phase_point.h"
#include "Lagrangian_points.h"

//GSL libraries
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <fstream>

//Boost libraries
#include <boost/numeric/odeint.hpp>

//Thrust libraries
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/remove.h>

#include <boost/numeric/odeint/external/thrust/thrust_algebra.hpp>
#include <boost/numeric/odeint/external/thrust/thrust_operations.hpp>
#include <boost/numeric/odeint/external/thrust/thrust_resize.hpp>

using namespace boost::numeric::odeint;

typedef double value_type;
typedef thrust::device_vector< value_type > state_type;
typedef thrust::device_vector< size_t > index_vector_type;
/*
typedef runge_kutta_fehlberg78< state_type , value_type ,
                                state_type , value_type ,
                                thrust_algebra , thrust_operations > stepper_type;
*/
typedef runge_kutta_dopri5< state_type , value_type ,
                            state_type , value_type ,
                            thrust_algebra , thrust_operations > stepper_type;

const int n=150;
double t0=2*M_PI/3;
double DT=2.0;

// Integrator precision
double abs_tol=1e-6;
double rel_tol=1e-6;
const value_type mu=0.1;
const value_type ecc=0.04;


struct er3bp_system
{

    struct er3bp_functor
    {
        double m_t;
        er3bp_functor( double t ) : m_t(t){}

        template< class T >
        __host__ __device__
        void operator()( T tpl ) const
        {
            value_type x = thrust::get< 0 >( tpl );
            value_type y = thrust::get< 1 >( tpl );
            value_type vx = thrust::get< 2 >( tpl );
            value_type vy = thrust::get< 3 >( tpl );
            thrust::get< 4 >(tpl)=vx; // set tuple's sixth element to vx
            thrust::get< 5 >(tpl)=vy;
            thrust::get< 6 >(tpl)=2*vy+(x-((1-mu)*(x+mu))/pow((x+mu)*(x+mu)+y*y,1.5)-
                                    (mu*(x-1+mu))/pow((x-1+mu)*(x-1+mu)+y*y,1.5))/(1+ecc*cos(m_t));
            thrust::get< 7 >(tpl)=-2*vx+(y-(1-mu)*y/pow((x+mu)*(x+mu)+y*y,1.5)-
                                    mu*y/pow((x-1+mu)*(x-1+mu)+y*y,1.5))/(1+ecc*cos(m_t));
        }
    };


    er3bp_system( size_t N ): m_N( N ) { }

    template< class State , class Deriv >
    void operator()( const State &x , Deriv &dxdt , value_type t ) const
    {
        thrust::for_each(
        // Create each tuple
                thrust::make_zip_iterator( thrust::make_tuple(
                        boost::begin( x ) ,             // x
                        boost::begin( x ) + m_N ,       // y
                        boost::begin( x ) + 2 * m_N ,   // vx
                        boost::begin( x ) + 3 * m_N,    // vy
                        boost::begin( dxdt ) ,
                        boost::begin( dxdt ) + m_N ,
                        boost::begin( dxdt ) + 2 * m_N,
                        boost::begin( dxdt ) + 3 * m_N ) ) ,
                thrust::make_zip_iterator( thrust::make_tuple(
                        boost::begin( x ) + m_N ,
                        boost::begin( x ) + 2 * m_N ,
                        boost::begin( x ) + 3 * m_N ,
                        boost::begin( x ) + 4 * m_N ,
                        boost::begin( dxdt ) + m_N ,
                        boost::begin( dxdt ) + 2 * m_N ,
                        boost::begin( dxdt ) + 3 * m_N ,
                        boost::begin( dxdt ) + 4 * m_N ) ) ,
                er3bp_functor(t) );
    }

size_t m_N;
};


struct Lagrangian_points* Lagr(double mu);
double Omega (double x, double y, double mu);


const int nx=n;
const int ny=3;
const int nvx=n;
const int ne=3;
const int n_tot=nx*ny*nvx*ne;

typedef  bool bool_matrix4d[nx][ny][nvx][ne];
typedef double matrix4d[nx][ny][nvx][ne];

bool_matrix4d filter;
bool_matrix4d filter_ftle;
phase_point pf[nx][ny][nvx][ne];
matrix4d ftle;

int main ()
{
double tempo;

  // Computation parameters
double Tf=t0+DT;
//double ecc=0.0;
double Ki;
Ki=1+ecc*cos(0.0);

  // Initial conditions ranges
  // Assigned initial conditions
double mu=0.10;
struct Lagrangian_points Lagrp;

Lagrp=*Lagr(mu);

double C_L1;
C_L1=2*Omega(Lagrp.xl1,Lagrp.yl1,mu)/Ki; //FIXME devo mettere /Ki?
//std::cout<<C_L1<<'\n';
double E_0=-C_L1/2.0+0.03715;
printf("C_L1 = %.5f E_0 = %.5f Ki = %.5f",C_L1,E_0,Ki);
//std::cout<<E_0; // FIXME generalizzare!

//std::cout<<E_0<<'\n';
double Y_0=0;

//int nx=n;
double x_0_min=-0.8;
double x_0_max=-0.15;
double dx=(x_0_max-x_0_min)/(nx-1);
double x_0[nx];
x_0[0]=x_0_min;
//std::cout<<"x_0[0]="<<x_0[0]<<'\t';
for (int i=1; i<nx ; i++)
{
    x_0[i]=x_0[i-1]+dx;
    //std::cout<<"x_0["<<i<<"]="<<x_0[i]<<'\t';
}

//int nvx=n;
double vx_0_min=-2;
double vx_0_max=2;
double dvx=(vx_0_max-vx_0_min)/(nvx-1);
double vx_0[nvx];
vx_0[0]=vx_0_min;
for (int i=1; i<nvx; i++)
{
    vx_0[i]=vx_0[i-1]+dvx;
}

//int ny=3;
double dy=(dx+dvx)/2; // FIXME assunzione sparata a caso
double y_0[3]={Y_0-dy, Y_0, Y_0+dy}; //FIXME da generalizzare!

//int ne=3;
double de=dy;  //FIXME assunzione
double e_0[3]={E_0-de, E_0, E_0+de}; //FIXME da generalizzare!

double vy_0;

thrust::device_vector< double > X_state(5*n_tot,0); //5: x,y,vx,vy,e
/*
value_type X[n_tot];
value_type Y[n_tot];
value_type VX[n_tot];
value_type VY[n_tot];
value_type E[n_tot];
*/
    //Construct vector for integration
int c=0;
for (int i=0; i<nx; i++)
{
    for (int j=0; j<ny; j++)
    {
        for (int k=0; k<nvx; k++)
        {
            for (int l=0; l<ne; l++)
            {
                vy_0=-sqrt(2*Omega(x_0[i],y_0[j],mu)/Ki+2*e_0[l]-vx_0[k]*vx_0[k]);
                if (!((j!=1) && (l!=1)) && (vy_0<0))
                {
                    X_state[c]=x_0[i];
                    X_state[c+n_tot]=y_0[j];
                    X_state[c+2*n_tot]=vx_0[k];
                    X_state[c+3*n_tot]=vy_0;
                    X_state[c+4*n_tot]=e_0[l];
                    filter[i][j][k][l]=1;
                    filter_ftle[i][j][k][l]=1;
                    c=c+1;
                }
                else
                {
                    filter[i][j][k][l]=0;
                    filter_ftle[i][j][k][l]=0;
                }
            }
        }
    }
}
int dim=c;

for(int i=0;i<dim;i++)
{
    X_state[i+dim]=X_state[i+n_tot];
}
for(int i=0;i<dim;i++)
{
    X_state[i+2*dim]=X_state[i+2*n_tot];
}
for(int i=0;i<dim;i++)
{
    X_state[i+3*dim]=X_state[i+3*n_tot];
}
for(int i=0;i<dim;i++)
{
    X_state[i+4*dim]=X_state[i+4*n_tot];
}

//X_state_host.resize(5*dim);

X_state.resize(5*dim);

/*
for(int i=0; i<dim;i++)
{
    X_state_host[i]=X[i];
    X_state_host[i+dim]=Y[i];
    X_state_host[i+2*dim]=VX[i];
    X_state_host[i+3*dim]=VY[i];
    X_state_host[i+4*dim]=E[i];
}

thrust::device_vector< value_type > X_state=X_state_host;
//thrust::remove(X.begin()+dim+1,X.begin()+n_tot,0);
*/

er3bp_system sistema (dim);
clock_t start=clock();
double steps;
steps = integrate_adaptive(
        make_controlled(abs_tol , rel_tol, stepper_type() ),
        sistema ,
        X_state ,
        t0 ,
        Tf ,
        1e-6);

std::cout<<steps<<'\n';
clock_t end=clock();
tempo=((double)(end-start))/CLOCKS_PER_SEC;
std::cout<<"integration time "<<tempo<<'\n';

    //Deconstruct integration results
    start=clock();

c=0;
for (int i=0; i<nx; i++)
{
    for (int j=0; j<ny; j++)
    {
        for (int k=0; k<nvx; k++)
        {
            for (int l=0; l<ne; l++)
            {
                pf[i][j][k][l].x=0;
                pf[i][j][k][l].y=0;
                pf[i][j][k][l].vx=0;
                pf[i][j][k][l].vy=0;
                pf[i][j][k][l].e=0;
                if (filter[i][j][k][l])
                {
                    pf[i][j][k][l].x=X_state[c];
                    pf[i][j][k][l].y=X_state[c+dim];
                    pf[i][j][k][l].vx=X_state[c+2*dim];
                    pf[i][j][k][l].vy=X_state[c+3*dim];
                    pf[i][j][k][l].e=0.5*((X_state[c+2*dim]*X_state[c+2*dim])+
                                        (X_state[c+3*dim]*X_state[c+3*dim]))-
                                        Omega(X_state[c],X_state[c+dim],mu)/(1+ecc*cos(Tf));
                    c=c+1;
                }
                if (i==0 || i==nx || j==0 || j==ny || k==0 || k==nvx || l==0 || l==ne)
                    {
                    filter_ftle[i][j][k][l]=0;
                    }
                else{
                    if (filter[i][j][k][l]==0)
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
end=clock();
tempo=((double)(end-start))/CLOCKS_PER_SEC;
std::cout<<"Deconstruction time "<<tempo<<'\n';
    //FTLE computation
//matrix4d ftle;
gsl_matrix *dphi=gsl_matrix_alloc (4, 4);
gsl_matrix *C=gsl_matrix_alloc (4,4);
double l_max;
//double alpha=1;
gsl_vector *eval = gsl_vector_alloc (4);
//gsl_matrix *evec = gsl_matrix_alloc (4, 4);
gsl_eigen_symm_workspace * w = gsl_eigen_symm_alloc (4);
for (int i=0; i<nx; i++)
{
    for (int j=0; j<ny; j++)
    {
        for (int k=0; k<nvx; k++)
        {
            for (int l=0; l<ne; l++)
            {
                ftle[i][j][k][l]=0;
                if (filter_ftle[i][j][k][l])
                {
                gsl_matrix_set(dphi,0,0,(pf[i+1][j][k][l].x-pf[i-1][j][k][l].x)/(2*dx));
                gsl_matrix_set(dphi,0,1,(pf[i][j+1][k][l].x-pf[i][j-1][k][l].x)/(2*dy));
                gsl_matrix_set(dphi,0,2,(pf[i][j][k+1][l].x-pf[i][j][k-1][l].x)/(2*dvx));
                gsl_matrix_set(dphi,0,3,(pf[i][j][k][l+1].x-pf[i][j][k][l-1].x)/(2*de));

                gsl_matrix_set(dphi,1,0,(pf[i+1][j][k][l].y-pf[i-1][j][k][l].y)/(2*dx));
                gsl_matrix_set(dphi,1,1,(pf[i][j+1][k][l].y-pf[i][j-1][k][l].y)/(2*dy));
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
FILE *ftle_stream;
char s[100];
sprintf(s,"ftle_ell_n=%.0i_t=%.2f.txt",n,t0);
ftle_stream=fopen(s,"w");
for (int i=0; i<nx; i++)
{   fprintf(ftle_stream,"\n");
    for (int k=0; k<nvx; k++)
    {
        fprintf(ftle_stream,"%.12f\t",ftle[i][1][k][1]);
    }
}

fclose(ftle_stream);
return 0;
}
