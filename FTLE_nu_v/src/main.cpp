#include <vector>
#include <boost/numeric/odeint.hpp>
#include <math.h>
#include "type_definitions.hpp"
#include "Configuration.hpp"
#include <fstream>
//GSL libraries
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>

using namespace boost::numeric::odeint;

typedef double value_type;
typedef std::vector< double > state_type;
//FIXME l'aggiunta di custom_error_checker sballa tutto quando non voglio usare il thrust
typedef runge_kutta_fehlberg78< state_type > error_stepper_type;

void er3bp( const state_type &x , state_type &dxdt , const double  t  );

int create_integration_points(double x, double y, double vx, double vy,
                            double dx, double dy, double dvx, double dvy,
                            double1d& X1, double1d& X2, double1d& X3, double1d& X4,
                            double1d& X5, double1d& X6, double1d& X7, double1d& X8);

int create_ic_vector (double c_max, double c_min, int nc, double1d& c0);

double calc_t0(double nu);

int main ()
{
    printf("n_nu=%i\tn_v=%i\n",n_nu,n_v);
    create_ic_vector(nu_max, nu_min, n_nu, nu_0);
    create_ic_vector(v_max, v_min, n_v, v_0);
    //printf("li ho creati\n");


    gsl_matrix *C=gsl_matrix_alloc (4,4);
    double l_max;
    gsl_vector *eval = gsl_vector_alloc (4);
    gsl_eigen_symm_workspace * w = gsl_eigen_symm_alloc (4);

    double2d ftle (n_nu, double1d(n_v,0) );

    for (int i=0; i<n_nu; i++)
    {   double nu;
        nu=nu_0[i];
        t0=calc_t0(nu);
        //printf("nu = %.6f \n",nu);
        //printf("t0 = %.6f\n",t0);
        printf("%i\n",i); //why it dows not work if I add \r?
        for (int j=0; j<n_v; j++)
        {
            double v,x,y,vx,vy,v_E_adim,v_tilde;
            v_tilde=sqrt(GM*(1+ecc*cos(t0))/(a_jup*(1-ecc*ecc)));
            v_E_adim=v_E/v_tilde;
            v=v_0[j]+v_E_adim;
            //printf("v_tilde  = %.6f \n",v_tilde);
            //printf("v_E_adim = %.6f \n",v_E_adim);
            //printf("v = %.6f\n",v);
            x=R*cos(nu);
            y=R*sin(nu);
            vx=-v*sin(nu);
            vy=v*cos(nu);

            create_integration_points(x, y, vx, vy, dx, dy, dvx, dvy, X1, X2, X3, X4, X5, X6, X7, X8);

            steps=integrate_adaptive(make_controlled(abs_tol, rel_tol, error_stepper_type()), er3bp, X1, t0, t0+T, in_step);
            steps=integrate_adaptive(make_controlled(abs_tol, rel_tol, error_stepper_type()), er3bp, X2, t0, t0+T, in_step);
            steps=integrate_adaptive(make_controlled(abs_tol, rel_tol, error_stepper_type()), er3bp, X3, t0, t0+T, in_step);
            steps=integrate_adaptive(make_controlled(abs_tol, rel_tol, error_stepper_type()), er3bp, X4, t0, t0+T, in_step);
            steps=integrate_adaptive(make_controlled(abs_tol, rel_tol, error_stepper_type()), er3bp, X5, t0, t0+T, in_step);
            steps=integrate_adaptive(make_controlled(abs_tol, rel_tol, error_stepper_type()), er3bp, X6, t0, t0+T, in_step);
            steps=integrate_adaptive(make_controlled(abs_tol, rel_tol, error_stepper_type()), er3bp, X7, t0, t0+T, in_step);
            steps=integrate_adaptive(make_controlled(abs_tol, rel_tol, error_stepper_type()), er3bp, X8, t0, t0+T, in_step);

            gsl_matrix *dphi=gsl_matrix_alloc (4,4);
            gsl_matrix_set(dphi,0,0,(X1[0]-X2[0])/(2*dx));
            gsl_matrix_set(dphi,0,1,(X3[0]-X4[0])/(2*dy));
            gsl_matrix_set(dphi,0,2,(X5[0]-X6[0])/(2*dvx));
            gsl_matrix_set(dphi,0,3,(X7[0]-X8[0])/(2*dvy));

            gsl_matrix_set(dphi,1,0,(X1[1]-X2[1])/(2*dx));
            gsl_matrix_set(dphi,1,1,(X3[1]-X4[1])/(2*dy));
            gsl_matrix_set(dphi,1,2,(X5[1]-X6[1])/(2*dvx));
            gsl_matrix_set(dphi,1,3,(X7[1]-X8[1])/(2*dvy));

            gsl_matrix_set(dphi,2,0,(X1[2]-X2[2])/(2*dx));
            gsl_matrix_set(dphi,2,1,(X3[2]-X4[2])/(2*dy));
            gsl_matrix_set(dphi,2,2,(X5[2]-X6[2])/(2*dvx));
            gsl_matrix_set(dphi,2,3,(X7[2]-X8[2])/(2*dvy));

            gsl_matrix_set(dphi,3,0,(X1[3]-X2[3])/(2*dx));
            gsl_matrix_set(dphi,3,1,(X3[3]-X4[3])/(2*dy));
            gsl_matrix_set(dphi,3,2,(X5[3]-X6[3])/(2*dvx));
            gsl_matrix_set(dphi,3,3,(X7[3]-X8[3])/(2*dvy));

            gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,dphi,dphi,0,C);

            gsl_eigen_symm (C, eval, w);

            gsl_vector_set(eval,0,fabs(gsl_vector_get(eval,0)));
            gsl_vector_set(eval,1,fabs(gsl_vector_get(eval,1)));
            gsl_vector_set(eval,2,fabs(gsl_vector_get(eval,2)));
            gsl_vector_set(eval,3,fabs(gsl_vector_get(eval,3)));
            l_max=gsl_vector_max(eval);
            ftle[i][j]=1/fabs(T)*log(sqrt(l_max));
/*
            done_perc=(double)(n_v*i+j)/(n_nu*n_v)*100;
            //printf("done_perc=%f\n",done_perc);
            if ( floor(done_perc) != floor(done_perc+1) )
                {
                printf("progress: %.1f %% \r",done_perc);
                }*/
        }
    }

    char file_name[150]=""; //String that will become the file name
    char s_temp[50]="";
    sprintf(file_name,"ftle_nuv_mu=%.4f_ecc=%.2f_n=%ix%i_DT=%f",mu,ecc,n_nu,n_v,T);
    char file_header[300]=""; //String that will be printed in the file
    sprintf(file_header,"mu=%.12f\necc=%.12f\nDT=%.2f\nt0=%.2f\nn_nu=%i\nn_v=%i\n",mu,ecc,T,t0,n_nu,n_v);
    sprintf(s_temp,"### Beginning of data ###\n");
    strcat(file_header,s_temp);
    strcat(file_name,".txt");
    ftle_stream=fopen(file_name,"w");
    fprintf(ftle_stream,"%s",file_header);
    for (int i=0; i<n_nu; i++)
        {
            for (int j=0; j<n_v; j++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", nu_0[i], v_0[j], ftle[i][j]);
            }
        }
    fclose(ftle_stream);
    return 0;
}
