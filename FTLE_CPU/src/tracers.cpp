#include "global_var.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;

typedef std::vector<double> double1d;
typedef std::vector<double1d> double2d;
typedef std::vector<double2d> double3d;
typedef std::vector<double3d> double4d;

typedef double value_type;
typedef std::vector< double > state_type;
typedef runge_kutta_cash_karp54< state_type > error_stepper_type;


int id_calc (int *arr);
int create_ic_vector (double c_max, double c_min, int nc, double1d& c0);
int create_missing_vector_t(int id, double1d& x_0, double1d& y_0, double1d& vx_0, double1d& vy_0, double1d& e_0, double4d& filter_integration, double Ki );
int print_tracers(FILE* stream, state_type X, int dim, double Tf, double4d filter);


int tracers (double Ki, std::vector<double> obs_times)
{
    double1d x0t(nxt);
    double1d y0t(nyt);
    double1d vx0t(nvxt);
    double1d vy0t(nvyt);
    double1d e0t(net);

    if (flags_t[0]==1){create_ic_vector(x_maxt, x_mint, nxt, x0t);}
    else {if (flags_t[0]==2) {x0t[0]=X_0t;} }
    if (flags_t[1]==1){create_ic_vector(y_maxt, y_mint, nyt, y0t);}
    else {if (flags_t[1]==2) {y0t[0]=Y_0t;} }
    if (flags_t[2]==1){create_ic_vector(vx_maxt, vx_mint, nvxt, vx0t);}
    else {if (flags_t[2]==2) {vx0t[0]=VX_0t;} }
    if (flags_t[3]==1){create_ic_vector(vy_maxt, vy_mint, nvyt, vy0t);}
    else {if (flags_t[3]==2) {vy0t[0]=VY_0t;} }
    if (flags_t[4]==1){create_ic_vector(e_maxt, e_mint, net, e0t);}
    else {if (flags_t[4]==2) {e0t[0]=E_0t;} }

    int id_tracers=id_calc(flags_t);
    int dim_tracers;

    double4d filter_tracers(nxt, double3d(nyt, double2d(n1t, double1d(n2t,0) ) ));
    dim_tracers=create_missing_vector_t(id_tracers, x0t, y0t, vx0t, vy0t, e0t, filter_tracers, Ki);

    double abs_tol=1e-6;
    double rel_tol=1e-6;
    char file_name[50]="tracers.txt";
    FILE *tracers_stream;
    tracers_stream=fopen(file_name,"w");
    //std::cout<<"Tracers integration\n";
    print_tracers(tracers_stream, X, dim_tracers,obs_times[0], filter_tracers);
    for (int i=0; i< n_frames-1; i++)
    {
        integrate_adaptive(make_controlled(abs_tol, rel_tol, stepper_type()), er3bp(),X_tracers, obs_times[i],obs_times[i+1],1e-6);
        print_tracers(tracers_stream, X_tracers, dim_tracers, obs_times[i+1], filter_tracers);
    }
    fclose(tracers_stream);
    return 0;
}

