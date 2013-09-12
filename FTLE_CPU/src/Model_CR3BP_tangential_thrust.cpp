#include "global_var.hpp"
#include <math.h>
#include <vector>

typedef std::vector< double > state_type;

/* The rhs of x' = f(x) */
void cr3bp_tt( const state_type &x , state_type &dxdt , const double  t  )
{
    dxdt[0]=x[2];
    dxdt[1]=x[3];
    dxdt[2]=2*x[3]+(x[0]-((1-mu)*(x[0]+mu))/(((x[0]+mu)*(x[0]+mu)+x[1]*x[1])*sqrt((x[0]+mu)*(x[0]+mu)+x[1]*x[1]))-(mu*(x[0]-1+mu))/(((x[0]-1+mu)*(x[0]-1+mu)+x[1]*x[1])*sqrt((x[0]-1+mu)*(x[0]-1+mu)+x[1]*x[1])))
            +1/(p-q*t)/(sqrt(x[2]*x[2]+x[3]*x[3]))*((x[2]-x[1])*cos(t)-(x[3]+x[0])*sin(t));
    dxdt[3]=-2*x[2]+(x[1]-(1-mu)*x[1]/(((x[0]+mu)*(x[0]+mu)+x[1]*x[1])*sqrt((x[0]+mu)*(x[0]+mu)+x[1]*x[1]))-mu*x[1]/(((x[0]-1+mu)*(x[0]-1+mu)+x[1]*x[1])*sqrt((x[0]-1+mu)*(x[0]-1+mu)+x[1]*x[1])))
            +1/(p-q*t)/(sqrt(x[2]*x[2]+x[3]*x[3]))*((x[2]-x[1])*sin(t)+(x[3]+x[0])*cos(t));
}
