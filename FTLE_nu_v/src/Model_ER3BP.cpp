//#include "global_var.hpp"
#include <math.h>
#include <vector>

typedef std::vector< double > state_type;
extern double mu, ecc;
/* The rhs of x' = f(x) */
void er3bp( const state_type &x , state_type &dxdt , const double  t  )
{
    dxdt[0]=x[2];
    dxdt[1]=x[3];
    dxdt[2]=2*x[3]+(x[0]-((1-mu)*(x[0]+mu))/(((x[0]+mu)*(x[0]+mu)+x[1]*x[1])*sqrt((x[0]+mu)*(x[0]+mu)+x[1]*x[1]))-(mu*(x[0]-1+mu))/(((x[0]-1+mu)*(x[0]-1+mu)+x[1]*x[1])*sqrt((x[0]-1+mu)*(x[0]-1+mu)+x[1]*x[1])))/(1+ecc*cos(t));
    dxdt[3]=-2*x[2]+(x[1]-(1-mu)*x[1]/(((x[0]+mu)*(x[0]+mu)+x[1]*x[1])*sqrt((x[0]+mu)*(x[0]+mu)+x[1]*x[1]))-mu*x[1]/(((x[0]-1+mu)*(x[0]-1+mu)+x[1]*x[1])*sqrt((x[0]-1+mu)*(x[0]-1+mu)+x[1]*x[1])))/(1+ecc*cos(t));
}
