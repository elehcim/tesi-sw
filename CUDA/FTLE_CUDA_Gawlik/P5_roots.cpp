/*
This function computes all the five roots of a 5th-order polynomial.
It returns a pointer to the arry of roots, which alternatively stores real and imaginary parts
*/
#include <iostream>
#include <stdio.h>
#include <gsl/gsl_poly.h>
const double mu=0.1;
double *P5_roots(double *a)
{
//int i;
double z[10];
gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc (6);
gsl_poly_complex_solve (a, 6, w, z);
gsl_poly_complex_workspace_free (w);
/*
for (i = 0; i < 5; i++)
{
    printf ("z%d = %+.18f %+.18f\n", i, z[2*i], z[2*i+1]);
}
*/
return &z[0];
}
