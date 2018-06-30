#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#ifndef HAVE_RUNGE_KUTTA
#define HAVE_RUNGE_KUTTA

void rkstep23(void f(double x, gsl_vector *y, gsl_vector *dydx),
              double x, gsl_vector* yx, double h, gsl_vector* yxh, gsl_vector* err);

int ode_driver(void f(double x, gsl_vector* y, gsl_vector* dydx),
                              gsl_vector *xlist, gsl_matrix *ylist,
                              double b, double h, double acc, double eps, int max);


#endif
