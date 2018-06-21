#ifndef HAVE_MONTECARLO
#define HAVE_MONTECARLO
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "adaptive.h"

void randomx(int dim, gsl_vector* a, gsl_vector* b, gsl_vector* x);
void plain_monte(int dim, gsl_vector* a, gsl_vector* b, double f(gsl_vector* x),
                  int N, double *result, double *err);
double Adaptive2D(double f(double x, double y),double a, double b, double c(double x),
                          double d(double x), double acc, double eps, double *err);


#endif
