#ifndef HAVE_SMP
#define HAVE_SMP
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include<omp.h>
#include<gsl/gsl_rng.h> //Used because RND() is npt thread safe


void randomx(int dim, gsl_vector* a, gsl_vector* b, gsl_vector* x, gsl_rng* r);

void plain_monte(int dim, gsl_vector* a, gsl_vector* b, double f(gsl_vector* x),
                  int N, double *result, double *err);

void plain_monte_OMP(int dim, gsl_vector* a, gsl_vector* b, double f(gsl_vector* x),
                  int N, double *result, double *err);

#endif
