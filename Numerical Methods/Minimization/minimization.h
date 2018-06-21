#ifndef HAVE_MINIMIZATION_H
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<stdio.h>
#include<assert.h>
#include<math.h>


void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R);
void qr_gs_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);
void qr_gs_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B);

int newton_barebones(double f(gsl_vector* x, gsl_vector* grad, gsl_matrix* H), gsl_vector* x0, double eps);
int quasinewton(double f(gsl_vector* x), gsl_vector* x0, double dx, double eps);
void numeric_grad(double f(gsl_vector* x), gsl_vector* x, gsl_vector* grad, double dx);

#define HAVE_MINIMIZATION_H
#endif
