#ifndef HAVE_ROOT_H
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<stdio.h>
#include<assert.h>
#include<math.h>

void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R);
void qr_gs_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);
void qr_gs_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B);

int Newton(void f(gsl_vector* x, gsl_vector* fx), gsl_vector* x0, double dx, double eps);
int Newton_Jacobian(void f(gsl_vector* x, gsl_vector* fx, gsl_matrix* J), gsl_vector* x0, double dx, double eps);
int Newton_Refined(void f(gsl_vector* x, gsl_vector* fx, gsl_matrix* J), gsl_vector* x0, double dx, double eps);

#define HAVE_ROOT_H
#endif
