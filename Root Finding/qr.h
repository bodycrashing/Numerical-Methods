#ifndef HAVE_QR
#define HAVE_QR
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>

void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R);
void qr_back_sub(gsl_matrix* A, gsl_vector* R);
void qr_gs_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b,gsl_vector* x);
void qr_gs_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix *B);

#endif
