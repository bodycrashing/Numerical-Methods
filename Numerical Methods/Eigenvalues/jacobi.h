#ifndef JACOBI_H
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_blas.h>

int jacobi(gsl_matrix* A, gsl_vector* e, gsl_matrix* V, int *num_rot);
int jacobi_eig_by_eig(gsl_matrix* A, gsl_vector* e, gsl_matrix* V, int row_number, int *num_rot, int sorting);




//int Jacobi_eig_by_eig(gsl_matrix* A ,gsl_vector* x, gsl_matrix* V,int order,int num_of_eig);

int jacobi_classic(gsl_matrix* A, gsl_vector* e, gsl_matrix* V, int row_number, int *num_rot);

#define JACOBI_H
#endif
