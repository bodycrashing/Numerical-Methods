#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include<stdio.h>
#include<stdlib.h>
#define FMT "%7.3f"
#define max_print 6

void matrix_print(gsl_matrix *A){
	for(int r=0;r<A->size1;r++){
		for(int c=0;c<A->size2;c++)printf(FMT,gsl_matrix_get(A,r,c));
		printf("\n");}}


void vector_print(gsl_vector *v){
	for(int i=0;i<v->size;i++) printf(FMT,gsl_vector_get(v,i));
	printf("\n");
	}


int main(int argc, char const *argv[]) {

  double A_row1[] = {6.13, -2.90, 5.86};
  double A_row2[] = {8.08, -6.31, -3.89};
  double A_row3[] = {4.13, 1.00, 0.19};
  double b_data[] = {6.23, 5.37, 2.29};
  int n = sizeof(b_data)/sizeof(b_data[0]);

  gsl_vector* b = gsl_vector_alloc(n);
  gsl_vector* y=gsl_vector_calloc(n);
  gsl_vector* x=gsl_vector_calloc(n);
  gsl_vector* M1 = gsl_vector_alloc(n);
  gsl_vector* M2 = gsl_vector_alloc(n);
  gsl_vector* M3 = gsl_vector_alloc(n);
  gsl_matrix* M = gsl_matrix_alloc(n,n);
  gsl_matrix* A=gsl_matrix_calloc(n,n);

  for(int i=0; i<n; i++){
    gsl_vector_set(b,i,b_data[i]);
    gsl_vector_set(M1,i,A_row1[i]);
    gsl_vector_set(M2,i,A_row2[i]);
    gsl_vector_set(M3,i,A_row3[i]);
  }
  gsl_matrix_set_row(M,0,M1);
  gsl_matrix_set_row(M,1,M2);
  gsl_matrix_set_row(M,2,M3);

	gsl_matrix_memcpy(A,M);
	printf("\nThe Matrix A:\n");
	matrix_print(M);

	printf("\nThe required solution b:\n");
	vector_print(b);

	/*Solving the linear equations*/
	gsl_linalg_HH_solve(M,b,x);
	printf("\nThe computed solution x is found to be:\n");
	vector_print(x);

	/*Testing that Ax=b*/
	gsl_blas_dgemv(CblasNoTrans,1,A,x,0,y);
	printf("\nTesting that Ax=b\n");
	vector_print(y);

gsl_matrix_free(A);
gsl_matrix_free(M);
gsl_vector_free(M1);
gsl_vector_free(M2);
gsl_vector_free(M3);
gsl_vector_free(b);
gsl_vector_free(x);
gsl_vector_free(y);

  return 0;
  }
