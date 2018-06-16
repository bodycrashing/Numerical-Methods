#include<assert.h>
#include <math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include"qr.h"
#include "least_squares.h"
#include "jacobi.h"

void lsfit(int m, double f(int i,double x),
	gsl_vector* x, gsl_vector* y, gsl_vector* dy,
	gsl_vector* c, gsl_matrix* S)
{
int n = (*x).size;

gsl_matrix *A    = gsl_matrix_alloc(n,m);
gsl_vector *b    = gsl_vector_alloc(n);
gsl_matrix *R    = gsl_matrix_alloc(m,m);
gsl_matrix *invR = gsl_matrix_alloc(m,m);
gsl_matrix *I    = gsl_matrix_alloc(m,m);

for(int i=0;i<n;i++){
	double xi  = gsl_vector_get(x,i);
	double yi  = gsl_vector_get(y,i);
	double dyi = gsl_vector_get(dy,i);
	assert(dyi>0);
	gsl_vector_set(b,i,yi/dyi);
	for(int k=0; k<m; k++)gsl_matrix_set(A,i,k,f(k,xi)/dyi);
	}
qr_gs_decomp(A,R);
qr_gs_solve(A,R,b,c);

gsl_matrix_set_identity(I);
qr_gs_inverse(I,R,invR);
gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,invR,invR,0,S);

gsl_matrix_free(A);
gsl_vector_free(b);
gsl_matrix_free(R);
gsl_matrix_free(invR);
gsl_matrix_free(I);
}
/*
void singular_val_decomp(gsl_matrix_set* A, gsl_matrix* V, gsl_vector* b, gsl_vector* c){
	int n = (*A).size1;
	int m = (*A).size2;
	//gsl_matrix* V = gsl_matrix_alloc(n,n); // The matrix used to store the eigenvectors
	gsl_vector* e = gsl_vector_alloc(n); // A vector used to store corresponding eigenvalues
	gsl_matrix*	D = gsl_matrix_alloc(n,n); // The diagonal matrix D = VTAV
	gsl_matrix* S = gsl_matrix_alloc(n,n);
	gsl_matrix* U = gsl_matrix_alloc(n,m);
	int num_rot = 0;
	jacobi(gsl_matrix* A, gsl_vector* e, gsl_matrix* V, int* num_rot)

	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,V,A,0,VTAV);

	for (int i = 0; i < (*VTAV).size1; i++) { //Takes the square-root of all diagonal entries of matrix D--> D^1/2 = S
		double D_ii = sqrt(gsl_matrix_get(VTAV,i,i));
		gsl_matrix_set(S,i,i,D_ii);
		gsl_matrix_set(D,i,i,1/D_ii);
	}

// We wish to compute the matrix U=AVD^-1/2. This is done by employing the BLAS-function "gsl_blas_dgemm" twice
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,A,V,0,U); // Computes the matrix product A*V--> U
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,U,D,0,U); // Computes the matrix product
qr_gs_solve(U,S,b,c);

}*/
