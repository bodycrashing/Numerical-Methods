#include<assert.h>
#include <math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include"qr.h"
#include "least_squares.h"
#include "jacobi.h"


#define FMT "%7.3f"
void printv(gsl_vector *A){
	for(int i=0;i<A->size;i++){
		printf(FMT,gsl_vector_get(A,i));
		printf("\n");
	}
}





void lsfit(int m, double f(int i,double x),
	gsl_vector* x, gsl_vector* y, gsl_vector* dy,
	gsl_vector* c, gsl_matrix* S)
{
int n = x->size;

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

gsl_matrix* singular_val_decomp(gsl_matrix* A, gsl_matrix* V, gsl_matrix* S){
	int n = (*A).size1; //A is tall matrix i.e. n>m.
	int m = (*A).size2;
	gsl_matrix* ATA = gsl_matrix_alloc(m,m);
	gsl_matrix*	D = gsl_matrix_alloc(m,m); // The diagonal matrix D = VTAV
	//gsl_matrix* S = gsl_matrix_alloc(m,m);
	gsl_matrix* U = gsl_matrix_alloc(n,m);
	gsl_vector* e = gsl_vector_alloc(n); // A vector used to store corresponding eigenvalues
	int num_rot = 0;

	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,A,A,0,ATA);
	jacobi(ATA, e, V, &num_rot);

	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,V,ATA,0,D);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,D,V,0,D);
	for (int i = 0; i < (*D).size1; i++) { //Takes the square-root of all diagonal entries of matrix D--> D^1/2 = S
		double D_ii = sqrt(gsl_matrix_get(D,i,i));
		gsl_matrix_set(S,i,i,D_ii);
		gsl_matrix_set(D,i,i,1/D_ii);
	}

// We wish to compute the matrix U=AVD^-1/2. This is done by employing the BLAS-function "gsl_blas_dgemm" twice
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,A,V,0,U); // Computes the matrix product A*V--> U
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,U,D,0,U); // Computes the matrix product

gsl_matrix_free(ATA);
gsl_matrix_free(D);
//gsl_matrix_free(S);
//gsl_matrix_free(U);
gsl_vector_free(e);

return U;
}
