#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "qr.h"
#define RND ((double)rand()/RAND_MAX)
#define FMT "%7.3f" // Formats prints such that it has a width of 7 "units" and 3 decimal places

void printm(gsl_matrix *A){
	for(int r=0;r<A->size1;r++){
		for(int c=0;c<A->size2;c++) printf(FMT,gsl_matrix_get(A,r,c));
		printf("\n");}
}

void printv(gsl_vector *A){
	for(int i=0;i<A->size;i++){
		printf(FMT,gsl_vector_get(A,i));
		printf("\n");
	}
}


int main(int argc, char** argv){
int n=4, m=3;
gsl_matrix *A = gsl_matrix_calloc(n,m);
gsl_matrix *R = gsl_matrix_calloc(m,m);

for(int i=0; i<n; i++){
	for(int j=0;j<m;j++){
		gsl_matrix_set(A,i,j,RND);
	}
}
printf("\n|-------------- Part A: QR decomposition --------------|\n");
printf("Random tall matrix A:\n"); printm(A);

qr_gs_decomp(A,R);
printf("\n Matrix Q:\n"); printm(A);
printf("\n Cheching that matrix R is upper-triangular:\n"); printm(R);

gsl_matrix *Q_TQ = gsl_matrix_calloc(m,m);
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,A,A,0.0,Q_TQ);
printf("\nChecking that Q is an orthogonal matrix by computing\nthe product Q^T*Q (should be 1):\n");
printm(Q_TQ);

gsl_matrix *qr = gsl_matrix_calloc(n,m);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A,R,0.0,qr);
printf("\nCheckhing the matrix product QR = A:\n"); printm(qr);



/*------------ Part B: Matrix inverse by Gram-Schmidt QR factorization ------------*/
printf("\n\n|-------------- Part B: Linear System --------------|\n");
m = 3;
gsl_matrix *AA = gsl_matrix_alloc(m,m);
gsl_matrix *QQ = gsl_matrix_alloc(m,m);
gsl_matrix *RR = gsl_matrix_alloc(m,m);

for(int i=0; i<m; i++){
	for(int j=0;j<m;j++){
		gsl_matrix_set(AA,i,j,RND);
	}
}
gsl_matrix_memcpy(QQ,AA); // Since the matrix AA is destroyed in the QR-decomp. process and we intend to eventually check that Ax=b
													// a copy of AA is made.
printf("Square random matrix A:\n"); printm(QQ);
gsl_vector *b = gsl_vector_alloc(m);
for(int i=0;i<m;i++)gsl_vector_set(b,i,RND); // Allocating, setting and printing the corresponding solution to the eqaution Ax=b.
printf("\nA random right-hand side b (Ax=b):\n"); printv(b);

qr_gs_decomp(QQ,RR);
gsl_vector *x = gsl_vector_alloc(m);
qr_gs_solve(QQ,RR,b,x);
printf("\nSolution x using QR-decomp:\n"); gsl_vector_fprintf(stdout,x,FMT);
gsl_blas_dgemv(CblasNoTrans,1.0,AA,x,0.0,b);
printf("\nAx, should be equal b:\n"); gsl_vector_fprintf(stdout,b,FMT);
gsl_matrix *AI = gsl_matrix_calloc(m,m);
qr_gs_inverse(QQ,RR,AI);
printf("\n Inverse matrix B=A^{-1}:\n"); printm(AI);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,AA,AI,0.0,RR);
printf("\n Matrix AA^-1, should be 1:\n"); printm(RR);



/*-------------- Part C: Givens Rotation ----------------*/
printf("\n\n|-------------- Part C:Givens Rotation --------------|\n");
gsl_matrix *Givens_Q = gsl_matrix_alloc(m,m);
gsl_matrix *B = gsl_matrix_alloc(m,m);
gsl_vector *Givens_b = gsl_vector_alloc(m);

gsl_matrix_memcpy(Givens_Q,AA); /* Again the decomposition functions as well as the solve function is built such that it destroyes
																	the original matrix and vectors used, hence we need to copy if we wish to used them for later comparison*/
gsl_vector_memcpy(Givens_b,b);
qr_givens(Givens_Q);
printf("Givens Q-matrix corresponding to rotation angles:\n"); printm(Givens_Q);
givens_solve(Givens_Q, Givens_b);
printf("\n Solution x found using Givens rotation:\n"); printv(Givens_b);

givens_inverse(Givens_Q, B);
printf("\nGivens inverse matrix B=A^{-1}:\n");printm(B);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,AA,B,0.0,Givens_Q);
printf("\n Matrix AA^-1, should be 1:\n"); printm(Givens_Q);


/*Freeing all allocated vectors and matrices*/
gsl_matrix_free(A);
gsl_matrix_free(R);
gsl_matrix_free(Q_TQ);
gsl_matrix_free(qr);
gsl_matrix_free(AA);
gsl_matrix_free(QQ);
gsl_matrix_free(RR);
gsl_matrix_free(Givens_Q);
gsl_matrix_free(B);

gsl_vector_free(x);
gsl_vector_free(b);
gsl_vector_free(Givens_b);

return 0;
}
