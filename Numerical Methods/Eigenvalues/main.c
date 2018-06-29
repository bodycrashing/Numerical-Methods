#include "jacobi.h"
#define RND ((double)rand()/RAND_MAX)
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


int main(int argc, char** argv)
{
int n = atoi(argv[1]);

gsl_matrix *A = gsl_matrix_alloc(n,n);
gsl_matrix *B = gsl_matrix_alloc(n,n);
for(int i=0;i<n;i++) for(int j=i;j<n;j++) {
	double x = RND;
	gsl_matrix_set(A,i,j,x);
	gsl_matrix_set(A,j,i,x);
	}
gsl_matrix_memcpy(B,A);

gsl_matrix *V = gsl_matrix_alloc(n,n);
gsl_vector *e = gsl_vector_alloc(n);
gsl_matrix*	VTAV = gsl_matrix_alloc(n,n);
int num_rot_cyclic = 0;


int sweeps_cyclic = jacobi(A, e, V, &num_rot_cyclic);
if(n<max_print){
printf("\n------ Part A: Cyclic Jacobi diagonalization of %ix%i matrix: -------\n",n,n);
printf("A random symmetric matrix A: \n"); matrix_print(B);
printf("\nSweeps\t = %d\n",sweeps_cyclic);
printf("Number of rotations = %i\n", num_rot_cyclic);
printf("Eigenvalues:\n");
vector_print(e);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,B,V,0,A);
gsl_blas_dgemm(CblasTrans  ,CblasNoTrans,1,V,A,0,VTAV);
printf("\nCheck: V^T*A*V should be diagonal with above eigenvalues:\n");
matrix_print(VTAV);
}

printf("\n------ Part B: Eig. by Eig Jacobi diagonalization of %ix%i matrix: -------\n",n,n);
int row_number = 5;
int sorting = 0;
int num_rot_eig_by_eig = 0;
gsl_matrix_set_identity(V);
gsl_vector_set_zero(e);
gsl_vector* sweeps_eig_by_eig = gsl_vector_alloc(n);
gsl_vector* rot_eig_by_eig = gsl_vector_alloc(n);
printf("\n1. Using the same matrix A as before for comparison\n");


if (n<max_print){
	printf("The first %i eigenvalue in Ascending order are: \n", row_number);
	for (int i = 0; i < row_number; i++) {
		gsl_matrix_memcpy(A,B);
		gsl_vector_set_zero(e);
		int sweeps = jacobi_eig_by_eig(A, e, V, i+1, &num_rot_eig_by_eig, sorting);
		gsl_vector_set(sweeps_eig_by_eig, i, sweeps); gsl_vector_set(rot_eig_by_eig, i, num_rot_eig_by_eig);
		printf("e_%i = %lg\n", i, gsl_vector_get(e,i));
		printf("Number of sweeps required: %i\nNumber of rotations required: %i\n\n", sweeps, num_rot_eig_by_eig);
	}
}


sorting = 1;
if (n<max_print){
	printf("The first %i eigenvalue in Descending order are: \n", row_number);
	for (int i = 0; i < row_number; i++) {
		gsl_matrix_memcpy(A,B);
		gsl_vector_set_zero(e);
		int sweeps = jacobi_eig_by_eig(A, e, V, i+1, &num_rot_eig_by_eig, sorting);
		printf("e_%i = %lg\n", i, gsl_vector_get(e,i));
	}
}

printf("\nComparision of the Cyclic and Value by Value method\n");
printf("\tSweeps\tRotations\n");
printf("Cyclic:\t%i\t%i\n", sweeps_cyclic, num_rot_cyclic);
for (double i = 0; i < row_number; i++) {
	double sweeps_i = gsl_vector_get(sweeps_eig_by_eig,i);
	double num_rot_i = gsl_vector_get(rot_eig_by_eig,i);
	printf("%lg Eig:\t%lg\t%lg\n", i+1, sweeps_i, num_rot_i);
}

// Freeing memory
gsl_matrix_free(A);
gsl_matrix_free(B);
gsl_matrix_free(V);
gsl_matrix_free(VTAV);
gsl_vector_free(e);

return 0;
}
