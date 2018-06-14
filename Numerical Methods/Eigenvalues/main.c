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

//printf("---------------\tExercise Part A:\t---------------\n");
/*int sweeps = jacobi(A,e,V,&num_rot); printf("n=%i, sweeps=%i\n",n,sweeps);
if(n<max_print){
fprintf(stderr,"\na random symmetric matrix A: \n"); matrix_print(B);
fprintf(stderr,"\nthe result of Jacobi diagonalization: \n");
fprintf(stderr,"sweeps\t = %d\n",sweeps);
fprintf(stderr, "eigenvalues:\n");
vector_print(e);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,B,V,0,A);
gsl_blas_dgemm(CblasTrans  ,CblasNoTrans,1,V,A,0,VTAV);
fprintf(stderr, "\ncheck: V^T*A*V should be diagonal with above eigenvalues:\n");
matrix_print(VTAV);
printf("number of rotations = %i\n", num_rot);
}
*/

int sweeps = jacobi(A, e, V, &num_rot_cyclic);
FILE* cyclic_data = fopen("cyclic_data.txt","w+");
fprintf(cyclic_data,"%i\t%i\t%i\n", n, sweeps, num_rot_cyclic);
fclose(cyclic_data);

if(n<max_print){
printf("\na random symmetric matrix A: \n"); matrix_print(B);
printf("\nthe result of Jacobi diagonalization: \n");
printf("sweeps\t = %d\n",sweeps);
printf("eigenvalues:\n");
vector_print(e);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,B,V,0,A);
gsl_blas_dgemm(CblasTrans  ,CblasNoTrans,1,V,A,0,VTAV);
printf("\ncheck: V^T*A*V should be diagonal with above eigenvalues:\n");
matrix_print(VTAV);
printf("number of rotations = %i\n", num_rot_cyclic);
}


//printf("---------------\tExercise Part B:\t---------------\n");
int row_number = 3;
int num_rot_eig_by_eig = 0;
gsl_matrix_set_identity(V);
gsl_vector_set_zero(e);
printf("1. Using the same matrix A as before for comparison\n");

sweeps = jacobi_eig_by_eig(A, e, V, row_number, &num_rot_eig_by_eig);
FILE* eig_by_eig_data = fopen("eig_by_eig_data.txt","w+");
fprintf(eig_by_eig_data,"%i\t%i\t%i\n", n, sweeps, num_rot_cyclic);
fclose(eig_by_eig_data);
/*
if (n<max_print) {
	printf("The first %i eigenvalue (in ascending order) are: \n", row_number);
	for (int i = 0; i < row_number; i++) {
		gsl_matrix_memcpy(A,B);
		int sweeps = jacobi_eig_by_eig(A,e,V,i);
		printf("lambda_%i = %lg (determined after n = %i sweeps)\n", i, gsl_vector_get(e, i), sweeps);
	}
}
*/



// Freeing memory
gsl_matrix_free(A);
gsl_matrix_free(B);
gsl_matrix_free(V);
gsl_matrix_free(VTAV);

gsl_vector_free(e);

return 0;
}
