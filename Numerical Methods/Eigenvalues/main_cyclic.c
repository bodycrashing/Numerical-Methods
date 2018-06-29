#include "jacobi.h"
#define RND ((double)rand()/RAND_MAX)


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


int sweeps = jacobi(A, e, V, &num_rot_cyclic);


// Freeing memory
gsl_matrix_free(A);
gsl_matrix_free(B);
gsl_matrix_free(V);
gsl_matrix_free(VTAV);

gsl_vector_free(e);

return 0;
}
