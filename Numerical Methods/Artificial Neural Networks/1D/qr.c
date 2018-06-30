#include"ann.h"

void vector_print(const gsl_vector* v){
	for(int j=0;j<v->size;j++){
		printf("%8.3g\n",gsl_vector_get(v,j));
	}
}

void matrix_print(const gsl_matrix* m){
	for(int j=0;j<m->size1;j++){
		for(int i=0;i<m->size2;i++){
			printf("%8.3g ",gsl_matrix_get(m,j,i));
		}
		printf("\n");
	}
}


void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R){
	int m=A->size2;
	assert(m==R->size1 && R->size1==R->size2);
	double R_ii, R_ij;

	for(int i=0;i<m;i++){
		gsl_vector_view a_i=gsl_matrix_column(A,i);
		R_ii = gsl_blas_dnrm2(&a_i.vector);
		gsl_matrix_set(R,i,i,R_ii);
		gsl_vector_scale(&a_i.vector,1.0/R_ii);

		for(int j = i+1;j<m;j++){
			gsl_vector_view a_j = gsl_matrix_column(A,j);
			gsl_blas_ddot (&a_i.vector,&a_j.vector,&R_ij);
			gsl_matrix_set(R,i,j,R_ij);
			gsl_blas_daxpy(-R_ij,&a_i.vector,&a_j.vector);
		}
	}
}

void qr_gs_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){
	gsl_blas_dgemv(CblasTrans,1,Q,b,0,x);
	for(int i=x->size-1;i>=0;i--){
		double x_i=gsl_vector_get(x,i);
		for(int j=i+1;j<x->size;j++){
			x_i -= gsl_matrix_get(R,i,j)*gsl_vector_get(x,j);
		}
		x_i/=gsl_matrix_get(R,i,i);
		gsl_vector_set(x,i,x_i);
	}
}

void qr_gs_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B){
	gsl_vector* temp = gsl_vector_alloc((*B).size2);
	gsl_matrix_set_identity(B);

	for(int i=0;i<(*B).size2;i++){
		gsl_vector_view b_i=gsl_matrix_column(B,i);
		gsl_blas_dcopy(&b_i.vector,temp);
		qr_gs_solve(Q,R,temp,&b_i.vector);
	}
	gsl_vector_free(temp);
}

void matrix_inverse(gsl_matrix* M){
	int n = M->size1;
	gsl_matrix* temp = gsl_matrix_alloc(n,n);
	for(int i = n-1;i>=0;i--){
		for(int j=0;j<n;j++){
			if(j==i){
				gsl_matrix_set(temp,i,j,1.0/gsl_matrix_get(M,i,j));
			}
			else{
				for(int k=i+1;k<n;k++){
					gsl_matrix_set(temp,i,j,
					gsl_matrix_get(temp,i,j)-
					gsl_matrix_get(M,i,k)*gsl_matrix_get(temp,k,j)/
					gsl_matrix_get(M,i,i));
				}
			}
		}
	}
	gsl_matrix_free(temp);
}
