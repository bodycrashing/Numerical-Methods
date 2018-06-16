
#include "qr.h"

/* The algortihm conducts QR-decomposition of matrix A, in which A is replaced with Q using the Gram-Schmidt approach.
The corresponding R matrix is built and filled seperately */
void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R) {

  int m = A->size2;
  for(int i=0; i<m; i++){
    gsl_vector_view col_i = gsl_matrix_column(A,i); // A so called "vector view" of the i'th column of matrix A is created
                                                // This allows one to directly work with the entries of the specific column without creating
                                                // an actual copy of the vector thus reducing the required memory and computational effort
    double r = gsl_blas_dnrm2(&col_i.vector);
    gsl_matrix_set(R,i,i,r);
    gsl_vector_scale(&col_i.vector, 1/r); //normalization

    for(int j=i+1; j<m; j++){
      gsl_vector_view q = gsl_matrix_column(A,j);
      double s = 0;
      gsl_blas_ddot(&col_i.vector, &q.vector, &s); // saves the dot-prodcut col_i*q -> s
      gsl_blas_daxpy(-s, &col_i.vector, &q.vector); // "gsl_blas_daxpy" is a GSL-Blas function returning -s*e+q -> q
      gsl_matrix_set(R,i,j,s);
      gsl_matrix_set(R,j,i,0);
    }
  }
}

void qr_back_sub(gsl_matrix *R, gsl_vector *x) {
int m = R->size1;
for(int i=m-1; i>=0; i--){
	double s = 0;
	for(int k=i+1; k<m; k++)
		s += gsl_matrix_get(R,i,k)*gsl_vector_get(x,k);
	  gsl_vector_set(x,i,(gsl_vector_get(x,i)-s)/gsl_matrix_get(R,i,i));
	}
}

void qr_gs_solve(gsl_matrix *Q, gsl_matrix *R, gsl_vector *b, gsl_vector *x) {
  gsl_blas_dgemv(CblasTrans,1.0,Q,b,0.0,x); // x = 1*Q^T*b + 0*x
  qr_back_sub(R,x); // backsustitution
}

void qr_gs_inverse(gsl_matrix *Q, gsl_matrix* R, gsl_matrix *B) {
  int n=R->size1;
  gsl_vector *b = gsl_vector_calloc(n);
  gsl_vector *x = gsl_vector_calloc(n);
  for(int i=0; i<n; i++){
    gsl_vector_set(b,i,1.0);
    qr_gs_solve(Q,R,b,x);
    gsl_vector_set(b,i,0.0);
    gsl_matrix_set_col(B,i,x);
  }
  gsl_vector_free(b);
  gsl_vector_free(x);
}




/*---------------- Givens Decomposition ----------------*/
void qr_givens(gsl_matrix* A){
  for(int q=0; q<A->size2; q++){
    for(int p=q+1; p<A->size1; p++){
      double theta = atan2(gsl_matrix_get(A,p,q),gsl_matrix_get(A,q,q));
      double s=sin(theta), c=cos(theta);
      for(int k=q; k<A->size2; k++){
        double xq=gsl_matrix_get(A,q,k), xp=gsl_matrix_get(A,p,k);
        gsl_matrix_set(A,q,k, xq*c+xp*s);
        gsl_matrix_set(A,p,k,-xq*s+xp*c);
      }
      gsl_matrix_set(A,p,q,theta);
    }
  }
}

void qivens_apply(gsl_matrix* QR, gsl_vector* v){
	for(int q=0; q<QR->size2; q++)for(int p=q+1; p<QR->size1; p++){
		double theta = gsl_matrix_get(QR,p,q);
		double vq=gsl_vector_get(v,q), vp=gsl_vector_get(v,p);
		gsl_vector_set(v,q, vq*cos(theta)+vp*sin(theta));
		gsl_vector_set(v,p,-vq*sin(theta)+vp*cos(theta));
		}
}



void givens_solve(gsl_matrix* QR, gsl_vector* b){
  qivens_apply(QR,b);
  //for(int i=0;i<b->size;i++)gsl_vector_set(x,i,gsl_vector_get(v,i));
  int m=QR->size2;
  for(int i=m-1;i>=0;i--){
    double s=0;
    for(int k=i+1;k<m;k++)
    s+=gsl_matrix_get(QR,i,k)*gsl_vector_get(b,k);
    gsl_vector_set(b,i,(gsl_vector_get(b,i)-s)/gsl_matrix_get(QR,i,i));
  }
}



void givens_inverse(gsl_matrix *QR, gsl_matrix *B) {
int n=QR->size1;
gsl_vector *b = gsl_vector_calloc(n);
gsl_matrix_set_identity(B);
for(int i=0;i<n;i++){
  gsl_matrix_get_col(b,B,i);
	givens_solve(QR,b);
	gsl_matrix_set_col(B,i,b);
	}
gsl_vector_free(b);
}
