#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
void max_finder(gsl_matrix* A){
  gsl_vector* v = gsl_vector_alloc(A.size2);
  double prev;
  double new;
  for (int i = 0; i < A.size1; i++) {
    for (int j = 0; j < A.size2-1; j++) {
      prev = fabs(gsl_matrix_get(A,i,j));
      new = fabs(gsl_matrix_get(A,i,j+1);
      if (prev<new) {
        gsl_vector_set(v,i);
      }
    }
  }
gsl_matrix_free(v);
