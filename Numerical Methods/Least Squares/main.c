#include "qr.h"
#include "least_squares.h"
#include "math.h"

double funs(int i, double x){
  switch(i){
    case 0: return log(x); break;
    case 1: return 1.0;   break;
    case 2: return x;     break;
    default: {fprintf(stderr,"funs: wrong i:%d",i); return NAN;}
  }
}

int main(int argc, char const *argv[]) {

double x[]  = {0.1, 1.33, 2.55, 3.78, 5, 6.22, 7.45, 8.68, 9.9};
double y[]  = {-15.3, 0.3, 2.45, 2.75, 2.27, 1.35, 0.157, -1.23, -2.75};
double dy[] = {1.04, 0.594, 0.983, 0.998, 1.11, 0.398, 0.535, 0.968, 0.478};
int n = sizeof(x)/sizeof(x[0]); // The command sizeof(x) returns the size (in bits) of variable x of data-type char, int etc. Hence The
                                // the following line of code effectively provides the length of the data-vector x.
  int m = 3;
  gsl_vector* c = gsl_vector_alloc(m);
  gsl_matrix* S = gsl_matrix_alloc(m,m);
  gsl_vector* xi = gsl_vector_alloc(n);
  gsl_vector* yi = gsl_vector_alloc(n);
  gsl_vector* dyi = gsl_vector_alloc(n);
  for(int i=0; i<n; i++){
    gsl_vector_set(xi,i,x[i]);
    gsl_vector_set(yi,i,y[i]);
    gsl_vector_set(dyi,i,dy[i]);
    printf("%g %g %g\n",x[i],y[i],dy[i]);
  } printf("\n\n");

  lsfit(m, funs, xi, yi, dyi, c, S); // This function creates n,m-matrix A where n=x->size and m="number of functions in the linear combination F_c"
  gsl_vector* dc = gsl_vector_alloc(m);
	for(int k=0;k<m;k++){
		double skk=gsl_matrix_get(S,k,k);
		gsl_vector_set(dc,k,sqrt(skk));
		}

  double fit(double x){
    double s=0;
    for(int k=0; k<m; k++){
      s += gsl_vector_get(c,k)*funs(k,x);
    }
    return s;
  }

	double fit_plus(int i, double x){
		return fit(x)+gsl_vector_get(dc,i)*funs(i,x);
		}

	double fit_minus(int i, double x){
		return fit(x)-gsl_vector_get(dc,i)*funs(i,x);
		}

	double z,dz=(x[n-1]-x[0])/90;
	for(int i=0;i<m;i++){
	z=x[0]-dz/2;
		do{
		printf("%g %g %g %g\n",z,fit(z),fit_plus(i,z),fit_minus(i,z));
		z+=dz;
		}while(z<x[n-1]+dz);
	printf("\n\n");
	}


gsl_matrix* A = gsl_matrix_alloc(n,m);
gsl_matrix* U = gsl_matrix_alloc(n,m);
gsl_matrix* V = gsl_matrix_alloc(m,m);
//singular_val_decomp(A, V, S);

  gsl_vector_free(c);
  gsl_vector_free(xi);
  gsl_vector_free(yi);
  gsl_vector_free(dyi);
  gsl_matrix_free(S);
  gsl_matrix_free(U);


  return 0;
}
