#include "qr.h"
#include "least_squares.h"

double funs(int i, double x){
  switch(i){
    case 0: return log(x); break;
    case 1: return 1.0;   break;
    case 2: return x;     break;
    default: {fprintf(stderr,"funs: wrong i:%d",i); return NAN;}
  }
}


int main(int argc, char const *argv[]) {


  double x[] = {0.100,0.145,0.211,0.307,0.447,0.649,0.944,1.372,1.995,2.900};
  double y[] = {12.644,9.235,7.377,6.460,5.555,5.896,5.673,6.964,8.896,11.355};
  double dy[] = {0.858,0.359,0.505,0.403,0.683,0.605,0.856,0.351,1.083,1.002};
  int n = sizeof(x)/sizeof(x[0]); // The command sizeof(x) returns the size (in bits) of variable x of data-type char, int etc. Hence The
                                // the following line of code effectively provides the length of the data-vector x.
  for(int i=0; i<n; i++)printf("%g %g %g\n",x[i],y[i],dy[i]);
  printf("\n\n");

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
  }

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
    return s;
  }
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














  return 0;
}
