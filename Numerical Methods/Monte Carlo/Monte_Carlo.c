#include "Monte_Carlo.h"
#define RND (double)rand()/RAND_MAX

void randomx(int dim, gsl_vector* a, gsl_vector* b, gsl_vector* x){
  for (int i = 0; i < dim; i++) {
    double ai = gsl_vector_get(a,i);
    double bi = gsl_vector_get(b,i);
    gsl_vector_set(x, i, ai + RND*(bi - ai)); // Setting "random" predefined nodes
  }
}

void plain_monte(int dim, gsl_vector* a, gsl_vector* b, double f(gsl_vector* x),
                  int N, double *result, double *err){
  //Determining the volume over which the integration is carried out
  double V = 1;
  for (int i = 0; i < dim; i++) {
    double ai = gsl_vector_get(a,i);
    double bi = gsl_vector_get(b,i);
    V *= (bi - ai);
  }

  double sum, sum2;
  double fx;
  gsl_vector* x = gsl_vector_alloc(dim);
  for (int i = 0; i < N; i++) {
    randomx(dim, a, b, x); // Creating N random nodes
    fx = f(x); // Calculating f and f^2 to be used in sigma^2=<f^2> - <f>^2
    sum += fx;
    sum2 += fx*fx;
  }
  double avr = sum/N; //Computing the average i.e. eq. (9.2)
  *result = avr*V;

  double var = sum2/N - avr*avr; //sigma^2=<f^2> - <f>^2
  *err = V*sqrt(var)*1/sqrt(N);

  gsl_vector_free(x);
}


double Adaptive2D(double f(double x, double y),double a, double b, double c(double x), double d(double x),
                            double acc, double eps, double *err){
double a_new, b_new;
double outer_int(double x){
	double inner_int(double y){
		return f(x,y);
	}
	a_new = c(x);
	b_new = d(x);
	return Adaptive_Infinity(inner_int,a_new,b_new,acc,eps,err);
}
	return Adaptive_Infinity(outer_int,a,b,acc,eps,err);
}
