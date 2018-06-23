#include "smp.h"
#include<omp.h>
void randomx(int dim, gsl_vector* a, gsl_vector* b, gsl_vector* x, gsl_rng* r){
  for (int i = 0; i < dim; i++) {
    double ai = gsl_vector_get(a,i);
    double bi = gsl_vector_get(b,i);
    gsl_vector_set(x, i, ai + gsl_rng_uniform(r)*(bi - ai)); // Setting "random" predefined nodes
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
  const gsl_rng_type *T;
  gsl_rng *r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  for (int i = 0; i < N; i++) {
    randomx(dim, a, b, x, r); // Creating N random nodes
    fx = f(x); // Calculating f and f^2 to be used in sigma^2=<f^2> - <f>^2
    sum += fx;
    sum2 += fx*fx;
  }
  double avr = sum/N; //Computing the average i.e. eq. (9.2)
  *result = avr*V;

  double var = sum2/N - avr*avr; //sigma^2=<f^2> - <f>^2
  *err = V*sqrt(var)*1/sqrt(N);

  gsl_vector_free(x);
  gsl_rng_free(r);
}



void plain_monte_OMP(int dim, gsl_vector* a, gsl_vector* b, double f(gsl_vector* x),
                  int N, double *result, double *err){
  //Determining the volume over which the integration is carried out
  double V = 1;
  for (int i = 0; i < dim; i++) {
    double ai = gsl_vector_get(a,i);
    double bi = gsl_vector_get(b,i);
    V *= (bi - ai);
  }

  double sum_f1 = 0;
  double sum_f1_2 = 0;
  double sum_f2 = 0;
  double sum_f2_2 = 0;
  double fx1;
  double fx2;
  gsl_vector* x1=gsl_vector_alloc(dim);
  gsl_vector* x2=gsl_vector_alloc(dim);

  	  const gsl_rng_type *T;
  	  gsl_rng *r1,*r2;
  	  gsl_rng_env_setup();
  	  T = gsl_rng_default;
  	  r1 = gsl_rng_alloc(T);
  	  r2 = gsl_rng_alloc(T);

  #pragma omp parallel sections
  {
  	#pragma omp section
    {
        for(int i = 0; i<N/2; i++){
          randomx(dim, a, b, x1, r1); // Creating N random nodes
          fx1 = f(x1); // Calculating f and f^2 to be used in sigma^2=<f^2> - <f>^2
          sum_f1 += fx1;
          sum_f1_2 += fx1*fx1;
        }
      }
  	#pragma omp section
    {
      for(int i = 0; i<N/2; i++){
        randomx(dim, a, b, x2, r2); // Creating N random nodes
        fx2 = f(x2); // Calculating f and f^2 to be used in sigma^2=<f^2> - <f>^2
        sum_f2 += fx2;
        sum_f2_2 += fx2*fx2;
      }
    }
  }

  double sum = sum_f1 + (sum_f2);
  double sum2 = sum_f1_2 + (sum_f2_2);

  double avr = sum/N;
  *result = avr*V;
  double var = sum2/N - avr*avr;
  *err = V*sqrt(var)*1/sqrt(N);

  gsl_rng_free(r1);
  gsl_rng_free(r2);
  gsl_vector_free(x1);
  gsl_vector_free(x2);

  return;
}
