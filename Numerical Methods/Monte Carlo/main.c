#include "Monte_Carlo.h"

double f_gamma(gsl_vector* x){
  double x0 = gsl_vector_get(x,0);
  double x1 = gsl_vector_get(x,1);
  double x2 = gsl_vector_get(x,2);
	return 1.0/(1-cos(x0)*cos(x1)*cos(x2))/(M_PI*M_PI*M_PI);
}

double f_3D_Sphere(gsl_vector* x){
  double x0 = gsl_vector_get(x,0);
  double x1 = gsl_vector_get(x,1);
  double x2 = gsl_vector_get(x,2);
	return x0*x0 + x1*x1 + x2*x2;
}

// Functions to be used as limits in Part C.
double c(double x){
  return x*x;
}

double d(double x){
  return sin(x);
}

double f_sin_cos(double x,double y){
  return sin(x)*cos(y);
}


int main(int argc, char const *argv[]) {
  /*-------------- Part A: Plain Monte Carlo integration: -------------*/
  printf("\n|--------- Part A: Plain Monte Carlo integration: --------|\n");
  int dim = 3; //Dimension of integration space
  int N = 10e+6; //Number of sampling points
  double err, res;

  gsl_vector* a = gsl_vector_alloc(dim);
  gsl_vector* b = gsl_vector_alloc(dim);
  gsl_vector_set_zero(a);
  gsl_vector_set_all(b,1);

  plain_monte(dim, a, b, f_3D_Sphere, N, &res, &err);
  printf("Integration of 3D unit-sphere from 0 to 1:\n");
  printf("Error = %lg\n", err);
  printf("Result: Q=%lg\nCorrect result: Q_exact=1\n",res);

  gsl_vector_set_zero(a);
  gsl_vector_set_all(b,M_PI);
  double exact = 1.3932039296856768591842462603255;
  plain_monte(dim, a, b, f_gamma, N, &res, &err);
  printf("\nIntegration of gamma function from 0 to Pi:\n");
  printf("Error = %lg\n", err);
  printf("Result: Q=%lg\nCorrect result: Q_exact=%.15g\n",res,exact);


  /*------------ Part B: Checking the O(1/sqrt(N)) behavior--------------*/
  int Nmax = 1000;
	FILE* data = fopen("data.txt","w+");
  for (int i = 100; i < Nmax; i+=10) {
    plain_monte(dim, a, b, f_3D_Sphere, i, &res, &err);
    fprintf(data, "%i\t%lg\n",i ,err);
  }
  fclose(data);


  /*------------ Part C: Adaptive 2D integrator --------------*/
  printf("\n|------------ Part C: Adaptive 2D integrator --------------|\n");
  double eps = 1e-5, acc = 1e-5;
  err = 0;
  res = Adaptive2D(f_sin_cos,0,M_PI,c,d,acc,eps, &err);
  printf("2D integration of f(x,y)=sin(x)*cos(y) using functions:\n\n\tc(x)=x^2  d(x)=sin(x)\n\nas curves limiting the integration area\n");
  printf("Result: Q=%lg\n\n", res);


  gsl_vector_free(a);
  gsl_vector_free(b);

  return 0;
}
