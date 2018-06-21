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


double c(double x){
  return cos(x);
}

double d(double x){
  return sin(x);
}

double f_sin_cos(double x,double y){
  return sin(x)*cos(y);
}







int main(int argc, char const *argv[]) {
  int dim = 3; //Dimension of integration space
  //int N = 10e+6; //Number of sampling points
  double err, res;

  gsl_vector* a = gsl_vector_alloc(dim);
  gsl_vector* b = gsl_vector_alloc(dim);
  gsl_vector_set_zero(a);
  gsl_vector_set_all(b,1);

/*
  plain_monte(dim, a, b, f_3D_Sphere, N, &res, &err);
  printf("Volume of 3D unit-sphere\n");
  printf("Result of the integration: res=%lg\nThe error of the integration: error=%lg\n", res, err);
*/

  /* Testing the function supplied in the exercise */
  /*
  gsl_vector_set_zero(a);
  gsl_vector_set_all(b,M_PI);
  plain_monte(dim, a, b, f_gamma, N, &res, &err);
  printf("Integration of gamma function\n");
  printf("Result of the integration: res=%lg\nThe error of the integration: error=%lg\n", res, err);
  */

  /*------------ Part B: Checking the O(1/sqrt(N)) behavior--------------*/
  int Nmax = 1000;

	FILE* data = fopen("data.txt","w+");
  for (int i = 100; i < Nmax; i+=10) {
    plain_monte(dim, a, b, f_3D_Sphere, i, &res, &err);
    fprintf(data, "%i\t%lg\n",i ,err);
  }
  fclose(data);


/*------------ Part C: 2D integrstion routine --------------*/
double eps = 1e-5, acc = 1e-5;
err = 0;
res = Adaptive2D(f_sin_cos,0,M_PI,c,d,acc,eps, &err);
printf("\n\nAdaptive2D result=%lg\n", res);


  return 0;
}
