#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>


int logistic_diff (double t, const double y[], double f[], void *params)
{
  (void)(t); /* avoid unused parameter warning */
  f[0] = y[0]*(1 - y[0]);
  return GSL_SUCCESS;
}

double logistic(double x)
{
  double result = 1/(1 + exp(-x));
  return result;
};


int orbit_eq(double phi, const double y[], double yPrime[], void *params) {
	double epsilon = *(double *) params;
	yPrime[0] = y[1];
	yPrime[1] = 1 - y[0] + epsilon*y[0]*y[0];
	return GSL_SUCCESS;
}



int main (void){

  FILE* logistic_stream = fopen("logistic.out.txt","w+");
  double  epsabs = 1e-8, epsrel = 1e-8;
  double hstart = 1e-3;

  gsl_odeiv2_system logistic_sys = {logistic_diff, NULL, 1, NULL};
  gsl_odeiv2_driver *logistic_driver = gsl_odeiv2_driver_alloc_y_new(&logistic_sys, gsl_odeiv2_step_rkf45, hstart, epsabs, epsrel);

  double t = 0.0;
  double t1 = 3.0;
  double y1[1] = {0.5};
  int logistic_result;
  for (double i=t; i < t1; i+=0.05) {
      logistic_result = gsl_odeiv2_driver_apply(logistic_driver, &t, i, y1);
      if (logistic_result != GSL_SUCCESS)
        {
          printf ("error: logistic driver returned %d\n", logistic_result);
          break;
        }
    fprintf(logistic_stream,"%lg\t%lg\t%lg\n", i, y1[0],logistic(i));
    }


  return 0;
}
