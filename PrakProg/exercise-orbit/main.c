#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
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


int orbit(double x, const double y[], double dydx[], void* params){
	double eps = *(double*)params;
	dydx[0] = y[1];
	dydx[1] = 1 + eps*y[0]*y[0] - y[0];
	return GSL_SUCCESS;
}

int main(int argc,char** argv){

	FILE* logistic_stream = fopen("logistic.out.txt","w+");
  double epsabs = 1e-8;
	double epsrel = 1e-8;
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


	double eps = 0,uprime=0;

	if (argc == 3){
		eps = atof(argv[1]), uprime = atof(argv[2]);
	}

	gsl_odeiv2_system orbit_sys;
	orbit_sys.function = orbit;
	orbit_sys.jacobian = NULL;
	orbit_sys.dimension = 2;
	orbit_sys.params = (void*)&eps;

	const gsl_odeiv2_step_type* steptype = gsl_odeiv2_step_rk8pd;
	gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&orbit_sys,steptype,hstart,epsabs,epsrel);

	double x, x0 = 0, x_max = 10*M_PI, x_step = 0.1;
	double y[2];
	y[0] = 1;
	y[1] = uprime;

	printf("x \t y\n");
	for(x = x0; x<x_max; x+=x_step){
		int status
		= gsl_odeiv2_driver_apply(driver, &x0, x, y);
		printf("%g \t %g\n", x, y[0]);

		if(status != GSL_SUCCESS)
			fprintf(stderr,"fun: status = %i",status);
	}
	gsl_odeiv2_driver_free(driver);

	return 0;
}
