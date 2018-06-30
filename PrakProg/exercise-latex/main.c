#include<stdio.h>
#include<math.h>
#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

int error_func(double x, const double y[], double dydx[], void* params){
	dydx[0] = 2/sqrt(M_PI)*exp(-x*x);
	return GSL_SUCCESS;
}


double solve_diff(double x){
	gsl_odeiv2_system sys;
	sys.function = error_func;
	sys.jacobian = NULL;
	sys.dimension = 1;
	sys.params = NULL;

	double hstart=copysign(0.1,x);
  fprintf(stderr, "%lg\n",hstart);
	double acc = 1e-8;
	double eps = 1e-8;
	gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, hstart,acc,eps);

	double x0 = 0;
  double y[1] = {0};
	gsl_odeiv2_driver_apply(driver,&x0,x,y);


	gsl_odeiv2_driver_free(driver);
	return y[0];
}






int main(int argc, char const *argv[]) {


  double a=atof(argv[1]), b=atof(argv[2]), dx=atof(argv[3]);
  	for(double x=a;x<b;x+=dx){
  		printf("%g %g\n ",x,solve_diff(x));

  }


  return 0;
}
