#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>
#include<math.h>
#include<assert.h>
#define STEPPER gsl_odeiv2_step_rkf45

int s_wave(double r, const double y[], double yp[], void* params){
	double e = *(double*)params;
	yp[0] = y[1];
	yp[1] = -2*(1/r+e)*y[0];
	return GSL_SUCCESS;
}

double F(double e, double r){
	assert(r>=0);
	double rmin = 1e-3;
	if(r<rmin) return r-r*r;

	gsl_odeiv2_system sys;
	sys.function = s_wave;
	sys.jacobian = NULL;
	sys.dimension = 2;
	sys.params = (void*)&e;

	gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new
		(&sys,STEPPER,1e-3,1e-6,1e-6);

	double t=rmin, y[] = {t-t*t, 1-2*t};
	int status = gsl_odeiv2_driver_apply(driver,&t,r,y);
	if(status!=GSL_SUCCESS) printf("Error in odeiv2 at r = %g\n",r);

	gsl_odeiv2_driver_free(driver);
	return y[0];
}
