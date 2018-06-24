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


int orbit_diff (double phi, const double u[], double dudphi[], void *params)
{
  double epsilon = *(double *)params;
  dudphi[0] = u[1];
  dudphi[1] = 1 - u[0] + epsilon*u[0]*u[0];
  return GSL_SUCCESS;
}



int main (void){
  FILE* logistic_stream = fopen("logistic.out.txt","w+");
  FILE* orbit_stream = fopen("orbit.out.txt","w+");

  double  epsilon, epsabs = 1e-8, epsrel = 1e-8;
  double hstart = 1e-3;

  gsl_odeiv2_system logistic_sys = {logistic_diff, NULL, 1, NULL};
  gsl_odeiv2_driver *logistic_driver = gsl_odeiv2_driver_alloc_y_new(&logistic_sys, gsl_odeiv2_step_rkf45, hstart, epsabs, epsrel);

  double t = 0.0;
  double t1 = 3.0;
  double y[1] = {0.5};
  int logistic_result;
  for (double i=t; i < t1; i+=0.05) {
      logistic_result = gsl_odeiv2_driver_apply(logistic_driver, &t, i, y);
      if (logistic_result != GSL_SUCCESS)
        {
          printf ("error: logistic driver returned %d\n", logistic_result);
          break;
        }
    fprintf(logistic_stream,"%lg\t%lg\t%lg\n", i, y[0],logistic(i));
    }




    gsl_odeiv2_system orbit;
  	orbit.function = orbit_diff;
  	orbit.jacobian = NULL;
  	orbit.dimension = 2;
  	orbit.params = (void *) &epsilon;

  	double hstart = 1e-3, epsabs = 1e-6, epsrel = 1e-6;
  	double phi_max = 39.5 * M_PI, delta_phi = 0.05;

  	gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&orbit, gsl_odeiv2_step_rk8pd, hstart, epsabs, epsrel);

  	double t = 0, y[2] = { 1, uprime };
  	for (double phi = 0; phi < phi_max; phi += delta_phi) {
  		int status = gsl_odeiv2_driver_apply (driver, &t, phi, y);
  		printf ("%g %g\n", phi, y[0]);
  		if (status != GSL_SUCCESS) fprintf (stderr, "fun: status=%i", status);
    }
  	gsl_odeiv2_driver_free (driver);










  /*
  for (double phi=0.0; phi < phi_max; phi+=0.05) {
      orbit_result = gsl_odeiv2_driver_apply(orbit_driver, &t, phi, u);

      if (orbit_result != GSL_SUCCESS)
        {
          printf ("error: orbit driver returned %d\n", orbit_result);
          break;
        }

      fprintf(orbit_stream, "%i\t%lg\n", i, u[0]);
    }
    fprintf(orbit_stream, "\n\n");
}

  //Freeing memory and closing data-streams
  fclose(logistic_stream);
  fclose(orbit_stream);
  gsl_odeiv2_driver_free (logistic_driver);

*/
  return 0;
}






	gsl_odeiv2_system orbit;
	orbit.function = orbital_equation;
	orbit.jacobian = NULL;
	orbit.dimension = 2;
	orbit.params = (void *) &epsilon;

	double hstart = 1e-3, epsabs = 1e-6, epsrel = 1e-6;
	double phi_max = 39.5 * M_PI, delta_phi = 0.05;

	gsl_odeiv2_driver *driver =
		gsl_odeiv2_driver_alloc_y_new
			(&orbit, gsl_odeiv2_step_rk8pd, hstart, epsabs, epsrel);

	double t = 0, y[2] = { 1, uprime };
	for (double phi = 0; phi < phi_max; phi += delta_phi) {
		int status = gsl_odeiv2_driver_apply (driver, &t, phi, y);
		printf ("%g %g\n", phi, y[0]);
		if (status != GSL_SUCCESS) fprintf (stderr, "fun: status=%i", status);
		}

	gsl_odeiv2_driver_free (driver);
return EXIT_SUCCESS;
}
