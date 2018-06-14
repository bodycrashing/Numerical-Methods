#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>

//--------------------------------------------------------------------------
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
//--------------------------------------------------------------------------


int orbit (double phi, const double u[], double dudphi[], void *params)
{
  (void)(t);
  double epsilon = *(double *)params;
  dudphi[0] = u[1];
  dudphi[1] = 1 - u[0] + epsilon*u[0]*u[0];
  return GSL_SUCCESS;
}



int main (void)
{
  //Defining designated file streams for each of the ODE's.
  FILE* logistic_stream = fopen("logistic.out.txt","w");
  FILE* orbit_stream = fopen("orbit.out.txt","w");

  double epsabs = 1e-8, epsrel = 1e-8;
  double hstart = 1e-3;

  gsl_odeiv2_system logistic_sys = {logistic_diff, NULL, 1, NULL};
  gsl_odeiv2_driver *logistic_driver = gsl_odeiv2_driver_alloc_y_new(&logistic_sys, gsl_odeiv2_step_rkf45, hstart, epsabs, epsrel);

  gsl_odeiv2_system orbit_sys = {orbit, NULL, 2, &epsilon};
  gsl_odeiv2_driver *orbit_driver = gsl_odeiv2_driver_alloc_y_new(&orbit_sys, gsl_odeiv2_step_rkf45, hstart, epsabs, epsrel);

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

    //  printf ("%.8e\t%.8e\t%.8e\n", i, y[0],logistic(i));
    fprintf(logistic_stream,"%.8e\t%.8e\t%.8e\n", i, y[0],logistic(i));
    }


  double uprime_initial[3] = {0.0,-0.5,-0.5};
  double u_initial = 1.0, t = 0.0;

  int orbit_result;
  for (int j = 0; j < 3; j++) {
    double u[2] = {u_initial,uprime_initial[j]};
    double phi_max = 39.5 * M_PI;
  for (double phi=0.0; phi < phi_max; phi+=0.05) {
      orbit_result = gsl_odeiv2_driver_apply(orbit_driver, &t, phi, u);

      if (orbit_result != GSL_SUCCESS)
        {
          printf ("error: orbit driver returned %d\n", orbit_result);
          break;
        }

      fprintf(orbit_stream, "%.8e\t%.8e\n", i, u[0];
    }
    fprintf(orbit_stream, "\n\n");
}

  //Clearing and closing
  fclose(logistic_stream);
  fclose(orbit_stream);
//  gsl_odeiv2_driver_free (logistic_result);
  gsl_odeiv2_driver_free (logistic_driver);

  return GSL_SUCCESS;
}
