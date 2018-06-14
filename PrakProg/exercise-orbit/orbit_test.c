int orbit (double phi, const double y[], double f[], void *params)
{
  (void)(t); /* avoid unused parameter warning */
  double epsilon = *(double *)params;
  f[0] = y[1];
  f[1] = 1 - y[0] + epsilon*y[0]*y[0];
  return GSL_SUCCESS;
}

gsl_odeiv2_system orbit_sys = {orbit, NULL, 2, &epsilon};
gsl_odeiv2_driver *orbit_driver = gsl_odeiv2_driver_alloc_y_new(&orbit_sys, gsl_odeiv2_step_rkf45, hstart, epsabs, epsrel);



int orbit (double phi, const double u[], double dudx[], void *params)
{
  (void)(t); /* avoid unused parameter warning */
  double epsilon = *(double *)params;
  f[0] = y[1];
  f[1] = 1 - y[0] + epsilon*y[0]*y[0];
  return GSL_SUCCESS;
}
