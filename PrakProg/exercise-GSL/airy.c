#include <stdio.h>
#include <gsl/gsl_sf_airy.h>
#include <math.h>

int main() {
  double range = 15.0;
  int n = 500; /*Number of points in the plotted interval*/
  for (double x = -range; x < 0; x+=(range/n)) {
    printf("%g\t%g\t%g\n", x, gsl_sf_airy_Ai(x,GSL_PREC_DOUBLE), gsl_sf_airy_Bi(x,GSL_PREC_DOUBLE));
  }

  return 0;
}
