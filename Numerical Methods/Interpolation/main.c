#include "interp.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>



int main(int argc, char const *argv[]) {

/*---------------------- Part A: Linear spline ------------------------*/
  int n = 20; //Number of known function values
  double vector_n[n];
  double* x = linspace(0, 5*M_PI, n, vector_n); //Creating a linspace array contaning the "n" x-values
  double y[n]; // Allocating vector for corresponding y-values (function values)

  //Intilizaing the allocated vectors
  printf("x\ty\n");
  for (int i = 0; i < n; i++) {
    y[i] = cos(x[i]);
    printf("%.5lg\t%.5lg\n", x[i], y[i]);
  }
  printf("\n\n");

  int m = 150; //Number of total interpolation points
  double vector_m[m];
  double* z = linspace(0, 5*M_PI, m, vector_m);
  printf("z\tspline\tint_spline\tcos(z)\n");
  for (int i = 0; i < m; i++) {
    double result = linterp(n, x, y, z[i]);
    double result_int = linterp_integral(n, x, y, z[i]);
    printf("%.5lg\t%.5lg\t%.5lg\t%.5lg\n", z[i], result, result_int, sin(z[i]));
  } printf("\n\n");



/* --------------Part B: Quadratic Interpolation--------------- */
qspline* s = qspline_alloc(n, x, y);
printf("z\tqspline\tqspline_int\tqspline_deriv\n");
for(int i=0; i<m; i++){
  double qspline = qspline_eval(s, z[i]);
  double qspline_int = qspline_integral(s, z[i]);
  double qspline_derivative = qspline_deriv(s, z[i]);
  printf("%.5lg\t%.5lg\t%.5lg\t%.5lg\n", z[i], qspline, qspline_int, qspline_derivative);
} printf("\n\n");

qspline_free(s);



/*  --------------Cubic Interpolation---------------  */
int nc=30;
double vector_c[nc];
double yc[nc];

double *xc = linspace(0, 5*M_PI, nc, vector_c);
for (int i = 0; i < nc; i++) {
  yc[i] = cos(2*xc[i])*sin(xc[i]);
}

cspline* c = cspline_alloc(nc, xc, yc);
printf("z\tcspline\tcspline_int\tcspline_deriv\n");
for(int i=0; i<m; i++){
  double cubic_spline = cspline_eval(c, z[i]);
  double cspline_int = cspline_integral(c, z[i]);
  double cspline_derivative = cspline_deriv(c, z[i]);
  printf("%.5lg\t%.5lg\t%.5lg\t%.5lg\n", z[i], cubic_spline, cspline_int, cspline_derivative);
} printf("\n\n");

/* Comparision with GSL Cubic Spline */
gsl_interp_accel *acc = gsl_interp_accel_alloc ();
gsl_spline *GSL_cspline = gsl_spline_alloc (gsl_interp_cspline, nc);
gsl_spline_init (GSL_cspline, xc, yc, nc);
printf("z\tGSL_cspline\n");
for (int i = 0; i < m; i++) {
  double GSL_cspline_res = gsl_spline_eval (GSL_cspline, z[i], acc);
  printf("%.5lg\t%.5lg\n", z[i], GSL_cspline_res);
}

cspline_free(c);
gsl_spline_free(GSL_cspline);
gsl_interp_accel_free(acc);

  return 0;
}
