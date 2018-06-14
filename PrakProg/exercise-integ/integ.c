#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

//----------------------------------------------------------
double f (double x) {
  double f = log(x) / sqrt(x);
  return f;
}

//----------------------------------------------------------
double norm_func (double x, void * params) {
  double alpha = *(double *) params;
  double func_val = exp(-alpha*x*x);

  return func_val;
}

double Hamilton_func (double x, void * params) {
  double alpha = *(double *) params;
  double func_val = (-alpha*x*x/2 + alpha/2 + x*x/2)*exp(-alpha*x*x);

  return func_val;
}
//----------------------------------------------------------


int
main (void)
{
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  double result_f, error_f;
  double expected = -4.0;

  gsl_function F;
  F.function = &f;
  F.params = NULL;

  gsl_integration_qags (&F, 0, 1, 0, 10e-7, 1000, w, &result_f, &error_f);
  printf ("result          = % .15f\n", result_f);
  printf ("estimated error = % .15f\n", error_f);
  printf ("actual error    = % .15f\n", result_f - expected);
  printf("\n\n");

  gsl_function Norm;
  Norm.function = &norm_func;

  gsl_function Hamilton;
  Hamilton.function = &Hamilton_func;

  double result_Norm, error_Norm;
  double result_Hamilton, error_Hamilton;

  for (double alpha = 0.01; alpha < 4.5; alpha+=0.05) {
    Norm.params = &alpha;
    Hamilton.params = &alpha;
    gsl_integration_qagi(&Norm, 1e-3, 1e-3, 1000, w, &result_Norm, &error_Norm);
    gsl_integration_qagi(&Hamilton, 1e-3, 1e-3, 1000, w, &result_Hamilton, &error_Hamilton);
    printf("%lg\t%lg\n", alpha, result_Hamilton/result_Norm);
  };

  gsl_integration_workspace_free (w);
  return 0;
}
