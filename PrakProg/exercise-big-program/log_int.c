#include "log_int.h"

double f(double x, void * params){
  double f = 1/x;
  return f;
}


double log_int(double x, int* calls){
	*calls += 1;
	assert(x > 0);

	// reducing x
	if(x < 1){return -log_int(1.0/x, calls);}
	if(x >= 2){return log(2.0) + log_int(x/2.0, calls);}

  double result, error;
  double epsrel = 1e-7;
  double epsabs = 0;
  double a=1.0, b=x;
  double limit = 1000;

  gsl_integration_workspace * w = gsl_integration_workspace_alloc(limit);
  gsl_function F;
  F.function = &f;
  F.params = NULL;

  gsl_integration_qags(&F, a, b, epsabs, epsrel, limit, w, &result, &error);
  gsl_integration_workspace_free(w);

  return result;
}
