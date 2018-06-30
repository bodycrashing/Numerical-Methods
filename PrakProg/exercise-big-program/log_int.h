#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_integration.h>
#ifndef HAVE_LOG_INT
#define HAVE_LOG_INT

double f(double x, void * params);

double log_int(double x, int* calls);

#endif
