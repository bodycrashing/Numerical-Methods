#ifndef HAVE_LEAST_SQUARES
#define HAVE_LEAST_SQUARES

#include "qr.h"
#include "jacobi.h"

void lsfit(int m, double f(int i,double x),
	gsl_vector* x, gsl_vector* y, gsl_vector* dy,
	gsl_vector* c, gsl_matrix* S);
#endif
