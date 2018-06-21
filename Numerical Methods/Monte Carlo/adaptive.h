#ifndef HAVE_ADAPTIVE_H
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include<gsl/gsl_integration.h>

double Adaptive(double f(double x),double a, double b, double acc, double eps, double *err);
double Adaptive24(double f(double x),double a, double b, double acc,
													 double eps, double f2, double f3, int nrec, double* err);
double Adaptive_Infinity(double f(double x),double a, double b, double acc, double eps,double *err);
double Clenshaw_Curtis(double f(double x),double a, double b, double acc, double eps, double *err);

#define HAVE_ADAPTIVE_H
#endif
