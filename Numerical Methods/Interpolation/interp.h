#ifndef INTERPOLATION_H /* For multiple includes	*/
#define INTERPOLATION_H
/* Libraries	*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

/* Linspace function 	*/
double* linspace(double a, double b, int n, double u[]);

/* Linear Interpolation 	*/
double linterp(int n, double *x, double *y, double z);
double linterp_integral(int n, double *x, double *y, double z);


/*	Quadratic Interpolation 	*/
typedef struct {int n; double *x, *y, *b, *c;} qspline;
qspline * qspline_alloc(int n, double* x, double *y);
double qspline_eval(qspline *s, double z);
double qspline_integral(qspline *s, double z);
double qspline_deriv(qspline *s, double z);
void qspline_free(qspline *s);


/*  Cubic Interpolation  */
typedef struct {int n; double *x,*y,*b,*c,*d;} cspline;
cspline * cspline_alloc (int n, double *x, double *y);
double cspline_eval (cspline * s, double z);
void cspline_free (cspline * s);
double cspline_integral(cspline *s, double z);
double cspline_deriv(cspline *s, double z);

#endif
