#ifndef HAVE_ARTI_H
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<stdio.h>
#include<assert.h>
#include<math.h>

typedef struct{int n; double (*f)(double);gsl_vector* data;} ann;
typedef struct{int n; double (*f)(double); gsl_vector* data;} ann2d;

void vector_print(const gsl_vector* V);
void matrix_print(const gsl_matrix* M);

void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R);
void qr_gs_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);
int quasinewton(double f(gsl_vector* x), gsl_vector* x, double dx, double eps);

ann* ann_alloc(int n, double(*act_fun)(double));
void ann_free(ann* network);
double ann_feed_forward(ann* network, double x);
void ann_train(ann* network, gsl_vector* xlist, gsl_vector* ylist,double eps,double step);

ann2d* ann2d_alloc(int n,double(*f)(double));
void ann2d_free(ann2d* network);
double ann2d_feed_forward(ann2d* nw, double x_1, double x_2);
void ann2d_train(ann2d* nw, gsl_matrix* xlist, gsl_vector* ylist, double eps, double step);

#define HAVE_ARTI_h
#endif
