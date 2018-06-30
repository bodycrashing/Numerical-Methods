#include<gsl/gsl_vector.h>
#include<stdio.h>
#include<math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#ifndef HAVE_NEURONS
#define HAVE_NEURONS

typedef struct {int n; double(*f)(double,double); gsl_vector* data;} ann;

ann* ann_alloc(int n,double(*f)(double,double));

void ann_free(ann* nw);

double ann_feed_forward(ann* nw,double x,double y);

void ann_train(ann* nw,gsl_vector* xlist,gsl_vector* ylist,gsl_vector* flist);

#endif
