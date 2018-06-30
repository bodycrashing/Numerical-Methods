#include<gsl/gsl_vector.h>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>
#include"ann.h"

ann* ann_alloc(int n,double(*f)(double, double)){
	ann* nw = malloc(sizeof(ann));
	nw->n=n;
	nw->f=f;
	nw->data=gsl_vector_alloc(5*n);
	return nw;
}

void ann_free(ann* nw){
	gsl_vector_free(nw->data);
	free(nw);
}

double ann_feed_forward(ann* nw,double x, double y){
	double s=0;
	for(int i=0;i<nw->n;i++){
		double a=gsl_vector_get(nw->data,0*nw->n+i);
		double b=gsl_vector_get(nw->data,1*nw->n+i);
		double c=gsl_vector_get(nw->data,2*nw->n+i);
		double d=gsl_vector_get(nw->data,3*nw->n+i);
		double w=gsl_vector_get(nw->data,4*nw->n+i);
		s+=nw->f((x+a)/b,(y+c)/d)*w;
	}
	return s;
}


void ann_train(ann* nw,gsl_vector* xlist,gsl_vector* ylist,gsl_vector* flist){
	double delta(const gsl_vector* p, void* params){
		gsl_vector_memcpy(nw->data,p);
		double s=0;
		for(int i=0;i<xlist->size;i++){
			double x=gsl_vector_get(xlist,i);
			double y=gsl_vector_get(ylist,i);
			double f=gsl_vector_get(flist,i);
			double yy=ann_feed_forward(nw,x,y);
			s+=fabs(yy-f);
		}
		return s/xlist->size;
	}
	gsl_vector* p=gsl_vector_alloc(nw->data->size);
	gsl_vector_memcpy(p,nw->data);

	gsl_vector* step_size = gsl_vector_alloc(nw->data->size);
	gsl_vector_set_all(step_size, 0.1);

	gsl_multimin_function F;
	F.n = p->size;
	F.f = delta;
	F.params = NULL;

	gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, F.n);
	gsl_multimin_fminimizer_set (s, &F, p, step_size);

	int iter = 0, status;
	do{
		iter++;
		int iteration_status = gsl_multimin_fminimizer_iterate(s);
			if(iteration_status != 0)
			{
			fprintf(stderr, "No Improvement\n");
			break;
			}

		double acc = 0.001;
		status = gsl_multimin_test_size(s->size, acc);

			}while( status == GSL_CONTINUE && iter < 1e+6);
/*
			for (int i = 0; i < nw->n; i++)
			{
			fprintf(stderr, "%g \t %g \t %g \t %g \t %g \n",
			gsl_vector_get(nw->data, 0*nw->n + i),
			gsl_vector_get(nw->data, 1*nw->n + i),
			gsl_vector_get(nw->data, 2*nw->n + i),
			gsl_vector_get(nw->data, 3*nw->n + i),
			gsl_vector_get(nw->data, 4*nw->n + i) );
			}
*/
	gsl_vector_memcpy(nw->data, s->x);
	gsl_vector_free(p);
	gsl_vector_free(step_size);
	gsl_multimin_fminimizer_free(s);

}
