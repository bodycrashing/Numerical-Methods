#include"ann.h"
#include "minimization.h"

ann* ann_alloc(int n, double(*act_fun)(double)){
	ann* a = malloc(sizeof(ann));
	a->n = n;
	a->f=act_fun;
	a->data = gsl_vector_alloc(3*n);
	return a;
}

void ann_free(ann* network){
	gsl_vector_free(network->data);
	free(network);
}




double ann_feed_forward(ann* network, double x){
	int n = network->n;
	double sum=0;
	for(int i=0; i<n; i++){
		double ai = gsl_vector_get(network->data,0*network->n+i);
		double bi = gsl_vector_get(network->data,1*network->n+i);
		double wi = gsl_vector_get(network->data,2*network->n+i);

		double arg = (x+ai)/bi;
		double result = network->f(arg)*wi;
		sum+=result;
	}
	return sum;
}





void ann_train(ann* network, gsl_vector* xlist, gsl_vector* ylist,double eps,double step){
	assert(xlist->size == ylist->size);

	int N=xlist->size;

	double func(gsl_vector* p){
		gsl_vector_memcpy(network->data,p);
		double sum=0;
		for(int k=0;k<N;k++){
			double x_k = gsl_vector_get(xlist,k);
			double y_k = gsl_vector_get(ylist,k);
			double F_k = ann_feed_forward(network,x_k);
			sum+=(F_k-y_k)*(F_k-y_k);
		}
		return sum/N;
	}

	gsl_vector* p = gsl_vector_alloc(network->data->size);
	gsl_vector_memcpy(p,network->data);
	quasinewton(func,p,step,eps);
	gsl_vector_memcpy(network->data,p);
	gsl_vector_free(p);
}



ann2d* ann2d_alloc(int n,double(*f)(double)){
	ann2d* a = malloc(sizeof(ann2d));
	a->n = n;
	a->f=	f;
	a->data=gsl_vector_alloc(5*n);
	return a;
}


void ann2d_free(ann2d* network){
	gsl_vector_free(network->data);
	free(network);
}


double ann2d_feed_forward(ann2d* nw, double x_1, double x_2){
	double sum=0;
	for(int i=0;i<nw->n;i++){
		double a1_i = gsl_vector_get(nw->data,i*5);
		double b1_i = gsl_vector_get(nw->data,i*5+1);
		double a2_i = gsl_vector_get(nw->data,i*5+2);
		double b2_i = gsl_vector_get(nw->data,i*5+3);
		double w_i = gsl_vector_get(nw->data,i*5+4);

		double arg1 = (x_1+a1_i)/b1_i;
		double arg2 = (x_2+a2_i)/b2_i;
		double result = nw->f(arg1)*nw->f(arg2)*w_i;

		sum+=result;
	}
	return sum;
}



void ann2d_train(ann2d* nw, gsl_matrix* xlist, gsl_vector* ylist, double eps, double step){
	int N = ylist->size;

	double func(gsl_vector* p){
		gsl_vector_memcpy(nw->data,p);
		double sum=0;
		for(int k=0;k<N;k++){
			double x1_k = gsl_matrix_get(xlist,k,0);
			double x2_k = gsl_matrix_get(xlist,k,1);
			double y_k = gsl_vector_get(ylist,k);

			double F_k = ann2d_feed_forward(nw,x1_k,x2_k);
			sum+=(F_k-y_k)*(F_k-y_k);
		}
		return sum/N;
	}

	gsl_vector* p = gsl_vector_alloc(5*nw->n);
	gsl_vector_memcpy(p,nw->data);

	quasinewton(func,p,step,eps);

	gsl_vector_memcpy(nw->data,p);
	gsl_vector_free(p);
}
