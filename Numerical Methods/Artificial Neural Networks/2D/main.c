#include<stdio.h>
#include<math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include"ann.h"

double train_func(double x, double y){return y*x*exp(-x*x)*exp(-y*y);}
//double train_func(double x){return x*x*exp(-x*x);}
double func_to_interpolate(double x,double y){return cos(5*x-1)*exp(-x*x)*exp(-y*y);}
//double func_to_interpolate(double x,double y){return x*x + y*y;}
int main(){

	int num_neurons=20;

	ann* nw=ann_alloc(num_neurons,train_func);
	double a=-1,b=1;
	int num_x=20;
	gsl_vector* vx=gsl_vector_alloc(num_x*num_x);
	gsl_vector* vy=gsl_vector_alloc(num_x*num_x);
	gsl_vector* vf=gsl_vector_alloc(num_x*num_x);
	for(int i=0;i<num_x;i++){
		for (int j = 0; j < num_x; j++) {
			double x=a+(b-a)*i/(num_x-1);
			double y=a+(b-a)*j/(num_x-1);
			double f=func_to_interpolate(x,y);
			gsl_vector_set(vx,i*num_x+j,x);
			gsl_vector_set(vy,i*num_x+j,y);
			gsl_vector_set(vf,i*num_x+j,f);
		}

	}


	for(int i=0;i<nw->n;i++){
		gsl_vector_set(nw->data,0*nw->n+i,a+(b-a)*i/(nw->n-1));
		gsl_vector_set(nw->data,1*nw->n+i,1);
		gsl_vector_set(nw->data,2*nw->n+i,a+(b-a)*i/(nw->n-1));
		gsl_vector_set(nw->data,3*nw->n+i,1);
		gsl_vector_set(nw->data,4*nw->n+i,1);
	}

	ann_train(nw,vx,vy,vf);

	for(int i=0;i<vx->size;i++){
		double x=gsl_vector_get(vx,i);
		double y=gsl_vector_get(vy,i);
		double f=gsl_vector_get(vf,i);
		printf("%g %g %g\n",x,y,f);
	}
	printf("\n\n");

  double x_min = -1.0, x_max = 1.0, y_min = -1.0, y_max = 1.0;
  for(double i = x_min; i < x_max; i+=fabs(x_max-x_min)/19) {
    for(double j=y_min; j<=y_max; j+=fabs(y_max-y_min)/19){
      double y=ann_feed_forward(nw,i,j);
      printf("%g %g %g\n",i,j,y);
    }
  }

ann_free(nw);
gsl_vector_free(vx);
gsl_vector_free(vy);
gsl_vector_free(vf);
return 0;
}
