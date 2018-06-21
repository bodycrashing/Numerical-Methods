#include"arti.h"
#include "minimization.h"
#define RND (double)rand()/RAND_MAX

int main(){
	//_____TAST 1_____
	printf("_____TASK 1_____\n\n");

	// Allocate space and define data
	double fit_fun(double x){
		return cos(x)*sin(2*x);
	}

	int n = 20;

	gsl_vector* x_data = gsl_vector_alloc(n);
	gsl_vector* y_data = gsl_vector_alloc(n);

	FILE* data = fopen("plot_data.txt","w+");

	double a=-2, b=2;
	for(int i=0;i<n;i++){
		double x = a+(b-a)*i/(n-1);
		double y = fit_fun(x);

		gsl_vector_set(x_data,i,x);
		gsl_vector_set(y_data,i,y);

		fprintf(data,"%lg\t%lg\n",x,y);
	}
	fprintf(data,"\n\n");

	// Define actication function
	double act_fun(double x){
		return x*exp(-x*x);
	}

	// Artificial Neuron Network
	int neurons = 10;
	double eps = 1e-6, step = 1e-6;

	ann* network = ann_alloc(neurons,&act_fun);
	for(int i=0;i<neurons;i++){
		gsl_vector_set(network->data,0*network->n+i,a+(b-a)*i/(network->n-1));
		gsl_vector_set(network->data,1*network->n+i,1);
		gsl_vector_set(network->data,2*network->n+i,1);
	}

	ann_train(network,x_data,y_data,eps,step);

	// Print data to file
	for(double i=a;i<b;i+=fabs(b-a)/1000){
		double result = ann_feed_forward(network,i);
		fprintf(data,"%lg\t%lg\n",i,result);
	}
	fprintf(data,"\n\n");

	printf("One dimensional ANN interpolation of cos(5x)sin(2x) from -2 to 2 can be seen in plot.svg. Interpolation is made with acitivation function x*exp(-x*x),  %i hidden neurons and eps = %lg.\n",neurons,eps);

	//_____TASK 2_____
	printf("\n_____TASK 2_____\n\n");

	// Allocate space and define data
	double d2_fun(double x,double y){
		//return (x*x+y-11)*(x*x+y-11)+(x+y*y-7)*(x+y*y-7);
		return cos(x+y)*sin(2*(x+y));
		//return (1-x)*(1-x)+pow(y-pow(x,2),2);
	}

	int n2 = 15;

	gsl_matrix* x_data2 = gsl_matrix_alloc(n2*n2,2);
	gsl_vector* y_data2 = gsl_vector_alloc(n2*n2);


	double x_min = 0.0, x_max = M_PI, y_min = 0.0, y_max = M_PI;
	double deltax = fabs(x_max-x_min)/(n2-1), deltay = fabs(y_max-y_min)/(n2-1);
	for(int i=0;i<n2;i++){
		double xx = x_min+i*deltax;
		for(int j=0;j<n2;j++){
			double yy = y_min+j*deltay;
			gsl_matrix_set(x_data2,i*n2+j,0,xx);
			gsl_matrix_set(x_data2,i*n2+j,1,yy);
			gsl_vector_set(y_data2,i*n2+j,d2_fun(xx,yy));
		fprintf(data,"%lg\t%lg\t%lg\n",xx,yy,d2_fun(xx,yy));
		}
	}
	fprintf(data,"\n\n");

	// 2 dim ANN
	neurons = 10;
	step = 1e-3;
	eps = 1e-6;

	ann2d* network2 = ann2d_alloc(neurons,&act_fun);
	for(int i=0;i<neurons;i++){
		gsl_vector_set(network2->data,i*5,RND*i/neurons);
		gsl_vector_set(network2->data,i*5+1,RND*10);
		gsl_vector_set(network2->data,i*5+2,RND*i/neurons);
		gsl_vector_set(network2->data,i*5+3,RND*8);
		gsl_vector_set(network2->data,i*5+4,RND*1);
	}

	ann2d_train(network2,x_data2,y_data2,eps,step);

	// print data to file
	for(double i=x_min;i<x_max;i+=fabs(x_max-x_min)/19){
		for(double j=y_min;j<y_max;j+=fabs(y_max-y_min)/19){
			double result = ann2d_feed_forward(network2,i,j);
			fprintf(data,"%lg\t%lg\t%lg\n",i,j,result);
		}
	}

        printf("Two dimensional ANN interpolation of cos(x+y)  function can be seen in plot.svg. Interpolation is made with acitivation function x*exp(-x*x), %i hidden neurons and eps = %lg. The interpolation doesn't seem to find the correct fucktion but more like the overall trend. With more calculating power, better results could possibly be achieved.\n",neurons,eps);

	fclose(data);


	// Free stuff
	gsl_vector_free(x_data);
	gsl_vector_free(y_data);
	ann_free(network);
	gsl_matrix_free(x_data2);
	gsl_vector_free(y_data2);
	ann2d_free(network2);

}
