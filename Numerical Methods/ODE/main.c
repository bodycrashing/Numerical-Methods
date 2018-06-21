#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<gsl/gsl_matrix.h>
#include"runge_kutta.h"
#include<math.h>
#define RND (double)rand()/RAND_MAX
#define FMT "%7.6f" //format of print "7 width, 3 digits after comma"


void system1(double x,gsl_vector* y, gsl_vector* dydx){
double y0=gsl_vector_get(y,0);
double y1=gsl_vector_get(y,1);
double f0=y1*x;
double f1=-y0;
gsl_vector_set(dydx,0,f0);
gsl_vector_set(dydx,1,f1);
}


void sine_diff_eq(double x, gsl_vector* y, gsl_vector* dydx){
	double y0 = gsl_vector_get(y,0);
	double y1 = gsl_vector_get(y,1);
	gsl_vector_set(dydx,0,y1);
	gsl_vector_set(dydx,1,-y0);
	//fprintf(stderr, "\n%lg",y0);
}


void int1(double x, gsl_vector* y, gsl_vector* dydx){
	gsl_vector_set(dydx,0,x);
}



int main(void) {
	// part A and B
int n=2, k, max=1e6; double a=0,b=1.0, h=0.001, acc=1e-6, eps=1e-6;
	gsl_vector* xlist = gsl_vector_alloc(max+1);
	gsl_matrix* ylist = gsl_matrix_alloc(max+1,n);

	// initialize
	gsl_vector_set(xlist,0,a);
	gsl_matrix_set(ylist,0,0,1.0);
	gsl_matrix_set(ylist,0,1,0.0);

k = ode_driver(sine_diff_eq,xlist,ylist,b,h,acc,eps,max);
printf("\n%i\n", k);

/*
FILE* data = fopen("out.txt","w+");
fclose(data);
*/

for (int i = 0; i < k; i++) {
	double x_i = gsl_vector_get(xlist,i);
	double y_i = gsl_matrix_get(ylist,i,0);
	printf("%g\t%g\t\n\n", x_i, y_i);
}

fprintf(stderr,"\nWe got here!\n");

n=1;
a = 0; b = 5; eps = 1e-6;

gsl_vector* xlist_int = gsl_vector_alloc(max+1);
gsl_matrix* ylist_int = gsl_matrix_alloc(max+1,n);
gsl_vector_set(xlist_int,0,a);
gsl_matrix_set(ylist_int,0,0,0);

double res=integrator(int1,xlist_int,ylist_int,b,h,acc,eps,max);

printf("\n%lg\n", res);
for (int i = max-20; i < max; i++) {
fprintf(stderr, "%lg\n",gsl_matrix_get(ylist_int,i,0));
}


gsl_vector_free(xlist);
gsl_matrix_free(ylist);
gsl_vector_free(xlist_int);
gsl_matrix_free(ylist_int);

	return 0;
}
