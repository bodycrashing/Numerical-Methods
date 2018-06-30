#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<gsl/gsl_matrix.h>
#include"runge_kutta.h"
#include<math.h>
#include <gsl/gsl_sf_erf.h>

void ode_logistic(double x, gsl_vector* f, gsl_vector* dfdx){
	double fx = gsl_vector_get(f,0);
	gsl_vector_set(dfdx,0, fx*(1.0 - fx));
}

void ode_sin(double x, gsl_vector* y, gsl_vector* dydx){
	double y0 = gsl_vector_get(y,0);
	double y1 = gsl_vector_get(y,1);
	gsl_vector_set(dydx,0,y1);
	gsl_vector_set(dydx,1,-y0);
}

void int_cos(double x, gsl_vector* y, gsl_vector* dydx){
	gsl_vector_set(dydx,0, cos(x));	/* dydx = f(x)	*/
}

void int_erf(double x, gsl_vector* y, gsl_vector* dydx){
	gsl_vector_set(dydx,0,1/sqrt(M_PI)*exp(-x*x));
}


int main(void){
	printf("\n------- Part A/B: Embedded Runge-Kutta ODE23 integrator with Path Storing -------\n");
	int n=2, k, max=1e6;
	double a=0, b=4*M_PI, h=0.001, acc=1e-6, eps=1e-6;
	gsl_vector* xlist = gsl_vector_alloc(max); // Corrected max+1
	gsl_matrix* ylist = gsl_matrix_alloc(max,n);

	// initialize
	gsl_vector_set(xlist,0,a);
	gsl_matrix_set(ylist,0,0,1.0);
	gsl_matrix_set(ylist,0,1,0.0);

	printf("\n|------- Harmonic Differential Equation -------|\n");
	k = ode_driver(ode_sin,xlist,ylist,b,h,acc,eps,max);
	printf("y''=-y\nIntegration limits: a=%lg  b=%lg\nNumber of steps=%i\n",a, b, k);
	printf("Result:	      y1=%lg		y2=%lg\n", gsl_matrix_get(ylist,k,1), gsl_matrix_get(ylist,k,0));
	printf("Exact result: sin(b)=%lg	cos(b)=%lg\n",sin(b), cos(b));

	FILE* data = fopen("data.txt","w+");
		for (int i = 0; i < k; i++) {
			double xi = gsl_vector_get(xlist,i);
			double yi = gsl_matrix_get(ylist,i,0);
			fprintf(data,"%g\t%g\t\n", xi, yi);
		} fprintf(data, "\n\n");



		n = 1; a = 0.0; b = 5.0;
		double y0 = 1.0/(1.0 + exp(-a));
		gsl_vector* xlist2 = gsl_vector_alloc(max);
		gsl_matrix* ylist2 = gsl_matrix_alloc(max,n);
		gsl_vector_set(xlist2,0,a);
		gsl_matrix_set(ylist2,0,0,y0);

	printf("\n|------- Logistic Differential Equation -------|\n");
	k = ode_driver(ode_logistic,xlist2,ylist2,b,h,acc,eps,max);
	printf("y'=y*(1-y)\nIntegration limits: a=%lg  b=%lg\nNumber of steps=%i\n",a, b, k);
	printf("Result:	      y=%lg\n", gsl_matrix_get(ylist2,k,0));
	printf("Exact result: 1/(1+exp(-b))=%lg\n",1.0/(1.0 + exp(-b)));
	for (int i = 0; i < k; i++) {
		double xi = gsl_vector_get(xlist2,i);
		double yi = gsl_matrix_get(ylist2,i,0);
		fprintf(data,"%g\t%g\t\n", xi, yi);
	} fprintf(data, "\n\n");


	printf("\n\n------------- Part C: A definite integral as an ODE --------------\n");
	printf("Integration of f(x)=cos(x):\n");
	a = 0; b = M_PI;
	gsl_vector* xlist_int = gsl_vector_alloc(max);
	gsl_matrix* ylist_int = gsl_matrix_alloc(max,n);
	gsl_vector_set(xlist_int,0,a);
	gsl_matrix_set(ylist_int,0,0,0.0);

	k = ode_driver(int_cos,xlist_int,ylist_int,b,h,acc,eps,max);
	printf("Integration limits: a=%lg  b=%lg\nNumber of steps=%i\n",a, b, k);
	printf("Result:	      I=%lg\n", gsl_matrix_get(ylist_int,k,0));
	printf("Exact result: I_exact=%lg\n",sin(b));

	printf("\nIntegration of Error function erf(x):\n");
	a = -M_PI; b = M_PI;
	gsl_vector* xlist_int2 = gsl_vector_alloc(max);
	gsl_matrix* ylist_int2 = gsl_matrix_alloc(max,n);
	gsl_vector_set(xlist_int2,0,a);
	gsl_matrix_set(ylist_int2,0,0,0.0);

	k = ode_driver(int_erf,xlist_int2,ylist_int2,b,h,acc,eps,max);
	printf("Integration limits: a=%lg  b=%lg\nNumber of steps=%i\n",a, b, k);
	printf("Result:	      I=%lg\n", gsl_matrix_get(ylist_int2,k,0));
	printf("Exact result: erf(Pi)=%lg\n",gsl_sf_erf(M_PI));


	// Free memory and closing data streams
	fclose(data);
	gsl_vector_free(xlist);
	gsl_matrix_free(ylist);
	gsl_vector_free(xlist2);
	gsl_matrix_free(ylist2);
	gsl_vector_free(xlist_int);
	gsl_matrix_free(ylist_int);
	gsl_vector_free(xlist_int2);
	gsl_matrix_free(ylist_int2);

	return 0;
}
