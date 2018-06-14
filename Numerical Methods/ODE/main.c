#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<gsl/gsl_matrix.h>
#include"runge_kutta.h"
#define RND (double)rand()/RAND_MAX
#define FMT "%7.6f" //format of print "7 width, 3 digits after comma"


void system1(int n,double x,gsl_vector* y,gsl_vector* dydx){
double y0=gsl_vector_get(y,0);
double y1=gsl_vector_get(y,1);
double f0=y1*x;
double f1=-y0;
gsl_vector_set(dydx,0,f0);
gsl_vector_set(dydx,1,f1);
}


int main(void) {
	// part A and B
int n=2, k, max=1e6; double b=7, h=b/20, acc=1e-26, eps=1e-2;
	gsl_vector* xlist = gsl_vector_alloc(max+1);
	gsl_matrix* ylist = gsl_matrix_alloc(max+1,n);

	// initialize
	gsl_vector_set(xlist,0,0.0);
	gsl_matrix_set(ylist,0,0,1.0);
	gsl_matrix_set(ylist,0,1,0.0);

//k = ode_driver(system1,xlist,ylist,b,h,acc,eps,max);
k=driver(system1,n,xlist,ylist,b,h,acc,eps,max);
/*
FILE* data = fopen("out.txt","w+");
for (int i = 0; i < k; i++) {
	double x_i = gsl_vector_get(xlist,i);
	double y_i = gsl_matrix_get(ylist,i,0);
	fprintf(data, "%g\t%g\t\n\n", x_i, y_i);
}
fclose(data);
*/
for (int i = 0; i < k; i++) {
	double x_i = gsl_vector_get(xlist,i);
	double y_i = gsl_matrix_get(ylist,i,0);
	printf("%g\t%g\t\n\n", x_i, y_i);
}
printf("\nWe got here!\n");



/*
	h=1e-3; acc=1e-6; eps=1e-6; b=M_PI*20;
	gsl_vector_set_zero(xlist);
	gsl_matrix_set_zero(ylist);
	gsl_vector_set(xlist,0,0.0);
	gsl_matrix_set(ylist,0,0,1.0);
	gsl_matrix_set(ylist,0,1,0.5);
	printf("Testing on the system from the orbital equation exercise in PP\n");
	printf("with eps=0.01 (relativistic correction)\n");
	printf("dydx_0= y1\n");
	printf("dydx_1= 1-y0+eps*y0*y0\n");
	printf("eps=%g, acc=%g, max=%i, b=%g, hstart=%g\n",eps,acc,max,b,h);
	printf("Startpoint set to:\n");
	printf("x=%g",gsl_vector_get(xlist,0));
	printf("y set to:\n");
	gsl_matrix_get_row(y,ylist,0);
	printv(y);
	k=ode_driver(orbital_equation,n,xlist,ylist,b,h,acc,eps,max);
	printf("And the solution at y(%g):\n",b);
	gsl_matrix_get_row(y,ylist,k);
	printv(y);


	printf("C ODE as integrator \n\n");

	double integral_function1(double x){
		double a = 2.0;
		int n = 1;
		return pow(x,2*n)*exp(-x*x/(a*a));
	}

	double integral_function2(double x){
		double a = 2.0;
		int n = 3;
		return pow(x,2*n)*exp(-x*x/(a*a));
	}

	double gaussian_theo(int n){
		double a=2.0;
		double fac=1;
		double fac2=1;
		for(int i=2;i<=n;i++){
			fac*=i;
		}
		for(int i=2;i<=2*n;i++){
			fac2*=i;
		}
		double integral = sqrt(M_PI)*fac2/fac*pow((a/2),2*n+1);
		return integral;
	}
	n=1;
	h=1e-3; acc=1e-6; eps=1e-6; b=100;
	gsl_vector* xlist_integ=gsl_vector_alloc(max);
	gsl_matrix* ylist_integ=gsl_matrix_alloc(max,n);
	printf("Performing the gaussian integral x^2*exp(-(x/a)^2) from 0 to 50 (should be infinity)\n");
	printf("eps=%g, acc=%g, max=%i, b=%g, hstart=%g\n",eps,acc,max,b,h);
	printf("Startpoint set to:\n");
	printf("x=%g\n",gsl_vector_get(xlist_integ,0));
	printf("y set to: 0\n");
	double result= integrator(integral_function1,n,xlist_integ,ylist_integ,b,h,acc,eps,max);
	printf("And the integral is %g\n",result);
	printf("Theoretically it should be if upper limit is infinity:\n");
	printf("%g\n",gaussian_theo(1));
	printf("\n");

	n=1;
	h=1e-3; acc=1e-6; eps=1e-6; b=1000;
	gsl_vector_set_zero(xlist);
	gsl_matrix_set_zero(ylist);
	printf("Performing the gaussian integral x^6*exp(-(x/a)^2) from 0 to 1000 (should be infinity)\n");
	printf("eps=%g, acc=%g, max=%i, b=%g, hstart=%g\n",eps,acc,max,b,h);
	printf("Startpoint set to:\n");
	printf("x=%g\n",gsl_vector_get(xlist_integ,0));
	printf("y set to: 0\n");
	result= integrator(integral_function2,n,xlist_integ,ylist_integ,b,h,acc,eps,max);
	printf("And the integral is %g\n",result);
	printf("Theoretically it should be if upper limit is infinity:\n");
	printf("%g\n",gaussian_theo(3));
	printf("\n");

	gsl_vector_free(xlist);
	gsl_matrix_free(ylist);
	gsl_vector_free(xlist_integ);
	gsl_matrix_free(ylist_integ);
	gsl_vector_free(y);*/

	return 0;
}
