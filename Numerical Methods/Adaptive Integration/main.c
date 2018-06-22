#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include "adaptive.h"
#include<gsl/gsl_integration.h>

int main(void){
printf("|--------- Part A: Adaptive Integration of listed functions: --------|\n");
printf("- f1(x)=sqrt(x)\n- f2(x)=1/sqrt(x)\n- f3(x)=ln(x)/sqrt(x)\n- f4(x)=4*sqrt(1-(1-x)^2)");

double a = 0.;
double b = 1.;
double acc = 1e-5;
double eps = 1e-5;

int Numcalls = 0;
double err = 0.;
double f_sqrt(double x){
  Numcalls++;
  return sqrt(x);
}
double Q = Adaptive_Infinity(f_sqrt,a,b,acc,eps,&err);
printf("\n\nIntegration of f1(x)=sqrt(x) from 0 to 1:\n");
printf("Number of Calls = %i\nError = %lg\n",Numcalls, err);
printf("Result: Q=%lg\nCorrect result: Q_exact=2/3\n",Q);

Numcalls = 0;
err = 0;
double f_invsqrt(double x){
  Numcalls++;
  return 1./sqrt(x);
}
Q = Adaptive(f_invsqrt,a,b,acc,eps,&err);
printf("\nIntegration of f2(x)=1/sqrt(x) from 0 to 1:\n");
printf("Number of Calls = %i\nError = %lg\n",Numcalls, err);
printf("Result: Q=%lg\nCorrect result: Q_exact=2\n",Q);


Numcalls = 0;
err = 0;
double f_lninvsqrt(double x){
  Numcalls++;
  return log(x)/sqrt(x);
}
Q = Adaptive(f_lninvsqrt,a,b,acc,eps,&err);
printf("\nIntegration of f3(x)=ln(x)/sqrt(x)from 0 to 1:\n");
printf("Number of Calls = %i\nError = %lg\n",Numcalls, err);
printf("Result: Q=%lg\nCorrect result: Q_exact=-4\n",Q);


Numcalls = 0;
err = 0;
double f_pi(double x){
  Numcalls++;
  return 4*sqrt(1-(1-x)*(1-x));
}

Q = Adaptive(f_pi,a,b,acc,eps,&err);
printf("\nIntegration of f4(x)=4*sqrt(1-(1-x)^2) from 0 to 1:\n");
printf("Number of Calls = %i\nError = %lg\n",Numcalls, err);
printf("Result: Q=%.20g\nCorrect result: Q_exact=%.20g\n",Q, M_PI);


/*-------------------- Part B:  Infinite limits ---------------------*/
printf("\n\n|-------------------- Part B:  Infinite limits ---------------------|\n");
Numcalls = 0;
err = 0;
b=INFINITY;
double f_exp(double x){
  Numcalls++;
  return exp(-x);
}
Q = Adaptive_Infinity(f_exp,a,b,acc,eps,&err);
printf("Integration of f(x)=exp(-x) from 0 to INFINITY:\n");
printf("Number of Calls = %i\nError = %lg\n",Numcalls, err);
printf("Result: Q=%lg\nCorrect result: Q_exact=1\n",Q);


Numcalls = 0;
err = 0;
a=-INFINITY;
b=INFINITY;
double f_gauss(double x){
  Numcalls++;
	return exp(-x*x);
}
Q = Adaptive_Infinity(f_gauss,a,b,acc,eps,&err);
printf("\nIntegration of gauss-function f(x)=exp(-x^2) from -INFINITY to INFINITY:\n");
printf("Number of Calls = %i\nError = %lg\n",Numcalls, err);
printf("Result: Q=%lg\nCorrect result: Q_exact=%lg <--sqrt(Pi)\n",Q, sqrt(M_PI));

//Comparing with the GSL gsl_integration_qagi integration routine
double gauss_gsl(double x, void* p){
		Numcalls++;
		return exp(-x*x);
	}
	gsl_function F;
  F.function = &gauss_gsl;
  F.params = NULL;
  double result;
  int limit = 1e6;
  gsl_integration_workspace* ws = gsl_integration_workspace_alloc(limit);
  int flag = gsl_integration_qagi(&F,eps,acc,limit,ws,&result,&err);
  printf("\nComparision with the GSL gsl_integration_qagi libary routine:\n");
  printf("f(x)=exp(-x^2) from -INFINITY to INFINITY:\n");
  printf("Number of Calls = %i\nError = %lg\n",Numcalls, err);
  printf("Result: QAGI=%lg\n",result);



/*--------------- Part C: Clenshaw–Curtis variable transformation ----------------*/
printf("\n\n|----------- Part C: Clenshaw–Curtis variable transformation -----------|\n");
Numcalls = 0;
err = 0;
a=0;
b=1;
Q = Clenshaw_Curtis(f_invsqrt,a,b,acc,eps,&err);
printf("Integration of f2(x)=1/sqrt(x) from 0 to 1 using Clenshaw–Curtis transform:\n");
printf("Number of Calls = %i\nError = %lg\n",Numcalls, err);
printf("Result: Q=%lg\nCorrect result: Q_exact=2\n",Q);
printf("\nHence we see a significant reduction in the number of function calls\nas compared to original adaptive integration approach\n");

return 0;
}
