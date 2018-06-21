#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include "adaptive.h"
int main(void){

printf("---------------------------------------------------------------------\n");
printf("Adaptive Integration of various functions\n");
printf("---------------------------------------------------------------------\n");

int Numcalls = 0;
double a = 0.;
double b = 1.;
double acc = 1e-5;
double eps = 1e-5;
double err = 0.;

////////////////////////////////////////////////////////
double f_sqrt(double x){
  Numcalls++;
  return sqrt(x);
}

double Q = Adaptive_Infinity(f_sqrt,a,b,acc,eps,&err);

printf("\n\nIntegration of sqrt(x) from 0 to 1: \n Q=%e \n Error=%lg\n",Q ,err);
printf("Number of Calls =%i\n",Numcalls );
printf("Correct result = 2/3\n");

////////////////////////////////////////////////////////

Numcalls = 0;
err = 0;

double f_invsqrt(double x){
  Numcalls++;
  return 1./sqrt(x);
}

Q = Adaptive_Infinity(f_invsqrt,a,b,acc,eps,&err);

printf("\n\nIntegration of 1/sqrt(x) from 0 to 1: \n Q=%e \n Error=%lg\n",Q ,err);
printf("Number of Calls =%i\n",Numcalls );
printf("Correct result = 2\n");


////////////////////////////////////////////////////////
Numcalls = 0;
err = 0;

double f_lninvsqrt(double x){
  Numcalls++;
  return log(x)/sqrt(x);
}

Q = Adaptive_Infinity(f_lninvsqrt,a,b,acc,eps,&err);

printf("\n\nIntegration of ln(x)/sqrt(x) from 0 to 1: \n Q=%e \n Error=%lg\n",Q ,err);
printf("Number of Calls =%i\n",Numcalls );
printf("Correct result = -4\n");


////////////////////////////////////////////////////////

Numcalls = 0;
err = 0;

double f_mix(double x){
  Numcalls++;
  return 4*sqrt(1-(1-x)*(1-x));
}

Q = Adaptive_Infinity(f_mix,a,b,acc,eps,&err);

printf("\n\nIntegration of 4*sqrt(1-(1-x)*(1-x)) from 0 to 1:\n Q= %e \n Error=%lg\n",Q ,err);
printf("Number of Calls =%i\n",Numcalls );
printf("Correct result = PI\n");


////////////////////////////////////////////////////////

Numcalls = 0;
err = 0;
b=INFINITY;

double f_exp(double x){
  Numcalls++;
  return exp(-x);
}

Q = Adaptive_Infinity(f_exp,a,b,acc,eps,&err);

printf("\n\nIntegration of exp(-x) from 0 to Infinity:\n Q= %e \n Error=%lg\n",Q ,err);
printf("Number of Calls =%i\n",Numcalls );
printf("Correct result = 1\n");


////////////////////////////////////////////////////////

Numcalls = 0;
err = 0;
a=-INFINITY;
b=INFINITY;

double f_fancy(double x){
  Numcalls++;
  return 1/sqrt(pow(x,12)+5);
}

Q = Adaptive_Infinity(f_fancy,a,b,acc,eps,&err);

printf("\n\nIntegration of 1/sqrt(pow(x,12)+5) from -Infinity to Infinity:\n Q=%e \n Error=%lg\n",Q ,err);
printf("Number of Calls =%i\n",Numcalls );
printf("Correct result = 1.1765...\n");

printf("---------------------------------------------------------------------\n");
printf("Adaptive Integration of various functions using Clenshaw-Curtis Method\n");
printf("---------------------------------------------------------------------\n");

Numcalls = 0;
err = 0;
a=0;
b=1;

Q = Clenshaw_Curtis(f_invsqrt,a,b,acc,eps,&err);

printf("\n\nIntegration of 1/sqrt(x) from 0 to 1: \n Q=%e \n Error=%lg\n",Q ,err);
printf("Number of Calls =%i\n",Numcalls );
printf("Correct result = 2\n");

return 0;
}
