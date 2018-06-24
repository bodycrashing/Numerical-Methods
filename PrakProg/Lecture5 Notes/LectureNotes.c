#include "stdio.h"
#include "math.h"

/*Nyt eksempel på en hjemmelavet funktion*/
double fun(double x){return x;}



int main() {
  double x=1;
  /*double* p=NULL;*/
double (*f)(double);
f=&sin;
/*Dette er den kanoniske måde at kalde funktioner på*/
  printf("f(%g) = %g\n",x,(*f)(x));
/*Man kan dog også gøre følgende*/
printf("sin(1)=%g\n",(&sin)(1));
  /*Med andre ord kan man kalde "pointers" til funktion og
  det vil fortsat virke*/




  return 0;
}
