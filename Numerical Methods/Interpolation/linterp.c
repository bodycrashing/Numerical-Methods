#include "interp.h"

double linterp (int n, double* x, double* y, double z){
  assert(n>1 && z>=x[0] && z<=x[n-1]);
  int i=0, j=n-1;
  while (j-i>1) {
    int m = (i+j)/2;
    if (z>x[m])
        i=m;
      else
        j=m;
  }
  return y[i] + (y[i+1]-y[i])/(x[i+1]-x[i])*(z-x[i]); /*Add last value of intraged linear interpolation function from x[i] to z*/
}


double linterp_integral(int n, double* x, double* y, double z){
assert(n > 1 && z >= x[0] && z <= x[n-1]);

int i;
double integral = 0.0;
double p;
for(i = 0; x[i+1]<z; i++){
  p = 1/2.0*(y[i+1]-y[i]) / (x[i+1]-x[i]);
  integral += y[i]*(x[i+1]-x[i]) + p*(x[i+1]-x[i])*(x[i+1]-x[i]);
  }
p = 1/2.0*(y[i+1]-y[i]) / (x[i+1]-x[i]);
integral += p*(z-x[i])*(z-x[i]) + y[i]*(z-x[i]);

return integral;
}
