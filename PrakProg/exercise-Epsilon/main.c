#include<stdio.h>
#include<limits.h>
#include<float.h>
#include <math.h>

void equal(double a, double b, double tau, double epsilon);

// Delopgave i)
int main(void) {
  printf("INT_MAX = %i\n", INT_MAX);
  int i=1;
  while (i+1>i) {
    i++;
  }
printf("my max int using WHILE-loop = %i\n",i);

int j;
  for (j = 0; j < j+1; j++) {
  }
printf("my max int using FOR-loop = %i\n",j);

int k=1;
  do {
    k++;
  } while(k+1>k);
printf("my max int using DO-WHILE-loop = %i\n",k);

//-------------------------------------------------------
// Delopgave ii)
printf("\nINT_MIN = %i\n", INT_MIN);
i=1;
while (i-1<i) {
  i--;
}
printf("my MIN int using WHILE-loop = %i\n",i);

j=1;
for (j = 1; j-1 < j; j--) {
}
printf("my MIN int using FOR-loop = %i\n",j);

k=1;
do {
  k--;
} while(k-1<k);
printf("my MIN int using DO-WHILE-loop = %i\n",k);

//---------------------------------------------------------
// Delopgave iii)
double x=1;
while(1+x!=1)
{x/=2;
            };
x*=2;
printf("\nx = %le\n", x);
printf("FLT_EPSILON = %le\n", FLT_EPSILON);
printf("DBL_EPSILON = %le\n", DBL_EPSILON);
printf("LDBL_EPSILON = %Lg\n", LDBL_EPSILON);


//-------------------------------------------------------

//Opgave 2
//-------------------------------------------------------
//Delopgave i)
int max=INT_MAX/2;
float sum_up_float = 1.0f;

for (int e = 0; e < max; e++) {
  sum_up_float/=e;
}
printf("\nsum_up_float = %g\n", sum_up_float);


float sum_down_float = 1.0f;
for (int e = 0; e < max; e++) {
  sum_down_float=sum_down_float/(max-e);
}
printf("sum_down_float = %g\n", sum_down_float);

//-------------------------------------------------------
// Delopgave iv)
double sum_up_double = 1.0;

for (int e = 0; e < max; e++) {
  sum_up_double/=e;
}
printf("sum_up_double = %g\n", sum_up_double);


double sum_down_double = 1.0;
for (int e = 0; e < max; e++) {
  sum_down_double=sum_down_double/(max-e);
}
printf("sum_down_double = %g\n", sum_down_double);
printf("\n1.0 delt med max = %g\n", 1.0f/max);

//-------------------------------------------------------
// Delopgave v)
printf("\nCalling function 'equal.c' to test whether to numbers of a given type are equal within tolerannce tau:\n");
equal(3.0, 2.0, 1.0, 1.0);



  return 0;
}
