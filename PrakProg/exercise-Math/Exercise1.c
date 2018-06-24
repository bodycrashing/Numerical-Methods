#include<math.h>
#include<stdio.h>
#include<stdlib.h>
int main () {
int bits = sizeof(void*)*8;
printf("\nThis is (probably) a %i bit system...\n\n",bits);
printf ("1.0f/9 float       : %.25g\n",  1.0f / 9);
printf ("1.0 /9 double      : %.25lg\n", 1.0  / 9);
printf ("1.0L/9 long double : %.25Lg\n", 1.0L / 9);
printf ("                      .123456789|123456789|123456789|\n");
double x = 7/2.0;
printf("double divided by int %g\n", x);


return 0;
}
