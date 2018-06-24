#include <stdio.h>
#include <complex.h>
#include <math.h>

int main () {
printf("\nComputing of values various functions and constants included in the 'math.h' and 'complex.h' libaries:\n");
/* Gamma-function and Bessel function for x=5 and x=0.5 respectively*/
printf("Gamma(5)=%lg\n", tgamma(5.0));
printf("Bessel_0(0.5)=%lg\n",  j1(0.5));

/* Complex square root of -2 */
double complex c0 = csqrt(-2);
printf("\nImaginary part of sqrt(-2)=%lg\n", cimag(c0));
printf("Real part of sqrt(-2)=%lg\n", creal(c0));

/* Calculation of complex expontial functions */
double complex c1=cpow(M_E,I); //M_E = value of "e"
printf("\nImaginary part of exp(I)=%lg\n", cimag(c1));
printf("Real part of sqrt exp(I)=%lg\n", creal(c1));
double complex c2=cpow(M_E,I*M_PI);
printf("IMAG exp(I*pi)=%lg\n", cimag(c2));
printf("REAL exp(I*pi)=%lg\n", creal(c2));

/* Calculating I^e (I imaginary unit)*/
double complex c3=cpow(I,M_E);
printf("\nImaginary part of e^I=%lf\n", cimag(c3));
printf("Real part of e^I=%lf\n", creal(c3));

/* Testing difference in variables ability to store significant digits i.e.
    differenece between float, double and long double */
long double test_longdouble = 0.1111111111111111111111111111L;
double test_double = 0.1111111111111111111111111111L;
float test_float = 0.1111111111111111111111111111L;

printf("\nDifferenece between float, double and long double:\n");
printf("Long double = %.25Lg\n",test_longdouble);
printf("Double = %.25g\n",test_double);
printf("Float = %.25f\n",test_float);

   return 0;
}
