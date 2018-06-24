#include"komplex.h"
#include<stdio.h>
#include <math.h>
#define TINY 1e-6

int main(){
  printf("--------- \nGhe following test the function libary 'komplex.h': ---------\n");

  komplex z1 = komplex_new(1,2); komplex_print("z1=",z1);
  komplex z2 = komplex_new(3,4); komplex_print("z2=",z2);

	komplex add = komplex_add(z1,z2);
  printf("\nTest addition-function:\n");
  komplex_print("z1+z2=",add);

  komplex sub = komplex_sub(z1,z2);
  printf("\nTest subtraction-function:\n");
  komplex_print("z1-z2=",sub);

  komplex divi = komplex_div(z1,z2);
  printf("\nTest division-function:\n");
  komplex_print("z1/z2=",divi);


return 0;
}
