#include"log_int.h"

int main(void){

	int calls; // number of calls of ode_ln
	printf("x\ty\tln(x)\tcalls\n");
  double a=1.0, b=10, dx = 1;
	for(double x=a; x<b; x+=dx){
		calls = 0;
		double y = log_int(x, &calls);
		printf("%8g \t%8g \t%8g \t%i\n", x, y, log(x), calls);
	}
	printf("\n\n");

	return 0;
}
