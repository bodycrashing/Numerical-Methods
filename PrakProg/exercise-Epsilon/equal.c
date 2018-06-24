#include <stdio.h>
#include <math.h>

void equal(double a, double b, double tau, double epsilon){

  if (fabs(a-b) < tau){
		printf("|a-b| < tau:\nreturning 1\n");

	}else if(fabs(a-b)/(fabs(a)+fabs(b)) < epsilon){
		printf("|a-b|/(|a|+|b|) < eps/2:\nreturning 1\n");
	}else{
		printf("Numbers are not equal:\nreturning 0\n");
	}

	return;
}
