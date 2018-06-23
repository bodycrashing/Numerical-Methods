#include "smp.h"
#include<time.h>

double f_gamma(gsl_vector* x){
  double x0 = gsl_vector_get(x,0);
  double x1 = gsl_vector_get(x,1);
  double x2 = gsl_vector_get(x,2);
	return 1.0/(1-cos(x0)*cos(x1)*cos(x2))/(M_PI*M_PI*M_PI);
}

int main(){

  int dim = 3; //Dimension of integration space
  int N = 10e+6; //Number of sampling points
  double err, res;

  gsl_vector* a = gsl_vector_alloc(dim);
  gsl_vector* b = gsl_vector_alloc(dim);
  double exact = 1.3932039296856768591842462603255;

	clock_t begin_time = clock();
  gsl_vector_set_zero(a);
  gsl_vector_set_all(b,M_PI);

  plain_monte(dim, a, b, f_gamma, N, &res, &err);
  printf("\nIntegration of gamma function from 0 to Pi:\n");
  printf("Error = %lg\n", err);
  printf("Result: Q=%lg\nCorrect result: Q_exact=%.15g\n",res,exact);


	clock_t end_time = clock();
	double time_spend = (double)(end_time-begin_time)/CLOCKS_PER_SEC;
	printf("Time spend without Multi threading is %g seconds \n",time_spend);

printf("\n---------------------------------------------------------------------\n");
printf("Test of Plain Monte Carlo Multi-Dimensional Integration With Multi Threading\n");
printf("---------------------------------------------------------------------\n");

	begin_time = clock();

  gsl_vector_set_zero(a);
  gsl_vector_set_all(b,M_PI);

	plain_monte_OMP(dim, a, b, f_gamma, N, &res, &err);

  end_time = clock();
	time_spend = (double)(end_time-begin_time)/CLOCKS_PER_SEC;
	printf("\nTime spend with Multi threading is %g seconds \n",time_spend);




gsl_vector_free(a);
gsl_vector_free(b);
return 0;
}
