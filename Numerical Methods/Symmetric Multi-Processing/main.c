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

  printf("----------- Plain Monte Carlo Integration WITHOUT Multi Threading -------------\n");
	clock_t begin_time = clock();
  gsl_vector_set_zero(a);
  gsl_vector_set_all(b,M_PI);

  plain_monte(dim, a, b, f_gamma, N, &res, &err);
  printf("\nIntegration of gamma function from 0 to Pi:\n");
  printf("Error = %lg\n", err);
  printf("Result: Q=%lg\nCorrect result: Q_exact=%.15g\n",res,exact);
	clock_t end_time = clock();
	double time_spend = (double)(end_time-begin_time)/CLOCKS_PER_SEC;
	printf("Time spend WITHOUT Multi threading: %g seconds \n",time_spend);


  printf("\n----------- Plain Monte Carlo Integration WITH Multi Threading -------------\n");
	begin_time = clock();
  gsl_vector_set_zero(a);
  gsl_vector_set_all(b,M_PI);

	plain_monte_OMP(dim, a, b, f_gamma, N, &res, &err);
  printf("\nIntegration of gamma function from 0 to Pi:\n");
  printf("Error = %lg\n", err);
  printf("Result: Q=%lg\nCorrect result: Q_exact=%.15g\n",res,exact);
  end_time = clock();
	time_spend = (double)(end_time-begin_time)/CLOCKS_PER_SEC;
	printf("Time spend WITH Multi threading: %g seconds \n",time_spend);


  printf("\nAs is evident from the results the multi-threaded version of the plain\nMonte Carlo integration actually takes long time to");
  printf("than the ordinary implementation.\nI have tested other students implentation on my Macbook and in all cases the multithreaded\n");
  printf("are slower. In other words it seems that the time it takes to allocate the 'workers'\nis greater than time gained by actual multithreading.\n");

  gsl_vector_free(a);
  gsl_vector_free(b);

return 0;
}
