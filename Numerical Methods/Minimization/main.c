#include "minimization.h"
#define FMT "%7.3f"
#define max_print 6

void matrix_print(gsl_matrix *A){
	for(int r=0;r<A->size1;r++){
		for(int c=0;c<A->size2;c++)printf(FMT,gsl_matrix_get(A,r,c));
		printf("\n");}}

void vector_print(gsl_vector *v){
	for(int i=0;i<v->size;i++) printf(FMT,gsl_vector_get(v,i));
	printf("\n");
	}


double Rosenbrock_Hessian(gsl_vector* x0, gsl_vector* grad, gsl_matrix* H){
  double x = gsl_vector_get(x0,0);
  double y = gsl_vector_get(x0,1);
  double fx = (1-x)*(1-x)+100*(y-x*x)*(y-x*x);

  gsl_vector_set(grad,0,2.0*x-2.0+400.0*(x*x*x-y*x));
  gsl_vector_set(grad,1,200.0*(y-x*x));

  gsl_matrix_set(H,0,0,2.0+100.0*(12*x*x-4*y));
  gsl_matrix_set(H,0,1,-400.0*x);
  gsl_matrix_set(H,1,0,-400*x);
  gsl_matrix_set(H,1,1,200.0);
  return fx;
}


double Rosenbrock(gsl_vector*v){
	double x=gsl_vector_get(v,0);
	double y=gsl_vector_get(v,1);
	return pow(1-x,2)+100*pow(y-x*x,2);
}


double Himmelblau_Hessian(gsl_vector* x0, gsl_vector* grad, gsl_matrix* H){
  double x=gsl_vector_get(x0,0);
  double y=gsl_vector_get(x0,1);
  double fx = (x*x+y-11)*(x*x+y-11)+(x+y*y-7)*(x+y*y-7);

  gsl_vector_set(grad,0,4*x*(x*x+y-11)+2*(x+y*y-7));
  gsl_vector_set(grad,1,2*(x*x+y-11)+4*y*(x+y*y-7));

  gsl_matrix_set(H,0,0,4*(3*x*x+y-11)+2);
  gsl_matrix_set(H,0,1,4*x+4*y);
  gsl_matrix_set(H,1,0,4*x+4*y);
  gsl_matrix_set(H,1,1,2+4*(x+3*y*y-7));
  return fx;
}


double Himmelblau(gsl_vector*v){
	double x=gsl_vector_get(v,0);
	double y=gsl_vector_get(v,1);
	return pow(x*x+y-11,2)+pow(x+y*y-7,2);
}


double fit_fun(gsl_vector* v){
	double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
	double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
	double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
	int N = sizeof(t)/sizeof(t[0]);

	double A = gsl_vector_get(v,0);
	double B = gsl_vector_get(v,1);
	double T = gsl_vector_get(v,2);
	double fx=0;
	for(int i=0;i<N;i++){
		fx += pow((A*exp(-t[i]/T)+B-y[i])/e[i],2);
	}
	return fx;
}


int main(int argc, char const *argv[]) {
	/*--------------- Part A: Newton's minimization method---------------------*/
	printf("\n|-------------- Part A: Newton's minimization method -------------------|\n");
	double eps = 10e-3;

	printf("Minimisation for Rosenbrock's Valley Function:\n");
  gsl_vector* x = gsl_vector_alloc(2);
  gsl_vector_set(x,0,2);
  gsl_vector_set(x,1,2);
  printf("Initial guess x0:\n");
  vector_print(x);

	int iter_ros = newton_barebones(Rosenbrock_Hessian, x, eps);
	printf("\n Number of steps = %i\n",iter_ros);
	printf("Computed roots of Rosenbrok's Valley Function:\n");
	vector_print(x);

	printf("\n\nMinimisation of Himmelblau's Function:\n");
	gsl_vector_set(x,0,-4);
	gsl_vector_set(x,1,-4);
	int iter_him = newton_barebones(Himmelblau_Hessian, x, eps);
	printf("\n Number of steps = %i\n",iter_him);
	printf("Computed roots of Himmelblau's Function:\n");
	vector_print(x);



	/*--------------- Part B: Quasi-Newton with Broyden's update -------------------*/
	printf("\n\n|-------------- Part B: Quasi-Newton with Broyden's update -------------------|\n");
	eps=10e-7;
	double dx = 1e-7;
	int nsteps;

	// part ii.)
	printf("\nQuasi-Newton: Rosenbrok's Function\n");
	gsl_vector_set(x,0,2);
	gsl_vector_set(x,1,2);
	printf("Initial guess x0:\n");
	vector_print(x);
	nsteps = quasinewton(Rosenbrock,x,dx,eps);
	printf("Steps=%i\n",nsteps);
	printf("Roots:\n"); vector_print(x);
	printf("f(x_root) = %lg\n", Rosenbrock(x));

	printf("\nQuasi-Newton: Himmelblau's Function\n");
	gsl_vector_set(x,0,18);
	gsl_vector_set(x,1,15);
	printf("Initial guess x0:\n");
	vector_print(x);
	nsteps=quasinewton(Himmelblau,x,dx,eps);
	printf("Steps=%i\n",nsteps);
	printf("Roots:\n"); vector_print(x);
	printf("f(x_root) = %lg\n", Himmelblau(x));
	printf("\n\n");

	// part iii.)
	double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
	double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
	double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
	int N = sizeof(t)/sizeof(t[0]);

	FILE* data = fopen("fit_data.txt","w+"); // Save data to file "data"
	for(int i=0;i<N;i++){
	fprintf(data,"%lg\t%lg\t%lg\n",t[i],y[i],e[i]);
	}
	fprintf(data,"\n");

	// Find constants of fit using Quasi_Newton minimization
	dx=1e-6;
	gsl_vector* x_fit = gsl_vector_alloc(3);
	gsl_vector_set(x_fit,0,0.5);
	gsl_vector_set(x_fit,1,1.0);
	gsl_vector_set(x_fit,2,2.5);

	printf("\n---- Expontial fit to supplied experimental data  ----\n");
	printf("Initial guess of fitting coefficients:\n");
	printf("A=%lg  B=%lg  T=%lg\n", gsl_vector_get(x_fit,0), gsl_vector_get(x_fit,1), gsl_vector_get(x_fit,2));
	int step = quasinewton(fit_fun,x_fit,dx,eps);
	printf("\nEstimate of fitting after Quasi-Newton Least Squares minimization:\n");
	printf("A=%lg  B=%lg  T=%lg\n", gsl_vector_get(x_fit,0), gsl_vector_get(x_fit,1), gsl_vector_get(x_fit,2));
	printf("Steps=%i \n\n",step);

	// Print fit to data file
	double t0=t[0], t_max = t[N-1], delta_t = t_max/200, A, B, T;
	 A = gsl_vector_get(x_fit,0);
	 B = gsl_vector_get(x_fit,1);
	 T = gsl_vector_get(x_fit,2);
	for(double t=t0; t<=t_max; t+=delta_t){
	fprintf(data,"%lg\t%lg\n",t,A*exp(-t/T)+B);
	}
	fclose(data);

// Free allocated memory
gsl_vector_free(x);
gsl_vector_free(x_fit);

  return 0;
}
