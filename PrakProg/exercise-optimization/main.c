#include<gsl/gsl_multimin.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_errno.h>
#include<stdio.h>
#include<math.h>

double rosenbrock(const gsl_vector* p0, void* params){
  double x = gsl_vector_get(p0,0);
  double y = gsl_vector_get(p0,1);
  double fx = (1-x)*(1-x)+100*(y-x*x)*(y-x*x);

  return fx;
}

struct experimental_data {int n; double *t,*y,*e;};
double function_to_minimize(const gsl_vector* x, void* params){
	double A=gsl_vector_get(x,0);
	double T = gsl_vector_get(x,1);
	double B = gsl_vector_get(x,2);

	struct experimental_data *p = (struct experimental_data*)params;
	int n=p->n;
	double *t=p->t;
	double *y=p->y;
	double *e=p->e;
	double sum = 0;
  #define model(t) A*exp(-(t)/T)+B
	for(int i=0;i<n;i++) sum+=pow((model(t[i])-y[i])/e[i],2);
	return sum;
}


int main(){
// Find minimum of Rosenbrack function by minimization tools i gsl
	gsl_multimin_function F;
	F.f = rosenbrock;
	F.n = 2;
	F.params = NULL;

	gsl_vector* start = gsl_vector_alloc(F.n);
	gsl_vector_set(start,0,0);
	gsl_vector_set(start,1,0);

	gsl_vector* step = gsl_vector_alloc(F.n);
	gsl_vector_set(step,0,0.1);
	gsl_vector_set(step,1,0.1);

	gsl_multimin_fminimizer* state = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,F.n);
	gsl_multimin_fminimizer_set(state,&F,start,step);

	gsl_vector* current = gsl_vector_alloc(2);
	int iter=0, status;
	double current_x, current_y;
	do{
		iter++;
		gsl_multimin_fminimizer_iterate(state);
		status = gsl_multimin_test_size(state->size,1e-10);
		current = gsl_multimin_fminimizer_x(state);
		current_x = gsl_vector_get(current,0);
		current_y = gsl_vector_get(current,1);
		printf("%g %g %g\n",current_x,current_y,
				gsl_multimin_fminimizer_minimum(state));
	if(status==GSL_SUCCESS)break;
	}while(1);

	printf("# Iterations needed for error of 1e-10: %i\n",iter);
	printf("# min at (x,y) = (%g,%g)\n\n\n"
			,gsl_vector_get(state->x,0),gsl_vector_get(state->x,1));

	gsl_vector_free(start);
	gsl_vector_free(step);
	gsl_vector_free(current);



// Fit to data
	double t[]= {0.47,1.41,2.36,3.30,4.24,5.18,6.13,7.07,8.01,8.95};
	double y[]= {5.49,4.08,3.54,2.61,2.09,1.91,1.55,1.47,1.45,1.25};
	double e[]= {0.26,0.12,0.27,0.10,0.15,0.11,0.13,0.07,0.15,0.09};
	int n = sizeof(t)/sizeof(t[0]);

	printf("# t[i], y[i], e[i]\n");
	for(int i=0;i<n;i++) printf("%g %g %g\n",t[i],y[i],e[i]);
	printf("\n\n");

	struct experimental_data params;
	params.n=n;
	params.t=t;
	params.y=y;
	params.e=e;

	gsl_multimin_function FUN;
	FUN.f = function_to_minimize;
	FUN.n = 3;
	FUN.params = (void*)&params;

	gsl_vector* START = gsl_vector_alloc(FUN.n);
	gsl_vector_set(START,0,0);
	gsl_vector_set(START,1,0);
	gsl_vector_set(START,2,0);

	gsl_vector* STEP = gsl_vector_alloc(FUN.n);
	gsl_vector_set(STEP,0,1);
	gsl_vector_set(STEP,1,1);
	gsl_vector_set(STEP,2,1);

	gsl_multimin_fminimizer* M  = gsl_multimin_fminimizer_alloc
		(gsl_multimin_fminimizer_nmsimplex2,FUN.n);
	gsl_multimin_fminimizer_set(M,&FUN,START,STEP);

	int ITER, stat;
	do{
		ITER++;
		gsl_multimin_fminimizer_iterate(M);
		stat = gsl_multimin_test_size(M->size,1e-10);
		if(stat==GSL_SUCCESS)break;
	}while(1);

	double A = gsl_vector_get(M->x,0);
	double T = gsl_vector_get(M->x,1);
	double B = gsl_vector_get(M->x,2);

	for(double x=t[0];x<=t[n-1];x+=0.01){
		printf("%g %g\n",x,A*exp(-x/T)+B);
	}
	printf("# A=%g, T=%g, B=%g\n",A,T,B);

	gsl_vector_free(START);
	gsl_vector_free(STEP);
	gsl_multimin_fminimizer_free(M);
	gsl_multimin_fminimizer_free(state);
	return 0;
}
