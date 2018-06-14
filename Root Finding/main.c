#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<stdio.h>

int newton (
	void f(gsl_vector* x,gsl_vector* fx),
	gsl_vector* x, double dx, double eps);

void vector_print(char* s,gsl_vector* v){
	printf("%s",s);
	for(int i=0;i<v->size;i++)printf("%10.4g ",gsl_vector_get(v,i));
	printf("\n");
	}

int main() {
	int ncalls=0;
	void f(gsl_vector* p,gsl_vector* fx){
		ncalls++;
		double x=gsl_vector_get(p,0), y=gsl_vector_get(p,1);
		gsl_vector_set(fx,0, 2*(1-x)*(-1)+100*2*(y-x*x)*(-1)*2*x);
		gsl_vector_set(fx,1, 100*2*(y-x*x));
		}
	gsl_vector* x=gsl_vector_alloc(2);
	gsl_vector_set(x,0,-2);
	gsl_vector_set(x,1,8);
	gsl_vector* fx=gsl_vector_alloc(2);
	printf("Root finding:\n");
	printf("extremum of the Rosenbrock's function:\n");
	vector_print("initial vector x: ",x);
	f(x,fx);
	vector_print("            f(x): ",fx);
	newton(f,x,1e-6,1e-3);
	printf("ncalls = %i\n",ncalls);
	gsl_vector_fprintf(stderr,x,"%g");
	vector_print("      solution x: ",x);
	f(x,fx);
	vector_print("            f(x): ",fx);
}
