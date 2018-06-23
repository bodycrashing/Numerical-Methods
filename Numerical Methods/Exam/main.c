
#include"rk5.h"
#define RND (double)rand()/RAND_MAX
#define FMT "%7.6f" //format of print "7 width, 3 digits after comma"


void sine_eq(double x, gsl_vector* y, gsl_vector* dydx){
	double y0 = gsl_vector_get(y,0);
	double y1 = gsl_vector_get(y,1);
	gsl_vector_set(dydx,0,y1);
	gsl_vector_set(dydx,1,-y0);
}


int main(void) {
	int n=2, k, max=1e6;
	double a=0, b=2*M_PI, h=0.001, acc=1e-6, eps=1e-6;
	gsl_vector* xlist = gsl_vector_alloc(max); // Corrected max+1
	gsl_matrix* ylist = gsl_matrix_alloc(max,n);

	// initialize
	gsl_vector_set(xlist,0,a);
	gsl_matrix_set(ylist,0,0,1.0);
	gsl_matrix_set(ylist,0,1,0.0);


	printf("\n|------- Solving Harmonic Differential Equation -------|\n");
	k = ode_driver(sine_eq,rkstep5,xlist,ylist,b,h,acc,eps,max);
	printf("The differential equation was integrated in:  steps=%i\n", k);

	FILE* data = fopen("data.txt","w+");
		for (int i = 0; i < k; i++) {
			double xi = gsl_vector_get(xlist,i);
			double yi = gsl_matrix_get(ylist,i,0);
			fprintf(data,"%g\t%g\t\n", xi, yi);
		} fprintf(data, "\n\n");

    gsl_vector_set(xlist,0,a);
    gsl_matrix_set(ylist,0,0,1.0);
    gsl_matrix_set(ylist,0,1,0.0);
    k = ode_driver(sine_eq,rkstep23,xlist,ylist,b,h,acc,eps,max);

    for (int i = 0; i < k; i++) {
			double xi = gsl_vector_get(xlist,i);
			double yi = gsl_matrix_get(ylist,i,0);
			fprintf(data,"%g\t%g\t\n", xi, yi);
		} fprintf(data, "\n\n");
fclose(data);


	gsl_vector_free(xlist);
	gsl_matrix_free(ylist);

	return 0;
}
