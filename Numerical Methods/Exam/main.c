
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
	int n=2, steps_rk5, steps_rk23, max=1e6;
	double a=0, b=4*M_PI, h=0.001, acc=1e-6, eps=1e-6;
	gsl_vector* xlist = gsl_vector_alloc(max); // Corrected max+1
	gsl_matrix* ylist = gsl_matrix_alloc(max,n);

	// initialize
	gsl_vector_set(xlist,0,a);
	gsl_matrix_set(ylist,0,0,1.0);
	gsl_matrix_set(ylist,0,1,0.0);



	steps_rk5 = ode_driver(sine_eq,rkstep5,xlist,ylist,b,h,acc,eps,max);
	FILE* data = fopen("data.txt","w+");
		for (int i = 0; i < steps_rk5; i++) {
			double xi = gsl_vector_get(xlist,i);
			double yi = gsl_matrix_get(ylist,i,0);
			fprintf(data,"%g\t%g\t\n", xi, yi);
		} fprintf(data, "\n\n");

    gsl_vector_set(xlist,0,a);
    gsl_matrix_set(ylist,0,0,1.0);
    gsl_matrix_set(ylist,0,1,0.0);
    steps_rk23 = ode_driver(sine_eq,rkstep23,xlist,ylist,b,h,acc,eps,max);
    for (int i = 0; i < steps_rk23; i++) {
			double xi = gsl_vector_get(xlist,i);
			double yi = gsl_matrix_get(ylist,i,0);
			fprintf(data,"%g\t%g\t\n", xi, yi);
		} fprintf(data, "\n\n");
	fclose(data);

	printf("\n|------- Solving Harmonic Differential Equation -------|\n");
	steps_rk5 = ode_driver(sine_eq,rkstep5,xlist,ylist,b,h,acc,eps,max);
	printf("Number of steps used by the two methods to complete the integration:\nrk5=%i\nk23=%i\n", steps_rk5, steps_rk23);


	gsl_vector_free(xlist);
	gsl_matrix_free(ylist);

	return 0;
}
