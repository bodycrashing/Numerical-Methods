#include"root.h"
#include<gsl/gsl_multiroots.h>
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

void f(gsl_vector* x, gsl_vector* fx){
	double A=10000;
	double X = gsl_vector_get(x,0);
	double Y = gsl_vector_get(x,1);
	gsl_vector_set(fx,0,A*Y*X-1.0);
	gsl_vector_set(fx,1,exp(-X)+exp(-Y)-1.0-1.0/A);
}

void f_Jacobian(gsl_vector* x, gsl_vector* fx, gsl_matrix* J){
	double A = 10000;
	double X = gsl_vector_get(x,0);
	double Y = gsl_vector_get(x,1);
	gsl_vector_set(fx,1,A*Y*X-1.0);
	gsl_vector_set(fx,1,exp(-X)+exp(-Y)-1.0-1.0/A);

	gsl_matrix_set(J,0,0,A*Y);
	gsl_matrix_set(J,0,1,A*X);
	gsl_matrix_set(J,1,0,-1.0*exp(-X));
	gsl_matrix_set(J,1,1,-1.0*exp(-Y));
}


void rosenbrock(gsl_vector* x, gsl_vector* fx){
	double X = gsl_vector_get(x,0);
	double Y = gsl_vector_get(x,1);
	gsl_vector_set(fx,0,2*(X-1)-400*X*(Y-X*X));
	gsl_vector_set(fx,1,200*(Y-X*X));
}

void rosenbrock_Jacobian(gsl_vector* x, gsl_vector* fx, gsl_matrix* J){
  double X = gsl_vector_get(x,0);
	double Y = gsl_vector_get(x,1);

        gsl_vector_set(fx,0,2.0*(X-1)-400.0*X*(Y-X*X));
        gsl_vector_set(fx,1,200.0*(Y-X*X));

	gsl_matrix_set(J,0,0, 2.0*(600.0*X*X - 200.0*Y + 1.0));
	gsl_matrix_set(J,0,1, -400.0*X);
	gsl_matrix_set(J,1,0, -400.0*X);
	gsl_matrix_set(J,1,1, 200.0);
}


void himmelblau(gsl_vector* x, gsl_vector* fx){
	double X = gsl_vector_get(x,0);
	double Y = gsl_vector_get(x,1);
	gsl_vector_set(fx, 0, 2.0*(2.0*X*(X*X + Y -11.0) + X + Y*Y - 7.0));
	gsl_vector_set(fx, 1, 2.0*(X*X + 2.0*Y*(X + Y*Y - 7.0) + Y - 11.0));
}

void himmelblau_Jacobian(gsl_vector* x, gsl_vector* fx, gsl_matrix* J){
        double X = gsl_vector_get(x,0);
        double Y = gsl_vector_get(x,1);
	gsl_vector_set(fx,0, 2.0*(2.0*X*(X*X + Y-11.0) + X + Y*Y-7.0));
	gsl_vector_set(fx,1, 2.0*(X*X + 2.0*Y*(X + Y*Y-7.0) + Y- 11.0));

	//gsl_matrix_set(J,0,0, 2*(4.0*X*X + 2.0*(X*X +Y -11.0) +1.0));
	gsl_matrix_set(J,0,0, 4*(3*X*X+Y-11) + 2);
	gsl_matrix_set(J,0,1, 2.0*(2.0*X+2.0*Y));
	gsl_matrix_set(J,1,0, 2.0*(2.0*X+2.0*Y));
	//gsl_matrix_set(J,1,1, 2.0*(4.0*Y*Y + 2.0*(Y*Y + X-7.0) + 1.0));
	gsl_matrix_set(J,1,1, 2 + 4*(X + 3*Y*Y-7));
}




int system_f_gsl(const gsl_vector* v, void *params, gsl_vector* fx){
double A = 10000.0;
double x = gsl_vector_get(v,0);
double y =	gsl_vector_get(v,1);
double fx0 = A*x*y -1;
double fx1 = exp(-x) + exp(-y) -(1.0 + 1.0/A);
gsl_vector_set(fx,0,fx0);
gsl_vector_set(fx,1,fx1);
return GSL_SUCCESS;
}

int Rosenbrock_grad_f_gsl(const gsl_vector* v, void *params, gsl_vector* fx){
double x = gsl_vector_get(v,0);
double y =	gsl_vector_get(v,1);
double fx0 = 2*x-2+400*(x*x*x-y*x);
double fx1 = 200*(y-x*x);
gsl_vector_set(fx,0,fx0);
gsl_vector_set(fx,1,fx1);
//function_calls++;
return GSL_SUCCESS;
}

int Himmelblau_grad_f_gsl(const gsl_vector* v, void *params, gsl_vector* fx){
double x = gsl_vector_get(v,0);
double y =	gsl_vector_get(v,1);
double fx0 = 4*x*(x*x+y-11)+2*(x+y*y-7);
double fx1 = 2*(x*x+y-11)+4*y*(x+y*y-7);
gsl_vector_set(fx,0,fx0);
gsl_vector_set(fx,1,fx1);
//function_calls++;
return GSL_SUCCESS;
}



int main(void){

	gsl_vector* x = gsl_vector_alloc(2);
	gsl_vector* fx = gsl_vector_alloc(2);
	double step = 1e-6, eps = 1e-6;

	printf("------ Part A: Back-tracking linesearch and numerical Jacobian ------\n\n");
	gsl_vector_set(x,0,1.0/10000);
	gsl_vector_set(x,1,1);
	printf("------------ Roots of System of Equations ------------\n");
	printf("dx=%lg  eps=%lg\n", step, eps);
	printf("Starting points: x=%lg  y=%lg \n", gsl_vector_get(x,0), gsl_vector_get(x,1));
	f(x,fx);
	int iter_f =  Newton(f,x,step,eps);
	printf("The computed roots are:\n");
	vector_print(x);
	printf("Determined in  %i iterations\n\n",iter_f);

	// Rosenbrock's valley function
	printf("------------ Minimum of Rosenbrock's Valley Function ------------\n");
	gsl_vector_set(x,0,2.0);
	gsl_vector_set(x,1,2.0);
	printf("Starting points: x=%lg  y=%lg \n", gsl_vector_get(x,0), gsl_vector_get(x,1));
	rosenbrock(x,fx);
	int iter_ros = Newton(rosenbrock,x,step,eps);
	printf("The computed roots are:\n");
	vector_print(x);
	printf("Determined in  %i iterations\n\n",iter_ros);

	// Himmelblau's function
	printf("------------ Minimum of Himmelblaus's Function ------------\n");
	gsl_vector_set(x,0,4.0);
	gsl_vector_set(x,1,4.0);
	printf("Starting points: x=%lg  y=%lg \n", gsl_vector_get(x,0), gsl_vector_get(x,1));
	himmelblau(x,fx);
	int iter_him = Newton(himmelblau,x,step,eps);
	printf("The computed roots are:\n");
	vector_print(x);
	printf("Determined in  %i iterations\n\n",iter_him);







	/*------ Part B: Newton's method with analytic Jacobian ------*/
	printf("\n------ Part B: Newton's method with analytic Jacobian ------\n\n");
	// System of equations using analytical jacobian
	printf("------------ System of Equations using Analytic Jacobian ------------\n");
	gsl_vector_set(x,0,1.0/10000);
	gsl_vector_set(x,1,1.0);
	printf("Starting points: x=%lg  y=%lg \n", gsl_vector_get(x,0), gsl_vector_get(x,1));
	int iter_f_analytical = Newton_Jacobian(f_Jacobian,x,step,eps);
	printf("The computed roots are:\n");
	vector_print(x);
	printf("Number of iterations analytical Jacobian: %i\nNumber of iterations numerical Jacobian: %i\n\n",iter_f_analytical,iter_f);

	// Rosenbrock with analytical jacobian
	gsl_vector_set(x,0,2.0);
	gsl_vector_set(x,1,2.0);
	printf("Starting points: x=%lg  y=%lg \n", gsl_vector_get(x,0), gsl_vector_get(x,1));
	int iter_ros_analytical = Newton_Jacobian(rosenbrock_Jacobian,x,step,eps);
	printf("The computed roots are:\n");
	vector_print(x);
	printf("Number of iterations analytical Jacobian: %i\nNumber of iterations numerical Jacobian: %i\n\n",iter_ros_analytical,iter_ros);

	// Himmelblau with analytical jacobian
	gsl_vector_set(x,0,4.0);
	gsl_vector_set(x,1,4.0);
	printf("Starting points: x=%lg  y=%lg \n", gsl_vector_get(x,0), gsl_vector_get(x,1));
	int iter_him_analytical = Newton_Jacobian(himmelblau_Jacobian,x,step,eps);
	printf("The computed roots are:\n");
	vector_print(x);
	printf("Number of iterations analytical Jacobian: %i\nNumber of iterations numerical Jacobian: %i\n\n",iter_him_analytical,iter_him);


	/*------- Comparing results to GSL Routine -------------*/
	int n=2, status;
	int iter=0;
	double epsilon_system=1e-3, epsilon=1e-6;
	const gsl_multiroot_fsolver_type* T =gsl_multiroot_fsolver_dnewton;
	gsl_multiroot_fsolver* s = gsl_multiroot_fsolver_alloc(T, n);

	gsl_multiroot_function gslsystem = {&system_f_gsl,n,NULL};
	gsl_multiroot_function gsl_Rosenbrock = {&Rosenbrock_grad_f_gsl,n,NULL};
	gsl_multiroot_function gsl_Himmelblau = {&Himmelblau_grad_f_gsl,n,NULL};

	gsl_multiroot_fsolver_set (s,&gslsystem,x);
	do{
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);
      if (status) break;
      status = gsl_multiroot_test_residual (s->f, epsilon_system);
    } while (status == GSL_CONTINUE && iter < 1000);
		printf("System of Equations using GSL routine:\n");
		vector_print(s->x);
		printf("in %i iterations\n",iter);

		iter=0;
		gsl_multiroot_fsolver_set (s,&gsl_Rosenbrock,x);
		do{
	      iter++;
	      status = gsl_multiroot_fsolver_iterate (s);
	      if (status) break;
	      status = gsl_multiroot_test_residual (s->f, epsilon);
	    }while (status == GSL_CONTINUE && iter < 1000);
			printf("Rosenbrock's Valley Function using GSL routine:\n");
			vector_print(s->x);
			printf("in %i iterations\n",iter);

		iter=0;
		gsl_multiroot_fsolver_set (s,&gsl_Himmelblau,x);
		do{
	      iter++;
	      status = gsl_multiroot_fsolver_iterate (s);
	      if (status) break;
	      status = gsl_multiroot_test_residual (s->f, epsilon);
	    }while (status == GSL_CONTINUE && iter < 1000);
			printf("Himmelblauss UFunction using GSL routine:\n");
			vector_print(s->x);
			printf("in %i iterations\n",iter);



/*--------- Part C: Newton's method with refined linesearch ----------*/
printf("\n\n");
printf("------ Part C: Newton's method with refined linesearch ------\n\n");

gsl_vector_set(x,0,1.0/10000);
gsl_vector_set(x,1,1.0);
printf("Starting points: x=%lg  y=%lg \n", gsl_vector_get(x,0), gsl_vector_get(x,1));
int iter_refined_f = Newton_Refined(f_Jacobian,x,step,eps);
printf("The computed roots are:\n");
vector_print(x);
printf("Number of iterations analytical Jacobian: %i\nNumber of iterations numerical Jacobian: %i\n\n",iter_f_analytical,iter_refined_f);

// Rosenbrock with analytical jacobian
gsl_vector_set(x,0,2.0);
gsl_vector_set(x,1,2.0);
printf("Minimum of Rosenbrock function with analytical jacobial given the starting points\n");
vector_print(x);
int iter_refined_ros = Newton_Jacobian(rosenbrock_Jacobian,x,step,eps);
printf("The computed roots are:\n");
vector_print(x);
printf("Number of iterations analytical Jacobian: %i\nNumber of iterations numerical Jacobian: %i\n\n",iter_ros_analytical,iter_refined_ros);

// Himmelblau with analytical jacobian
gsl_vector_set(x,0,4.0);
gsl_vector_set(x,1,4.0);
printf("Minimum of Himmalblau function with analytical jacobial given the starting points\n");
vector_print(x);
int iter_refined_him = Newton_Jacobian(himmelblau_Jacobian,x,step,eps);
printf("The computed roots are:\n");
vector_print(x);
printf("Number of iterations analytical Jacobian: %i\nNumber of iterations numerical Jacobian: %i\n\n",iter_him_analytical,iter_refined_him);



	gsl_multiroot_fsolver_free(s);
	gsl_vector_free(x);
	gsl_vector_free(fx);

	return 0;
}
