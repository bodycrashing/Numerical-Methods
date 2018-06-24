#include<stdio.h>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_multiroots.h>
#include<assert.h>

double F(double e, double r);

int aux(const gsl_vector* current_point, void* params, gsl_vector* f){
	double e = gsl_vector_get(current_point,0);
	assert(e<0);
	double rmax = *(double*)params;
	double fval = F(e,rmax);
	gsl_vector_set(f,0,fval);
	return GSL_SUCCESS;
}

int rosenbrock(const gsl_vector* point, void* params, gsl_vector* f){
  double x = gsl_vector_get(point,0);
  double y = gsl_vector_get(point,1);
  double dx = -2*(1-x)-400*x*(y-x*x);
  double dy = 200*(y-x*x);
  gsl_vector_set(f,0,dx);
  gsl_vector_set(f,1,dy);
  return GSL_SUCCESS;
}


int main(){
    /*--------------- Exercise 1: Gradient of Rosenbrock function --------------*/
  const gsl_multiroot_fsolver_type* T_rosenbrock = gsl_multiroot_fsolver_hybrids;
  gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc(T_rosenbrock,2);

  gsl_multiroot_function F_rosen;
  F_rosen.f = rosenbrock;
  F_rosen.n = 2;
  F_rosen.params = NULL;

  gsl_vector* p0 = gsl_vector_alloc(2);
  gsl_vector_set(p0,0,10);
  gsl_vector_set(p0,1,10);
  gsl_multiroot_fsolver_set(solver,&F_rosen,p0);

  int flag, iter=0;
  do{
    iter++;
    gsl_multiroot_fsolver_iterate(solver);
    flag = gsl_multiroot_test_residual(solver->f,1e-12);
  }
  while(flag==GSL_CONTINUE);

  double xmin = gsl_vector_get(solver->x,0);
  double ymin = gsl_vector_get(solver->x,1);
  printf("# The minimum of the Rosenbrock function is (x,y) = (%g,%g)\n",xmin,ymin);
  printf("# It took %i interations.\n\n",iter);


  /*--------------- Exercise 2: Bound states of hydrogen atom --------------*/
	double rmax = 8.0;
	const gsl_multiroot_fsolver_type* T = gsl_multiroot_fsolver_hybrids;
	gsl_multiroot_fsolver* S = gsl_multiroot_fsolver_alloc(T,1);

	gsl_multiroot_function wave_func;
	wave_func.f = aux;
	wave_func.n = 1;
	wave_func.params = (void*)&rmax;

	gsl_vector* p00 = gsl_vector_alloc(1);
	gsl_vector_set(p00,0,-2.0);
	gsl_multiroot_fsolver_set(S,&wave_func,p00);

	do{
		gsl_multiroot_fsolver_iterate(S);
		flag = gsl_multiroot_test_residual(S->f,1e-12);
	}
	while(flag==GSL_CONTINUE);

	double e_0 = gsl_vector_get(S->x,0);
	for(double r = 1e-3;r<rmax;r+=0.1){
		printf("%g %g %g\n",r,F(e_0,r),r*exp(-r));
	}


  gsl_vector_free(p0);
  gsl_multiroot_fsolver_free(solver);
  gsl_vector_free(p00);
  gsl_multiroot_fsolver_free(S);

	return 0;
}
