#include "runge_kutta.h"
/*---------------- Stepper function --------------------------------*/
void rkstep23(void f(double x, gsl_vector *y, gsl_vector *dydx),
              double x, gsl_vector* yx, double h, gsl_vector* yxh, gsl_vector* err){
// x is the independet variable of the ODE
// yx is the desired function y(x)
// yxh is the function at the vaule x+h i.e. yxh = y(x+h)

int i;
int n = yx->size;
gsl_vector* y_temp = gsl_vector_alloc(n);

gsl_vector* k1 = gsl_vector_alloc(n);
gsl_vector* k2 = gsl_vector_alloc(n);
gsl_vector* k3 = gsl_vector_alloc(n);
gsl_vector* k4 = gsl_vector_alloc(n);

f(x,yx,k1);            for(i=0;i<n;i++) gsl_vector_set(y_temp,i,gsl_vector_get(yx,i) + 1./2*gsl_vector_get(k1,i)*h);
f(x+1./2*h,y_temp,k2); for(i=0;i<n;i++) gsl_vector_set(y_temp,i,gsl_vector_get(yx,i) + 3./4*gsl_vector_get(k2,i)*h);
f(x+3./4*h,y_temp,k3); for(i=0;i<n;i++){
                        gsl_vector_set(yxh,i,
				gsl_vector_get(yx,i) + (2./9*gsl_vector_get(k1,i) +
                        1./3*gsl_vector_get(k2,i) + 4./9*gsl_vector_get(k3,i))*h);
                      }
f(x+h,yxh,k4);    for(i=0;i<n;i++){
                  gsl_vector_set(y_temp,i, gsl_vector_get(yx,i) + (7./24*gsl_vector_get(k1,i) +
                  1./4*gsl_vector_get(k2,i) + 1./3*gsl_vector_get(k3,i) +
                  1./8*gsl_vector_get(k4,i))*h);

                  gsl_vector_set(err,i,gsl_vector_get(yxh,i) - gsl_vector_get(y_temp,i));
//                  fprintf(stderr, "\n%lg\n",gsl_vector_get(yxh,i) - gsl_vector_get(y_temp,i));
                }

  gsl_vector_free(y_temp);
  gsl_vector_free(k1);
  gsl_vector_free(k2);
  gsl_vector_free(k3);
  gsl_vector_free(k4);
}


int ode_driver(void f(double x, gsl_vector *y, gsl_vector *dydx),
                gsl_vector *xlist, gsl_matrix *ylist,
                double b, double h, double acc, double eps, int max){

  int i, k=0;
  int n = ylist->size2;
  double x, err, norm_y, tol, a=gsl_vector_get(xlist,0);

  gsl_vector* y = gsl_vector_alloc(n);
  gsl_vector* yxh = gsl_vector_alloc(n);
  gsl_vector* err_vec = gsl_vector_alloc(n);

  while(gsl_vector_get(xlist,k)<b){
    x = gsl_vector_get(xlist,k);
    gsl_matrix_get_row(y,ylist,k);

    if(x+h>b) h=b-x;
    rkstep23(f,x,y,h,yxh,err_vec); // En error vektor????

    err = gsl_blas_dnrm2(err_vec);
    //fprintf(stderr, "x=%g, h=%g, err=%g\n",x,h,err);
    norm_y = gsl_blas_dnrm2(yxh);
    tol=(norm_y*eps+acc)*sqrt(h/(b-a));

    if(err<tol){
      k++;
      if(k>max-1){
        fprintf(stderr, "\n\nMaximum number of steps reached!\n\n");
        break;
      }
      gsl_vector_set(xlist,k,x+h);
      for(i=0; i<n; i++) gsl_matrix_set(ylist,k,i,gsl_vector_get(yxh,i));
    }

    if(err>0){
      h*=pow(tol/err,0.25)*0.95;
    }
    else{
      h*=2;
    }
  }
  return k+1;

  gsl_vector_free(y);
  gsl_vector_free(yxh);
  gsl_vector_free(err_vec);
}



double integrator(void f(double x, gsl_vector* y, gsl_vector* dydx),gsl_vector *xlist, gsl_matrix *ylist,
double b, double h, double acc, double eps, int max){

	int k = ode_driver(f,xlist,ylist,b,h,acc,eps,max);
	printf("Number of steps %i\n",k);
	return gsl_matrix_get(ylist,k-1,0);
}
