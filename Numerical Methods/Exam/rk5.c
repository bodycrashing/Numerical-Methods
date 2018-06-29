
#include "rk5.h"

void rkstep5(void f(double t, gsl_vector* y, gsl_vector* dydt), double t, gsl_vector* y, double h, gsl_vector* yh, gsl_vector* err){

  int n = y->size;
	int i;
	int p = 5; //<-- order of method
	gsl_vector* k1 = gsl_vector_alloc(n);
	gsl_vector* k2 = gsl_vector_alloc(n);
	gsl_vector* k3 = gsl_vector_alloc(n);
	gsl_vector* k4 = gsl_vector_alloc(n);
	gsl_vector* k5 = gsl_vector_alloc(n);
	gsl_vector* k6 = gsl_vector_alloc(n);
	gsl_vector* yt = gsl_vector_alloc(n);

	f(t,y,k1);
	for(i=0; i<n; i++){
		gsl_vector_set(yt,i,gsl_vector_get(y,i) + 1./4*gsl_vector_get(k1,i)*h);
	}

	f(t+1./4*h,yt,k2);
	for(i=0; i<n; i++){
		gsl_vector_set(yt,i,gsl_vector_get(y,i) + 3./32*gsl_vector_get(k1,i)*h + 9./32*gsl_vector_get(k2,i)*h);
	}

	f(t+3./8*h,yt,k3);
	for(i=0; i<n; i++){
		gsl_vector_set(yt,i,gsl_vector_get(y,i) + 1932./2197*gsl_vector_get(k1,i)*h - 7200./2197*gsl_vector_get(k2,i)*h + 7296./2197*gsl_vector_get(k3,i)*h);
	}

	f(t+12./13*h,yt,k4);
	for(i=0; i<n; i++){
		gsl_vector_set(yt,i,gsl_vector_get(y,i) + 439./216*gsl_vector_get(k1,i)*h - 8*gsl_vector_get(k2,i)*h + 3680./513*gsl_vector_get(k3,i)*h - 845./4104*gsl_vector_get(k4,i)*h);
	}

	f(t+h,yt,k5);
	for(i=0; i<n; i++){
		gsl_vector_set(yt,i,gsl_vector_get(y,i) - 8./27*gsl_vector_get(k1,i)*h + 2*gsl_vector_get(k2,i)*h - 3544./2565*gsl_vector_get(k3,i)*h + 1859./4104*gsl_vector_get(k4,i)*h - 11./40*gsl_vector_get(k5,i)*h);
	}

	f(t+h/2,yt,k6);
	for(i=0; i<n; i++){
		gsl_vector_set(yh,i,gsl_vector_get(y,i) + 16./135*gsl_vector_get(k1,i)*h + 6656./12825*gsl_vector_get(k3,i)*h + 28561./56430*gsl_vector_get(k4,i)*h - 9./50*gsl_vector_get(k5,i)*h + 2./55*gsl_vector_get(k6,i)*h);

	gsl_vector_set(err,i,(gsl_vector_get(k1,i)*h - gsl_vector_get(k1,i)*h/2 - gsl_vector_get(k6,i)*h/2)/(pow(2,p)-1));
	}

  gsl_vector_free(yt);
  gsl_vector_free(k1);
  gsl_vector_free(k2);
  gsl_vector_free(k3);
  gsl_vector_free(k4);
  gsl_vector_free(k5);
  gsl_vector_free(k6);

	return;
}



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
f(x+3./4*h,y_temp,k3); for(i=0;i<n;i++){gsl_vector_set(yxh,i, gsl_vector_get(yx,i) + (2./9*gsl_vector_get(k1,i) +
                                        1./3*gsl_vector_get(k2,i) + 4./9*gsl_vector_get(k3,i))*h);
                }
f(x+h,yxh,k4);    for(i=0;i<n;i++){
                  gsl_vector_set(y_temp,i, gsl_vector_get(yx,i) + (7./24*gsl_vector_get(k1,i) +
                  1./4*gsl_vector_get(k2,i) + 1./3*gsl_vector_get(k3,i) +
                  1./8*gsl_vector_get(k4,i))*h);
                  gsl_vector_set(err,i,gsl_vector_get(yxh,i) - gsl_vector_get(y_temp,i));
                }

  gsl_vector_free(y_temp);
  gsl_vector_free(k1);
  gsl_vector_free(k2);
  gsl_vector_free(k3);
  gsl_vector_free(k4);

  return;
}


int ode_driver(void f(double x, gsl_vector *y, gsl_vector *dydx),
                void stepper(void f(double x, gsl_vector* y, gsl_vector* dydt), double x, gsl_vector* y, double h, gsl_vector* yh, gsl_vector* err),
                gsl_vector *xlist, gsl_matrix *ylist, double b, double h, double acc, double eps, int max){

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
    //rkstep5(f,x,y,h,yxh,err_vec);
    stepper(f,x,y,h,yxh,err_vec);

    err = gsl_blas_dnrm2(err_vec);
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
