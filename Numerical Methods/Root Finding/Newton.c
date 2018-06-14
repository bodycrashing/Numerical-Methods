#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include"qr.h"

void newton(void f(gsl_vector* x, gsl_vector* fx), gsl_vector* x, double dx, double eps){
	int n=x->size;
	gsl_matrix* J = gsl_matrix_alloc(n,n);
	gsl_matrix* R = gsl_matrix_alloc(n,n);
	gsl_vector* fx = gsl_vector_alloc(n);
	gsl_vector* z  = gsl_vector_alloc(n);
	gsl_vector* fz = gsl_vector_alloc(n);
	gsl_vector* df = gsl_vector_alloc(n);
	gsl_vector* Dx = gsl_vector_alloc(n);

	while(1){
		f(x,fx);
		for (int j=0;j<n;j++){
			gsl_vector_set(x,j,gsl_vector_get(x,j)+dx);
			f(x,df);
			gsl_vector_sub(df,fx); /* df=f(x+dx)-f(x) */
			for(int i=0;i<n;i++) gsl_matrix_set(J,i,j,gsl_vector_get(df,i)/dx);
			gsl_vector_set(x,j,gsl_vector_get(x,j)-dx);
			}
		qr_gs_decomp(J,R);
		qr_gs_solve(J,R,fx,Dx);
		gsl_vector_scale(Dx,-1);

		double s=1;
		while(1){
			gsl_vector_memcpy(z,x);
			gsl_vector_add(z,Dx);
			f(z,fz);
			if( gsl_blas_dnrm2(fz)<(1-s/2)*gsl_blas_dnrm2(fx) || s<0.02 ) break;
			s*=0.5;
			gsl_vector_scale(Dx,0.5);
			}

		gsl_vector_memcpy(x,z);
		gsl_vector_memcpy(fx,fz);
		if( gsl_blas_dnrm2(Dx)<dx || gsl_blas_dnrm2(fx)<eps ) break;
		}

	gsl_matrix_free(J);
	gsl_matrix_free(R);
	gsl_vector_free(fx);
	gsl_vector_free(z);
	gsl_vector_free(fz);
	gsl_vector_free(df);
	gsl_vector_free(Dx);
}
/*
def newton_jacobian(f:"function", jacobian:"jacobian", x_start:vector, eps:float=1e-3):
	x=x_start.copy()
	while True :
		fx=f(x)
		J=jacobian(x)
		givens.qr(J)
		Dx = givens.solve(J,-fx)
		s=2
		while True :
			s/=2
			y=x+Dx*s
			fy=f(y)
			if fy.norm()<(1-s/2)*fx.norm() or s<0.02 : break
		x=y; fx=fy;
		if fx.norm()<eps : break
	return x;
*/
