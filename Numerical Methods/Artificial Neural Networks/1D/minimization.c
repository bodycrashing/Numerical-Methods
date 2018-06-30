
#include "minimization.h"
int newton_barebones(double f(gsl_vector* x, gsl_vector* grad, gsl_matrix* H), gsl_vector* x0, double eps){
	int n = x0->size;
	gsl_vector* grad = gsl_vector_alloc(n);
	gsl_vector* z = gsl_vector_alloc(n);
	//gsl_vector* fz = gsl_vector_alloc(n);
	gsl_vector* Dx = gsl_vector_alloc(n);
  gsl_matrix* R = gsl_matrix_alloc(n,n);
  gsl_matrix* H = gsl_matrix_alloc(n,n);
	int iter = 0;
	while (1) {
		iter++;

		double fx = f(x0,grad,H); // the function "f" automatically fills the Hessian H and the
                  // the gradient of the given function when supllied with a vector
                  // of size(n) and a matirx of size(n,n).

		//Dx*H = -Grad(f(x))
		qr_gs_decomp(H,R);
		gsl_vector_scale(grad,-1.0);
		qr_gs_solve(H,R,grad,Dx);
		gsl_vector_scale(grad,-1.0);

		// Backtracking
		double lambda = 1.0;
    double alpha = 10e-3;
		double res;

		while(1){
			iter++;
			gsl_vector_memcpy(z,x0);
			gsl_vector_add(z,Dx); // Trying out full step z = lambda*Dx
			double fz = f(z,grad,H);
			gsl_blas_ddot(z,grad,&res);
			double armijo = fx + alpha*res;
			if (fz <  armijo || lambda < 1.0/64.0) break;
			lambda/=2.0;
			gsl_vector_scale(Dx,lambda);

		}
	gsl_vector_memcpy(x0,z);
	//gsl_vector_memcpy(fx,fz);
	if(gsl_blas_dnrm2(grad)<eps) break;
	}

	// Free allocated memory
	gsl_matrix_free(H);
	gsl_matrix_free(R);
	//gsl_vector_free(fx);
	gsl_vector_free(z);
	gsl_vector_free(grad);
	gsl_vector_free(Dx);
	//gsl_vector_free(fz);
	return iter;
}




void numeric_grad(double f(gsl_vector* x), gsl_vector* x, gsl_vector* grad, double dx){
	int n = x->size;
	for(int i=0; i<n; i++){
		double fx = f(x);
		double xi = gsl_vector_get(x,i);
		gsl_vector_set(x,i,xi+dx);
		double fxi = f(x);
		double df =(fxi-fx)/dx;
		gsl_vector_set(grad,i,df);
		gsl_vector_set(x,i,xi);
	}
}





int quasinewton(double f(gsl_vector* x),
													gsl_vector* x0, double dx, double eps){
	int iter = 0 ;
	int n = x0->size;

	gsl_matrix* invH = gsl_matrix_alloc(n,n);
	gsl_vector* grad = gsl_vector_alloc(n);
	gsl_vector* Dx = gsl_vector_alloc(n);
	gsl_vector* z=gsl_vector_alloc(n);
	gsl_vector* gradz=gsl_vector_alloc(n);
	gsl_vector* y=gsl_vector_alloc(n);
	gsl_vector* u=gsl_vector_alloc(n);
	gsl_matrix_set_identity(invH);

	numeric_grad(f,x0,grad,dx);
	double fx = f(x0);
	while (iter<1000) {
	iter++;

	if(gsl_blas_dnrm2(grad)<eps){fprintf(stderr,"qnewton: |grad|<acc\n"); break;}

	gsl_blas_dgemv(CblasNoTrans,-1,invH,grad,0,Dx);
	double lambda = 1.0;
	double alpha = 10e-5;
	double sTgrad;
	double fz;
		while(1){
			iter++;
			gsl_vector_memcpy(z,x0);
			gsl_vector_add(z,Dx); // Trying out full step z = lambda*Dx
			fz = f(z);
			gsl_blas_ddot(z,grad,&sTgrad);
			if (fz <  fx + alpha*sTgrad || lambda < 1.0/64.0) break;
			if(gsl_blas_dnrm2(Dx) < dx){gsl_matrix_set_identity(invH); break;}
			lambda/=2.0;
			gsl_vector_scale(Dx,lambda);
		}

		numeric_grad(f,z,gradz,dx);
		gsl_vector_memcpy(y,gradz);
		gsl_blas_daxpy(-1,grad,y); /* y = grad(z)-grad(x) */
		gsl_vector_memcpy(u,Dx); /* u=s */
		gsl_blas_dgemv(CblasNoTrans,-1,invH,y,1,u); /* u=s-By */
		double sTy,uTy;
		gsl_blas_ddot(Dx,y,&sTy); // The denominator of eq.(14)
		if(fabs(sTy)>1e-12){
			gsl_blas_ddot(u,y,&uTy); //
			double gamma=uTy/2/sTy;
			gsl_blas_daxpy(-gamma,Dx,u); /* u=u-gamma*s */
			// Add an vector to explicitly show what you are doing...
			gsl_blas_dger(1.0/sTy,u,Dx,invH);
			gsl_blas_dger(1.0/sTy,Dx,u,invH);
		}
		gsl_vector_memcpy(x0,z);
		gsl_vector_memcpy(grad,gradz);
		fx=fz;
}

gsl_matrix_free(invH);
gsl_vector_free(grad);
gsl_vector_free(Dx);
gsl_vector_free(z);
gsl_vector_free(gradz);
gsl_vector_free(y);
gsl_vector_free(u);

	return iter;

}
