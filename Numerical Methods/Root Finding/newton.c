#include"root.h"

int Newton(void f(gsl_vector* x, gsl_vector* fx), gsl_vector* x0, double dx, double eps){
	int n = x0->size;
	gsl_matrix* J = gsl_matrix_alloc(n,n);
	gsl_matrix* R = gsl_matrix_alloc(n,n);
	gsl_vector* fx = gsl_vector_alloc(n);
	gsl_vector* z = gsl_vector_alloc(n);
	gsl_vector* fz = gsl_vector_alloc(n);
	gsl_vector* df = gsl_vector_alloc(n);
	gsl_vector* Dx = gsl_vector_alloc(n);

	int iter = 0;
	while (1) {
		iter++;
		f(x0,fx);
		// Calculating the numeric Jacobian
		for (int j=0;j<n;j++){
			gsl_vector_set(x0,j,gsl_vector_get(x0,j)+dx); // Adding a small increment dx to x
			f(x0,df);
			gsl_vector_sub(df,fx); // Calculating "upper-half" of finite differnece df=f(x+dx)-f(x)
			gsl_vector_scale(df,1.0/dx);
			gsl_matrix_set_col(J,j,df);
			gsl_vector_set(x0,j,gsl_vector_get(x0,j)-dx); // Subtrackting dx again to rest x0 in order not to "forget" the starting point
		}																							// 	of the previous step

		// Solve J*Dx = -f(x)
		qr_gs_decomp(J,R);
		gsl_vector_scale(fx,-1.0);
		qr_gs_solve(J,R,fx,Dx);
		gsl_vector_scale(fx,-1.0);

		// Backtracking
		double lambda = 1.0;
		while(1){
			iter++;
			gsl_vector_memcpy(z,x0);
			gsl_vector_add(z,Dx); // Trying out full step
			f(z,fz);
<<<<<<< HEAD
			if (gsl_blas_dnrm2(fz) < (1-lambda/2.0)*gsl_blas_dnrm2(fx) || lambda < 1.0/64.0) break;
=======
			if (gsl_blas_dnrm2(fz) < (1-lambda/2.0)*gsl_blas_dnrm2(fx) || lambda > 1.0/64.0) break;
>>>>>>> 881c86d8c4bc3472f0396274d7f1cb1df2f5f285
			lambda/=2.0;
			gsl_vector_scale(Dx,0.5);
		}
	gsl_vector_memcpy(x0,z);
	gsl_vector_memcpy(fx,fz);
	if( gsl_blas_dnrm2(Dx)<dx || gsl_blas_dnrm2(fx)<eps ) break;
	}

	// Free allocated memory
	gsl_matrix_free(J);
	gsl_matrix_free(R);
	gsl_vector_free(fx);
	gsl_vector_free(z);
	gsl_vector_free(df);
	gsl_vector_free(Dx);
	gsl_vector_free(fz);
	return iter;
}


int Newton_Jacobian(void f(gsl_vector* x, gsl_vector* fx, gsl_matrix* J), gsl_vector* x0, double dx, double eps){
	int n = x0->size;
	gsl_matrix* J = gsl_matrix_alloc(n,n);
	gsl_matrix* R = gsl_matrix_alloc(n,n);
	gsl_vector* fx = gsl_vector_alloc(n);
	gsl_vector* z = gsl_vector_alloc(n);
	gsl_vector* fz = gsl_vector_alloc(n);
	gsl_vector* df = gsl_vector_alloc(n);
	gsl_vector* Dx = gsl_vector_alloc(n);

	int iter = 0;
	while (1) {
		iter++;
		// Calculating the numeric Jacobian
		f(x0,fx,J);

		// Solve J*Dx = -f(x)
		qr_gs_decomp(J,R);
		gsl_vector_scale(fx,-1.0);
		qr_gs_solve(J,R,fx,Dx);
		gsl_vector_scale(fx,-1.0);

		// Backtracking
		double lambda = 1.0;
		while(1){
			iter++;
			gsl_vector_memcpy(z,x0);
			gsl_vector_add(z,Dx); // Trying out full step
			f(z,fz,J);
<<<<<<< HEAD
			if (gsl_blas_dnrm2(fz) < (1-lambda/2.0)*gsl_blas_dnrm2(fx) || lambda < 1.0/64.0) break;
=======
			if (gsl_blas_dnrm2(fz) < (1-lambda/2.0)*gsl_blas_dnrm2(fx) || lambda > 1.0/64.0) break;
>>>>>>> 881c86d8c4bc3472f0396274d7f1cb1df2f5f285
			lambda/=2.0;
			gsl_vector_scale(Dx,0.5);
		}
	gsl_vector_memcpy(x0,z);
	gsl_vector_memcpy(fx,fz);
	if( gsl_blas_dnrm2(Dx)<dx || gsl_blas_dnrm2(fx)<eps ) break;
	}

	// Free allocated memory
	gsl_matrix_free(J);
	gsl_matrix_free(R);
	gsl_vector_free(fx);
	gsl_vector_free(z);
	gsl_vector_free(df);
	gsl_vector_free(Dx);
	gsl_vector_free(fz);
	return iter;
}






int Newton_Refined(void f(gsl_vector* x, gsl_vector* fx, gsl_matrix* J), gsl_vector* x0, double dx, double eps){
	int n = x0->size;
	gsl_matrix* J = gsl_matrix_alloc(n,n);
	gsl_matrix* R = gsl_matrix_alloc(n,n);
	gsl_vector* fx = gsl_vector_alloc(n);
	gsl_vector* z = gsl_vector_alloc(n);
	gsl_vector* fz = gsl_vector_alloc(n);
	gsl_vector* df = gsl_vector_alloc(n);
	gsl_vector* Dx = gsl_vector_alloc(n);

	int iter = 0;
	while (1) {
		iter++;
		// Calculating the numeric Jacobian
		f(x0,fx,J);

		// Solve J*Dx = -f(x)
		qr_gs_decomp(J,R);
		gsl_vector_scale(fx,-1.0);
		qr_gs_solve(J,R,fx,Dx);
		gsl_vector_scale(fx,-1.0);

		// Backtracking
		double lambda=1;
		double temp=0;
		gsl_blas_ddot(fx,fx,&temp);
		double g0 = 0.5*temp;
		double gprime0 = -temp;
		double glambda = 0;
		while(1){

			iter++;
			gsl_vector_memcpy(z,x0);
			gsl_vector_add(z,Dx); // Trying out full step
			f(z,fz,J);
<<<<<<< HEAD
			if (gsl_blas_dnrm2(fz) < (1-lambda/2.0)*gsl_blas_dnrm2(fx) || lambda < 1.0/64.0) break;
=======
			if (gsl_blas_dnrm2(fz) < (1-lambda/2.0)*gsl_blas_dnrm2(fx) || lambda > 1.0/64.0) break;
>>>>>>> 881c86d8c4bc3472f0396274d7f1cb1df2f5f285
			double c = (glambda-g0-gprime0*lambda)/(lambda*lambda);
			glambda = g0+gprime0*lambda+c*lambda*lambda;
			lambda = gprime0/(2*c);
			gsl_vector_scale(Dx,lambda);
		}

	gsl_vector_memcpy(x0,z);
	gsl_vector_memcpy(fx,fz);
	if( gsl_blas_dnrm2(Dx)<dx || gsl_blas_dnrm2(fx)<eps ) break;
	}

	// Free allocated memory
	gsl_matrix_free(J);
	gsl_matrix_free(R);
	gsl_vector_free(fx);
	gsl_vector_free(z);
	gsl_vector_free(df);
	gsl_vector_free(Dx);
	gsl_vector_free(fz);
	return iter;

}
