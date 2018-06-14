int driverPathStoring(
	gsl_vector *tPath,
	double b,
	double *h,
	gsl_matrix *yPath,
	double abs,
	double eps,
	void stepper(
		double t, double h, gsl_vector *y,
		void f(double t, gsl_vector *y, gsl_vector *dydt),
		gsl_vector *yh, gsl_vector *err
		),
	void f(double t, gsl_vector *y, gsl_vector *dydt)
)
{
	int n = yPath->size1;
	int maxSteps = (*tPath).size, step = 0, DRIVER_FAIL=0;
	double a = gsl_vector_get(tPath,0), tol, normErr, t;
	gsl_vector *yh = gsl_vector_alloc(n);
	gsl_vector *err = gsl_vector_alloc(n);
	gsl_vector_view y, yNext;

	while( gsl_vector_get(tPath,step) < b )
	{
		t = gsl_vector_get(tPath,step);
		y = gsl_matrix_column(yPath,step);

		if(t+*h>b)
		{
			*h = b-t;
		}
		stepper(t,*h,&y.vector,f,yh,err);

		tol = ( eps*gsl_blas_dnrm2(yh) + abs ) * sqrt(*h/(b-a));
		normErr = gsl_blas_dnrm2(err);
		//printf("\n%g\n", normErr);

		if( normErr < tol )
		{
			step++;
			if( step+1 > maxSteps )
			{
				return DRIVER_FAIL;
			}
			gsl_vector_set(tPath,step,t+*h);
			yNext = gsl_matrix_column(yPath,step);
			gsl_vector_memcpy(&yNext.vector,yh);
		}

		*h *= pow(tol/normErr,0.25)*0.95;
	}

	gsl_vector_free(yh);
	gsl_vector_free(err);

	return step+1;
}
