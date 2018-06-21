
#include"interp.h"

cspline *cspline_alloc(int n, double *x, double *y)
{				// builds natural cubic spline
	cspline *s = (cspline *) malloc(sizeof(cspline));
	s->x = (double *)malloc(n * sizeof(double));
	s->y = (double *)malloc(n * sizeof(double));
	s->b = (double *)malloc(n * sizeof(double));
	s->c = (double *)malloc((n-1) * sizeof(double));
	s->d = (double *)malloc((n-1) * sizeof(double));
	s->n = n;
	for (int i = 0; i<n; i++) {
		s->x[i] = x[i];
		s->y[i] = y[i];
	}
	double h[n-1], p[n-1];
	for (int i = 0; i<n-1; i++) {
		h[i] = x[i+1] - x[i];
		assert(h[i] > 0);
	}
	for (int i=0; i<n-1; i++)
		p[i] = (y[i+1] - y[i]) / h[i];
	double D[n], Q[n-1], B[n];
	D[0] = 2;
	for (int i = 0; i<n-2; i++)
		D[i + 1] = 2*h[i] / h[i+1] + 2;
	D[n - 1] = 2;
	Q[0] = 1;
	for (int i = 0; i < n-2; i++)
		Q[i+1] = h[i] / h[i+1];
	for (int i=0; i < n-2; i++)
		B[i+1] = 3 * (p[i] + p[i+1] * h[i] / h[i+1]);
	B[0] = 3 * p[0];
	B[n-1] = 3 * p[n-2];
	for (int i=1; i<n; i++) {
		D[i] -= Q[i-1] / D[i-1];
		B[i] -= B[i-1] / D[i-1];
	}
	s->b[n-1] = B[n-1]/D[n-1];	//back-substitution :
	for (int i = n-2; i>=0; i--)
		s->b[i] = (B[i] - Q[i] * s->b[i+1]) / D[i];
	for (int i = 0; i<n-1; i++) {
		s->c[i] = (-2*s->b[i] - s->b[i+1] + 3*p[i]) / h[i];
		s->d[i] = (s->b[i] + s->b[i+1] - 2*p[i]) / h[i] / h[i];
	}
	return s;
}


double cspline_eval(cspline *s, double z)
{
	assert(z >= s->x[0] && z <= s->x[s->n-1]);
	int i = 0, j = s->n-1;	// binary search for the interval for z :
	while (j-i > 1) {
		int m = (i+j)/2;
		if (z > s->x[m])
			i = m;
		else
			j = m;
	}
	double h = z-s->x[i];	// calculate the inerpolating spline :
	return s->y[i] + h*(s->b[i] + h*(s->c[i] + h*s->d[i]));
}


double cspline_integral(cspline *s, double z){
	assert(z>=s->x[0] && z<=s->x[s->n-1]);

	int i=0, j=s->n-1, k=0;
	while(j-i>1){
		int m=(j+i)/2;
		if(z>s->x[m])
			i=m;
		else
			j=m;
	}
	double integral=0.0, h;
	while(k<i){
		h = s->x[k+1]-s->x[k];
		integral += s->y[k]*h + 1./2*s->b[k]*h*h+1./3*s->c[k]*h*h*h + 1./4*s->d[k]*h*h*h*h;
		k++;
	}
	h = z-s->x[i];
	integral += s->y[i]*h + 1./2*s->b[i]*h*h+1./3*s->c[i]*h*h*h + 1./4*s->d[i]*h*h*h*h;
	return integral;
}

double cspline_deriv(cspline *s, double z){
	assert(z>=s->x[0] && z<=s->x[s->n-1]);

	int i=0, j=s->n-1;
	while(j-i>1){
		int m=(j+i)/2;
		if(z>s->x[m])
			i=m;
		else
			j=m;
	}
	double h=z-s->x[i];
	double deriv = s->b[i] + 2.0*s->c[i]*h + 3.0*s->d[i]*h*h;
	return deriv;
}








void cspline_free(cspline * s)
{				//free the allocated memory
	free(s->x);
	free(s->y);
	free(s->b);
	free(s->c);
	free(s->d);
	free(s);
}
