#include "adaptive.h"

double Adaptive24(double f(double x),double a, double b, double acc,
															double eps, double f2, double f3, int nrec, double* err){
															//Notice that in comparision to Dimitris implementation a dereffereneced "double* err"
															// is also supplied to the function storing the final absolute error of the integration
assert(nrec<1e+6);

	double f1 = f(a+(b-a)/6);
	double f4 = f(a+5*(b-a)/6);
	double Q = (2*f1+f2+f3+2*f4)/6*(b-a);
	double q = (f1+f2+f3+f4)/4*(b-a);
	double tol = acc+eps*fabs(Q);
	*err = fabs(Q-q);

	if(*err<tol)return Q;

	else{
		double Q1 = Adaptive24(f,a,(a+b)/2,acc/sqrt(2.),eps,f1,f2,nrec+1,err);
		double Q2 = Adaptive24(f,(a+b)/2,b,acc/sqrt(2.),eps,f3,f4,nrec+1,err);
		return Q1+Q2;
	}
}

double Adaptive(double f(double x),double a, double b, double acc, double eps, double *err){
	int nrec = 0;
	double f2 = f(a+2*(b-a)/6);
	double f3 = f(a+4*(b-a)/6);
	return Adaptive24(f,a,b,acc,eps,f2,f3,nrec,err);
}


double Adaptive_Infinity(double f(double x),double a, double b, double acc, double eps,double *err){

int a_test = isinf(-a);
int b_test = isinf(b);

	if(a_test && b_test){
		double var_trans(double t){
			return f(t/(1-t*t))*(1+t*t)/((1-t*t)*(1-t*t)); //eq. (8.57)
		};
		int nrec = 0;
		double a =-1, b=1;
		double f2 = f(a+2*(b-a)/6);
		double f3 = f(a+4*(b-a)/6);
		return Adaptive24(var_trans,a,b,acc,eps,f2,f3,nrec,err);
	}

	else if(a_test){
		double var_trans(double t){
			return f(b-(1-t)/t)*1/(t*t); //eq. (8.62)
		};
		double a = 0, b = 1;
		double f2 = f(a+2*(b-a)/6);
		double f3 = f(a+4*(b-a)/6);
		int nrec = 0;
		return Adaptive24(var_trans,a,b,acc,eps,f2,f3,nrec,err);
	}

	else if(b_test){
		double var_trans(double t){
			return f(a+(1-t)/t)*1/(t*t); //eq. (8.60)
		};
		double a = 0, b = 1;
		double f2 = f(a+2*(b-a)/6);
		double f3 = f(a+4*(b-a)/6);
		int nrec = 0;
		return Adaptive24(var_trans,a,b,acc,eps,f2,f3,nrec,err);
	}

	else {
		int nrec = 0;
		double f2 = f(a+2*(b-a)/6); //no rescalling needed 
		double f3 = f(a+4*(b-a)/6);
		return Adaptive24(f,a,b,acc,eps,f2,f3,nrec,err);
	}
}


double Clenshaw_Curtis(double f(double x),double a, double b, double acc, double eps, double *err){
double a_new=0;
double b_new=M_PI;

double f_clenshaw_curtis(double theta){
	double x = 0.5*(a+b)+0.5*(a-b)*cos(theta);
	return f(x)*sin(theta)*(b-a)/2;
	}
return Adaptive_Infinity(f_clenshaw_curtis,a_new,b_new,acc,eps,err);
}
