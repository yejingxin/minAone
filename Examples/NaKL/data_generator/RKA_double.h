/*
Adaptive RK4 - based on NR

Chris Knowlton 7/23/12
*/

#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif


#ifndef _RKAdaptive
#define _RKAdaptive

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void rk4(int n, double *x, double t, double tau, void func(double*,double,double*,double*), double *param, double *xout);

void rkA(int n, double *x, double t, double tau, void func(double*,double,double*,double*), double *param, double *taus, double *xout);
/*
n = number of variables - because c does not have a function that returns the number of
	bytes or elements allocated to a pointer
x = state vector
t = current time
tau = fixed time step
func = vector field - declared elsewhere, but needs form of func(x,t,p,dxdt)
param = fixed parameter values
		- you can include stimuli as a parameter whose value is updated each iteration
taus = currently used integration timestep - passing as variable allows value to be saved
		for the next iteration
xout = x(t+tau)
*/

void rkA(int n, double *x, double t, double tau, void func(double*,double,double*,double*), double *param, double *taus, double *xout){

	const double Safe1 = 0.9;
	const double Safe2 = 4.;
	const double eps = 1e-8;
	const int maxtry = 100;	
	const double err = 0.95;

	int i,j;
	double tsave = t;
	double tstop = t + tau;

	double tau_old = taus[0];
	double tau_new = taus[0];
	//double *xtemp = malloc(n*sizeof(double));
	//double *xsave = malloc(n*sizeof(double));

	double xtemp[n];
	double xsave[n];
	
	for(i=0;i<n;i++)
		xsave[i]=x[i];

	double half_tau;
	double scale;
	double errmax;

	double locerr;

	int off = 0;

	//printf("tau in = %e\n", tau_new);

   //printf("x in = %e\n", xsave[0]);		

	while(tsave < tstop){
	//int iter;
	//for(iter = 0;iter<5; iter++){

		//for(i=0;i<n;i++)
			//xsave[i]=x[i];
		
		if(off) break;

		for(i=0;i<maxtry;i++){

			//printf("time  is %e\n", t);

   		taus[0] = tau_new;
			tau_new = min(tau_new,tstop-t);
			//make sure we are not going past end of fixed iterations

			half_tau = 0.5*tau_new;
			rk4(n,xsave,tsave, half_tau, func, param, xtemp);
			t = tsave + half_tau;
			rk4(n,xtemp, t, half_tau, func, param, xout);
			//go forward 2 half steps
         //printf("x out = %e\n", xout[0]);
			
			t = tsave + tau_new;
			rk4(n,xsave, tsave, tau_new, func, param, xtemp);
			//go forward whole step

         //printf("x in prime = %e\n", xsave[0]);		
			errmax = 0;

			//get normalized error for each element 
			for(j=0; j<n; j++){
				scale = err*(fabs(xtemp[j])+fabs(xout[j]))/2.;
				locerr = fabs(xout[j]-xtemp[j]);
				locerr /= (scale+eps);
				errmax = max(errmax,locerr);
			}
			//printf("errmax = %e\n", errmax);

			if(errmax!=errmax){
				off = 1;
				break;
			}

			tau_old = tau_new;
			tau_new = Safe1*tau_old*pow(errmax,-0.2);
			//update tau based on largest error

			tau_new = max(tau_new,tau_old/Safe2);
			tau_new = min(tau_new,Safe2*tau_old);
			//limit size of change

			//printf("tau = %e\n",tau_new);

			if(errmax < 1){
				tsave = t;
				if(i>0) printf("used %d steps\n", i+1);
				for(j=0;j<n;j++){
					xsave[j]=xout[j];
				}
				//oops that line above is kind of important
				break;
			}

			else t = tsave;

			if(i==maxtry-1) printf("max attempts used\n");
		}
	}
	
	//printf("x out final = %e\n" , xout[0]);
	//copy values to output
	// while it may appear that t -> t+tau, because t is passed as a double and not double*
	// it does not change the value passed as an argument - only a copy of that variable
	//free(xsave);
	//free(xtemp);
}

void rk4(int n, double *x, double t, double tau, void func(double*,double,double*,double*), double *param, double *xout){

	double k1[n];
	double k2[n];
	double k3[n];
	double k4[n];


	double xtemp1[n];
	
	int i;

	func(x,t,param,k1);

	for(i=0;i<n;i++)
		xtemp1[i] = x[i]+0.5*k1[i]*tau;

	func(xtemp1,t+.5*tau,param,k2);

	for(i=0;i<n;i++)
		xtemp1[i] = x[i]+0.5*k2[i]*tau;

	func(xtemp1,t+0.5*tau,param,k3);

	for(i=0;i<n;i++)
		xtemp1[i] = x[i]+k3[i]*tau;

	func(xtemp1,t+tau,param,k4);

	for(i=0;i<n;i++)
		xout[i] = x[i]+tau/6.0*(k1[i]+k4[i]+2*(k2[i]+k3[i]));


}


#endif
