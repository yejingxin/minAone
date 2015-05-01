// simple_nakl.cpp
// Nonlinear Ipopt program

// Author: Bryan A. Toth
// btoth@physics.ucsd.edu

#include "simple_naklminAone_nlp.hpp"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <cstring>
#ifndef HAVE_CSTDIO
#define HAVE_CSTDIO
# include <cstdio>
# include <iostream>
# include "myfunctions.hpp"
# include <fstream>
# include <string>
# include <stdlib.h>
#else
# ifndef HAVE_STDIO_H
#define HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif

using namespace Ipopt;

using namespace std;

double ghk(double);
double ghkjac(double, int);
double ghkhes(double, int, int);
// constructor
SIMPLE_NAKL_NLP::SIMPLE_NAKL_NLP(int id)
{
  nU=0;
  nP=0;
  nY=7;
  nM=1;
  nI=1;

     K11val = new double[nU];
     K11val2 = new double[nU];
     K11valp1 = new double[nU];
     dK11val = new double[nU];
     dK11val2 = new double[nU];
     dK11valp1 = new double[nU];
     Xdval = new double[nM];
     Xdval2 = new double[nM];
     Xdvalp1 = new double[nM];
     Xval = new double[nY];
     Xval2 = new double[nY];
     Xvalp1 = new double[nY];
     Pval = new double[nP];
     Ival = new double[nI];
     Ival2 = new double[nI];
     Ivalp1 = new double[nI];
     Rf0 = new double[nY];
     taskid = id;

     string buffer;
     specs = new string[6+nP+nY+nI+2*nU+nM+1];
     
     int count;
     count = 0;
     
     ifstream fin ("specs.txt");
     if (fin.is_open())
     {
       while (! fin.eof())
       {
         getline (fin,buffer);
	 if (buffer[0] !='#')
	   {
	   specs[count] = buffer;
	   count++;
	   }
       }
       fin.close();
     }
     
     else cout << "Unable to open file";
    Time = atoi(specs[0].c_str());
    skip = atoi(specs[1].c_str());
    hstep = atof(specs[2].c_str());

    string filename;
    int ret;
    VDATA0 = new double[2*Time+1];
    VDATA0dummy = new double[skip];

    Ntotal = (2*Time+1)*nY+nP;
    solution = new double[Ntotal];
    FILE *pFile0;
    filename = specs[3];
    pFile0 = fopen(filename.c_str(),"r");

    for(Index jt=0;jt<skip;jt++)
	{
	ret = fscanf (pFile0, "%lf", &VDATA0dummy[jt]);
	if (ret == EOF) break;
	}
    for(Index jt=0;jt<2*Time+1;jt++)
	{
	ret = fscanf (pFile0, "%lf", &VDATA0[jt]);
	if (ret == EOF) break;
	}
    fclose (pFile0);
    Iinj = new double[2*Time+1];
    Iinjdummy = new double[skip];

    FILE *qFile0;
    filename = specs[4];
    qFile0 = fopen(filename.c_str(),"r");

    for(Index jt=0;jt<skip;jt++)
        {
	ret = fscanf (qFile0, "%lf", &Iinjdummy[jt]);
	if (ret == EOF) break;
        }
    for(Index jt=0;jt<2*Time+1;jt++)
        {
	ret = fscanf (qFile0, "%lf", &Iinj[jt]);
	if (ret == EOF) break;
	}
    fclose (qFile0);
    int rows = nY+2*nU+nP+1;
    bounds = new double*[rows];
    for (Index i=0;i<rows;i++) bounds[i] = new double[4];
    int toggle=0;
    if (specs[3+nM+nI] == "1") toggle = 1;
    int counter;
    for(Index k=0;k<rows;k++)
       {
       counter=0;
       char* tmp = new char[specs[4+nM+nI+toggle+k].size()+1];
       strcpy( tmp, specs[4+nM+nI+toggle+k].c_str() );
       char *ptr = strtok(tmp,",");
       bounds[k][3] = 0.0;
       while(ptr != 0) {
          if(counter<3) {
	     bounds[k][counter] = atof(ptr);
	     }
          if(counter==3) {
             bounds[k][counter] = atof(ptr);
             }
	  ptr = strtok(0,",");
	  counter++;
          }
    }

        for (Index i=0;i<nY;i++){
		Rf0[i]=bounds[i][2];
		bounds[i][3]=bounds[i][2];
	}
    beta=0;
    alpha = bounds[nY+2*nU+nP][0];
    delta_beta=(int) bounds[nY+2*nU+nP][1];
    max_beta=(int) bounds[nY+2*nU+nP][2];
    if (specs[3+nM+nI] == "1")
       {
       filename = specs[4+nM+nI];
       }

}

// destructor
SIMPLE_NAKL_NLP::~SIMPLE_NAKL_NLP()
{
  delete [] K11val;
  delete [] K11val2;
  delete [] K11valp1;
  delete [] dK11val;
  delete [] dK11val2;
  delete [] dK11valp1;
  delete [] Xdval;
  delete [] Xdval2;
  delete [] Xdvalp1;
  delete [] Xval;
  delete [] Xval2;
  delete [] Xvalp1;
  delete [] Pval;
  delete [] Ival;
  delete [] Ival2;
  delete [] Ivalp1;
  delete [] specs;
    delete [] VDATA0;
    delete [] VDATA0dummy;
    delete [] Iinj;
    delete [] Iinjdummy;
  int rows = nY+2*nU+nP;
  for (Index i=0;i<rows;i++) delete [] bounds[i];
  delete [] bounds;

}


bool SIMPLE_NAKL_NLP::changeRf(){
	if((beta+delta_beta)>(max_beta-1))
		return false;
	else
		beta = beta + delta_beta;
	for (Index i=0;i<nY;i++) {
		bounds[i][3]=pow(alpha,beta)*Rf0[i];
	}
	printf("beta=%d\n",beta);
	return true;
}

// returns the size of the problem
bool SIMPLE_NAKL_NLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
				Index& nnz_h_lag, IndexStyleEnum& index_style)

{
  // Number of variables
  n = 14*Time+7;

  // Number of equality constraints
  m = 0;

  // Number of Jacobian nonzero entries
  nnz_jac_g = 0;

  // Number of Hessian nonzero entries
  nnz_h_lag = 22*(Time+1)+133*Time+0;

  // use the C style indexing (0-based)
  index_style = TNLP::C_STYLE;

  return true;
}


// returns the variable bounds
bool SIMPLE_NAKL_NLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{
  // Here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.
  assert(n == 14*Time+7);
  assert(m == 0);

  for(Index jt=0;jt<Time+1;jt++) {
     for(Index var=0;var<nY;var++) {
        // Bounds for x
        x_l[(Time+1)*var+jt]=bounds[var][0];
        x_u[(Time+1)*var+jt]=bounds[var][1];
        // Bounds for midpoints
        if(jt<Time) {
       x_l[(Time+1)*(nY+2*nU)+Time*var+jt]=bounds[var][0];
       x_u[(Time+1)*(nY+2*nU)+Time*var+jt]=bounds[var][1];
       }
    }
     for(Index con=0;con<2*nU;con++) {
       // Bounds for k
       x_l[(Time+1)*(nY+con)+jt]=bounds[nY+con][0];
       x_u[(Time+1)*(nY+con)+jt]=bounds[nY+con][1];
       // Bounds for midpoints
       if(jt<Time) {
          x_l[(Time+1)*(nY+2*nU)+Time*(nY+con)+jt]=bounds[nY+con][0];
	  x_u[(Time+1)*(nY+2*nU)+Time*(nY+con)+jt]=bounds[nY+con][1];
	  }
     }

  } // End for loop

     for(Index par=0;par<nP;par++) {
        // Bounds for parameters
        x_l[2*Time*(nY+2*nU)+nY+2*nU+par]=bounds[nY+2*nU+par][0];
        x_u[2*Time*(nY+2*nU)+nY+2*nU+par]=bounds[nY+2*nU+par][1];
              }

  return true;
}


// returns the initial point for the problem
bool SIMPLE_NAKL_NLP::get_starting_point(Index n, bool init_x, Number* x,
                                   bool init_z, Number* z_L, Number* z_U,
                                   Index m, bool init_lambda,
                                   Number* lambda)
{
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  for (Index i=0; i<n; i++) {
        x[i] = 0.0;
      }

    char filename[20];
    FILE *initFILE;
	sprintf(filename,"D%d_M%d_PATH%d.dat", nY,nM,taskid);
    if(specs[3+nM+nI] =="1" && (initFILE = fopen(filename,"r") ) ){
    
        specs[3+nM+nI] ="2"; // "2" indicates it is not the initial value any more.
        int tmp1;
        double tmp2;
        while (!feof(initFILE)){
        
        fscanf(initFILE, "%d %d %lf ", &beta, &tmp1, &tmp2);
        printf("changed beta in start: %d\n\n", beta);
        for (Index i=0;i<Time;i++) {
     	    for (Index j=0;j<nY;j++) {
        	    fscanf(initFILE,"%lf ", &x[j*(Time+1)+i]);
            }
      	    for (Index j=0;j<nY;j++) {
        	    fscanf(initFILE,"%lf ", &x[(nY+2*nU)*(Time+1)+j*Time+i]);
		    }
        }
  	    for (Index j=0;j<nY;j++) {
     	    fscanf(initFILE,"%lf ", &x[j*(Time+1)+Time]);
        }
  	    for (Index j=0;j<nP;j++) {
     	    fscanf(initFILE,"%lf ", &x[(2*Time+1)*(nY+2*nU)+j]);
        }
        
        }//endwhile
	
        beta+=delta_beta;
        for (Index i=0;i<nY;i++) {
		    bounds[i][3]=pow(alpha,beta)*Rf0[i];
	    }
        fclose (initFILE);
        for(Index i=0;i<Ntotal;i++) solution[i] = x[i];
    
    }else if(specs[3+nM+nI] =="2"){
    
        for(Index i=0;i<Ntotal;i++) x[i] = solution[i];
        
	}else{
	
	    specs[3+nM+nI] ="2"; // "2" indicates it is not the initial value any more.
        for(Index jt=0;jt<Time+1;jt++) {
            for(Index var=0;var<nY;var++) {
                // Initial conditions for x
                for(int i=0; i<taskid+2;i++) x[(Time+1)*var+jt] = rand()*1.0/RAND_MAX*(bounds[var][1]-bounds[var][0])+bounds[var][0];
                // Initial conditions for midpoints
                if(jt<Time) {
                    for(int i=0; i<taskid+2;i++) x[(Time+1)*(nY+2*nU)+Time*var+jt] = rand()*1.0/RAND_MAX*(bounds[var][1]-bounds[var][0])+bounds[var][0];  
		        }
		    }
        } // End for loop
        
        for(Index par=0;par<nP;par++) {
            // Initial conditions for p5
            for(int i=0; i<taskid+2;i++) x[2*Time*(nY+2*nU)+nY+2*nU+par]=rand()*1.0/RAND_MAX*(bounds[nY+2*nU+par][1]-bounds[nY+2*nU+par][0])+bounds[nY+2*nU+par][0];
        }
	}
    return true;
}   




// returns the value of the objective function
bool SIMPLE_NAKL_NLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  assert(n == 14*Time+7);
  obj_value = 0;

  for(Index jt=0;jt<Time;jt++) {

     for(Index i=0;i<nY;i++) {
        Xval[i] = x[jt + i*(Time+1)];
        Xvalp1[i] = x[jt + i*(Time+1) + 1];
        Xval2[i] = x[(Time+1)*(nY+2*nU) + jt + i*(Time)];
     } //end for loop


     for(Index i=0;i<nU;i++) {
        K11val[i] = x[jt + nY*(Time+1) + 2*i*(Time+1)];
        K11valp1[i] = x[jt + nY*(Time+1) + 2*i*(Time+1) + 1];
        K11val2[i] = x[(Time+1)*(nY+2*nU) + (nY+2*i)*Time + jt];
        dK11val[i] = x[jt + (nY+2*i+1)*(Time+1)];
        dK11valp1[i] = x[jt + (nY+2*i+1)*(Time+1)+1];
        dK11val2[i] = x[(Time+1)*(nY+2*nU) + (nY+2*i+1)*Time + jt];
     } //end for loop

     Xdval[0] = VDATA0[2*jt];
     Xdval2[0] = VDATA0[2*jt+1];
     Xdvalp1[0] = VDATA0[2*jt+2];
     Ival[0] = Iinj[2*jt];
     Ival2[0] = Iinj[2*jt+1];
     Ivalp1[0] = Iinj[2*jt+2];

     for(Index i=0;i<nP;i++) {
        Pval[i] = x[(2*Time+1)*(nY+2*nU)+i];
     } //end for loop


    obj_value += pow(Xdval[0] - Xval[0], 2) + pow(Xdval2[0] - Xval2[0], 2); 

    obj_value += bounds[0][3]*(pow(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48), 2) + pow(0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0] + 0.125*hstep*((1.0L/10.0L)*Ival[0] - 1.0L/10.0L*Ivalp1[0] + (9.0L/10.0L)*Xval[1] - 9.0L/10.0L*Xvalp1[1] - Xval[0] + Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) - 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) - 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90)), 2)); 

    obj_value += bounds[1][3]*(pow(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48), 2) + pow(0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1] + 0.125*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) - 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] + Xvalp1[1] + (9.0L/10.0L)*Xval[0] - 9.0L/10.0L*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])), 2)); 

    obj_value += bounds[2][3]*(pow(0.166666666666667*hstep*((-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-4*Xval2[2] + 2.0*tanh(0.05*Xval2[0] + 1.75) + 2.0)/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) + (-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35)) + Xval[2] - Xvalp1[2], 2) + pow(0.125*hstep*((Xvalp1[2] - 0.5*tanh(0.05*Xvalp1[0] + 1.75) - 0.5)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35)) + 0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2], 2)); 

    obj_value += bounds[3][3]*(pow(0.166666666666667*hstep*(-100.0*Xval[3] - 400.0*Xval2[3] - 100.0*Xvalp1[3] + 50.0*tanh(0.0526*Xval[0] + 1.578) + 200.0*tanh(0.0526*Xval2[0] + 1.578) + 50.0*tanh(0.0526*Xvalp1[0] + 1.578) + 300.0) + Xval[3] - Xvalp1[3], 2) + pow(0.125*hstep*(-100.0*Xval[3] + 100.0*Xvalp1[3] + 50.0*tanh(0.0526*Xval[0] + 1.578) - 50.0*tanh(0.0526*Xvalp1[0] + 1.578)) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3], 2)); 

    obj_value += bounds[4][3]*(pow(Xval[4] - Xvalp1[4] + 0.166666666666667*hstep*((-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-4*Xval2[4] + 2.0*tanh(-0.071429*Xval2[0] - 3.214305) + 2.0)/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) + (-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475)), 2) + pow(0.5*Xval[4] - Xval2[4] + 0.5*Xvalp1[4] + 0.125*hstep*((Xvalp1[4] - 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) - 0.5)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475)), 2)); 

    obj_value += bounds[5][3]*(pow(0.166666666666667*hstep*(-Xval[5] - 4*Xval2[5] - Xvalp1[5] + 0.5*tanh(0.05*Xval[1] + 2.0) + 2.0*tanh(0.05*Xval2[1] + 2.0) + 0.5*tanh(0.05*Xvalp1[1] + 2.0) + 3.0) + Xval[5] - Xvalp1[5], 2) + pow(0.125*hstep*(-Xval[5] + Xvalp1[5] + 0.5*tanh(0.05*Xval[1] + 2.0) - 0.5*tanh(0.05*Xvalp1[1] + 2.0)) + 0.5*Xval[5] - Xval2[5] + 0.5*Xvalp1[5], 2)); 

    obj_value += bounds[6][3]*(pow(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1), 2) + pow(0.5*Xval[6] - Xval2[6] + 0.5*Xvalp1[6] + 0.125*hstep*(-0.0303030303030303*Xval[6] + 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])), 2)); 

  } //end for loop

// Add last element
     for(Index i=0;i<nY;i++) {
        Xval[i] = x[Time + i*(Time+1)];
        Xvalp1[i] = 0;
        Xval2[i] = 0;
     } //end for loop


     for(Index i=0;i<nU;i++) {
        K11val[i] = x[Time + nY*(Time+1) + 2*i*(Time+1)];
        K11valp1[i] = 0;
        K11val2[i] = 0;
        dK11val[i] = x[Time + (nY+2*i+1)*(Time+1)];
        dK11valp1[i] = 0;
        dK11val2[i] = 0;
     } //end for loop

     Xdval[0] = VDATA0[2*Time];
     Xdval2[0] = 0;
     Xdvalp1[0] = 0;
     Ival[0] = Iinj[2*Time];
     Ival2[0] = 0;
     Ivalp1[0] = 0;

     for(Index i=0;i<nP;i++) {
        Pval[i] = x[(2*Time+1)*(nY+2*nU)+i];
     } //end for loop


  obj_value += pow(Xdval[0] - Xval[0], 2) + pow(Xdval2[0] - Xval2[0], 2);

  obj_value = obj_value/(2*Time+1);

  return true;
}


// return the gradient of the objective function grad_{x} f(x)
bool SIMPLE_NAKL_NLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  assert(n == 14*Time+7);

  for(Index i=0;i<n;i++) {
     grad_f[i] = 0;
  }

  for(Index jt=0;jt<Time;jt++) {

     for(Index i=0;i<nY;i++) {
        Xval[i] = x[jt + i*(Time+1)];
        Xvalp1[i] = x[jt + i*(Time+1) + 1];
        Xval2[i] = x[(Time+1)*(nY+2*nU) + jt + i*(Time)];
     } //end for loop

     for(Index i=0;i<nU;i++) {
        K11val[i] = x[jt + nY*(Time+1) + 2*i*(Time+1)];
        K11valp1[i] = x[jt + nY*(Time+1) + 2*i*(Time+1) + 1];
        K11val2[i] = x[(Time+1)*(nY+2*nU) + (nY+2*i)*Time + jt];
        dK11val[i] = x[jt + (nY+2*i+1)*(Time+1)];
        dK11valp1[i] = x[jt + (nY+2*i+1)*(Time+1)+1];
        dK11val2[i] = x[(Time+1)*(nY+2*nU) + (nY+2*i+1)*Time + jt];
     } //end for loop

     Xdval[0] = VDATA0[2*jt];
     Xdval2[0] = VDATA0[2*jt+1];
     Xdvalp1[0] = VDATA0[2*jt+2];
     Ival[0] = Iinj[2*jt];
     Ival2[0] = Iinj[2*jt+1];
     Ivalp1[0] = Iinj[2*jt+2];

     for(Index i=0;i<nP;i++) {
        Pval[i] = x[(2*Time+1)*(nY+2*nU)+i];
     } //end for loop

    grad_f[jt+0*(Time+1)] += (-2*Xdval[0] + 2*Xval[0])/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 0*Time + jt] += (-2*Xdval2[0] + 2*Xval2[0])/(2*Time+1);

    grad_f[jt+0*(Time+1)] += bounds[0][3]*((0.25*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 1.0)*(0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0] + 0.125*hstep*((1.0L/10.0L)*Ival[0] - 1.0L/10.0L*Ivalp1[0] + (9.0L/10.0L)*Xval[1] - 9.0L/10.0L*Xvalp1[1] - Xval[0] + Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) - 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) - 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90))) + (0.333333333333333*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 2)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)))/(2*Time+1);
    grad_f[jt+1*(Time+1)] += bounds[0][3]*(0.3*hstep*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)) + 0.225*hstep*(0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0] + 0.125*hstep*((1.0L/10.0L)*Ival[0] - 1.0L/10.0L*Ivalp1[0] + (9.0L/10.0L)*Xval[1] - 9.0L/10.0L*Xvalp1[1] - Xval[0] + Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) - 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) - 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90))))/(2*Time+1);
    grad_f[jt+2*(Time+1)] += bounds[0][3]*(10.6666666666667*hstep*pow(Xval[2], 3)*(-Xval[0] - 90)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)) + 8.0*hstep*pow(Xval[2], 3)*(-Xval[0] - 90)*(0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0] + 0.125*hstep*((1.0L/10.0L)*Ival[0] - 1.0L/10.0L*Ivalp1[0] + (9.0L/10.0L)*Xval[1] - 9.0L/10.0L*Xvalp1[1] - Xval[0] + Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) - 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) - 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90))))/(2*Time+1);
    grad_f[jt+3*(Time+1)] += bounds[0][3]*(70.0*Xval[4]*hstep*pow(Xval[3], 2)*(-Xval[0] + 55)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)) + 52.5*Xval[4]*hstep*pow(Xval[3], 2)*(-Xval[0] + 55)*(0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0] + 0.125*hstep*((1.0L/10.0L)*Ival[0] - 1.0L/10.0L*Ivalp1[0] + (9.0L/10.0L)*Xval[1] - 9.0L/10.0L*Xvalp1[1] - Xval[0] + Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) - 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) - 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90))))/(2*Time+1);
    grad_f[jt+4*(Time+1)] += bounds[0][3]*(23.3333333333333*hstep*pow(Xval[3], 3)*(-Xval[0] + 55)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)) + 17.5*hstep*pow(Xval[3], 3)*(-Xval[0] + 55)*(0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0] + 0.125*hstep*((1.0L/10.0L)*Ival[0] - 1.0L/10.0L*Ivalp1[0] + (9.0L/10.0L)*Xval[1] - 9.0L/10.0L*Xvalp1[1] - Xval[0] + Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) - 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) - 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90))))/(2*Time+1);
    grad_f[jt+1+0*(Time+1)] += bounds[0][3]*((0.333333333333333*hstep*(-70.0*Xvalp1[4]*pow(Xvalp1[3], 3) - 8*pow(Xvalp1[2], 4) - 1) - 2)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)) + (0.25*hstep*(70.0*Xvalp1[4]*pow(Xvalp1[3], 3) + 8*pow(Xvalp1[2], 4) + 1) + 1.0)*(0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0] + 0.125*hstep*((1.0L/10.0L)*Ival[0] - 1.0L/10.0L*Ivalp1[0] + (9.0L/10.0L)*Xval[1] - 9.0L/10.0L*Xvalp1[1] - Xval[0] + Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) - 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) - 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90))))/(2*Time+1);
    grad_f[jt+1+1*(Time+1)] += bounds[0][3]*(0.3*hstep*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)) - 0.225*hstep*(0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0] + 0.125*hstep*((1.0L/10.0L)*Ival[0] - 1.0L/10.0L*Ivalp1[0] + (9.0L/10.0L)*Xval[1] - 9.0L/10.0L*Xvalp1[1] - Xval[0] + Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) - 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) - 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90))))/(2*Time+1);
    grad_f[jt+1+2*(Time+1)] += bounds[0][3]*(10.6666666666667*hstep*pow(Xvalp1[2], 3)*(-Xvalp1[0] - 90)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)) - 8.0*hstep*pow(Xvalp1[2], 3)*(-Xvalp1[0] - 90)*(0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0] + 0.125*hstep*((1.0L/10.0L)*Ival[0] - 1.0L/10.0L*Ivalp1[0] + (9.0L/10.0L)*Xval[1] - 9.0L/10.0L*Xvalp1[1] - Xval[0] + Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) - 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) - 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90))))/(2*Time+1);
    grad_f[jt+1+3*(Time+1)] += bounds[0][3]*(70.0*Xvalp1[4]*hstep*pow(Xvalp1[3], 2)*(-Xvalp1[0] + 55)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)) - 52.5*Xvalp1[4]*hstep*pow(Xvalp1[3], 2)*(-Xvalp1[0] + 55)*(0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0] + 0.125*hstep*((1.0L/10.0L)*Ival[0] - 1.0L/10.0L*Ivalp1[0] + (9.0L/10.0L)*Xval[1] - 9.0L/10.0L*Xvalp1[1] - Xval[0] + Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) - 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) - 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90))))/(2*Time+1);
    grad_f[jt+1+4*(Time+1)] += bounds[0][3]*(23.3333333333333*hstep*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)) - 17.5*hstep*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55)*(0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0] + 0.125*hstep*((1.0L/10.0L)*Ival[0] - 1.0L/10.0L*Ivalp1[0] + (9.0L/10.0L)*Xval[1] - 9.0L/10.0L*Xvalp1[1] - Xval[0] + Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) - 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) - 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90))))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 0*Time + jt] += bounds[0][3]*(-1.0*Xval[0] + 2*Xval2[0] - 1.0*Xvalp1[0] + 0.333333333333333*hstep*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48))*(-280.0*Xval2[4]*pow(Xval2[3], 3) - 32*pow(Xval2[2], 4) - 4) - 0.25*hstep*((1.0L/10.0L)*Ival[0] - 1.0L/10.0L*Ivalp1[0] + (9.0L/10.0L)*Xval[1] - 9.0L/10.0L*Xvalp1[1] - Xval[0] + Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) - 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) - 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90)))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 1*Time + jt] += bounds[0][3]*(1.2*hstep*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 2*Time + jt] += bounds[0][3]*(42.6666666666667*hstep*pow(Xval2[2], 3)*(-Xval2[0] - 90)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 3*Time + jt] += bounds[0][3]*(280.0*Xval2[4]*hstep*pow(Xval2[3], 2)*(-Xval2[0] + 55)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 4*Time + jt] += bounds[0][3]*(93.3333333333333*hstep*pow(Xval2[3], 3)*(-Xval2[0] + 55)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)))/(2*Time+1);
    grad_f[jt+0*(Time+1)] += bounds[1][3]*(0.3*hstep*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48)) + 0.225*hstep*(0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1] + 0.125*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) - 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] + Xvalp1[1] + (9.0L/10.0L)*Xval[0] - 9.0L/10.0L*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]))))/(2*Time+1);
    grad_f[jt+1*(Time+1)] += bounds[1][3]*(// Not C:
// ghkjac(Xval[1],1)
(0.25*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) + 1.0)*(0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1] + 0.125*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) - 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] + Xvalp1[1] + (9.0L/10.0L)*Xval[0] - 9.0L/10.0L*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]))) + (0.333333333333333*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) + 2)*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48)))/(2*Time+1);
    grad_f[jt+5*(Time+1)] += bounds[1][3]*(-0.000666666666666667*hstep*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*ghk(Xval[1]) - 0.0005*hstep*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1] + 0.125*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) - 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] + Xvalp1[1] + (9.0L/10.0L)*Xval[0] - 9.0L/10.0L*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])))*ghk(Xval[1]))/(2*Time+1);
    grad_f[jt+6*(Time+1)] += bounds[1][3]*(0.333333333333333*hstep*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*(-150*pow(Xval[6], 3)*(-Xval[1] - 90)/pow(pow(Xval[6], 2) + 100, 2) + 150*Xval[6]*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 0.001*pow(Xval[5], 2)*ghk(Xval[1])) + 0.25*hstep*(-150*pow(Xval[6], 3)*(-Xval[1] - 90)/pow(pow(Xval[6], 2) + 100, 2) + 150*Xval[6]*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 0.001*pow(Xval[5], 2)*ghk(Xval[1]))*(0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1] + 0.125*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) - 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] + Xvalp1[1] + (9.0L/10.0L)*Xval[0] - 9.0L/10.0L*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]))))/(2*Time+1);
    grad_f[jt+1+0*(Time+1)] += bounds[1][3]*(0.3*hstep*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48)) - 0.225*hstep*(0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1] + 0.125*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) - 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] + Xvalp1[1] + (9.0L/10.0L)*Xval[0] - 9.0L/10.0L*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]))))/(2*Time+1);
    grad_f[jt+1+1*(Time+1)] += bounds[1][3]*(// Not C:
// ghkjac(Xvalp1[1],1)
(0.333333333333333*hstep*(-75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) + 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) - 1) - 2)*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48)) + (0.25*hstep*(75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) - 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) + 1) + 1.0)*(0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1] + 0.125*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) - 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] + Xvalp1[1] + (9.0L/10.0L)*Xval[0] - 9.0L/10.0L*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]))))/(2*Time+1);
    grad_f[jt+1+5*(Time+1)] += bounds[1][3]*(-0.000666666666666667*hstep*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*ghk(Xvalp1[1]) + 0.0005*hstep*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1] + 0.125*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) - 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] + Xvalp1[1] + (9.0L/10.0L)*Xval[0] - 9.0L/10.0L*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])))*ghk(Xvalp1[1]))/(2*Time+1);
    grad_f[jt+1+6*(Time+1)] += bounds[1][3]*(0.333333333333333*hstep*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*(-150*pow(Xvalp1[6], 3)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) + 150*Xvalp1[6]*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])) + 0.25*hstep*(150*pow(Xvalp1[6], 3)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) - 150*Xvalp1[6]*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]))*(0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1] + 0.125*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) - 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] + Xvalp1[1] + (9.0L/10.0L)*Xval[0] - 9.0L/10.0L*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]))))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 0*Time + jt] += bounds[1][3]*(1.2*hstep*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48)))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 1*Time + jt] += bounds[1][3]*(// Not C:
// ghkjac(Xval2[1],1)
-1.0*Xval[1] + 2*Xval2[1] - 1.0*Xvalp1[1] + 0.333333333333333*hstep*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*(-300*pow(Xval2[6], 2)/(pow(Xval2[6], 2) + 100) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkjac(Xval2[1],1) + 0.740740740740741*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1]) - 4) - 0.25*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) - 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] + Xvalp1[1] + (9.0L/10.0L)*Xval[0] - 9.0L/10.0L*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 5*Time + jt] += bounds[1][3]*(-0.00266666666666667*hstep*Xval2[5]*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*ghk(Xval2[1]))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 6*Time + jt] += bounds[1][3]*(0.333333333333333*hstep*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*(-600*pow(Xval2[6], 3)*(-Xval2[1] - 90)/pow(pow(Xval2[6], 2) + 100, 2) + 600*Xval2[6]*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 0.004*pow(Xval2[5], 2)*ghk(Xval2[1])))/(2*Time+1);
    grad_f[jt+0*(Time+1)] += bounds[2][3]*(0.333333333333333*hstep*((-0.025*pow(tanh(0.05*Xval[0] + 1.75), 2) + 0.025)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + (-0.008325*pow(tanh(-0.0333*Xval[0] - 0.8991), 2) + 0.008325)*(-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35, 2))*(0.166666666666667*hstep*((-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-4*Xval2[2] + 2.0*tanh(0.05*Xval2[0] + 1.75) + 2.0)/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) + (-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35)) + Xval[2] - Xvalp1[2]) + 0.25*hstep*((-0.025*pow(tanh(0.05*Xval[0] + 1.75), 2) + 0.025)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + (-0.008325*pow(tanh(-0.0333*Xval[0] - 0.8991), 2) + 0.008325)*(-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35, 2))*(0.125*hstep*((Xvalp1[2] - 0.5*tanh(0.05*Xvalp1[0] + 1.75) - 0.5)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35)) + 0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2]))/(2*Time+1);
    grad_f[jt+2*(Time+1)] += bounds[2][3]*((-0.333333333333333*hstep/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + 2)*(0.166666666666667*hstep*((-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-4*Xval2[2] + 2.0*tanh(0.05*Xval2[0] + 1.75) + 2.0)/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) + (-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35)) + Xval[2] - Xvalp1[2]) + (-0.25*hstep/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + 1.0)*(0.125*hstep*((Xvalp1[2] - 0.5*tanh(0.05*Xvalp1[0] + 1.75) - 0.5)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35)) + 0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2]))/(2*Time+1);
    grad_f[jt+1+0*(Time+1)] += bounds[2][3]*(0.333333333333333*hstep*((-0.025*pow(tanh(0.05*Xvalp1[0] + 1.75), 2) + 0.025)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-0.008325*pow(tanh(-0.0333*Xvalp1[0] - 0.8991), 2) + 0.008325)*(-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35, 2))*(0.166666666666667*hstep*((-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-4*Xval2[2] + 2.0*tanh(0.05*Xval2[0] + 1.75) + 2.0)/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) + (-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35)) + Xval[2] - Xvalp1[2]) + 0.25*hstep*((0.025*pow(tanh(0.05*Xvalp1[0] + 1.75), 2) - 0.025)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) - (-0.008325*pow(tanh(-0.0333*Xvalp1[0] - 0.8991), 2) + 0.008325)*(-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35, 2))*(0.125*hstep*((Xvalp1[2] - 0.5*tanh(0.05*Xvalp1[0] + 1.75) - 0.5)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35)) + 0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2]))/(2*Time+1);
    grad_f[jt+1+2*(Time+1)] += bounds[2][3]*((-0.333333333333333*hstep/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) - 2)*(0.166666666666667*hstep*((-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-4*Xval2[2] + 2.0*tanh(0.05*Xval2[0] + 1.75) + 2.0)/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) + (-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35)) + Xval[2] - Xvalp1[2]) + (0.25*hstep/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + 1.0)*(0.125*hstep*((Xvalp1[2] - 0.5*tanh(0.05*Xvalp1[0] + 1.75) - 0.5)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35)) + 0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2]))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 0*Time + jt] += bounds[2][3]*(0.333333333333333*hstep*((-0.1*pow(tanh(0.05*Xval2[0] + 1.75), 2) + 0.1)/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) + (-0.0333*pow(tanh(-0.0333*Xval2[0] - 0.8991), 2) + 0.0333)*(-Xval2[2] + 0.5*tanh(0.05*Xval2[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35, 2))*(0.166666666666667*hstep*((-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-4*Xval2[2] + 2.0*tanh(0.05*Xval2[0] + 1.75) + 2.0)/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) + (-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35)) + Xval[2] - Xvalp1[2]))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 2*Time + jt] += bounds[2][3]*(-0.25*hstep*((Xvalp1[2] - 0.5*tanh(0.05*Xvalp1[0] + 1.75) - 0.5)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35)) - 1.33333333333333*hstep*(0.166666666666667*hstep*((-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-4*Xval2[2] + 2.0*tanh(0.05*Xval2[0] + 1.75) + 2.0)/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) + (-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35)) + Xval[2] - Xvalp1[2])/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) - 1.0*Xval[2] + 2*Xval2[2] - 1.0*Xvalp1[2])/(2*Time+1);
    grad_f[jt+0*(Time+1)] += bounds[3][3]*(0.333333333333333*hstep*(-2.63*pow(tanh(0.0526*Xval[0] + 1.578), 2) + 2.63)*(0.166666666666667*hstep*(-100.0*Xval[3] - 400.0*Xval2[3] - 100.0*Xvalp1[3] + 50.0*tanh(0.0526*Xval[0] + 1.578) + 200.0*tanh(0.0526*Xval2[0] + 1.578) + 50.0*tanh(0.0526*Xvalp1[0] + 1.578) + 300.0) + Xval[3] - Xvalp1[3]) + 0.25*hstep*(-2.63*pow(tanh(0.0526*Xval[0] + 1.578), 2) + 2.63)*(0.125*hstep*(-100.0*Xval[3] + 100.0*Xvalp1[3] + 50.0*tanh(0.0526*Xval[0] + 1.578) - 50.0*tanh(0.0526*Xvalp1[0] + 1.578)) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3]))/(2*Time+1);
    grad_f[jt+3*(Time+1)] += bounds[3][3]*((-33.3333333333333*hstep + 2)*(0.166666666666667*hstep*(-100.0*Xval[3] - 400.0*Xval2[3] - 100.0*Xvalp1[3] + 50.0*tanh(0.0526*Xval[0] + 1.578) + 200.0*tanh(0.0526*Xval2[0] + 1.578) + 50.0*tanh(0.0526*Xvalp1[0] + 1.578) + 300.0) + Xval[3] - Xvalp1[3]) + (-25.0*hstep + 1.0)*(0.125*hstep*(-100.0*Xval[3] + 100.0*Xvalp1[3] + 50.0*tanh(0.0526*Xval[0] + 1.578) - 50.0*tanh(0.0526*Xvalp1[0] + 1.578)) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3]))/(2*Time+1);
    grad_f[jt+1+0*(Time+1)] += bounds[3][3]*(0.333333333333333*hstep*(-2.63*pow(tanh(0.0526*Xvalp1[0] + 1.578), 2) + 2.63)*(0.166666666666667*hstep*(-100.0*Xval[3] - 400.0*Xval2[3] - 100.0*Xvalp1[3] + 50.0*tanh(0.0526*Xval[0] + 1.578) + 200.0*tanh(0.0526*Xval2[0] + 1.578) + 50.0*tanh(0.0526*Xvalp1[0] + 1.578) + 300.0) + Xval[3] - Xvalp1[3]) + 0.25*hstep*(2.63*pow(tanh(0.0526*Xvalp1[0] + 1.578), 2) - 2.63)*(0.125*hstep*(-100.0*Xval[3] + 100.0*Xvalp1[3] + 50.0*tanh(0.0526*Xval[0] + 1.578) - 50.0*tanh(0.0526*Xvalp1[0] + 1.578)) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3]))/(2*Time+1);
    grad_f[jt+1+3*(Time+1)] += bounds[3][3]*((-33.3333333333333*hstep - 2)*(0.166666666666667*hstep*(-100.0*Xval[3] - 400.0*Xval2[3] - 100.0*Xvalp1[3] + 50.0*tanh(0.0526*Xval[0] + 1.578) + 200.0*tanh(0.0526*Xval2[0] + 1.578) + 50.0*tanh(0.0526*Xvalp1[0] + 1.578) + 300.0) + Xval[3] - Xvalp1[3]) + (25.0*hstep + 1.0)*(0.125*hstep*(-100.0*Xval[3] + 100.0*Xvalp1[3] + 50.0*tanh(0.0526*Xval[0] + 1.578) - 50.0*tanh(0.0526*Xvalp1[0] + 1.578)) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3]))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 0*Time + jt] += bounds[3][3]*(0.333333333333333*hstep*(-10.52*pow(tanh(0.0526*Xval2[0] + 1.578), 2) + 10.52)*(0.166666666666667*hstep*(-100.0*Xval[3] - 400.0*Xval2[3] - 100.0*Xvalp1[3] + 50.0*tanh(0.0526*Xval[0] + 1.578) + 200.0*tanh(0.0526*Xval2[0] + 1.578) + 50.0*tanh(0.0526*Xvalp1[0] + 1.578) + 300.0) + Xval[3] - Xvalp1[3]))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 3*Time + jt] += bounds[3][3]*(-133.333333333333*hstep*(0.166666666666667*hstep*(-100.0*Xval[3] - 400.0*Xval2[3] - 100.0*Xvalp1[3] + 50.0*tanh(0.0526*Xval[0] + 1.578) + 200.0*tanh(0.0526*Xval2[0] + 1.578) + 50.0*tanh(0.0526*Xvalp1[0] + 1.578) + 300.0) + Xval[3] - Xvalp1[3]) - 0.25*hstep*(-100.0*Xval[3] + 100.0*Xvalp1[3] + 50.0*tanh(0.0526*Xval[0] + 1.578) - 50.0*tanh(0.0526*Xvalp1[0] + 1.578)) - 1.0*Xval[3] + 2*Xval2[3] - 1.0*Xvalp1[3])/(2*Time+1);
    grad_f[jt+0*(Time+1)] += bounds[4][3]*(0.333333333333333*hstep*((0.0357145*pow(tanh(-0.071429*Xval[0] - 3.214305), 2) - 0.0357145)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + (-0.031249875*pow(tanh(-0.083333*Xval[0] - 3.3749865), 2) + 0.031249875)*(-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475, 2))*(Xval[4] - Xvalp1[4] + 0.166666666666667*hstep*((-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-4*Xval2[4] + 2.0*tanh(-0.071429*Xval2[0] - 3.214305) + 2.0)/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) + (-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475))) + 0.25*hstep*((0.0357145*pow(tanh(-0.071429*Xval[0] - 3.214305), 2) - 0.0357145)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + (-0.031249875*pow(tanh(-0.083333*Xval[0] - 3.3749865), 2) + 0.031249875)*(-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475, 2))*(0.5*Xval[4] - Xval2[4] + 0.5*Xvalp1[4] + 0.125*hstep*((Xvalp1[4] - 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) - 0.5)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475))))/(2*Time+1);
    grad_f[jt+4*(Time+1)] += bounds[4][3]*((-0.333333333333333*hstep/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + 2)*(Xval[4] - Xvalp1[4] + 0.166666666666667*hstep*((-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-4*Xval2[4] + 2.0*tanh(-0.071429*Xval2[0] - 3.214305) + 2.0)/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) + (-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475))) + (-0.25*hstep/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + 1.0)*(0.5*Xval[4] - Xval2[4] + 0.5*Xvalp1[4] + 0.125*hstep*((Xvalp1[4] - 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) - 0.5)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475))))/(2*Time+1);
    grad_f[jt+1+0*(Time+1)] += bounds[4][3]*(0.25*hstep*((-0.0357145*pow(tanh(-0.071429*Xvalp1[0] - 3.214305), 2) + 0.0357145)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) - (-0.031249875*pow(tanh(-0.083333*Xvalp1[0] - 3.3749865), 2) + 0.031249875)*(-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475, 2))*(0.5*Xval[4] - Xval2[4] + 0.5*Xvalp1[4] + 0.125*hstep*((Xvalp1[4] - 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) - 0.5)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475))) + 0.333333333333333*hstep*((0.0357145*pow(tanh(-0.071429*Xvalp1[0] - 3.214305), 2) - 0.0357145)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-0.031249875*pow(tanh(-0.083333*Xvalp1[0] - 3.3749865), 2) + 0.031249875)*(-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475, 2))*(Xval[4] - Xvalp1[4] + 0.166666666666667*hstep*((-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-4*Xval2[4] + 2.0*tanh(-0.071429*Xval2[0] - 3.214305) + 2.0)/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) + (-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475))))/(2*Time+1);
    grad_f[jt+1+4*(Time+1)] += bounds[4][3]*((-0.333333333333333*hstep/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) - 2)*(Xval[4] - Xvalp1[4] + 0.166666666666667*hstep*((-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-4*Xval2[4] + 2.0*tanh(-0.071429*Xval2[0] - 3.214305) + 2.0)/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) + (-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475))) + (0.25*hstep/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + 1.0)*(0.5*Xval[4] - Xval2[4] + 0.5*Xvalp1[4] + 0.125*hstep*((Xvalp1[4] - 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) - 0.5)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475))))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 0*Time + jt] += bounds[4][3]*(0.333333333333333*hstep*((0.142858*pow(tanh(-0.071429*Xval2[0] - 3.214305), 2) - 0.142858)/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) + (-0.1249995*pow(tanh(-0.083333*Xval2[0] - 3.3749865), 2) + 0.1249995)*(-Xval2[4] + 0.5*tanh(-0.071429*Xval2[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475, 2))*(Xval[4] - Xvalp1[4] + 0.166666666666667*hstep*((-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-4*Xval2[4] + 2.0*tanh(-0.071429*Xval2[0] - 3.214305) + 2.0)/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) + (-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475))))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 4*Time + jt] += bounds[4][3]*(-1.0*Xval[4] + 2*Xval2[4] - 1.0*Xvalp1[4] - 0.25*hstep*((Xvalp1[4] - 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) - 0.5)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475)) - 1.33333333333333*hstep*(Xval[4] - Xvalp1[4] + 0.166666666666667*hstep*((-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-4*Xval2[4] + 2.0*tanh(-0.071429*Xval2[0] - 3.214305) + 2.0)/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) + (-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475)))/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475))/(2*Time+1);
    grad_f[jt+1*(Time+1)] += bounds[5][3]*(0.333333333333333*hstep*(-0.025*pow(tanh(0.05*Xval[1] + 2.0), 2) + 0.025)*(0.166666666666667*hstep*(-Xval[5] - 4*Xval2[5] - Xvalp1[5] + 0.5*tanh(0.05*Xval[1] + 2.0) + 2.0*tanh(0.05*Xval2[1] + 2.0) + 0.5*tanh(0.05*Xvalp1[1] + 2.0) + 3.0) + Xval[5] - Xvalp1[5]) + 0.25*hstep*(-0.025*pow(tanh(0.05*Xval[1] + 2.0), 2) + 0.025)*(0.125*hstep*(-Xval[5] + Xvalp1[5] + 0.5*tanh(0.05*Xval[1] + 2.0) - 0.5*tanh(0.05*Xvalp1[1] + 2.0)) + 0.5*Xval[5] - Xval2[5] + 0.5*Xvalp1[5]))/(2*Time+1);
    grad_f[jt+5*(Time+1)] += bounds[5][3]*((-0.333333333333333*hstep + 2)*(0.166666666666667*hstep*(-Xval[5] - 4*Xval2[5] - Xvalp1[5] + 0.5*tanh(0.05*Xval[1] + 2.0) + 2.0*tanh(0.05*Xval2[1] + 2.0) + 0.5*tanh(0.05*Xvalp1[1] + 2.0) + 3.0) + Xval[5] - Xvalp1[5]) + (-0.25*hstep + 1.0)*(0.125*hstep*(-Xval[5] + Xvalp1[5] + 0.5*tanh(0.05*Xval[1] + 2.0) - 0.5*tanh(0.05*Xvalp1[1] + 2.0)) + 0.5*Xval[5] - Xval2[5] + 0.5*Xvalp1[5]))/(2*Time+1);
    grad_f[jt+1+1*(Time+1)] += bounds[5][3]*(0.333333333333333*hstep*(-0.025*pow(tanh(0.05*Xvalp1[1] + 2.0), 2) + 0.025)*(0.166666666666667*hstep*(-Xval[5] - 4*Xval2[5] - Xvalp1[5] + 0.5*tanh(0.05*Xval[1] + 2.0) + 2.0*tanh(0.05*Xval2[1] + 2.0) + 0.5*tanh(0.05*Xvalp1[1] + 2.0) + 3.0) + Xval[5] - Xvalp1[5]) + 0.25*hstep*(0.025*pow(tanh(0.05*Xvalp1[1] + 2.0), 2) - 0.025)*(0.125*hstep*(-Xval[5] + Xvalp1[5] + 0.5*tanh(0.05*Xval[1] + 2.0) - 0.5*tanh(0.05*Xvalp1[1] + 2.0)) + 0.5*Xval[5] - Xval2[5] + 0.5*Xvalp1[5]))/(2*Time+1);
    grad_f[jt+1+5*(Time+1)] += bounds[5][3]*((-0.333333333333333*hstep - 2)*(0.166666666666667*hstep*(-Xval[5] - 4*Xval2[5] - Xvalp1[5] + 0.5*tanh(0.05*Xval[1] + 2.0) + 2.0*tanh(0.05*Xval2[1] + 2.0) + 0.5*tanh(0.05*Xvalp1[1] + 2.0) + 3.0) + Xval[5] - Xvalp1[5]) + (0.25*hstep + 1.0)*(0.125*hstep*(-Xval[5] + Xvalp1[5] + 0.5*tanh(0.05*Xval[1] + 2.0) - 0.5*tanh(0.05*Xvalp1[1] + 2.0)) + 0.5*Xval[5] - Xval2[5] + 0.5*Xvalp1[5]))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 1*Time + jt] += bounds[5][3]*(0.333333333333333*hstep*(-0.1*pow(tanh(0.05*Xval2[1] + 2.0), 2) + 0.1)*(0.166666666666667*hstep*(-Xval[5] - 4*Xval2[5] - Xvalp1[5] + 0.5*tanh(0.05*Xval[1] + 2.0) + 2.0*tanh(0.05*Xval2[1] + 2.0) + 0.5*tanh(0.05*Xvalp1[1] + 2.0) + 3.0) + Xval[5] - Xvalp1[5]))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 5*Time + jt] += bounds[5][3]*(-1.33333333333333*hstep*(0.166666666666667*hstep*(-Xval[5] - 4*Xval2[5] - Xvalp1[5] + 0.5*tanh(0.05*Xval[1] + 2.0) + 2.0*tanh(0.05*Xval2[1] + 2.0) + 0.5*tanh(0.05*Xvalp1[1] + 2.0) + 3.0) + Xval[5] - Xvalp1[5]) - 0.25*hstep*(-Xval[5] + Xvalp1[5] + 0.5*tanh(0.05*Xval[1] + 2.0) - 0.5*tanh(0.05*Xvalp1[1] + 2.0)) - 1.0*Xval[5] + 2*Xval2[5] - 1.0*Xvalp1[5])/(2*Time+1);
    grad_f[jt+1*(Time+1)] += bounds[6][3]*(// Not C:
// ghkjac(Xval[1],1)
0.333333333333333*hstep*(-1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.000222222222222222*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]))*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1)) + 0.25*hstep*(-1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.000222222222222222*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]))*(0.5*Xval[6] - Xval2[6] + 0.5*Xvalp1[6] + 0.125*hstep*(-0.0303030303030303*Xval[6] + 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]))))/(2*Time+1);
    grad_f[jt+5*(Time+1)] += bounds[6][3]*(-8.0e-7*hstep*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1))*ghk(Xval[1]) - 6.0e-7*hstep*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(0.5*Xval[6] - Xval2[6] + 0.5*Xvalp1[6] + 0.125*hstep*(-0.0303030303030303*Xval[6] + 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])))*ghk(Xval[1]))/(2*Time+1);
    grad_f[jt+6*(Time+1)] += bounds[6][3]*((0.25*hstep*(1.2e-6*pow(Xval[5], 2)*ghk(Xval[1]) - 0.0303030303030303) + 1.0)*(0.5*Xval[6] - Xval2[6] + 0.5*Xvalp1[6] + 0.125*hstep*(-0.0303030303030303*Xval[6] + 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]))) + (0.333333333333333*hstep*(1.2e-6*pow(Xval[5], 2)*ghk(Xval[1]) - 0.0303030303030303) + 2)*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1)))/(2*Time+1);
    grad_f[jt+1+1*(Time+1)] += bounds[6][3]*(// Not C:
// ghkjac(Xvalp1[1],1)
0.333333333333333*hstep*(-1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) + 0.000222222222222222*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]))*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1)) + 0.25*hstep*(1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) - 0.000222222222222222*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]))*(0.5*Xval[6] - Xval2[6] + 0.5*Xvalp1[6] + 0.125*hstep*(-0.0303030303030303*Xval[6] + 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]))))/(2*Time+1);
    grad_f[jt+1+5*(Time+1)] += bounds[6][3]*(-8.0e-7*hstep*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1))*ghk(Xvalp1[1]) + 6.0e-7*hstep*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(0.5*Xval[6] - Xval2[6] + 0.5*Xvalp1[6] + 0.125*hstep*(-0.0303030303030303*Xval[6] + 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])))*ghk(Xvalp1[1]))/(2*Time+1);
    grad_f[jt+1+6*(Time+1)] += bounds[6][3]*((0.25*hstep*(-1.2e-6*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]) + 0.0303030303030303) + 1.0)*(0.5*Xval[6] - Xval2[6] + 0.5*Xvalp1[6] + 0.125*hstep*(-0.0303030303030303*Xval[6] + 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]))) + (0.333333333333333*hstep*(1.2e-6*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]) - 0.0303030303030303) - 2)*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1)))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 1*Time + jt] += bounds[6][3]*(// Not C:
// ghkjac(Xval2[1],1)
0.333333333333333*hstep*(-4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkjac(Xval2[1],1) + 0.000888888888888889*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1]))*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1)))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 5*Time + jt] += bounds[6][3]*(-3.2e-6*hstep*Xval2[5]*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1))*ghk(Xval2[1]))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 6*Time + jt] += bounds[6][3]*(-1.0*Xval[6] + 2*Xval2[6] - 1.0*Xvalp1[6] + 0.333333333333333*hstep*(4.8e-6*pow(Xval2[5], 2)*ghk(Xval2[1]) - 0.121212121212121)*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1)) - 0.25*hstep*(-0.0303030303030303*Xval[6] + 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])))/(2*Time+1);

  } //end for loop

// Add last element
     for(Index i=0;i<nY;i++) {
        Xval[i] = x[Time + i*(Time+1)];
        Xvalp1[i] = 0;
        Xval2[i] = 0;
     } //end for loop

     for(Index i=0;i<nU;i++) {
        K11val[i] = x[Time + nY*(Time+1) + 2*i*(Time+1)];
        K11valp1[i] = 0;
        K11val2[i] = 0;
        dK11val[i] = x[Time + (nY+2*i+1)*(Time+1)];
        dK11valp1[i] = 0;
        dK11val2[i] = 0;
     } //end for loop

     Xdval[0] = VDATA0[2*Time];
     Xdval2[0] = 0;
     Xdvalp1[0] = 0;
     Ival[0] = Iinj[2*Time];
     Ival2[0] = 0;
     Ivalp1[0] = 0;

     for(Index i=0;i<nP;i++) {
        Pval[i] = x[(2*Time+1)*(nY+2*nU)+i];
     } //end for loop


    grad_f[Time+0*(Time+1)] += (-2*Xdval[0] + 2*Xval[0])/(2*Time+1);

  return true;
}


// return the value of the constraints: g(x)
bool SIMPLE_NAKL_NLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  assert(n == 14*Time+7);
  assert(m == 0);

  return true;
}


// return the structure or values of the jacobian
bool SIMPLE_NAKL_NLP::eval_jac_g(Index n, const Number* x, bool new_x,
                            Index m, Index nele_jac, Index* iRow, Index* jCol,
                            Number* values)
{

return true;
}


// return the structure or values of the hessian
bool SIMPLE_NAKL_NLP::eval_h(Index n, const Number* x, bool new_x,
                       Number obj_factor, Index m, const Number* lambda,
                       bool new_lambda, Index nele_hess, Index* iRow,
                       Index* jCol, Number* values)
{

if (values == NULL) {
   // return the structure.  This is a symmetric matrix, fill in the lower left
   // triangle only.

   // Each non-one Hessian element has its own explicit loop
   // since each element needs a different number of matrix elements


   for(Index jt=0;jt<Time+1;jt++) {
     iRow[0*(Time+1)+0*(Time)+0+jt] = (Time+1)*0+jt;
     jCol[0*(Time+1)+0*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[1*(Time+1)+0*(Time)+0+jt] = (Time+1)*0+jt+1;
     jCol[1*(Time+1)+0*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[1*(Time+1)+1*(Time)+0+jt] = (Time+1)*1+jt;
     jCol[1*(Time+1)+1*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[2*(Time+1)+1*(Time)+0+jt] = (Time+1)*1+jt;
     jCol[2*(Time+1)+1*(Time)+0+jt] = (Time+1)*0+jt+1;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[2*(Time+1)+2*(Time)+0+jt] = (Time+1)*1+jt;
     jCol[2*(Time+1)+2*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[3*(Time+1)+2*(Time)+0+jt] = (Time+1)*1+jt+1;
     jCol[3*(Time+1)+2*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[3*(Time+1)+3*(Time)+0+jt] = (Time+1)*1+jt+1;
     jCol[3*(Time+1)+3*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[3*(Time+1)+4*(Time)+0+jt] = (Time+1)*2+jt;
     jCol[3*(Time+1)+4*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[4*(Time+1)+4*(Time)+0+jt] = (Time+1)*2+jt;
     jCol[4*(Time+1)+4*(Time)+0+jt] = (Time+1)*0+jt+1;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[4*(Time+1)+5*(Time)+0+jt] = (Time+1)*2+jt;
     jCol[4*(Time+1)+5*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[5*(Time+1)+5*(Time)+0+jt] = (Time+1)*2+jt;
     jCol[5*(Time+1)+5*(Time)+0+jt] = (Time+1)*1+jt+1;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[5*(Time+1)+6*(Time)+0+jt] = (Time+1)*2+jt;
     jCol[5*(Time+1)+6*(Time)+0+jt] = (Time+1)*2+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[6*(Time+1)+6*(Time)+0+jt] = (Time+1)*2+jt+1;
     jCol[6*(Time+1)+6*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[6*(Time+1)+7*(Time)+0+jt] = (Time+1)*2+jt+1;
     jCol[6*(Time+1)+7*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[6*(Time+1)+8*(Time)+0+jt] = (Time+1)*2+jt+1;
     jCol[6*(Time+1)+8*(Time)+0+jt] = (Time+1)*2+jt;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[6*(Time+1)+9*(Time)+0+jt] = (Time+1)*3+jt;
     jCol[6*(Time+1)+9*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[7*(Time+1)+9*(Time)+0+jt] = (Time+1)*3+jt;
     jCol[7*(Time+1)+9*(Time)+0+jt] = (Time+1)*0+jt+1;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[7*(Time+1)+10*(Time)+0+jt] = (Time+1)*3+jt;
     jCol[7*(Time+1)+10*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[8*(Time+1)+10*(Time)+0+jt] = (Time+1)*3+jt;
     jCol[8*(Time+1)+10*(Time)+0+jt] = (Time+1)*1+jt+1;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[8*(Time+1)+11*(Time)+0+jt] = (Time+1)*3+jt;
     jCol[8*(Time+1)+11*(Time)+0+jt] = (Time+1)*2+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[9*(Time+1)+11*(Time)+0+jt] = (Time+1)*3+jt;
     jCol[9*(Time+1)+11*(Time)+0+jt] = (Time+1)*2+jt+1;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[9*(Time+1)+12*(Time)+0+jt] = (Time+1)*3+jt;
     jCol[9*(Time+1)+12*(Time)+0+jt] = (Time+1)*3+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+12*(Time)+0+jt] = (Time+1)*3+jt+1;
     jCol[10*(Time+1)+12*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+13*(Time)+0+jt] = (Time+1)*3+jt+1;
     jCol[10*(Time+1)+13*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+14*(Time)+0+jt] = (Time+1)*3+jt+1;
     jCol[10*(Time+1)+14*(Time)+0+jt] = (Time+1)*2+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+15*(Time)+0+jt] = (Time+1)*3+jt+1;
     jCol[10*(Time+1)+15*(Time)+0+jt] = (Time+1)*3+jt;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[10*(Time+1)+16*(Time)+0+jt] = (Time+1)*4+jt;
     jCol[10*(Time+1)+16*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[11*(Time+1)+16*(Time)+0+jt] = (Time+1)*4+jt;
     jCol[11*(Time+1)+16*(Time)+0+jt] = (Time+1)*0+jt+1;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[11*(Time+1)+17*(Time)+0+jt] = (Time+1)*4+jt;
     jCol[11*(Time+1)+17*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[12*(Time+1)+17*(Time)+0+jt] = (Time+1)*4+jt;
     jCol[12*(Time+1)+17*(Time)+0+jt] = (Time+1)*1+jt+1;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[12*(Time+1)+18*(Time)+0+jt] = (Time+1)*4+jt;
     jCol[12*(Time+1)+18*(Time)+0+jt] = (Time+1)*2+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[13*(Time+1)+18*(Time)+0+jt] = (Time+1)*4+jt;
     jCol[13*(Time+1)+18*(Time)+0+jt] = (Time+1)*2+jt+1;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[13*(Time+1)+19*(Time)+0+jt] = (Time+1)*4+jt;
     jCol[13*(Time+1)+19*(Time)+0+jt] = (Time+1)*3+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[14*(Time+1)+19*(Time)+0+jt] = (Time+1)*4+jt;
     jCol[14*(Time+1)+19*(Time)+0+jt] = (Time+1)*3+jt+1;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[14*(Time+1)+20*(Time)+0+jt] = (Time+1)*4+jt;
     jCol[14*(Time+1)+20*(Time)+0+jt] = (Time+1)*4+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+20*(Time)+0+jt] = (Time+1)*4+jt+1;
     jCol[15*(Time+1)+20*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+21*(Time)+0+jt] = (Time+1)*4+jt+1;
     jCol[15*(Time+1)+21*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+22*(Time)+0+jt] = (Time+1)*4+jt+1;
     jCol[15*(Time+1)+22*(Time)+0+jt] = (Time+1)*2+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+23*(Time)+0+jt] = (Time+1)*4+jt+1;
     jCol[15*(Time+1)+23*(Time)+0+jt] = (Time+1)*3+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+24*(Time)+0+jt] = (Time+1)*4+jt+1;
     jCol[15*(Time+1)+24*(Time)+0+jt] = (Time+1)*4+jt;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[15*(Time+1)+25*(Time)+0+jt] = (Time+1)*5+jt;
     jCol[15*(Time+1)+25*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[16*(Time+1)+25*(Time)+0+jt] = (Time+1)*5+jt;
     jCol[16*(Time+1)+25*(Time)+0+jt] = (Time+1)*0+jt+1;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[16*(Time+1)+26*(Time)+0+jt] = (Time+1)*5+jt;
     jCol[16*(Time+1)+26*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[17*(Time+1)+26*(Time)+0+jt] = (Time+1)*5+jt;
     jCol[17*(Time+1)+26*(Time)+0+jt] = (Time+1)*1+jt+1;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[17*(Time+1)+27*(Time)+0+jt] = (Time+1)*5+jt;
     jCol[17*(Time+1)+27*(Time)+0+jt] = (Time+1)*5+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[18*(Time+1)+27*(Time)+0+jt] = (Time+1)*5+jt+1;
     jCol[18*(Time+1)+27*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[18*(Time+1)+28*(Time)+0+jt] = (Time+1)*5+jt+1;
     jCol[18*(Time+1)+28*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[18*(Time+1)+29*(Time)+0+jt] = (Time+1)*5+jt+1;
     jCol[18*(Time+1)+29*(Time)+0+jt] = (Time+1)*5+jt;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[18*(Time+1)+30*(Time)+0+jt] = (Time+1)*6+jt;
     jCol[18*(Time+1)+30*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[19*(Time+1)+30*(Time)+0+jt] = (Time+1)*6+jt;
     jCol[19*(Time+1)+30*(Time)+0+jt] = (Time+1)*0+jt+1;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[19*(Time+1)+31*(Time)+0+jt] = (Time+1)*6+jt;
     jCol[19*(Time+1)+31*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[20*(Time+1)+31*(Time)+0+jt] = (Time+1)*6+jt;
     jCol[20*(Time+1)+31*(Time)+0+jt] = (Time+1)*1+jt+1;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[20*(Time+1)+32*(Time)+0+jt] = (Time+1)*6+jt;
     jCol[20*(Time+1)+32*(Time)+0+jt] = (Time+1)*5+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[21*(Time+1)+32*(Time)+0+jt] = (Time+1)*6+jt;
     jCol[21*(Time+1)+32*(Time)+0+jt] = (Time+1)*5+jt+1;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[21*(Time+1)+33*(Time)+0+jt] = (Time+1)*6+jt;
     jCol[21*(Time+1)+33*(Time)+0+jt] = (Time+1)*6+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+33*(Time)+0+jt] = (Time+1)*6+jt+1;
     jCol[22*(Time+1)+33*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+34*(Time)+0+jt] = (Time+1)*6+jt+1;
     jCol[22*(Time+1)+34*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+35*(Time)+0+jt] = (Time+1)*6+jt+1;
     jCol[22*(Time+1)+35*(Time)+0+jt] = (Time+1)*5+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+36*(Time)+0+jt] = (Time+1)*6+jt+1;
     jCol[22*(Time+1)+36*(Time)+0+jt] = (Time+1)*6+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+37*(Time)+0+jt] = (Time+1)*7+Time*0+jt;
     jCol[22*(Time+1)+37*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+38*(Time)+0+jt] = (Time+1)*7+Time*0+jt;
     jCol[22*(Time+1)+38*(Time)+0+jt] = (Time+1)*0+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+39*(Time)+0+jt] = (Time+1)*7+Time*0+jt;
     jCol[22*(Time+1)+39*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+40*(Time)+0+jt] = (Time+1)*7+Time*0+jt;
     jCol[22*(Time+1)+40*(Time)+0+jt] = (Time+1)*1+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+41*(Time)+0+jt] = (Time+1)*7+Time*0+jt;
     jCol[22*(Time+1)+41*(Time)+0+jt] = (Time+1)*2+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+42*(Time)+0+jt] = (Time+1)*7+Time*0+jt;
     jCol[22*(Time+1)+42*(Time)+0+jt] = (Time+1)*2+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+43*(Time)+0+jt] = (Time+1)*7+Time*0+jt;
     jCol[22*(Time+1)+43*(Time)+0+jt] = (Time+1)*3+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+44*(Time)+0+jt] = (Time+1)*7+Time*0+jt;
     jCol[22*(Time+1)+44*(Time)+0+jt] = (Time+1)*3+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+45*(Time)+0+jt] = (Time+1)*7+Time*0+jt;
     jCol[22*(Time+1)+45*(Time)+0+jt] = (Time+1)*4+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+46*(Time)+0+jt] = (Time+1)*7+Time*0+jt;
     jCol[22*(Time+1)+46*(Time)+0+jt] = (Time+1)*4+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+47*(Time)+0+jt] = (Time+1)*7+Time*0+jt;
     jCol[22*(Time+1)+47*(Time)+0+jt] = (Time+1)*5+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+48*(Time)+0+jt] = (Time+1)*7+Time*0+jt;
     jCol[22*(Time+1)+48*(Time)+0+jt] = (Time+1)*5+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+49*(Time)+0+jt] = (Time+1)*7+Time*0+jt;
     jCol[22*(Time+1)+49*(Time)+0+jt] = (Time+1)*6+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+50*(Time)+0+jt] = (Time+1)*7+Time*0+jt;
     jCol[22*(Time+1)+50*(Time)+0+jt] = (Time+1)*6+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+51*(Time)+0+jt] = (Time+1)*7+Time*0+jt;
     jCol[22*(Time+1)+51*(Time)+0+jt] = (Time+1)*7+Time*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+52*(Time)+0+jt] = (Time+1)*7+Time*1+jt;
     jCol[22*(Time+1)+52*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+53*(Time)+0+jt] = (Time+1)*7+Time*1+jt;
     jCol[22*(Time+1)+53*(Time)+0+jt] = (Time+1)*0+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+54*(Time)+0+jt] = (Time+1)*7+Time*1+jt;
     jCol[22*(Time+1)+54*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+55*(Time)+0+jt] = (Time+1)*7+Time*1+jt;
     jCol[22*(Time+1)+55*(Time)+0+jt] = (Time+1)*1+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+56*(Time)+0+jt] = (Time+1)*7+Time*1+jt;
     jCol[22*(Time+1)+56*(Time)+0+jt] = (Time+1)*2+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+57*(Time)+0+jt] = (Time+1)*7+Time*1+jt;
     jCol[22*(Time+1)+57*(Time)+0+jt] = (Time+1)*2+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+58*(Time)+0+jt] = (Time+1)*7+Time*1+jt;
     jCol[22*(Time+1)+58*(Time)+0+jt] = (Time+1)*3+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+59*(Time)+0+jt] = (Time+1)*7+Time*1+jt;
     jCol[22*(Time+1)+59*(Time)+0+jt] = (Time+1)*3+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+60*(Time)+0+jt] = (Time+1)*7+Time*1+jt;
     jCol[22*(Time+1)+60*(Time)+0+jt] = (Time+1)*4+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+61*(Time)+0+jt] = (Time+1)*7+Time*1+jt;
     jCol[22*(Time+1)+61*(Time)+0+jt] = (Time+1)*4+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+62*(Time)+0+jt] = (Time+1)*7+Time*1+jt;
     jCol[22*(Time+1)+62*(Time)+0+jt] = (Time+1)*5+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+63*(Time)+0+jt] = (Time+1)*7+Time*1+jt;
     jCol[22*(Time+1)+63*(Time)+0+jt] = (Time+1)*5+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+64*(Time)+0+jt] = (Time+1)*7+Time*1+jt;
     jCol[22*(Time+1)+64*(Time)+0+jt] = (Time+1)*6+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+65*(Time)+0+jt] = (Time+1)*7+Time*1+jt;
     jCol[22*(Time+1)+65*(Time)+0+jt] = (Time+1)*6+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+66*(Time)+0+jt] = (Time+1)*7+Time*1+jt;
     jCol[22*(Time+1)+66*(Time)+0+jt] = (Time+1)*7+Time*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+67*(Time)+0+jt] = (Time+1)*7+Time*1+jt;
     jCol[22*(Time+1)+67*(Time)+0+jt] = (Time+1)*7+Time*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+68*(Time)+0+jt] = (Time+1)*7+Time*2+jt;
     jCol[22*(Time+1)+68*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+69*(Time)+0+jt] = (Time+1)*7+Time*2+jt;
     jCol[22*(Time+1)+69*(Time)+0+jt] = (Time+1)*0+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+70*(Time)+0+jt] = (Time+1)*7+Time*2+jt;
     jCol[22*(Time+1)+70*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+71*(Time)+0+jt] = (Time+1)*7+Time*2+jt;
     jCol[22*(Time+1)+71*(Time)+0+jt] = (Time+1)*1+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+72*(Time)+0+jt] = (Time+1)*7+Time*2+jt;
     jCol[22*(Time+1)+72*(Time)+0+jt] = (Time+1)*2+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+73*(Time)+0+jt] = (Time+1)*7+Time*2+jt;
     jCol[22*(Time+1)+73*(Time)+0+jt] = (Time+1)*2+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+74*(Time)+0+jt] = (Time+1)*7+Time*2+jt;
     jCol[22*(Time+1)+74*(Time)+0+jt] = (Time+1)*3+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+75*(Time)+0+jt] = (Time+1)*7+Time*2+jt;
     jCol[22*(Time+1)+75*(Time)+0+jt] = (Time+1)*3+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+76*(Time)+0+jt] = (Time+1)*7+Time*2+jt;
     jCol[22*(Time+1)+76*(Time)+0+jt] = (Time+1)*4+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+77*(Time)+0+jt] = (Time+1)*7+Time*2+jt;
     jCol[22*(Time+1)+77*(Time)+0+jt] = (Time+1)*4+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+78*(Time)+0+jt] = (Time+1)*7+Time*2+jt;
     jCol[22*(Time+1)+78*(Time)+0+jt] = (Time+1)*7+Time*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+79*(Time)+0+jt] = (Time+1)*7+Time*2+jt;
     jCol[22*(Time+1)+79*(Time)+0+jt] = (Time+1)*7+Time*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+80*(Time)+0+jt] = (Time+1)*7+Time*2+jt;
     jCol[22*(Time+1)+80*(Time)+0+jt] = (Time+1)*7+Time*2+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+81*(Time)+0+jt] = (Time+1)*7+Time*3+jt;
     jCol[22*(Time+1)+81*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+82*(Time)+0+jt] = (Time+1)*7+Time*3+jt;
     jCol[22*(Time+1)+82*(Time)+0+jt] = (Time+1)*0+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+83*(Time)+0+jt] = (Time+1)*7+Time*3+jt;
     jCol[22*(Time+1)+83*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+84*(Time)+0+jt] = (Time+1)*7+Time*3+jt;
     jCol[22*(Time+1)+84*(Time)+0+jt] = (Time+1)*1+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+85*(Time)+0+jt] = (Time+1)*7+Time*3+jt;
     jCol[22*(Time+1)+85*(Time)+0+jt] = (Time+1)*2+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+86*(Time)+0+jt] = (Time+1)*7+Time*3+jt;
     jCol[22*(Time+1)+86*(Time)+0+jt] = (Time+1)*2+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+87*(Time)+0+jt] = (Time+1)*7+Time*3+jt;
     jCol[22*(Time+1)+87*(Time)+0+jt] = (Time+1)*3+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+88*(Time)+0+jt] = (Time+1)*7+Time*3+jt;
     jCol[22*(Time+1)+88*(Time)+0+jt] = (Time+1)*3+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+89*(Time)+0+jt] = (Time+1)*7+Time*3+jt;
     jCol[22*(Time+1)+89*(Time)+0+jt] = (Time+1)*4+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+90*(Time)+0+jt] = (Time+1)*7+Time*3+jt;
     jCol[22*(Time+1)+90*(Time)+0+jt] = (Time+1)*4+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+91*(Time)+0+jt] = (Time+1)*7+Time*3+jt;
     jCol[22*(Time+1)+91*(Time)+0+jt] = (Time+1)*7+Time*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+92*(Time)+0+jt] = (Time+1)*7+Time*3+jt;
     jCol[22*(Time+1)+92*(Time)+0+jt] = (Time+1)*7+Time*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+93*(Time)+0+jt] = (Time+1)*7+Time*3+jt;
     jCol[22*(Time+1)+93*(Time)+0+jt] = (Time+1)*7+Time*2+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+94*(Time)+0+jt] = (Time+1)*7+Time*3+jt;
     jCol[22*(Time+1)+94*(Time)+0+jt] = (Time+1)*7+Time*3+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+95*(Time)+0+jt] = (Time+1)*7+Time*4+jt;
     jCol[22*(Time+1)+95*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+96*(Time)+0+jt] = (Time+1)*7+Time*4+jt;
     jCol[22*(Time+1)+96*(Time)+0+jt] = (Time+1)*0+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+97*(Time)+0+jt] = (Time+1)*7+Time*4+jt;
     jCol[22*(Time+1)+97*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+98*(Time)+0+jt] = (Time+1)*7+Time*4+jt;
     jCol[22*(Time+1)+98*(Time)+0+jt] = (Time+1)*1+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+99*(Time)+0+jt] = (Time+1)*7+Time*4+jt;
     jCol[22*(Time+1)+99*(Time)+0+jt] = (Time+1)*2+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+100*(Time)+0+jt] = (Time+1)*7+Time*4+jt;
     jCol[22*(Time+1)+100*(Time)+0+jt] = (Time+1)*2+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+101*(Time)+0+jt] = (Time+1)*7+Time*4+jt;
     jCol[22*(Time+1)+101*(Time)+0+jt] = (Time+1)*3+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+102*(Time)+0+jt] = (Time+1)*7+Time*4+jt;
     jCol[22*(Time+1)+102*(Time)+0+jt] = (Time+1)*3+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+103*(Time)+0+jt] = (Time+1)*7+Time*4+jt;
     jCol[22*(Time+1)+103*(Time)+0+jt] = (Time+1)*4+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+104*(Time)+0+jt] = (Time+1)*7+Time*4+jt;
     jCol[22*(Time+1)+104*(Time)+0+jt] = (Time+1)*4+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+105*(Time)+0+jt] = (Time+1)*7+Time*4+jt;
     jCol[22*(Time+1)+105*(Time)+0+jt] = (Time+1)*7+Time*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+106*(Time)+0+jt] = (Time+1)*7+Time*4+jt;
     jCol[22*(Time+1)+106*(Time)+0+jt] = (Time+1)*7+Time*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+107*(Time)+0+jt] = (Time+1)*7+Time*4+jt;
     jCol[22*(Time+1)+107*(Time)+0+jt] = (Time+1)*7+Time*2+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+108*(Time)+0+jt] = (Time+1)*7+Time*4+jt;
     jCol[22*(Time+1)+108*(Time)+0+jt] = (Time+1)*7+Time*3+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+109*(Time)+0+jt] = (Time+1)*7+Time*4+jt;
     jCol[22*(Time+1)+109*(Time)+0+jt] = (Time+1)*7+Time*4+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+110*(Time)+0+jt] = (Time+1)*7+Time*5+jt;
     jCol[22*(Time+1)+110*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+111*(Time)+0+jt] = (Time+1)*7+Time*5+jt;
     jCol[22*(Time+1)+111*(Time)+0+jt] = (Time+1)*0+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+112*(Time)+0+jt] = (Time+1)*7+Time*5+jt;
     jCol[22*(Time+1)+112*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+113*(Time)+0+jt] = (Time+1)*7+Time*5+jt;
     jCol[22*(Time+1)+113*(Time)+0+jt] = (Time+1)*1+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+114*(Time)+0+jt] = (Time+1)*7+Time*5+jt;
     jCol[22*(Time+1)+114*(Time)+0+jt] = (Time+1)*5+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+115*(Time)+0+jt] = (Time+1)*7+Time*5+jt;
     jCol[22*(Time+1)+115*(Time)+0+jt] = (Time+1)*5+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+116*(Time)+0+jt] = (Time+1)*7+Time*5+jt;
     jCol[22*(Time+1)+116*(Time)+0+jt] = (Time+1)*6+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+117*(Time)+0+jt] = (Time+1)*7+Time*5+jt;
     jCol[22*(Time+1)+117*(Time)+0+jt] = (Time+1)*6+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+118*(Time)+0+jt] = (Time+1)*7+Time*5+jt;
     jCol[22*(Time+1)+118*(Time)+0+jt] = (Time+1)*7+Time*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+119*(Time)+0+jt] = (Time+1)*7+Time*5+jt;
     jCol[22*(Time+1)+119*(Time)+0+jt] = (Time+1)*7+Time*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+120*(Time)+0+jt] = (Time+1)*7+Time*5+jt;
     jCol[22*(Time+1)+120*(Time)+0+jt] = (Time+1)*7+Time*5+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+121*(Time)+0+jt] = (Time+1)*7+Time*6+jt;
     jCol[22*(Time+1)+121*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+122*(Time)+0+jt] = (Time+1)*7+Time*6+jt;
     jCol[22*(Time+1)+122*(Time)+0+jt] = (Time+1)*0+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+123*(Time)+0+jt] = (Time+1)*7+Time*6+jt;
     jCol[22*(Time+1)+123*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+124*(Time)+0+jt] = (Time+1)*7+Time*6+jt;
     jCol[22*(Time+1)+124*(Time)+0+jt] = (Time+1)*1+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+125*(Time)+0+jt] = (Time+1)*7+Time*6+jt;
     jCol[22*(Time+1)+125*(Time)+0+jt] = (Time+1)*5+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+126*(Time)+0+jt] = (Time+1)*7+Time*6+jt;
     jCol[22*(Time+1)+126*(Time)+0+jt] = (Time+1)*5+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+127*(Time)+0+jt] = (Time+1)*7+Time*6+jt;
     jCol[22*(Time+1)+127*(Time)+0+jt] = (Time+1)*6+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+128*(Time)+0+jt] = (Time+1)*7+Time*6+jt;
     jCol[22*(Time+1)+128*(Time)+0+jt] = (Time+1)*6+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+129*(Time)+0+jt] = (Time+1)*7+Time*6+jt;
     jCol[22*(Time+1)+129*(Time)+0+jt] = (Time+1)*7+Time*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+130*(Time)+0+jt] = (Time+1)*7+Time*6+jt;
     jCol[22*(Time+1)+130*(Time)+0+jt] = (Time+1)*7+Time*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+131*(Time)+0+jt] = (Time+1)*7+Time*6+jt;
     jCol[22*(Time+1)+131*(Time)+0+jt] = (Time+1)*7+Time*5+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+132*(Time)+0+jt] = (Time+1)*7+Time*6+jt;
     jCol[22*(Time+1)+132*(Time)+0+jt] = (Time+1)*7+Time*6+jt;
   }
}

else {
  // return the values.  This is a symmetric matrix, fill the lower left
  // triangle only
  // initialize the values array
  // Point to the initial starting spot for the Hessian elements

  for(Index jt=0;jt<22*(Time+1)+133*Time+0;jt++) values[jt] = 0.; // Initialize matrix

   // fill the objective portion

  for(Index jt=0;jt<Time;jt++) {

     for(Index i=0;i<nY;i++) {
        Xval[i] = x[jt + i*(Time+1)];
        Xvalp1[i] = x[jt + i*(Time+1) + 1];
        Xval2[i] = x[(Time+1)*(nY+2*nU) + jt + i*(Time)];
     } //end for loop

     for(Index i=0;i<nU;i++) {
        K11val[i] = x[jt + nY*(Time+1) + 2*i*(Time+1)];
        K11valp1[i] = x[jt + nY*(Time+1) + 2*i*(Time+1) + 1];
        K11val2[i] = x[(Time+1)*(nY+2*nU) + (nY+2*i)*Time + jt];
        dK11val[i] = x[jt + (nY+2*i+1)*(Time+1)];
        dK11valp1[i] = x[jt + (nY+2*i+1)*(Time+1)+1];
        dK11val2[i] = x[(Time+1)*(nY+2*nU) + (nY+2*i+1)*Time + jt];
     } //end for loop

     Xdval[0] = VDATA0[2*jt];
     Xdval2[0] = VDATA0[2*jt+1];
     Xdvalp1[0] = VDATA0[2*jt+2];
     Ival[0] = Iinj[2*jt];
     Ival2[0] = Iinj[2*jt+1];
     Ivalp1[0] = Iinj[2*jt+2];

     for(Index i=0;i<nP;i++) {
        Pval[i] = x[(2*Time+1)*(nY+2*nU)+i];
     } //end for loop

   values[0*(Time+1)+0*(Time)+0+jt] += obj_factor*1*(2)/(2*Time+1);

   values[22*(Time+1)+51*(Time)+0+jt] += obj_factor*1*(2)/(2*Time+1);

   values[0*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[0][3]*((0.125*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 0.5)*(0.25*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 1.0) + (0.166666666666667*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 1)*(0.333333333333333*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 2))/(2*Time+1);

   values[1*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[0][3]*((0.125*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 0.5)*(0.25*hstep*(70.0*Xvalp1[4]*pow(Xvalp1[3], 3) + 8*pow(Xvalp1[2], 4) + 1) + 1.0) + (0.166666666666667*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 1)*(0.333333333333333*hstep*(-70.0*Xvalp1[4]*pow(Xvalp1[3], 3) - 8*pow(Xvalp1[2], 4) - 1) - 2))/(2*Time+1);

   values[0*(Time+1)+0*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*((0.166666666666667*hstep*(-70.0*Xvalp1[4]*pow(Xvalp1[3], 3) - 8*pow(Xvalp1[2], 4) - 1) - 1)*(0.333333333333333*hstep*(-70.0*Xvalp1[4]*pow(Xvalp1[3], 3) - 8*pow(Xvalp1[2], 4) - 1) - 2) + (0.125*hstep*(70.0*Xvalp1[4]*pow(Xvalp1[3], 3) + 8*pow(Xvalp1[2], 4) + 1) + 0.5)*(0.25*hstep*(70.0*Xvalp1[4]*pow(Xvalp1[3], 3) + 8*pow(Xvalp1[2], 4) + 1) + 1.0))/(2*Time+1);

   values[1*(Time+1)+1*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.225*hstep*(0.125*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 0.5) + 0.3*hstep*(0.166666666666667*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 1))/(2*Time+1);

   values[2*(Time+1)+1*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.3*hstep*(0.166666666666667*hstep*(-70.0*Xvalp1[4]*pow(Xvalp1[3], 3) - 8*pow(Xvalp1[2], 4) - 1) - 1) + 0.225*hstep*(0.125*hstep*(70.0*Xvalp1[4]*pow(Xvalp1[3], 3) + 8*pow(Xvalp1[2], 4) + 1) + 0.5))/(2*Time+1);

   values[2*(Time+1)+2*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0703125*pow(hstep, 2))/(2*Time+1);

   values[3*(Time+1)+2*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.225*hstep*(0.125*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 0.5) + 0.3*hstep*(0.166666666666667*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 1))/(2*Time+1);

   values[1*(Time+1)+1*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(0.3*hstep*(0.166666666666667*hstep*(-70.0*Xvalp1[4]*pow(Xvalp1[3], 3) - 8*pow(Xvalp1[2], 4) - 1) - 1) - 0.225*hstep*(0.125*hstep*(70.0*Xvalp1[4]*pow(Xvalp1[3], 3) + 8*pow(Xvalp1[2], 4) + 1) + 0.5))/(2*Time+1);

   values[3*(Time+1)+3*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0196875*pow(hstep, 2))/(2*Time+1);

   values[2*(Time+1)+2*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(0.0703125*pow(hstep, 2))/(2*Time+1);

   values[3*(Time+1)+4*(Time)+0+jt] += obj_factor*bounds[0][3]*(8.0*hstep*pow(Xval[2], 3)*(-Xval[0] - 90)*(0.125*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 0.5) + 10.6666666666667*hstep*pow(Xval[2], 3)*(-Xval[0] - 90)*(0.166666666666667*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 1) - 10.6666666666667*hstep*pow(Xval[2], 3)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)) - 8.0*hstep*pow(Xval[2], 3)*(0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0] + 0.125*hstep*((1.0L/10.0L)*Ival[0] - 1.0L/10.0L*Ivalp1[0] + (9.0L/10.0L)*Xval[1] - 9.0L/10.0L*Xvalp1[1] - Xval[0] + Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) - 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) - 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90))))/(2*Time+1);

   values[4*(Time+1)+4*(Time)+0+jt] += obj_factor*bounds[0][3]*(10.6666666666667*hstep*pow(Xval[2], 3)*(-Xval[0] - 90)*(0.166666666666667*hstep*(-70.0*Xvalp1[4]*pow(Xvalp1[3], 3) - 8*pow(Xvalp1[2], 4) - 1) - 1) + 8.0*hstep*pow(Xval[2], 3)*(-Xval[0] - 90)*(0.125*hstep*(70.0*Xvalp1[4]*pow(Xvalp1[3], 3) + 8*pow(Xvalp1[2], 4) + 1) + 0.5))/(2*Time+1);

   values[4*(Time+1)+5*(Time)+0+jt] += obj_factor*bounds[0][3]*(2.5*pow(hstep, 2)*pow(Xval[2], 3)*(-Xval[0] - 90))/(2*Time+1);

   values[5*(Time+1)+5*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.7*pow(hstep, 2)*pow(Xval[2], 3)*(-Xval[0] - 90))/(2*Time+1);

   values[5*(Time+1)+6*(Time)+0+jt] += obj_factor*bounds[0][3]*(88.8888888888889*pow(hstep, 2)*pow(Xval[2], 6)*pow(-Xval[0] - 90, 2) + 32.0*hstep*pow(Xval[2], 2)*(-Xval[0] - 90)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)) + 24.0*hstep*pow(Xval[2], 2)*(-Xval[0] - 90)*(0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0] + 0.125*hstep*((1.0L/10.0L)*Ival[0] - 1.0L/10.0L*Ivalp1[0] + (9.0L/10.0L)*Xval[1] - 9.0L/10.0L*Xvalp1[1] - Xval[0] + Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) - 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) - 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90))))/(2*Time+1);

   values[6*(Time+1)+6*(Time)+0+jt] += obj_factor*bounds[0][3]*(-8.0*hstep*pow(Xvalp1[2], 3)*(-Xvalp1[0] - 90)*(0.125*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 0.5) + 10.6666666666667*hstep*pow(Xvalp1[2], 3)*(-Xvalp1[0] - 90)*(0.166666666666667*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 1))/(2*Time+1);

   values[3*(Time+1)+4*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(10.6666666666667*hstep*pow(Xvalp1[2], 3)*(-Xvalp1[0] - 90)*(0.166666666666667*hstep*(-70.0*Xvalp1[4]*pow(Xvalp1[3], 3) - 8*pow(Xvalp1[2], 4) - 1) - 1) - 8.0*hstep*pow(Xvalp1[2], 3)*(-Xvalp1[0] - 90)*(0.125*hstep*(70.0*Xvalp1[4]*pow(Xvalp1[3], 3) + 8*pow(Xvalp1[2], 4) + 1) + 0.5) - 10.6666666666667*hstep*pow(Xvalp1[2], 3)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)) + 8.0*hstep*pow(Xvalp1[2], 3)*(0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0] + 0.125*hstep*((1.0L/10.0L)*Ival[0] - 1.0L/10.0L*Ivalp1[0] + (9.0L/10.0L)*Xval[1] - 9.0L/10.0L*Xvalp1[1] - Xval[0] + Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) - 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) - 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90))))/(2*Time+1);

   values[6*(Time+1)+7*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.7*pow(hstep, 2)*pow(Xvalp1[2], 3)*(-Xvalp1[0] - 90))/(2*Time+1);

   values[4*(Time+1)+5*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(2.5*pow(hstep, 2)*pow(Xvalp1[2], 3)*(-Xvalp1[0] - 90))/(2*Time+1);

   values[6*(Time+1)+8*(Time)+0+jt] += obj_factor*bounds[0][3]*(24.8888888888889*pow(hstep, 2)*pow(Xval[2], 3)*pow(Xvalp1[2], 3)*(-Xval[0] - 90)*(-Xvalp1[0] - 90))/(2*Time+1);

   values[5*(Time+1)+6*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(88.8888888888889*pow(hstep, 2)*pow(Xvalp1[2], 6)*pow(-Xvalp1[0] - 90, 2) + 32.0*hstep*pow(Xvalp1[2], 2)*(-Xvalp1[0] - 90)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)) - 24.0*hstep*pow(Xvalp1[2], 2)*(-Xvalp1[0] - 90)*(0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0] + 0.125*hstep*((1.0L/10.0L)*Ival[0] - 1.0L/10.0L*Ivalp1[0] + (9.0L/10.0L)*Xval[1] - 9.0L/10.0L*Xvalp1[1] - Xval[0] + Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) - 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) - 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90))))/(2*Time+1);

   values[6*(Time+1)+9*(Time)+0+jt] += obj_factor*bounds[0][3]*(52.5*Xval[4]*hstep*pow(Xval[3], 2)*(-Xval[0] + 55)*(0.125*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 0.5) + 70.0*Xval[4]*hstep*pow(Xval[3], 2)*(-Xval[0] + 55)*(0.166666666666667*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 1) - 70.0*Xval[4]*hstep*pow(Xval[3], 2)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)) - 52.5*Xval[4]*hstep*pow(Xval[3], 2)*(0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0] + 0.125*hstep*((1.0L/10.0L)*Ival[0] - 1.0L/10.0L*Ivalp1[0] + (9.0L/10.0L)*Xval[1] - 9.0L/10.0L*Xvalp1[1] - Xval[0] + Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) - 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) - 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90))))/(2*Time+1);

   values[7*(Time+1)+9*(Time)+0+jt] += obj_factor*bounds[0][3]*(70.0*Xval[4]*hstep*pow(Xval[3], 2)*(-Xval[0] + 55)*(0.166666666666667*hstep*(-70.0*Xvalp1[4]*pow(Xvalp1[3], 3) - 8*pow(Xvalp1[2], 4) - 1) - 1) + 52.5*Xval[4]*hstep*pow(Xval[3], 2)*(-Xval[0] + 55)*(0.125*hstep*(70.0*Xvalp1[4]*pow(Xvalp1[3], 3) + 8*pow(Xvalp1[2], 4) + 1) + 0.5))/(2*Time+1);

   values[7*(Time+1)+10*(Time)+0+jt] += obj_factor*bounds[0][3]*(16.40625*Xval[4]*pow(hstep, 2)*pow(Xval[3], 2)*(-Xval[0] + 55))/(2*Time+1);

   values[8*(Time+1)+10*(Time)+0+jt] += obj_factor*bounds[0][3]*(4.59375*Xval[4]*pow(hstep, 2)*pow(Xval[3], 2)*(-Xval[0] + 55))/(2*Time+1);

   values[8*(Time+1)+11*(Time)+0+jt] += obj_factor*bounds[0][3]*(583.333333333333*Xval[4]*pow(hstep, 2)*pow(Xval[3], 2)*pow(Xval[2], 3)*(-Xval[0] - 90)*(-Xval[0] + 55))/(2*Time+1);

   values[9*(Time+1)+11*(Time)+0+jt] += obj_factor*bounds[0][3]*(163.333333333333*Xval[4]*pow(hstep, 2)*pow(Xval[3], 2)*pow(Xvalp1[2], 3)*(-Xval[0] + 55)*(-Xvalp1[0] - 90))/(2*Time+1);

   values[9*(Time+1)+12*(Time)+0+jt] += obj_factor*bounds[0][3]*(3828.125*pow(Xval[4], 2)*pow(hstep, 2)*pow(Xval[3], 4)*pow(-Xval[0] + 55, 2) + 140.0*Xval[4]*hstep*Xval[3]*(-Xval[0] + 55)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)) + 105.0*Xval[4]*hstep*Xval[3]*(-Xval[0] + 55)*(0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0] + 0.125*hstep*((1.0L/10.0L)*Ival[0] - 1.0L/10.0L*Ivalp1[0] + (9.0L/10.0L)*Xval[1] - 9.0L/10.0L*Xvalp1[1] - Xval[0] + Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) - 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) - 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90))))/(2*Time+1);

   values[10*(Time+1)+12*(Time)+0+jt] += obj_factor*bounds[0][3]*(-52.5*Xvalp1[4]*hstep*pow(Xvalp1[3], 2)*(-Xvalp1[0] + 55)*(0.125*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 0.5) + 70.0*Xvalp1[4]*hstep*pow(Xvalp1[3], 2)*(-Xvalp1[0] + 55)*(0.166666666666667*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 1))/(2*Time+1);

   values[6*(Time+1)+9*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(70.0*Xvalp1[4]*hstep*pow(Xvalp1[3], 2)*(-Xvalp1[0] + 55)*(0.166666666666667*hstep*(-70.0*Xvalp1[4]*pow(Xvalp1[3], 3) - 8*pow(Xvalp1[2], 4) - 1) - 1) - 52.5*Xvalp1[4]*hstep*pow(Xvalp1[3], 2)*(-Xvalp1[0] + 55)*(0.125*hstep*(70.0*Xvalp1[4]*pow(Xvalp1[3], 3) + 8*pow(Xvalp1[2], 4) + 1) + 0.5) - 70.0*Xvalp1[4]*hstep*pow(Xvalp1[3], 2)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)) + 52.5*Xvalp1[4]*hstep*pow(Xvalp1[3], 2)*(0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0] + 0.125*hstep*((1.0L/10.0L)*Ival[0] - 1.0L/10.0L*Ivalp1[0] + (9.0L/10.0L)*Xval[1] - 9.0L/10.0L*Xvalp1[1] - Xval[0] + Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) - 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) - 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90))))/(2*Time+1);

   values[10*(Time+1)+13*(Time)+0+jt] += obj_factor*bounds[0][3]*(4.59375*Xvalp1[4]*pow(hstep, 2)*pow(Xvalp1[3], 2)*(-Xvalp1[0] + 55))/(2*Time+1);

   values[7*(Time+1)+10*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(16.40625*Xvalp1[4]*pow(hstep, 2)*pow(Xvalp1[3], 2)*(-Xvalp1[0] + 55))/(2*Time+1);

   values[10*(Time+1)+14*(Time)+0+jt] += obj_factor*bounds[0][3]*(163.333333333333*Xvalp1[4]*pow(hstep, 2)*pow(Xvalp1[3], 2)*pow(Xval[2], 3)*(-Xval[0] - 90)*(-Xvalp1[0] + 55))/(2*Time+1);

   values[8*(Time+1)+11*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(583.333333333333*Xvalp1[4]*pow(hstep, 2)*pow(Xvalp1[3], 2)*pow(Xvalp1[2], 3)*(-Xvalp1[0] - 90)*(-Xvalp1[0] + 55))/(2*Time+1);

   values[10*(Time+1)+15*(Time)+0+jt] += obj_factor*bounds[0][3]*(1071.875*Xval[4]*Xvalp1[4]*pow(hstep, 2)*pow(Xval[3], 2)*pow(Xvalp1[3], 2)*(-Xval[0] + 55)*(-Xvalp1[0] + 55))/(2*Time+1);

   values[9*(Time+1)+12*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(3828.125*pow(Xvalp1[4], 2)*pow(hstep, 2)*pow(Xvalp1[3], 4)*pow(-Xvalp1[0] + 55, 2) + 140.0*Xvalp1[4]*hstep*Xvalp1[3]*(-Xvalp1[0] + 55)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)) - 105.0*Xvalp1[4]*hstep*Xvalp1[3]*(-Xvalp1[0] + 55)*(0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0] + 0.125*hstep*((1.0L/10.0L)*Ival[0] - 1.0L/10.0L*Ivalp1[0] + (9.0L/10.0L)*Xval[1] - 9.0L/10.0L*Xvalp1[1] - Xval[0] + Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) - 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) - 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90))))/(2*Time+1);

   values[10*(Time+1)+16*(Time)+0+jt] += obj_factor*bounds[0][3]*(17.5*hstep*pow(Xval[3], 3)*(-Xval[0] + 55)*(0.125*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 0.5) + 23.3333333333333*hstep*pow(Xval[3], 3)*(-Xval[0] + 55)*(0.166666666666667*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 1) - 23.3333333333333*hstep*pow(Xval[3], 3)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)) - 17.5*hstep*pow(Xval[3], 3)*(0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0] + 0.125*hstep*((1.0L/10.0L)*Ival[0] - 1.0L/10.0L*Ivalp1[0] + (9.0L/10.0L)*Xval[1] - 9.0L/10.0L*Xvalp1[1] - Xval[0] + Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) - 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) - 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90))))/(2*Time+1);

   values[11*(Time+1)+16*(Time)+0+jt] += obj_factor*bounds[0][3]*(23.3333333333333*hstep*pow(Xval[3], 3)*(-Xval[0] + 55)*(0.166666666666667*hstep*(-70.0*Xvalp1[4]*pow(Xvalp1[3], 3) - 8*pow(Xvalp1[2], 4) - 1) - 1) + 17.5*hstep*pow(Xval[3], 3)*(-Xval[0] + 55)*(0.125*hstep*(70.0*Xvalp1[4]*pow(Xvalp1[3], 3) + 8*pow(Xvalp1[2], 4) + 1) + 0.5))/(2*Time+1);

   values[11*(Time+1)+17*(Time)+0+jt] += obj_factor*bounds[0][3]*(5.46875*pow(hstep, 2)*pow(Xval[3], 3)*(-Xval[0] + 55))/(2*Time+1);

   values[12*(Time+1)+17*(Time)+0+jt] += obj_factor*bounds[0][3]*(1.53125*pow(hstep, 2)*pow(Xval[3], 3)*(-Xval[0] + 55))/(2*Time+1);

   values[12*(Time+1)+18*(Time)+0+jt] += obj_factor*bounds[0][3]*(194.444444444444*pow(hstep, 2)*pow(Xval[3], 3)*pow(Xval[2], 3)*(-Xval[0] - 90)*(-Xval[0] + 55))/(2*Time+1);

   values[13*(Time+1)+18*(Time)+0+jt] += obj_factor*bounds[0][3]*(54.4444444444444*pow(hstep, 2)*pow(Xval[3], 3)*pow(Xvalp1[2], 3)*(-Xval[0] + 55)*(-Xvalp1[0] - 90))/(2*Time+1);

   values[13*(Time+1)+19*(Time)+0+jt] += obj_factor*bounds[0][3]*(1276.04166666667*Xval[4]*pow(hstep, 2)*pow(Xval[3], 5)*pow(-Xval[0] + 55, 2) + 70.0*hstep*pow(Xval[3], 2)*(-Xval[0] + 55)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)) + 52.5*hstep*pow(Xval[3], 2)*(-Xval[0] + 55)*(0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0] + 0.125*hstep*((1.0L/10.0L)*Ival[0] - 1.0L/10.0L*Ivalp1[0] + (9.0L/10.0L)*Xval[1] - 9.0L/10.0L*Xvalp1[1] - Xval[0] + Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) - 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) - 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90))))/(2*Time+1);

   values[14*(Time+1)+19*(Time)+0+jt] += obj_factor*bounds[0][3]*(357.291666666667*Xvalp1[4]*pow(hstep, 2)*pow(Xval[3], 3)*pow(Xvalp1[3], 2)*(-Xval[0] + 55)*(-Xvalp1[0] + 55))/(2*Time+1);

   values[14*(Time+1)+20*(Time)+0+jt] += obj_factor*bounds[0][3]*(425.347222222222*pow(hstep, 2)*pow(Xval[3], 6)*pow(-Xval[0] + 55, 2))/(2*Time+1);

   values[15*(Time+1)+20*(Time)+0+jt] += obj_factor*bounds[0][3]*(-17.5*hstep*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55)*(0.125*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 0.5) + 23.3333333333333*hstep*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55)*(0.166666666666667*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 1))/(2*Time+1);

   values[10*(Time+1)+16*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(23.3333333333333*hstep*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55)*(0.166666666666667*hstep*(-70.0*Xvalp1[4]*pow(Xvalp1[3], 3) - 8*pow(Xvalp1[2], 4) - 1) - 1) - 17.5*hstep*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55)*(0.125*hstep*(70.0*Xvalp1[4]*pow(Xvalp1[3], 3) + 8*pow(Xvalp1[2], 4) + 1) + 0.5) - 23.3333333333333*hstep*pow(Xvalp1[3], 3)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)) + 17.5*hstep*pow(Xvalp1[3], 3)*(0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0] + 0.125*hstep*((1.0L/10.0L)*Ival[0] - 1.0L/10.0L*Ivalp1[0] + (9.0L/10.0L)*Xval[1] - 9.0L/10.0L*Xvalp1[1] - Xval[0] + Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) - 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) - 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90))))/(2*Time+1);

   values[15*(Time+1)+21*(Time)+0+jt] += obj_factor*bounds[0][3]*(1.53125*pow(hstep, 2)*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55))/(2*Time+1);

   values[11*(Time+1)+17*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(5.46875*pow(hstep, 2)*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55))/(2*Time+1);

   values[15*(Time+1)+22*(Time)+0+jt] += obj_factor*bounds[0][3]*(54.4444444444444*pow(hstep, 2)*pow(Xvalp1[3], 3)*pow(Xval[2], 3)*(-Xval[0] - 90)*(-Xvalp1[0] + 55))/(2*Time+1);

   values[12*(Time+1)+18*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(194.444444444444*pow(hstep, 2)*pow(Xvalp1[3], 3)*pow(Xvalp1[2], 3)*(-Xvalp1[0] - 90)*(-Xvalp1[0] + 55))/(2*Time+1);

   values[15*(Time+1)+23*(Time)+0+jt] += obj_factor*bounds[0][3]*(357.291666666667*Xval[4]*pow(hstep, 2)*pow(Xval[3], 2)*pow(Xvalp1[3], 3)*(-Xval[0] + 55)*(-Xvalp1[0] + 55))/(2*Time+1);

   values[13*(Time+1)+19*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(1276.04166666667*Xvalp1[4]*pow(hstep, 2)*pow(Xvalp1[3], 5)*pow(-Xvalp1[0] + 55, 2) + 70.0*hstep*pow(Xvalp1[3], 2)*(-Xvalp1[0] + 55)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)) - 52.5*hstep*pow(Xvalp1[3], 2)*(-Xvalp1[0] + 55)*(0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0] + 0.125*hstep*((1.0L/10.0L)*Ival[0] - 1.0L/10.0L*Ivalp1[0] + (9.0L/10.0L)*Xval[1] - 9.0L/10.0L*Xvalp1[1] - Xval[0] + Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) - 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) - 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90))))/(2*Time+1);

   values[15*(Time+1)+24*(Time)+0+jt] += obj_factor*bounds[0][3]*(119.097222222222*pow(hstep, 2)*pow(Xval[3], 3)*pow(Xvalp1[3], 3)*(-Xval[0] + 55)*(-Xvalp1[0] + 55))/(2*Time+1);

   values[14*(Time+1)+20*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(425.347222222222*pow(hstep, 2)*pow(Xvalp1[3], 6)*pow(-Xvalp1[0] + 55, 2))/(2*Time+1);

   values[22*(Time+1)+37*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.333333333333333*hstep*(0.166666666666667*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 1)*(-280.0*Xval2[4]*pow(Xval2[3], 3) - 32*pow(Xval2[2], 4) - 4) - 0.25*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) - 1.0)/(2*Time+1);

   values[22*(Time+1)+38*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.333333333333333*hstep*(0.166666666666667*hstep*(-70.0*Xvalp1[4]*pow(Xvalp1[3], 3) - 8*pow(Xvalp1[2], 4) - 1) - 1)*(-280.0*Xval2[4]*pow(Xval2[3], 3) - 32*pow(Xval2[2], 4) - 4) - 0.25*hstep*(70.0*Xvalp1[4]*pow(Xvalp1[3], 3) + 8*pow(Xvalp1[2], 4) + 1) - 1.0)/(2*Time+1);

   values[22*(Time+1)+39*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.05*pow(hstep, 2)*(-280.0*Xval2[4]*pow(Xval2[3], 3) - 32*pow(Xval2[2], 4) - 4) - 0.225*hstep)/(2*Time+1);

   values[22*(Time+1)+40*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.05*pow(hstep, 2)*(-280.0*Xval2[4]*pow(Xval2[3], 3) - 32*pow(Xval2[2], 4) - 4) + 0.225*hstep)/(2*Time+1);

   values[22*(Time+1)+41*(Time)+0+jt] += obj_factor*bounds[0][3]*(1.77777777777778*pow(hstep, 2)*pow(Xval[2], 3)*(-Xval[0] - 90)*(-280.0*Xval2[4]*pow(Xval2[3], 3) - 32*pow(Xval2[2], 4) - 4) - 8.0*hstep*pow(Xval[2], 3)*(-Xval[0] - 90))/(2*Time+1);

   values[22*(Time+1)+42*(Time)+0+jt] += obj_factor*bounds[0][3]*(1.77777777777778*pow(hstep, 2)*pow(Xvalp1[2], 3)*(-Xvalp1[0] - 90)*(-280.0*Xval2[4]*pow(Xval2[3], 3) - 32*pow(Xval2[2], 4) - 4) + 8.0*hstep*pow(Xvalp1[2], 3)*(-Xvalp1[0] - 90))/(2*Time+1);

   values[22*(Time+1)+43*(Time)+0+jt] += obj_factor*bounds[0][3]*(11.6666666666667*Xval[4]*pow(hstep, 2)*pow(Xval[3], 2)*(-Xval[0] + 55)*(-280.0*Xval2[4]*pow(Xval2[3], 3) - 32*pow(Xval2[2], 4) - 4) - 52.5*Xval[4]*hstep*pow(Xval[3], 2)*(-Xval[0] + 55))/(2*Time+1);

   values[22*(Time+1)+44*(Time)+0+jt] += obj_factor*bounds[0][3]*(11.6666666666667*Xvalp1[4]*pow(hstep, 2)*pow(Xvalp1[3], 2)*(-Xvalp1[0] + 55)*(-280.0*Xval2[4]*pow(Xval2[3], 3) - 32*pow(Xval2[2], 4) - 4) + 52.5*Xvalp1[4]*hstep*pow(Xvalp1[3], 2)*(-Xvalp1[0] + 55))/(2*Time+1);

   values[22*(Time+1)+45*(Time)+0+jt] += obj_factor*bounds[0][3]*(3.88888888888889*pow(hstep, 2)*pow(Xval[3], 3)*(-Xval[0] + 55)*(-280.0*Xval2[4]*pow(Xval2[3], 3) - 32*pow(Xval2[2], 4) - 4) - 17.5*hstep*pow(Xval[3], 3)*(-Xval[0] + 55))/(2*Time+1);

   values[22*(Time+1)+46*(Time)+0+jt] += obj_factor*bounds[0][3]*(3.88888888888889*pow(hstep, 2)*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55)*(-280.0*Xval2[4]*pow(Xval2[3], 3) - 32*pow(Xval2[2], 4) - 4) + 17.5*hstep*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55))/(2*Time+1);

   values[22*(Time+1)+51*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0555555555555556*pow(hstep, 2)*pow(-280.0*Xval2[4]*pow(Xval2[3], 3) - 32*pow(Xval2[2], 4) - 4, 2) + 2)/(2*Time+1);

   values[22*(Time+1)+52*(Time)+0+jt] += obj_factor*bounds[0][3]*(1.2*hstep*(0.166666666666667*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 1))/(2*Time+1);

   values[22*(Time+1)+53*(Time)+0+jt] += obj_factor*bounds[0][3]*(1.2*hstep*(0.166666666666667*hstep*(-70.0*Xvalp1[4]*pow(Xvalp1[3], 3) - 8*pow(Xvalp1[2], 4) - 1) - 1))/(2*Time+1);

   values[22*(Time+1)+54*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.18*pow(hstep, 2))/(2*Time+1);

   values[22*(Time+1)+55*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.18*pow(hstep, 2))/(2*Time+1);

   values[22*(Time+1)+56*(Time)+0+jt] += obj_factor*bounds[0][3]*(6.4*pow(hstep, 2)*pow(Xval[2], 3)*(-Xval[0] - 90))/(2*Time+1);

   values[22*(Time+1)+57*(Time)+0+jt] += obj_factor*bounds[0][3]*(6.4*pow(hstep, 2)*pow(Xvalp1[2], 3)*(-Xvalp1[0] - 90))/(2*Time+1);

   values[22*(Time+1)+58*(Time)+0+jt] += obj_factor*bounds[0][3]*(42.0*Xval[4]*pow(hstep, 2)*pow(Xval[3], 2)*(-Xval[0] + 55))/(2*Time+1);

   values[22*(Time+1)+59*(Time)+0+jt] += obj_factor*bounds[0][3]*(42.0*Xvalp1[4]*pow(hstep, 2)*pow(Xvalp1[3], 2)*(-Xvalp1[0] + 55))/(2*Time+1);

   values[22*(Time+1)+60*(Time)+0+jt] += obj_factor*bounds[0][3]*(14.0*pow(hstep, 2)*pow(Xval[3], 3)*(-Xval[0] + 55))/(2*Time+1);

   values[22*(Time+1)+61*(Time)+0+jt] += obj_factor*bounds[0][3]*(14.0*pow(hstep, 2)*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55))/(2*Time+1);

   values[22*(Time+1)+66*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.2*pow(hstep, 2)*(-280.0*Xval2[4]*pow(Xval2[3], 3) - 32*pow(Xval2[2], 4) - 4))/(2*Time+1);

   values[22*(Time+1)+67*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.72*pow(hstep, 2))/(2*Time+1);

   values[22*(Time+1)+68*(Time)+0+jt] += obj_factor*bounds[0][3]*(42.6666666666667*hstep*pow(Xval2[2], 3)*(-Xval2[0] - 90)*(0.166666666666667*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 1))/(2*Time+1);

   values[22*(Time+1)+69*(Time)+0+jt] += obj_factor*bounds[0][3]*(42.6666666666667*hstep*pow(Xval2[2], 3)*(-Xval2[0] - 90)*(0.166666666666667*hstep*(-70.0*Xvalp1[4]*pow(Xvalp1[3], 3) - 8*pow(Xvalp1[2], 4) - 1) - 1))/(2*Time+1);

   values[22*(Time+1)+70*(Time)+0+jt] += obj_factor*bounds[0][3]*(6.4*pow(hstep, 2)*pow(Xval2[2], 3)*(-Xval2[0] - 90))/(2*Time+1);

   values[22*(Time+1)+71*(Time)+0+jt] += obj_factor*bounds[0][3]*(6.4*pow(hstep, 2)*pow(Xval2[2], 3)*(-Xval2[0] - 90))/(2*Time+1);

   values[22*(Time+1)+72*(Time)+0+jt] += obj_factor*bounds[0][3]*(227.555555555556*pow(hstep, 2)*pow(Xval[2], 3)*pow(Xval2[2], 3)*(-Xval[0] - 90)*(-Xval2[0] - 90))/(2*Time+1);

   values[22*(Time+1)+73*(Time)+0+jt] += obj_factor*bounds[0][3]*(227.555555555556*pow(hstep, 2)*pow(Xval2[2], 3)*pow(Xvalp1[2], 3)*(-Xval2[0] - 90)*(-Xvalp1[0] - 90))/(2*Time+1);

   values[22*(Time+1)+74*(Time)+0+jt] += obj_factor*bounds[0][3]*(1493.33333333333*Xval[4]*pow(hstep, 2)*pow(Xval[3], 2)*pow(Xval2[2], 3)*(-Xval[0] + 55)*(-Xval2[0] - 90))/(2*Time+1);

   values[22*(Time+1)+75*(Time)+0+jt] += obj_factor*bounds[0][3]*(1493.33333333333*Xvalp1[4]*pow(hstep, 2)*pow(Xvalp1[3], 2)*pow(Xval2[2], 3)*(-Xval2[0] - 90)*(-Xvalp1[0] + 55))/(2*Time+1);

   values[22*(Time+1)+76*(Time)+0+jt] += obj_factor*bounds[0][3]*(497.777777777778*pow(hstep, 2)*pow(Xval[3], 3)*pow(Xval2[2], 3)*(-Xval[0] + 55)*(-Xval2[0] - 90))/(2*Time+1);

   values[22*(Time+1)+77*(Time)+0+jt] += obj_factor*bounds[0][3]*(497.777777777778*pow(hstep, 2)*pow(Xvalp1[3], 3)*pow(Xval2[2], 3)*(-Xval2[0] - 90)*(-Xvalp1[0] + 55))/(2*Time+1);

   values[22*(Time+1)+78*(Time)+0+jt] += obj_factor*bounds[0][3]*(7.11111111111111*pow(hstep, 2)*pow(Xval2[2], 3)*(-Xval2[0] - 90)*(-280.0*Xval2[4]*pow(Xval2[3], 3) - 32*pow(Xval2[2], 4) - 4) - 42.6666666666667*hstep*pow(Xval2[2], 3)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)))/(2*Time+1);

   values[22*(Time+1)+79*(Time)+0+jt] += obj_factor*bounds[0][3]*(25.6*pow(hstep, 2)*pow(Xval2[2], 3)*(-Xval2[0] - 90))/(2*Time+1);

   values[22*(Time+1)+80*(Time)+0+jt] += obj_factor*bounds[0][3]*(910.222222222222*pow(hstep, 2)*pow(Xval2[2], 6)*pow(-Xval2[0] - 90, 2) + 128.0*hstep*pow(Xval2[2], 2)*(-Xval2[0] - 90)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)))/(2*Time+1);

   values[22*(Time+1)+81*(Time)+0+jt] += obj_factor*bounds[0][3]*(280.0*Xval2[4]*hstep*pow(Xval2[3], 2)*(-Xval2[0] + 55)*(0.166666666666667*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 1))/(2*Time+1);

   values[22*(Time+1)+82*(Time)+0+jt] += obj_factor*bounds[0][3]*(280.0*Xval2[4]*hstep*pow(Xval2[3], 2)*(-Xval2[0] + 55)*(0.166666666666667*hstep*(-70.0*Xvalp1[4]*pow(Xvalp1[3], 3) - 8*pow(Xvalp1[2], 4) - 1) - 1))/(2*Time+1);

   values[22*(Time+1)+83*(Time)+0+jt] += obj_factor*bounds[0][3]*(42.0*Xval2[4]*pow(hstep, 2)*pow(Xval2[3], 2)*(-Xval2[0] + 55))/(2*Time+1);

   values[22*(Time+1)+84*(Time)+0+jt] += obj_factor*bounds[0][3]*(42.0*Xval2[4]*pow(hstep, 2)*pow(Xval2[3], 2)*(-Xval2[0] + 55))/(2*Time+1);

   values[22*(Time+1)+85*(Time)+0+jt] += obj_factor*bounds[0][3]*(1493.33333333333*Xval2[4]*pow(hstep, 2)*pow(Xval2[3], 2)*pow(Xval[2], 3)*(-Xval[0] - 90)*(-Xval2[0] + 55))/(2*Time+1);

   values[22*(Time+1)+86*(Time)+0+jt] += obj_factor*bounds[0][3]*(1493.33333333333*Xval2[4]*pow(hstep, 2)*pow(Xval2[3], 2)*pow(Xvalp1[2], 3)*(-Xval2[0] + 55)*(-Xvalp1[0] - 90))/(2*Time+1);

   values[22*(Time+1)+87*(Time)+0+jt] += obj_factor*bounds[0][3]*(9800.0*Xval[4]*Xval2[4]*pow(hstep, 2)*pow(Xval[3], 2)*pow(Xval2[3], 2)*(-Xval[0] + 55)*(-Xval2[0] + 55))/(2*Time+1);

   values[22*(Time+1)+88*(Time)+0+jt] += obj_factor*bounds[0][3]*(9800.0*Xval2[4]*Xvalp1[4]*pow(hstep, 2)*pow(Xval2[3], 2)*pow(Xvalp1[3], 2)*(-Xval2[0] + 55)*(-Xvalp1[0] + 55))/(2*Time+1);

   values[22*(Time+1)+89*(Time)+0+jt] += obj_factor*bounds[0][3]*(3266.66666666667*Xval2[4]*pow(hstep, 2)*pow(Xval[3], 3)*pow(Xval2[3], 2)*(-Xval[0] + 55)*(-Xval2[0] + 55))/(2*Time+1);

   values[22*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[0][3]*(3266.66666666667*Xval2[4]*pow(hstep, 2)*pow(Xval2[3], 2)*pow(Xvalp1[3], 3)*(-Xval2[0] + 55)*(-Xvalp1[0] + 55))/(2*Time+1);

   values[22*(Time+1)+91*(Time)+0+jt] += obj_factor*bounds[0][3]*(46.6666666666667*Xval2[4]*pow(hstep, 2)*pow(Xval2[3], 2)*(-Xval2[0] + 55)*(-280.0*Xval2[4]*pow(Xval2[3], 3) - 32*pow(Xval2[2], 4) - 4) - 280.0*Xval2[4]*hstep*pow(Xval2[3], 2)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)))/(2*Time+1);

   values[22*(Time+1)+92*(Time)+0+jt] += obj_factor*bounds[0][3]*(168.0*Xval2[4]*pow(hstep, 2)*pow(Xval2[3], 2)*(-Xval2[0] + 55))/(2*Time+1);

   values[22*(Time+1)+93*(Time)+0+jt] += obj_factor*bounds[0][3]*(5973.33333333333*Xval2[4]*pow(hstep, 2)*pow(Xval2[3], 2)*pow(Xval2[2], 3)*(-Xval2[0] - 90)*(-Xval2[0] + 55))/(2*Time+1);

   values[22*(Time+1)+94*(Time)+0+jt] += obj_factor*bounds[0][3]*(39200.0*pow(Xval2[4], 2)*pow(hstep, 2)*pow(Xval2[3], 4)*pow(-Xval2[0] + 55, 2) + 560.0*Xval2[4]*hstep*Xval2[3]*(-Xval2[0] + 55)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)))/(2*Time+1);

   values[22*(Time+1)+95*(Time)+0+jt] += obj_factor*bounds[0][3]*(93.3333333333333*hstep*pow(Xval2[3], 3)*(-Xval2[0] + 55)*(0.166666666666667*hstep*(-70.0*Xval[4]*pow(Xval[3], 3) - 8*pow(Xval[2], 4) - 1) + 1))/(2*Time+1);

   values[22*(Time+1)+96*(Time)+0+jt] += obj_factor*bounds[0][3]*(93.3333333333333*hstep*pow(Xval2[3], 3)*(-Xval2[0] + 55)*(0.166666666666667*hstep*(-70.0*Xvalp1[4]*pow(Xvalp1[3], 3) - 8*pow(Xvalp1[2], 4) - 1) - 1))/(2*Time+1);

   values[22*(Time+1)+97*(Time)+0+jt] += obj_factor*bounds[0][3]*(14.0*pow(hstep, 2)*pow(Xval2[3], 3)*(-Xval2[0] + 55))/(2*Time+1);

   values[22*(Time+1)+98*(Time)+0+jt] += obj_factor*bounds[0][3]*(14.0*pow(hstep, 2)*pow(Xval2[3], 3)*(-Xval2[0] + 55))/(2*Time+1);

   values[22*(Time+1)+99*(Time)+0+jt] += obj_factor*bounds[0][3]*(497.777777777778*pow(hstep, 2)*pow(Xval2[3], 3)*pow(Xval[2], 3)*(-Xval[0] - 90)*(-Xval2[0] + 55))/(2*Time+1);

   values[22*(Time+1)+100*(Time)+0+jt] += obj_factor*bounds[0][3]*(497.777777777778*pow(hstep, 2)*pow(Xval2[3], 3)*pow(Xvalp1[2], 3)*(-Xval2[0] + 55)*(-Xvalp1[0] - 90))/(2*Time+1);

   values[22*(Time+1)+101*(Time)+0+jt] += obj_factor*bounds[0][3]*(3266.66666666667*Xval[4]*pow(hstep, 2)*pow(Xval[3], 2)*pow(Xval2[3], 3)*(-Xval[0] + 55)*(-Xval2[0] + 55))/(2*Time+1);

   values[22*(Time+1)+102*(Time)+0+jt] += obj_factor*bounds[0][3]*(3266.66666666667*Xvalp1[4]*pow(hstep, 2)*pow(Xval2[3], 3)*pow(Xvalp1[3], 2)*(-Xval2[0] + 55)*(-Xvalp1[0] + 55))/(2*Time+1);

   values[22*(Time+1)+103*(Time)+0+jt] += obj_factor*bounds[0][3]*(1088.88888888889*pow(hstep, 2)*pow(Xval[3], 3)*pow(Xval2[3], 3)*(-Xval[0] + 55)*(-Xval2[0] + 55))/(2*Time+1);

   values[22*(Time+1)+104*(Time)+0+jt] += obj_factor*bounds[0][3]*(1088.88888888889*pow(hstep, 2)*pow(Xval2[3], 3)*pow(Xvalp1[3], 3)*(-Xval2[0] + 55)*(-Xvalp1[0] + 55))/(2*Time+1);

   values[22*(Time+1)+105*(Time)+0+jt] += obj_factor*bounds[0][3]*(15.5555555555556*pow(hstep, 2)*pow(Xval2[3], 3)*(-Xval2[0] + 55)*(-280.0*Xval2[4]*pow(Xval2[3], 3) - 32*pow(Xval2[2], 4) - 4) - 93.3333333333333*hstep*pow(Xval2[3], 3)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)))/(2*Time+1);

   values[22*(Time+1)+106*(Time)+0+jt] += obj_factor*bounds[0][3]*(56.0*pow(hstep, 2)*pow(Xval2[3], 3)*(-Xval2[0] + 55))/(2*Time+1);

   values[22*(Time+1)+107*(Time)+0+jt] += obj_factor*bounds[0][3]*(1991.11111111111*pow(hstep, 2)*pow(Xval2[3], 3)*pow(Xval2[2], 3)*(-Xval2[0] - 90)*(-Xval2[0] + 55))/(2*Time+1);

   values[22*(Time+1)+108*(Time)+0+jt] += obj_factor*bounds[0][3]*(13066.6666666667*Xval2[4]*pow(hstep, 2)*pow(Xval2[3], 5)*pow(-Xval2[0] + 55, 2) + 280.0*hstep*pow(Xval2[3], 2)*(-Xval2[0] + 55)*(Xval[0] - Xvalp1[0] + 0.166666666666667*hstep*((1.0L/10.0L)*Ival[0] + (2.0L/5.0L)*Ival2[0] + (1.0L/10.0L)*Ivalp1[0] + (9.0L/10.0L)*Xval[1] + (18.0L/5.0L)*Xval2[1] + (9.0L/10.0L)*Xvalp1[1] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + 70.0*Xval[4]*pow(Xval[3], 3)*(-Xval[0] + 55) + 280.0*Xval2[4]*pow(Xval2[3], 3)*(-Xval2[0] + 55) + 70.0*Xvalp1[4]*pow(Xvalp1[3], 3)*(-Xvalp1[0] + 55) + 8*pow(Xval[2], 4)*(-Xval[0] - 90) + 32*pow(Xval2[2], 4)*(-Xval2[0] - 90) + 8*pow(Xvalp1[2], 4)*(-Xvalp1[0] - 90) - 48)))/(2*Time+1);

   values[22*(Time+1)+109*(Time)+0+jt] += obj_factor*bounds[0][3]*(4355.55555555555*pow(hstep, 2)*pow(Xval2[3], 6)*pow(-Xval2[0] + 55, 2))/(2*Time+1);

   values[0*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.0703125*pow(hstep, 2))/(2*Time+1);

   values[1*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.0196875*pow(hstep, 2))/(2*Time+1);

   values[0*(Time+1)+0*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(0.0703125*pow(hstep, 2))/(2*Time+1);

   values[1*(Time+1)+1*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xval[1],1)
0.1125*hstep*(0.25*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) + 1.0) + 0.15*hstep*(0.333333333333333*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) + 2))/(2*Time+1);

   values[2*(Time+1)+1*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xval[1],1)
-0.1125*hstep*(0.25*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) + 1.0) + 0.15*hstep*(0.333333333333333*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) + 2))/(2*Time+1);

   values[2*(Time+1)+2*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xval[1],1)
// ghkhes(Xval[1],1,0)
0.333333333333333*hstep*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*(-0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkhes(Xval[1],1,0) - 0.0137174211248285*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) + 0.37037037037037*pow(Xval[5], 2)*exp(-0.0740740740740741*Xval[1])*ghkjac(Xval[1],1)) + 0.25*hstep*(-0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkhes(Xval[1],1,0) - 0.0137174211248285*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) + 0.37037037037037*pow(Xval[5], 2)*exp(-0.0740740740740741*Xval[1])*ghkjac(Xval[1],1))*(0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1] + 0.125*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) - 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] + Xvalp1[1] + (9.0L/10.0L)*Xval[0] - 9.0L/10.0L*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]))) + (0.125*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) + 0.5)*(0.25*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) + 1.0) + (0.166666666666667*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) + 1)*(0.333333333333333*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) + 2))/(2*Time+1);

   values[3*(Time+1)+2*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xvalp1[1],1)
0.15*hstep*(0.333333333333333*hstep*(-75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) + 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) - 1) - 2) + 0.1125*hstep*(0.25*hstep*(75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) - 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) + 1) + 1.0))/(2*Time+1);

   values[1*(Time+1)+1*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xvalp1[1],1)
0.15*hstep*(0.333333333333333*hstep*(-75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) + 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) - 1) - 2) - 0.1125*hstep*(0.25*hstep*(75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) - 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) + 1) + 1.0))/(2*Time+1);

   values[3*(Time+1)+3*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xval[1],1)
// ghkjac(Xvalp1[1],1)
(0.125*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) + 0.5)*(0.25*hstep*(75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) - 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) + 1) + 1.0) + (0.166666666666667*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) + 1)*(0.333333333333333*hstep*(-75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) + 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) - 1) - 2))/(2*Time+1);

   values[2*(Time+1)+2*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xvalp1[1],1)
// ghkhes(Xvalp1[1],1,0)
0.333333333333333*hstep*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*(-0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkhes(Xvalp1[1],1,0) - 0.0137174211248285*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) + 0.37037037037037*pow(Xvalp1[5], 2)*exp(-0.0740740740740741*Xvalp1[1])*ghkjac(Xvalp1[1],1)) + 0.25*hstep*(0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkhes(Xvalp1[1],1,0) + 0.0137174211248285*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) - 0.37037037037037*pow(Xvalp1[5], 2)*exp(-0.0740740740740741*Xvalp1[1])*ghkjac(Xvalp1[1],1))*(0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1] + 0.125*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) - 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] + Xvalp1[1] + (9.0L/10.0L)*Xval[0] - 9.0L/10.0L*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]))) + (0.166666666666667*hstep*(-75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) + 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) - 1) - 1)*(0.333333333333333*hstep*(-75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) + 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) - 1) - 2) + (0.125*hstep*(75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) - 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) + 1) + 0.5)*(0.25*hstep*(75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) - 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) + 1) + 1.0))/(2*Time+1);

   values[15*(Time+1)+25*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.00015625*pow(hstep, 2)*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]))/(2*Time+1);

   values[16*(Time+1)+25*(Time)+0+jt] += obj_factor*bounds[1][3]*(-4.375e-5*pow(hstep, 2)*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]))/(2*Time+1);

   values[16*(Time+1)+26*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xval[1],1)
-0.0005*hstep*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(0.125*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) + 0.5)*ghk(Xval[1]) - 0.000666666666666667*hstep*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(0.166666666666667*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) + 1)*ghk(Xval[1]) - 0.000666666666666667*hstep*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*ghkjac(Xval[1],1) - 0.0005*hstep*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1] + 0.125*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) - 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] + Xvalp1[1] + (9.0L/10.0L)*Xval[0] - 9.0L/10.0L*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])))*ghkjac(Xval[1],1) + 0.123456790123457*hstep*Xval[5]*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) + 0.0925925925925926*hstep*Xval[5]*(0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1] + 0.125*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) - 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] + Xvalp1[1] + (9.0L/10.0L)*Xval[0] - 9.0L/10.0L*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])))*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]))/(2*Time+1);

   values[17*(Time+1)+26*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xvalp1[1],1)
-0.000666666666666667*hstep*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(0.166666666666667*hstep*(-75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) + 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) - 1) - 1)*ghk(Xval[1]) - 0.0005*hstep*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(0.125*hstep*(75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) - 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) + 1) + 0.5)*ghk(Xval[1]))/(2*Time+1);

   values[17*(Time+1)+27*(Time)+0+jt] += obj_factor*bounds[1][3]*(3.47222222222222e-7*pow(hstep, 2)*pow(Xval[5], 2)*pow(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]), 2)*pow(ghk(Xval[1]), 2) - 0.000666666666666667*hstep*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*ghk(Xval[1]) - 0.0005*hstep*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1] + 0.125*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) - 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] + Xvalp1[1] + (9.0L/10.0L)*Xval[0] - 9.0L/10.0L*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])))*ghk(Xval[1]))/(2*Time+1);

   values[18*(Time+1)+27*(Time)+0+jt] += obj_factor*bounds[1][3]*(-4.375e-5*pow(hstep, 2)*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]))/(2*Time+1);

   values[15*(Time+1)+25*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(-0.00015625*pow(hstep, 2)*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]))/(2*Time+1);

   values[18*(Time+1)+28*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xval[1],1)
0.0005*hstep*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(0.125*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) + 0.5)*ghk(Xvalp1[1]) - 0.000666666666666667*hstep*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(0.166666666666667*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) + 1)*ghk(Xvalp1[1]))/(2*Time+1);

   values[16*(Time+1)+26*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xvalp1[1],1)
-0.000666666666666667*hstep*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(0.166666666666667*hstep*(-75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) + 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) - 1) - 1)*ghk(Xvalp1[1]) + 0.0005*hstep*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(0.125*hstep*(75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) - 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) + 1) + 0.5)*ghk(Xvalp1[1]) - 0.000666666666666667*hstep*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*ghkjac(Xvalp1[1],1) + 0.0005*hstep*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1] + 0.125*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) - 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] + Xvalp1[1] + (9.0L/10.0L)*Xval[0] - 9.0L/10.0L*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])))*ghkjac(Xvalp1[1],1) + 0.123456790123457*hstep*Xvalp1[5]*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) - 0.0925925925925926*hstep*Xvalp1[5]*(0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1] + 0.125*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) - 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] + Xvalp1[1] + (9.0L/10.0L)*Xval[0] - 9.0L/10.0L*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])))*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]))/(2*Time+1);

   values[18*(Time+1)+29*(Time)+0+jt] += obj_factor*bounds[1][3]*(9.72222222222222e-8*pow(hstep, 2)*Xval[5]*Xvalp1[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xval[1])*ghk(Xvalp1[1]))/(2*Time+1);

   values[17*(Time+1)+27*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(3.47222222222222e-7*pow(hstep, 2)*pow(Xvalp1[5], 2)*pow(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]), 2)*pow(ghk(Xvalp1[1]), 2) - 0.000666666666666667*hstep*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*ghk(Xvalp1[1]) + 0.0005*hstep*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1] + 0.125*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) - 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] + Xvalp1[1] + (9.0L/10.0L)*Xval[0] - 9.0L/10.0L*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])))*ghk(Xvalp1[1]))/(2*Time+1);

   values[18*(Time+1)+30*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.078125*pow(hstep, 2)*(-150*pow(Xval[6], 3)*(-Xval[1] - 90)/pow(pow(Xval[6], 2) + 100, 2) + 150*Xval[6]*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 0.001*pow(Xval[5], 2)*ghk(Xval[1])))/(2*Time+1);

   values[19*(Time+1)+30*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.021875*pow(hstep, 2)*(-150*pow(Xval[6], 3)*(-Xval[1] - 90)/pow(pow(Xval[6], 2) + 100, 2) + 150*Xval[6]*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 0.001*pow(Xval[5], 2)*ghk(Xval[1])))/(2*Time+1);

   values[19*(Time+1)+31*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xval[1],1)
0.25*hstep*(0.125*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) + 0.5)*(-150*pow(Xval[6], 3)*(-Xval[1] - 90)/pow(pow(Xval[6], 2) + 100, 2) + 150*Xval[6]*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 0.001*pow(Xval[5], 2)*ghk(Xval[1])) + 0.333333333333333*hstep*(0.166666666666667*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) + 1)*(-150*pow(Xval[6], 3)*(-Xval[1] - 90)/pow(pow(Xval[6], 2) + 100, 2) + 150*Xval[6]*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 0.001*pow(Xval[5], 2)*ghk(Xval[1])) + 0.333333333333333*hstep*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*(150*pow(Xval[6], 3)/pow(pow(Xval[6], 2) + 100, 2) - 150*Xval[6]/(pow(Xval[6], 2) + 100) + 0.001*pow(Xval[5], 2)*ghkjac(Xval[1],1)) + 0.25*hstep*(150*pow(Xval[6], 3)/pow(pow(Xval[6], 2) + 100, 2) - 150*Xval[6]/(pow(Xval[6], 2) + 100) + 0.001*pow(Xval[5], 2)*ghkjac(Xval[1],1))*(0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1] + 0.125*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) - 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] + Xvalp1[1] + (9.0L/10.0L)*Xval[0] - 9.0L/10.0L*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]))))/(2*Time+1);

   values[20*(Time+1)+31*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xvalp1[1],1)
0.333333333333333*hstep*(0.166666666666667*hstep*(-75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) + 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) - 1) - 1)*(-150*pow(Xval[6], 3)*(-Xval[1] - 90)/pow(pow(Xval[6], 2) + 100, 2) + 150*Xval[6]*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 0.001*pow(Xval[5], 2)*ghk(Xval[1])) + 0.25*hstep*(0.125*hstep*(75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) - 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) + 1) + 0.5)*(-150*pow(Xval[6], 3)*(-Xval[1] - 90)/pow(pow(Xval[6], 2) + 100, 2) + 150*Xval[6]*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 0.001*pow(Xval[5], 2)*ghk(Xval[1])))/(2*Time+1);

   values[20*(Time+1)+32*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.000173611111111111*pow(hstep, 2)*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(-150*pow(Xval[6], 3)*(-Xval[1] - 90)/pow(pow(Xval[6], 2) + 100, 2) + 150*Xval[6]*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 0.001*pow(Xval[5], 2)*ghk(Xval[1]))*ghk(Xval[1]) + 0.000666666666666667*hstep*Xval[5]*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*ghk(Xval[1]) + 0.0005*hstep*Xval[5]*(0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1] + 0.125*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) - 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] + Xvalp1[1] + (9.0L/10.0L)*Xval[0] - 9.0L/10.0L*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])))*ghk(Xval[1]))/(2*Time+1);

   values[21*(Time+1)+32*(Time)+0+jt] += obj_factor*bounds[1][3]*(-4.86111111111111e-5*pow(hstep, 2)*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(-150*pow(Xval[6], 3)*(-Xval[1] - 90)/pow(pow(Xval[6], 2) + 100, 2) + 150*Xval[6]*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 0.001*pow(Xval[5], 2)*ghk(Xval[1]))*ghk(Xvalp1[1]))/(2*Time+1);

   values[21*(Time+1)+33*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.0868055555555556*pow(hstep, 2)*pow(-150*pow(Xval[6], 3)*(-Xval[1] - 90)/pow(pow(Xval[6], 2) + 100, 2) + 150*Xval[6]*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 0.001*pow(Xval[5], 2)*ghk(Xval[1]), 2) + 0.333333333333333*hstep*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*(600*pow(Xval[6], 4)*(-Xval[1] - 90)/pow(pow(Xval[6], 2) + 100, 3) - 750*pow(Xval[6], 2)*(-Xval[1] - 90)/pow(pow(Xval[6], 2) + 100, 2) + (-150*Xval[1] - 13500)/(pow(Xval[6], 2) + 100)) + 0.25*hstep*(600*pow(Xval[6], 4)*(-Xval[1] - 90)/pow(pow(Xval[6], 2) + 100, 3) - 750*pow(Xval[6], 2)*(-Xval[1] - 90)/pow(pow(Xval[6], 2) + 100, 2) + (-150*Xval[1] - 13500)/(pow(Xval[6], 2) + 100))*(0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1] + 0.125*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) - 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] + Xvalp1[1] + (9.0L/10.0L)*Xval[0] - 9.0L/10.0L*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]))))/(2*Time+1);

   values[22*(Time+1)+33*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.05*pow(hstep, 2)*(-150*pow(Xvalp1[6], 3)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) + 150*Xvalp1[6]*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])) + 0.028125*pow(hstep, 2)*(150*pow(Xvalp1[6], 3)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) - 150*Xvalp1[6]*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])))/(2*Time+1);

   values[18*(Time+1)+30*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(0.05*pow(hstep, 2)*(-150*pow(Xvalp1[6], 3)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) + 150*Xvalp1[6]*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])) - 0.028125*pow(hstep, 2)*(150*pow(Xvalp1[6], 3)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) - 150*Xvalp1[6]*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])))/(2*Time+1);

   values[22*(Time+1)+34*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xval[1],1)
0.25*hstep*(0.125*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) + 0.5)*(150*pow(Xvalp1[6], 3)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) - 150*Xvalp1[6]*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])) + 0.333333333333333*hstep*(0.166666666666667*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) + 1)*(-150*pow(Xvalp1[6], 3)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) + 150*Xvalp1[6]*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])))/(2*Time+1);

   values[19*(Time+1)+31*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xvalp1[1],1)
0.333333333333333*hstep*(0.166666666666667*hstep*(-75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) + 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) - 1) - 1)*(-150*pow(Xvalp1[6], 3)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) + 150*Xvalp1[6]*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])) + 0.25*hstep*(0.125*hstep*(75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) - 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) + 1) + 0.5)*(150*pow(Xvalp1[6], 3)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) - 150*Xvalp1[6]*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])) + 0.333333333333333*hstep*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*(150*pow(Xvalp1[6], 3)/pow(pow(Xvalp1[6], 2) + 100, 2) - 150*Xvalp1[6]/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*ghkjac(Xvalp1[1],1)) + 0.25*hstep*(-150*pow(Xvalp1[6], 3)/pow(pow(Xvalp1[6], 2) + 100, 2) + 150*Xvalp1[6]/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*ghkjac(Xvalp1[1],1))*(0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1] + 0.125*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) - 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] + Xvalp1[1] + (9.0L/10.0L)*Xval[0] - 9.0L/10.0L*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]))))/(2*Time+1);

   values[22*(Time+1)+35*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.000111111111111111*pow(hstep, 2)*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(-150*pow(Xvalp1[6], 3)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) + 150*Xvalp1[6]*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]))*ghk(Xval[1]) - 6.25e-5*pow(hstep, 2)*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(150*pow(Xvalp1[6], 3)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) - 150*Xvalp1[6]*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]))*ghk(Xval[1]))/(2*Time+1);

   values[20*(Time+1)+32*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(-0.000111111111111111*pow(hstep, 2)*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(-150*pow(Xvalp1[6], 3)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) + 150*Xvalp1[6]*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]))*ghk(Xvalp1[1]) + 6.25e-5*pow(hstep, 2)*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(150*pow(Xvalp1[6], 3)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) - 150*Xvalp1[6]*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]))*ghk(Xvalp1[1]) + 0.000666666666666667*hstep*Xvalp1[5]*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*ghk(Xvalp1[1]) - 0.0005*hstep*Xvalp1[5]*(0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1] + 0.125*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) - 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] + Xvalp1[1] + (9.0L/10.0L)*Xval[0] - 9.0L/10.0L*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])))*ghk(Xvalp1[1]))/(2*Time+1);

   values[22*(Time+1)+36*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.0555555555555556*pow(hstep, 2)*(-150*pow(Xval[6], 3)*(-Xval[1] - 90)/pow(pow(Xval[6], 2) + 100, 2) + 150*Xval[6]*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 0.001*pow(Xval[5], 2)*ghk(Xval[1]))*(-150*pow(Xvalp1[6], 3)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) + 150*Xvalp1[6]*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])) + 0.03125*pow(hstep, 2)*(-150*pow(Xval[6], 3)*(-Xval[1] - 90)/pow(pow(Xval[6], 2) + 100, 2) + 150*Xval[6]*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 0.001*pow(Xval[5], 2)*ghk(Xval[1]))*(150*pow(Xvalp1[6], 3)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) - 150*Xvalp1[6]*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])))/(2*Time+1);

   values[21*(Time+1)+33*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(0.0555555555555556*pow(hstep, 2)*pow(-150*pow(Xvalp1[6], 3)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) + 150*Xvalp1[6]*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]), 2) + 0.03125*pow(hstep, 2)*pow(150*pow(Xvalp1[6], 3)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) - 150*Xvalp1[6]*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]), 2) + 0.333333333333333*hstep*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*(600*pow(Xvalp1[6], 4)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 3) - 750*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) + (-150*Xvalp1[1] - 13500)/(pow(Xvalp1[6], 2) + 100)) + 0.25*hstep*(-600*pow(Xvalp1[6], 4)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 3) + 750*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) - (-150*Xvalp1[1] - 13500)/(pow(Xvalp1[6], 2) + 100))*(0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1] + 0.125*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) - 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] + Xvalp1[1] + (9.0L/10.0L)*Xval[0] - 9.0L/10.0L*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]))))/(2*Time+1);

   values[22*(Time+1)+37*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.18*pow(hstep, 2))/(2*Time+1);

   values[22*(Time+1)+38*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.18*pow(hstep, 2))/(2*Time+1);

   values[22*(Time+1)+39*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xval[1],1)
1.2*hstep*(0.166666666666667*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) + 1))/(2*Time+1);

   values[22*(Time+1)+40*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xvalp1[1],1)
1.2*hstep*(0.166666666666667*hstep*(-75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) + 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) - 1) - 1))/(2*Time+1);

   values[22*(Time+1)+47*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.0004*pow(hstep, 2)*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]))/(2*Time+1);

   values[22*(Time+1)+48*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.0004*pow(hstep, 2)*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]))/(2*Time+1);

   values[22*(Time+1)+49*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.2*pow(hstep, 2)*(-150*pow(Xval[6], 3)*(-Xval[1] - 90)/pow(pow(Xval[6], 2) + 100, 2) + 150*Xval[6]*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 0.001*pow(Xval[5], 2)*ghk(Xval[1])))/(2*Time+1);

   values[22*(Time+1)+50*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.2*pow(hstep, 2)*(-150*pow(Xvalp1[6], 3)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) + 150*Xvalp1[6]*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])))/(2*Time+1);

   values[22*(Time+1)+51*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.72*pow(hstep, 2))/(2*Time+1);

   values[22*(Time+1)+52*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xval2[1],1)
0.05*pow(hstep, 2)*(-300*pow(Xval2[6], 2)/(pow(Xval2[6], 2) + 100) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkjac(Xval2[1],1) + 0.740740740740741*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1]) - 4) - 0.225*hstep)/(2*Time+1);

   values[22*(Time+1)+53*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xval2[1],1)
0.05*pow(hstep, 2)*(-300*pow(Xval2[6], 2)/(pow(Xval2[6], 2) + 100) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkjac(Xval2[1],1) + 0.740740740740741*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1]) - 4) + 0.225*hstep)/(2*Time+1);

   values[22*(Time+1)+54*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xval[1],1)
// ghkjac(Xval2[1],1)
0.333333333333333*hstep*(0.166666666666667*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) + 1)*(-300*pow(Xval2[6], 2)/(pow(Xval2[6], 2) + 100) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkjac(Xval2[1],1) + 0.740740740740741*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1]) - 4) - 0.25*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) - 1.0)/(2*Time+1);

   values[22*(Time+1)+55*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xval2[1],1)
// ghkjac(Xvalp1[1],1)
0.333333333333333*hstep*(0.166666666666667*hstep*(-75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) + 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) - 1) - 1)*(-300*pow(Xval2[6], 2)/(pow(Xval2[6], 2) + 100) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkjac(Xval2[1],1) + 0.740740740740741*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1]) - 4) - 0.25*hstep*(75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) - 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) + 1) - 1.0)/(2*Time+1);

   values[22*(Time+1)+62*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xval2[1],1)
-0.000111111111111111*pow(hstep, 2)*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(-300*pow(Xval2[6], 2)/(pow(Xval2[6], 2) + 100) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkjac(Xval2[1],1) + 0.740740740740741*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1]) - 4)*ghk(Xval[1]) + 0.0005*hstep*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]))/(2*Time+1);

   values[22*(Time+1)+63*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xval2[1],1)
-0.000111111111111111*pow(hstep, 2)*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(-300*pow(Xval2[6], 2)/(pow(Xval2[6], 2) + 100) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkjac(Xval2[1],1) + 0.740740740740741*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1]) - 4)*ghk(Xvalp1[1]) - 0.0005*hstep*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]))/(2*Time+1);

   values[22*(Time+1)+64*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xval2[1],1)
0.0555555555555556*pow(hstep, 2)*(-150*pow(Xval[6], 3)*(-Xval[1] - 90)/pow(pow(Xval[6], 2) + 100, 2) + 150*Xval[6]*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 0.001*pow(Xval[5], 2)*ghk(Xval[1]))*(-300*pow(Xval2[6], 2)/(pow(Xval2[6], 2) + 100) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkjac(Xval2[1],1) + 0.740740740740741*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1]) - 4) - 0.25*hstep*(-150*pow(Xval[6], 3)*(-Xval[1] - 90)/pow(pow(Xval[6], 2) + 100, 2) + 150*Xval[6]*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 0.001*pow(Xval[5], 2)*ghk(Xval[1])))/(2*Time+1);

   values[22*(Time+1)+65*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xval2[1],1)
0.0555555555555556*pow(hstep, 2)*(-150*pow(Xvalp1[6], 3)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) + 150*Xvalp1[6]*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]))*(-300*pow(Xval2[6], 2)/(pow(Xval2[6], 2) + 100) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkjac(Xval2[1],1) + 0.740740740740741*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1]) - 4) - 0.25*hstep*(150*pow(Xvalp1[6], 3)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) - 150*Xvalp1[6]*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])))/(2*Time+1);

   values[22*(Time+1)+66*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xval2[1],1)
0.2*pow(hstep, 2)*(-300*pow(Xval2[6], 2)/(pow(Xval2[6], 2) + 100) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkjac(Xval2[1],1) + 0.740740740740741*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1]) - 4))/(2*Time+1);

   values[22*(Time+1)+67*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xval2[1],1)
// ghkhes(Xval2[1],1,0)
0.0555555555555556*pow(hstep, 2)*pow(-300*pow(Xval2[6], 2)/(pow(Xval2[6], 2) + 100) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkjac(Xval2[1],1) + 0.740740740740741*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1]) - 4, 2) + 0.333333333333333*hstep*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*(-0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkhes(Xval2[1],1,0) - 0.0548696844993141*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1]) + 1.48148148148148*pow(Xval2[5], 2)*exp(-0.0740740740740741*Xval2[1])*ghkjac(Xval2[1],1)) + 2)/(2*Time+1);

   values[22*(Time+1)+110*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.0004*pow(hstep, 2)*Xval2[5]*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]))/(2*Time+1);

   values[22*(Time+1)+111*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.0004*pow(hstep, 2)*Xval2[5]*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]))/(2*Time+1);

   values[22*(Time+1)+112*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xval[1],1)
-0.00266666666666667*hstep*Xval2[5]*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*(0.166666666666667*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) + 1)*ghk(Xval2[1]))/(2*Time+1);

   values[22*(Time+1)+113*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xvalp1[1],1)
-0.00266666666666667*hstep*Xval2[5]*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*(0.166666666666667*hstep*(-75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) + 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) - 1) - 1)*ghk(Xval2[1]))/(2*Time+1);

   values[22*(Time+1)+114*(Time)+0+jt] += obj_factor*bounds[1][3]*(8.88888888888889e-7*pow(hstep, 2)*Xval[5]*Xval2[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval[1])*ghk(Xval2[1]))/(2*Time+1);

   values[22*(Time+1)+115*(Time)+0+jt] += obj_factor*bounds[1][3]*(8.88888888888889e-7*pow(hstep, 2)*Xval2[5]*Xvalp1[5]*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xval2[1])*ghk(Xvalp1[1]))/(2*Time+1);

   values[22*(Time+1)+116*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.000444444444444444*pow(hstep, 2)*Xval2[5]*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*(-150*pow(Xval[6], 3)*(-Xval[1] - 90)/pow(pow(Xval[6], 2) + 100, 2) + 150*Xval[6]*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 0.001*pow(Xval[5], 2)*ghk(Xval[1]))*ghk(Xval2[1]))/(2*Time+1);

   values[22*(Time+1)+117*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.000444444444444444*pow(hstep, 2)*Xval2[5]*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*(-150*pow(Xvalp1[6], 3)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) + 150*Xvalp1[6]*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]))*ghk(Xval2[1]))/(2*Time+1);

   values[22*(Time+1)+118*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.0016*pow(hstep, 2)*Xval2[5]*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]))/(2*Time+1);

   values[22*(Time+1)+119*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xval2[1],1)
-0.000444444444444444*pow(hstep, 2)*Xval2[5]*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*(-300*pow(Xval2[6], 2)/(pow(Xval2[6], 2) + 100) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkjac(Xval2[1],1) + 0.740740740740741*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1]) - 4)*ghk(Xval2[1]) - 0.00266666666666667*hstep*Xval2[5]*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*ghkjac(Xval2[1],1) + 0.493827160493827*hstep*Xval2[5]*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1]))/(2*Time+1);

   values[22*(Time+1)+120*(Time)+0+jt] += obj_factor*bounds[1][3]*(3.55555555555556e-6*pow(hstep, 2)*pow(Xval2[5], 2)*pow(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]), 2)*pow(ghk(Xval2[1]), 2) - 0.00266666666666667*hstep*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*ghk(Xval2[1]))/(2*Time+1);

   values[22*(Time+1)+121*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.05*pow(hstep, 2)*(-600*pow(Xval2[6], 3)*(-Xval2[1] - 90)/pow(pow(Xval2[6], 2) + 100, 2) + 600*Xval2[6]*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 0.004*pow(Xval2[5], 2)*ghk(Xval2[1])))/(2*Time+1);

   values[22*(Time+1)+122*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.05*pow(hstep, 2)*(-600*pow(Xval2[6], 3)*(-Xval2[1] - 90)/pow(pow(Xval2[6], 2) + 100, 2) + 600*Xval2[6]*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 0.004*pow(Xval2[5], 2)*ghk(Xval2[1])))/(2*Time+1);

   values[22*(Time+1)+123*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xval[1],1)
0.333333333333333*hstep*(0.166666666666667*hstep*(-75*pow(Xval[6], 2)/(pow(Xval[6], 2) + 100) - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.185185185185185*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) - 1) + 1)*(-600*pow(Xval2[6], 3)*(-Xval2[1] - 90)/pow(pow(Xval2[6], 2) + 100, 2) + 600*Xval2[6]*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 0.004*pow(Xval2[5], 2)*ghk(Xval2[1])))/(2*Time+1);

   values[22*(Time+1)+124*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xvalp1[1],1)
0.333333333333333*hstep*(0.166666666666667*hstep*(-75*pow(Xvalp1[6], 2)/(pow(Xvalp1[6], 2) + 100) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) + 0.185185185185185*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) - 1) - 1)*(-600*pow(Xval2[6], 3)*(-Xval2[1] - 90)/pow(pow(Xval2[6], 2) + 100, 2) + 600*Xval2[6]*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 0.004*pow(Xval2[5], 2)*ghk(Xval2[1])))/(2*Time+1);

   values[22*(Time+1)+125*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.000111111111111111*pow(hstep, 2)*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(-600*pow(Xval2[6], 3)*(-Xval2[1] - 90)/pow(pow(Xval2[6], 2) + 100, 2) + 600*Xval2[6]*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 0.004*pow(Xval2[5], 2)*ghk(Xval2[1]))*ghk(Xval[1]))/(2*Time+1);

   values[22*(Time+1)+126*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.000111111111111111*pow(hstep, 2)*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(-600*pow(Xval2[6], 3)*(-Xval2[1] - 90)/pow(pow(Xval2[6], 2) + 100, 2) + 600*Xval2[6]*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 0.004*pow(Xval2[5], 2)*ghk(Xval2[1]))*ghk(Xvalp1[1]))/(2*Time+1);

   values[22*(Time+1)+127*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.0555555555555556*pow(hstep, 2)*(-150*pow(Xval[6], 3)*(-Xval[1] - 90)/pow(pow(Xval[6], 2) + 100, 2) + 150*Xval[6]*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 0.001*pow(Xval[5], 2)*ghk(Xval[1]))*(-600*pow(Xval2[6], 3)*(-Xval2[1] - 90)/pow(pow(Xval2[6], 2) + 100, 2) + 600*Xval2[6]*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 0.004*pow(Xval2[5], 2)*ghk(Xval2[1])))/(2*Time+1);

   values[22*(Time+1)+128*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.0555555555555556*pow(hstep, 2)*(-600*pow(Xval2[6], 3)*(-Xval2[1] - 90)/pow(pow(Xval2[6], 2) + 100, 2) + 600*Xval2[6]*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 0.004*pow(Xval2[5], 2)*ghk(Xval2[1]))*(-150*pow(Xvalp1[6], 3)*(-Xvalp1[1] - 90)/pow(pow(Xvalp1[6], 2) + 100, 2) + 150*Xvalp1[6]*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) + 0.001*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])))/(2*Time+1);

   values[22*(Time+1)+129*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.2*pow(hstep, 2)*(-600*pow(Xval2[6], 3)*(-Xval2[1] - 90)/pow(pow(Xval2[6], 2) + 100, 2) + 600*Xval2[6]*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 0.004*pow(Xval2[5], 2)*ghk(Xval2[1])))/(2*Time+1);

   values[22*(Time+1)+130*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not C:
// ghkjac(Xval2[1],1)
0.0555555555555556*pow(hstep, 2)*(-600*pow(Xval2[6], 3)*(-Xval2[1] - 90)/pow(pow(Xval2[6], 2) + 100, 2) + 600*Xval2[6]*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 0.004*pow(Xval2[5], 2)*ghk(Xval2[1]))*(-300*pow(Xval2[6], 2)/(pow(Xval2[6], 2) + 100) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkjac(Xval2[1],1) + 0.740740740740741*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1]) - 4) + 0.333333333333333*hstep*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*(600*pow(Xval2[6], 3)/pow(pow(Xval2[6], 2) + 100, 2) - 600*Xval2[6]/(pow(Xval2[6], 2) + 100) + 0.004*pow(Xval2[5], 2)*ghkjac(Xval2[1],1)))/(2*Time+1);

   values[22*(Time+1)+131*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.000444444444444444*pow(hstep, 2)*Xval2[5]*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*(-600*pow(Xval2[6], 3)*(-Xval2[1] - 90)/pow(pow(Xval2[6], 2) + 100, 2) + 600*Xval2[6]*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 0.004*pow(Xval2[5], 2)*ghk(Xval2[1]))*ghk(Xval2[1]) + 0.00266666666666667*hstep*Xval2[5]*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*ghk(Xval2[1]))/(2*Time+1);

   values[22*(Time+1)+132*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.0555555555555556*pow(hstep, 2)*pow(-600*pow(Xval2[6], 3)*(-Xval2[1] - 90)/pow(pow(Xval2[6], 2) + 100, 2) + 600*Xval2[6]*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 0.004*pow(Xval2[5], 2)*ghk(Xval2[1]), 2) + 0.333333333333333*hstep*(Xval[1] - Xvalp1[1] + 0.166666666666667*hstep*(75*pow(Xval[6], 2)*(-Xval[1] - 90)/(pow(Xval[6], 2) + 100) + 300*pow(Xval2[6], 2)*(-Xval2[1] - 90)/(pow(Xval2[6], 2) + 100) + 75*pow(Xvalp1[6], 2)*(-Xvalp1[1] - 90)/(pow(Xvalp1[6], 2) + 100) - Xval[1] - 4*Xval2[1] - Xvalp1[1] + (9.0L/10.0L)*Xval[0] + (18.0L/5.0L)*Xval2[0] + (9.0L/10.0L)*Xvalp1[0] - 0.001*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 0.004*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 0.001*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 48))*(2400*pow(Xval2[6], 4)*(-Xval2[1] - 90)/pow(pow(Xval2[6], 2) + 100, 3) - 3000*pow(Xval2[6], 2)*(-Xval2[1] - 90)/pow(pow(Xval2[6], 2) + 100, 2) + (-600*Xval2[1] - 54000)/(pow(Xval2[6], 2) + 100)))/(2*Time+1);

   values[0*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.0868055555555556*pow(hstep, 2)*pow((-0.025*pow(tanh(0.05*Xval[0] + 1.75), 2) + 0.025)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + (-0.008325*pow(tanh(-0.0333*Xval[0] - 0.8991), 2) + 0.008325)*(-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35, 2), 2) + 0.333333333333333*hstep*(0.166666666666667*hstep*((-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-4*Xval2[2] + 2.0*tanh(0.05*Xval2[0] + 1.75) + 2.0)/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) + (-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35)) + Xval[2] - Xvalp1[2])*((0.0025*pow(tanh(0.05*Xval[0] + 1.75), 2) - 0.0025)*tanh(0.05*Xval[0] + 1.75)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + (-0.01665*pow(tanh(-0.0333*Xval[0] - 0.8991), 2) + 0.01665)*(-0.025*pow(tanh(0.05*Xval[0] + 1.75), 2) + 0.025)/pow(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35, 2) - (0.000554445*pow(tanh(-0.0333*Xval[0] - 0.8991), 2) - 0.000554445)*(-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)*tanh(-0.0333*Xval[0] - 0.8991)/pow(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35, 2) + (-0.01665*pow(tanh(-0.0333*Xval[0] - 0.8991), 2) + 0.01665)*(-0.008325*pow(tanh(-0.0333*Xval[0] - 0.8991), 2) + 0.008325)*(-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35, 3)) + 0.25*hstep*(0.125*hstep*((Xvalp1[2] - 0.5*tanh(0.05*Xvalp1[0] + 1.75) - 0.5)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35)) + 0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2])*((0.0025*pow(tanh(0.05*Xval[0] + 1.75), 2) - 0.0025)*tanh(0.05*Xval[0] + 1.75)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + (-0.01665*pow(tanh(-0.0333*Xval[0] - 0.8991), 2) + 0.01665)*(-0.025*pow(tanh(0.05*Xval[0] + 1.75), 2) + 0.025)/pow(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35, 2) - (0.000554445*pow(tanh(-0.0333*Xval[0] - 0.8991), 2) - 0.000554445)*(-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)*tanh(-0.0333*Xval[0] - 0.8991)/pow(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35, 2) + (-0.01665*pow(tanh(-0.0333*Xval[0] - 0.8991), 2) + 0.01665)*(-0.008325*pow(tanh(-0.0333*Xval[0] - 0.8991), 2) + 0.008325)*(-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35, 3)))/(2*Time+1);

   values[1*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.0555555555555556*pow(hstep, 2)*((-0.025*pow(tanh(0.05*Xval[0] + 1.75), 2) + 0.025)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + (-0.008325*pow(tanh(-0.0333*Xval[0] - 0.8991), 2) + 0.008325)*(-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35, 2))*((-0.025*pow(tanh(0.05*Xvalp1[0] + 1.75), 2) + 0.025)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-0.008325*pow(tanh(-0.0333*Xvalp1[0] - 0.8991), 2) + 0.008325)*(-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35, 2)) + 0.03125*pow(hstep, 2)*((-0.025*pow(tanh(0.05*Xval[0] + 1.75), 2) + 0.025)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + (-0.008325*pow(tanh(-0.0333*Xval[0] - 0.8991), 2) + 0.008325)*(-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35, 2))*((0.025*pow(tanh(0.05*Xvalp1[0] + 1.75), 2) - 0.025)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) - (-0.008325*pow(tanh(-0.0333*Xvalp1[0] - 0.8991), 2) + 0.008325)*(-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35, 2)))/(2*Time+1);

   values[0*(Time+1)+0*(Time)+0+ 1+jt] += obj_factor*bounds[2][3]*(0.0555555555555556*pow(hstep, 2)*pow((-0.025*pow(tanh(0.05*Xvalp1[0] + 1.75), 2) + 0.025)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-0.008325*pow(tanh(-0.0333*Xvalp1[0] - 0.8991), 2) + 0.008325)*(-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35, 2), 2) + 0.03125*pow(hstep, 2)*pow((0.025*pow(tanh(0.05*Xvalp1[0] + 1.75), 2) - 0.025)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) - (-0.008325*pow(tanh(-0.0333*Xvalp1[0] - 0.8991), 2) + 0.008325)*(-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35, 2), 2) + 0.333333333333333*hstep*(0.166666666666667*hstep*((-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-4*Xval2[2] + 2.0*tanh(0.05*Xval2[0] + 1.75) + 2.0)/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) + (-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35)) + Xval[2] - Xvalp1[2])*((0.0025*pow(tanh(0.05*Xvalp1[0] + 1.75), 2) - 0.0025)*tanh(0.05*Xvalp1[0] + 1.75)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-0.01665*pow(tanh(-0.0333*Xvalp1[0] - 0.8991), 2) + 0.01665)*(-0.025*pow(tanh(0.05*Xvalp1[0] + 1.75), 2) + 0.025)/pow(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35, 2) - (0.000554445*pow(tanh(-0.0333*Xvalp1[0] - 0.8991), 2) - 0.000554445)*(-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)*tanh(-0.0333*Xvalp1[0] - 0.8991)/pow(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35, 2) + (-0.01665*pow(tanh(-0.0333*Xvalp1[0] - 0.8991), 2) + 0.01665)*(-0.008325*pow(tanh(-0.0333*Xvalp1[0] - 0.8991), 2) + 0.008325)*(-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35, 3)) + 0.25*hstep*(0.125*hstep*((Xvalp1[2] - 0.5*tanh(0.05*Xvalp1[0] + 1.75) - 0.5)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35)) + 0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2])*((-0.0025*pow(tanh(0.05*Xvalp1[0] + 1.75), 2) + 0.0025)*tanh(0.05*Xvalp1[0] + 1.75)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) - (-0.01665*pow(tanh(-0.0333*Xvalp1[0] - 0.8991), 2) + 0.01665)*(-0.025*pow(tanh(0.05*Xvalp1[0] + 1.75), 2) + 0.025)/pow(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35, 2) + (0.000554445*pow(tanh(-0.0333*Xvalp1[0] - 0.8991), 2) - 0.000554445)*(-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)*tanh(-0.0333*Xvalp1[0] - 0.8991)/pow(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35, 2) - (-0.01665*pow(tanh(-0.0333*Xvalp1[0] - 0.8991), 2) + 0.01665)*(-0.008325*pow(tanh(-0.0333*Xvalp1[0] - 0.8991), 2) + 0.008325)*(-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35, 3)))/(2*Time+1);

   values[3*(Time+1)+4*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + 2)*((-0.025*pow(tanh(0.05*Xval[0] + 1.75), 2) + 0.025)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + (-0.008325*pow(tanh(-0.0333*Xval[0] - 0.8991), 2) + 0.008325)*(-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35, 2)) + 0.125*hstep*(-0.25*hstep/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + 1.0)*((-0.025*pow(tanh(0.05*Xval[0] + 1.75), 2) + 0.025)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + (-0.008325*pow(tanh(-0.0333*Xval[0] - 0.8991), 2) + 0.008325)*(-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35, 2)) - 0.333333333333333*hstep*(-0.008325*pow(tanh(-0.0333*Xval[0] - 0.8991), 2) + 0.008325)*(0.166666666666667*hstep*((-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-4*Xval2[2] + 2.0*tanh(0.05*Xval2[0] + 1.75) + 2.0)/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) + (-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35)) + Xval[2] - Xvalp1[2])/pow(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35, 2) - 0.25*hstep*(-0.008325*pow(tanh(-0.0333*Xval[0] - 0.8991), 2) + 0.008325)*(0.125*hstep*((Xvalp1[2] - 0.5*tanh(0.05*Xvalp1[0] + 1.75) - 0.5)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35)) + 0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2])/pow(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35, 2))/(2*Time+1);

   values[4*(Time+1)+4*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + 2)*((-0.025*pow(tanh(0.05*Xvalp1[0] + 1.75), 2) + 0.025)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-0.008325*pow(tanh(-0.0333*Xvalp1[0] - 0.8991), 2) + 0.008325)*(-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35, 2)) + 0.125*hstep*(-0.25*hstep/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + 1.0)*((0.025*pow(tanh(0.05*Xvalp1[0] + 1.75), 2) - 0.025)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) - (-0.008325*pow(tanh(-0.0333*Xvalp1[0] - 0.8991), 2) + 0.008325)*(-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35, 2)))/(2*Time+1);

   values[5*(Time+1)+6*(Time)+0+jt] += obj_factor*bounds[2][3]*((-0.333333333333333*hstep/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + 2)*(-0.166666666666667*hstep/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + 1) + (-0.25*hstep/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + 1.0)*(-0.125*hstep/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + 0.5))/(2*Time+1);

   values[6*(Time+1)+6*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) - 2)*((-0.025*pow(tanh(0.05*Xval[0] + 1.75), 2) + 0.025)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + (-0.008325*pow(tanh(-0.0333*Xval[0] - 0.8991), 2) + 0.008325)*(-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35, 2)) + 0.125*hstep*(0.25*hstep/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + 1.0)*((-0.025*pow(tanh(0.05*Xval[0] + 1.75), 2) + 0.025)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + (-0.008325*pow(tanh(-0.0333*Xval[0] - 0.8991), 2) + 0.008325)*(-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35, 2)))/(2*Time+1);

   values[3*(Time+1)+4*(Time)+0+ 1+jt] += obj_factor*bounds[2][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) - 2)*((-0.025*pow(tanh(0.05*Xvalp1[0] + 1.75), 2) + 0.025)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-0.008325*pow(tanh(-0.0333*Xvalp1[0] - 0.8991), 2) + 0.008325)*(-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35, 2)) + 0.125*hstep*(0.25*hstep/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + 1.0)*((0.025*pow(tanh(0.05*Xvalp1[0] + 1.75), 2) - 0.025)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) - (-0.008325*pow(tanh(-0.0333*Xvalp1[0] - 0.8991), 2) + 0.008325)*(-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35, 2)) - 0.333333333333333*hstep*(-0.008325*pow(tanh(-0.0333*Xvalp1[0] - 0.8991), 2) + 0.008325)*(0.166666666666667*hstep*((-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-4*Xval2[2] + 2.0*tanh(0.05*Xval2[0] + 1.75) + 2.0)/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) + (-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35)) + Xval[2] - Xvalp1[2])/pow(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35, 2) + 0.25*hstep*(-0.008325*pow(tanh(-0.0333*Xvalp1[0] - 0.8991), 2) + 0.008325)*(0.125*hstep*((Xvalp1[2] - 0.5*tanh(0.05*Xvalp1[0] + 1.75) - 0.5)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35)) + 0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2])/pow(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35, 2))/(2*Time+1);

   values[6*(Time+1)+8*(Time)+0+jt] += obj_factor*bounds[2][3]*((-0.166666666666667*hstep/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + 1)*(-0.333333333333333*hstep/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) - 2) + (-0.125*hstep/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + 0.5)*(0.25*hstep/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + 1.0))/(2*Time+1);

   values[5*(Time+1)+6*(Time)+0+ 1+jt] += obj_factor*bounds[2][3]*((-0.333333333333333*hstep/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) - 2)*(-0.166666666666667*hstep/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) - 1) + (0.125*hstep/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + 0.5)*(0.25*hstep/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + 1.0))/(2*Time+1);

   values[22*(Time+1)+37*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.0555555555555556*pow(hstep, 2)*((-0.025*pow(tanh(0.05*Xval[0] + 1.75), 2) + 0.025)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + (-0.008325*pow(tanh(-0.0333*Xval[0] - 0.8991), 2) + 0.008325)*(-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35, 2))*((-0.1*pow(tanh(0.05*Xval2[0] + 1.75), 2) + 0.1)/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) + (-0.0333*pow(tanh(-0.0333*Xval2[0] - 0.8991), 2) + 0.0333)*(-Xval2[2] + 0.5*tanh(0.05*Xval2[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35, 2)))/(2*Time+1);

   values[22*(Time+1)+38*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.0555555555555556*pow(hstep, 2)*((-0.1*pow(tanh(0.05*Xval2[0] + 1.75), 2) + 0.1)/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) + (-0.0333*pow(tanh(-0.0333*Xval2[0] - 0.8991), 2) + 0.0333)*(-Xval2[2] + 0.5*tanh(0.05*Xval2[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35, 2))*((-0.025*pow(tanh(0.05*Xvalp1[0] + 1.75), 2) + 0.025)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-0.008325*pow(tanh(-0.0333*Xvalp1[0] - 0.8991), 2) + 0.008325)*(-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35, 2)))/(2*Time+1);

   values[22*(Time+1)+41*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.333333333333333*hstep*(-0.166666666666667*hstep/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + 1)*((-0.1*pow(tanh(0.05*Xval2[0] + 1.75), 2) + 0.1)/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) + (-0.0333*pow(tanh(-0.0333*Xval2[0] - 0.8991), 2) + 0.0333)*(-Xval2[2] + 0.5*tanh(0.05*Xval2[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35, 2)))/(2*Time+1);

   values[22*(Time+1)+42*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.333333333333333*hstep*(-0.166666666666667*hstep/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) - 1)*((-0.1*pow(tanh(0.05*Xval2[0] + 1.75), 2) + 0.1)/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) + (-0.0333*pow(tanh(-0.0333*Xval2[0] - 0.8991), 2) + 0.0333)*(-Xval2[2] + 0.5*tanh(0.05*Xval2[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35, 2)))/(2*Time+1);

   values[22*(Time+1)+51*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.0555555555555556*pow(hstep, 2)*pow((-0.1*pow(tanh(0.05*Xval2[0] + 1.75), 2) + 0.1)/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) + (-0.0333*pow(tanh(-0.0333*Xval2[0] - 0.8991), 2) + 0.0333)*(-Xval2[2] + 0.5*tanh(0.05*Xval2[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35, 2), 2) + 0.333333333333333*hstep*(0.166666666666667*hstep*((-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-4*Xval2[2] + 2.0*tanh(0.05*Xval2[0] + 1.75) + 2.0)/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) + (-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35)) + Xval[2] - Xvalp1[2])*((0.01*pow(tanh(0.05*Xval2[0] + 1.75), 2) - 0.01)*tanh(0.05*Xval2[0] + 1.75)/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) + (-0.0666*pow(tanh(-0.0333*Xval2[0] - 0.8991), 2) + 0.0666)*(-0.025*pow(tanh(0.05*Xval2[0] + 1.75), 2) + 0.025)/pow(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35, 2) - (0.00221778*pow(tanh(-0.0333*Xval2[0] - 0.8991), 2) - 0.00221778)*(-Xval2[2] + 0.5*tanh(0.05*Xval2[0] + 1.75) + 0.5)*tanh(-0.0333*Xval2[0] - 0.8991)/pow(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35, 2) + (-0.0666*pow(tanh(-0.0333*Xval2[0] - 0.8991), 2) + 0.0666)*(-0.008325*pow(tanh(-0.0333*Xval2[0] - 0.8991), 2) + 0.008325)*(-Xval2[2] + 0.5*tanh(0.05*Xval2[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35, 3)))/(2*Time+1);

   values[22*(Time+1)+68*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.222222222222222*pow(hstep, 2)*((-0.025*pow(tanh(0.05*Xval[0] + 1.75), 2) + 0.025)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + (-0.008325*pow(tanh(-0.0333*Xval[0] - 0.8991), 2) + 0.008325)*(-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35, 2))/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) - 0.25*hstep*((-0.025*pow(tanh(0.05*Xval[0] + 1.75), 2) + 0.025)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + (-0.008325*pow(tanh(-0.0333*Xval[0] - 0.8991), 2) + 0.008325)*(-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35, 2)))/(2*Time+1);

   values[22*(Time+1)+69*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.222222222222222*pow(hstep, 2)*((-0.025*pow(tanh(0.05*Xvalp1[0] + 1.75), 2) + 0.025)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-0.008325*pow(tanh(-0.0333*Xvalp1[0] - 0.8991), 2) + 0.008325)*(-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35, 2))/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) - 0.25*hstep*((0.025*pow(tanh(0.05*Xvalp1[0] + 1.75), 2) - 0.025)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) - (-0.008325*pow(tanh(-0.0333*Xvalp1[0] - 0.8991), 2) + 0.008325)*(-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35, 2)))/(2*Time+1);

   values[22*(Time+1)+72*(Time)+0+jt] += obj_factor*bounds[2][3]*(-1.33333333333333*hstep*(-0.166666666666667*hstep/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) + 1)/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) + 0.25*hstep/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35) - 1.0)/(2*Time+1);

   values[22*(Time+1)+73*(Time)+0+jt] += obj_factor*bounds[2][3]*(-1.33333333333333*hstep*(-0.166666666666667*hstep/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) - 1)/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) - 0.25*hstep/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) - 1.0)/(2*Time+1);

   values[22*(Time+1)+78*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.222222222222222*pow(hstep, 2)*((-0.1*pow(tanh(0.05*Xval2[0] + 1.75), 2) + 0.1)/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) + (-0.0333*pow(tanh(-0.0333*Xval2[0] - 0.8991), 2) + 0.0333)*(-Xval2[2] + 0.5*tanh(0.05*Xval2[0] + 1.75) + 0.5)/pow(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35, 2))/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) - 1.33333333333333*hstep*(-0.008325*pow(tanh(-0.0333*Xval2[0] - 0.8991), 2) + 0.008325)*(0.166666666666667*hstep*((-Xvalp1[2] + 0.5*tanh(0.05*Xvalp1[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xvalp1[0] - 0.8991) + 0.35) + (-4*Xval2[2] + 2.0*tanh(0.05*Xval2[0] + 1.75) + 2.0)/(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35) + (-Xval[2] + 0.5*tanh(0.05*Xval[0] + 1.75) + 0.5)/(0.25*tanh(-0.0333*Xval[0] - 0.8991) + 0.35)) + Xval[2] - Xvalp1[2])/pow(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35, 2))/(2*Time+1);

   values[22*(Time+1)+80*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.888888888888889*pow(hstep, 2)/pow(0.25*tanh(-0.0333*Xval2[0] - 0.8991) + 0.35, 2) + 2)/(2*Time+1);

   values[0*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.0868055555555556*pow(hstep, 2)*pow(-2.63*pow(tanh(0.0526*Xval[0] + 1.578), 2) + 2.63, 2) - 0.876666666666667*hstep*(-0.1052*pow(tanh(0.0526*Xval[0] + 1.578), 2) + 0.1052)*(0.166666666666667*hstep*(-100.0*Xval[3] - 400.0*Xval2[3] - 100.0*Xvalp1[3] + 50.0*tanh(0.0526*Xval[0] + 1.578) + 200.0*tanh(0.0526*Xval2[0] + 1.578) + 50.0*tanh(0.0526*Xvalp1[0] + 1.578) + 300.0) + Xval[3] - Xvalp1[3])*tanh(0.0526*Xval[0] + 1.578) - 0.6575*hstep*(-0.1052*pow(tanh(0.0526*Xval[0] + 1.578), 2) + 0.1052)*(0.125*hstep*(-100.0*Xval[3] + 100.0*Xvalp1[3] + 50.0*tanh(0.0526*Xval[0] + 1.578) - 50.0*tanh(0.0526*Xvalp1[0] + 1.578)) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3])*tanh(0.0526*Xval[0] + 1.578))/(2*Time+1);

   values[1*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.0555555555555556*pow(hstep, 2)*(-2.63*pow(tanh(0.0526*Xval[0] + 1.578), 2) + 2.63)*(-2.63*pow(tanh(0.0526*Xvalp1[0] + 1.578), 2) + 2.63) + 0.03125*pow(hstep, 2)*(-2.63*pow(tanh(0.0526*Xval[0] + 1.578), 2) + 2.63)*(2.63*pow(tanh(0.0526*Xvalp1[0] + 1.578), 2) - 2.63))/(2*Time+1);

   values[0*(Time+1)+0*(Time)+0+ 1+jt] += obj_factor*bounds[3][3]*(0.0555555555555556*pow(hstep, 2)*pow(-2.63*pow(tanh(0.0526*Xvalp1[0] + 1.578), 2) + 2.63, 2) + 0.03125*pow(hstep, 2)*pow(2.63*pow(tanh(0.0526*Xvalp1[0] + 1.578), 2) - 2.63, 2) - 0.876666666666667*hstep*(-0.1052*pow(tanh(0.0526*Xvalp1[0] + 1.578), 2) + 0.1052)*(0.166666666666667*hstep*(-100.0*Xval[3] - 400.0*Xval2[3] - 100.0*Xvalp1[3] + 50.0*tanh(0.0526*Xval[0] + 1.578) + 200.0*tanh(0.0526*Xval2[0] + 1.578) + 50.0*tanh(0.0526*Xvalp1[0] + 1.578) + 300.0) + Xval[3] - Xvalp1[3])*tanh(0.0526*Xvalp1[0] + 1.578) + 0.6575*hstep*(-0.1052*pow(tanh(0.0526*Xvalp1[0] + 1.578), 2) + 0.1052)*(0.125*hstep*(-100.0*Xval[3] + 100.0*Xvalp1[3] + 50.0*tanh(0.0526*Xval[0] + 1.578) - 50.0*tanh(0.0526*Xvalp1[0] + 1.578)) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3])*tanh(0.0526*Xvalp1[0] + 1.578))/(2*Time+1);

   values[6*(Time+1)+9*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.166666666666667*hstep*(-33.3333333333333*hstep + 2)*(-2.63*pow(tanh(0.0526*Xval[0] + 1.578), 2) + 2.63) + 0.125*hstep*(-25.0*hstep + 1.0)*(-2.63*pow(tanh(0.0526*Xval[0] + 1.578), 2) + 2.63))/(2*Time+1);

   values[7*(Time+1)+9*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.166666666666667*hstep*(-33.3333333333333*hstep + 2)*(-2.63*pow(tanh(0.0526*Xvalp1[0] + 1.578), 2) + 2.63) + 0.125*hstep*(-25.0*hstep + 1.0)*(2.63*pow(tanh(0.0526*Xvalp1[0] + 1.578), 2) - 2.63))/(2*Time+1);

   values[9*(Time+1)+12*(Time)+0+jt] += obj_factor*bounds[3][3]*((-33.3333333333333*hstep + 2)*(-16.6666666666667*hstep + 1) + (-25.0*hstep + 1.0)*(-12.5*hstep + 0.5))/(2*Time+1);

   values[10*(Time+1)+12*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.166666666666667*hstep*(-33.3333333333333*hstep - 2)*(-2.63*pow(tanh(0.0526*Xval[0] + 1.578), 2) + 2.63) + 0.125*hstep*(25.0*hstep + 1.0)*(-2.63*pow(tanh(0.0526*Xval[0] + 1.578), 2) + 2.63))/(2*Time+1);

   values[6*(Time+1)+9*(Time)+0+ 1+jt] += obj_factor*bounds[3][3]*(0.166666666666667*hstep*(-33.3333333333333*hstep - 2)*(-2.63*pow(tanh(0.0526*Xvalp1[0] + 1.578), 2) + 2.63) + 0.125*hstep*(25.0*hstep + 1.0)*(2.63*pow(tanh(0.0526*Xvalp1[0] + 1.578), 2) - 2.63))/(2*Time+1);

   values[10*(Time+1)+15*(Time)+0+jt] += obj_factor*bounds[3][3]*((-33.3333333333333*hstep - 2)*(-16.6666666666667*hstep + 1) + (-12.5*hstep + 0.5)*(25.0*hstep + 1.0))/(2*Time+1);

   values[9*(Time+1)+12*(Time)+0+ 1+jt] += obj_factor*bounds[3][3]*((-33.3333333333333*hstep - 2)*(-16.6666666666667*hstep - 1) + (12.5*hstep + 0.5)*(25.0*hstep + 1.0))/(2*Time+1);

   values[22*(Time+1)+37*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.0555555555555556*pow(hstep, 2)*(-2.63*pow(tanh(0.0526*Xval[0] + 1.578), 2) + 2.63)*(-10.52*pow(tanh(0.0526*Xval2[0] + 1.578), 2) + 10.52))/(2*Time+1);

   values[22*(Time+1)+38*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.0555555555555556*pow(hstep, 2)*(-10.52*pow(tanh(0.0526*Xval2[0] + 1.578), 2) + 10.52)*(-2.63*pow(tanh(0.0526*Xvalp1[0] + 1.578), 2) + 2.63))/(2*Time+1);

   values[22*(Time+1)+43*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.333333333333333*hstep*(-16.6666666666667*hstep + 1)*(-10.52*pow(tanh(0.0526*Xval2[0] + 1.578), 2) + 10.52))/(2*Time+1);

   values[22*(Time+1)+44*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.333333333333333*hstep*(-16.6666666666667*hstep - 1)*(-10.52*pow(tanh(0.0526*Xval2[0] + 1.578), 2) + 10.52))/(2*Time+1);

   values[22*(Time+1)+51*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.0555555555555556*pow(hstep, 2)*pow(-10.52*pow(tanh(0.0526*Xval2[0] + 1.578), 2) + 10.52, 2) - 3.50666666666667*hstep*(-0.1052*pow(tanh(0.0526*Xval2[0] + 1.578), 2) + 0.1052)*(0.166666666666667*hstep*(-100.0*Xval[3] - 400.0*Xval2[3] - 100.0*Xvalp1[3] + 50.0*tanh(0.0526*Xval[0] + 1.578) + 200.0*tanh(0.0526*Xval2[0] + 1.578) + 50.0*tanh(0.0526*Xvalp1[0] + 1.578) + 300.0) + Xval[3] - Xvalp1[3])*tanh(0.0526*Xval2[0] + 1.578))/(2*Time+1);

   values[22*(Time+1)+81*(Time)+0+jt] += obj_factor*bounds[3][3]*(-22.2222222222222*pow(hstep, 2)*(-2.63*pow(tanh(0.0526*Xval[0] + 1.578), 2) + 2.63) - 0.25*hstep*(-2.63*pow(tanh(0.0526*Xval[0] + 1.578), 2) + 2.63))/(2*Time+1);

   values[22*(Time+1)+82*(Time)+0+jt] += obj_factor*bounds[3][3]*(-22.2222222222222*pow(hstep, 2)*(-2.63*pow(tanh(0.0526*Xvalp1[0] + 1.578), 2) + 2.63) - 0.25*hstep*(2.63*pow(tanh(0.0526*Xvalp1[0] + 1.578), 2) - 2.63))/(2*Time+1);

   values[22*(Time+1)+87*(Time)+0+jt] += obj_factor*bounds[3][3]*(-133.333333333333*hstep*(-16.6666666666667*hstep + 1) + 25.0*hstep - 1.0)/(2*Time+1);

   values[22*(Time+1)+88*(Time)+0+jt] += obj_factor*bounds[3][3]*(-133.333333333333*hstep*(-16.6666666666667*hstep - 1) - 25.0*hstep - 1.0)/(2*Time+1);

   values[22*(Time+1)+91*(Time)+0+jt] += obj_factor*bounds[3][3]*(-22.2222222222222*pow(hstep, 2)*(-10.52*pow(tanh(0.0526*Xval2[0] + 1.578), 2) + 10.52))/(2*Time+1);

   values[22*(Time+1)+94*(Time)+0+jt] += obj_factor*bounds[3][3]*(8888.88888888889*pow(hstep, 2) + 2)/(2*Time+1);

   values[0*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.0868055555555556*pow(hstep, 2)*pow((0.0357145*pow(tanh(-0.071429*Xval[0] - 3.214305), 2) - 0.0357145)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + (-0.031249875*pow(tanh(-0.083333*Xval[0] - 3.3749865), 2) + 0.031249875)*(-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475, 2), 2) + 0.333333333333333*hstep*(Xval[4] - Xvalp1[4] + 0.166666666666667*hstep*((-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-4*Xval2[4] + 2.0*tanh(-0.071429*Xval2[0] - 3.214305) + 2.0)/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) + (-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475)))*((0.005102102041*pow(tanh(-0.071429*Xval[0] - 3.214305), 2) - 0.005102102041)*tanh(-0.071429*Xval[0] - 3.214305)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + (-0.06249975*pow(tanh(-0.083333*Xval[0] - 3.3749865), 2) + 0.06249975)*(0.0357145*pow(tanh(-0.071429*Xval[0] - 3.214305), 2) - 0.0357145)/pow(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475, 2) - (0.00520829166675*pow(tanh(-0.083333*Xval[0] - 3.3749865), 2) - 0.00520829166675)*(-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)*tanh(-0.083333*Xval[0] - 3.3749865)/pow(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475, 2) + (-0.06249975*pow(tanh(-0.083333*Xval[0] - 3.3749865), 2) + 0.06249975)*(-0.031249875*pow(tanh(-0.083333*Xval[0] - 3.3749865), 2) + 0.031249875)*(-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475, 3)) + 0.25*hstep*(0.5*Xval[4] - Xval2[4] + 0.5*Xvalp1[4] + 0.125*hstep*((Xvalp1[4] - 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) - 0.5)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475)))*((0.005102102041*pow(tanh(-0.071429*Xval[0] - 3.214305), 2) - 0.005102102041)*tanh(-0.071429*Xval[0] - 3.214305)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + (-0.06249975*pow(tanh(-0.083333*Xval[0] - 3.3749865), 2) + 0.06249975)*(0.0357145*pow(tanh(-0.071429*Xval[0] - 3.214305), 2) - 0.0357145)/pow(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475, 2) - (0.00520829166675*pow(tanh(-0.083333*Xval[0] - 3.3749865), 2) - 0.00520829166675)*(-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)*tanh(-0.083333*Xval[0] - 3.3749865)/pow(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475, 2) + (-0.06249975*pow(tanh(-0.083333*Xval[0] - 3.3749865), 2) + 0.06249975)*(-0.031249875*pow(tanh(-0.083333*Xval[0] - 3.3749865), 2) + 0.031249875)*(-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475, 3)))/(2*Time+1);

   values[1*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.03125*pow(hstep, 2)*((0.0357145*pow(tanh(-0.071429*Xval[0] - 3.214305), 2) - 0.0357145)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + (-0.031249875*pow(tanh(-0.083333*Xval[0] - 3.3749865), 2) + 0.031249875)*(-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475, 2))*((-0.0357145*pow(tanh(-0.071429*Xvalp1[0] - 3.214305), 2) + 0.0357145)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) - (-0.031249875*pow(tanh(-0.083333*Xvalp1[0] - 3.3749865), 2) + 0.031249875)*(-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475, 2)) + 0.0555555555555556*pow(hstep, 2)*((0.0357145*pow(tanh(-0.071429*Xval[0] - 3.214305), 2) - 0.0357145)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + (-0.031249875*pow(tanh(-0.083333*Xval[0] - 3.3749865), 2) + 0.031249875)*(-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475, 2))*((0.0357145*pow(tanh(-0.071429*Xvalp1[0] - 3.214305), 2) - 0.0357145)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-0.031249875*pow(tanh(-0.083333*Xvalp1[0] - 3.3749865), 2) + 0.031249875)*(-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475, 2)))/(2*Time+1);

   values[0*(Time+1)+0*(Time)+0+ 1+jt] += obj_factor*bounds[4][3]*(0.03125*pow(hstep, 2)*pow((-0.0357145*pow(tanh(-0.071429*Xvalp1[0] - 3.214305), 2) + 0.0357145)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) - (-0.031249875*pow(tanh(-0.083333*Xvalp1[0] - 3.3749865), 2) + 0.031249875)*(-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475, 2), 2) + 0.0555555555555556*pow(hstep, 2)*pow((0.0357145*pow(tanh(-0.071429*Xvalp1[0] - 3.214305), 2) - 0.0357145)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-0.031249875*pow(tanh(-0.083333*Xvalp1[0] - 3.3749865), 2) + 0.031249875)*(-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475, 2), 2) + 0.333333333333333*hstep*(Xval[4] - Xvalp1[4] + 0.166666666666667*hstep*((-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-4*Xval2[4] + 2.0*tanh(-0.071429*Xval2[0] - 3.214305) + 2.0)/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) + (-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475)))*((0.005102102041*pow(tanh(-0.071429*Xvalp1[0] - 3.214305), 2) - 0.005102102041)*tanh(-0.071429*Xvalp1[0] - 3.214305)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-0.06249975*pow(tanh(-0.083333*Xvalp1[0] - 3.3749865), 2) + 0.06249975)*(0.0357145*pow(tanh(-0.071429*Xvalp1[0] - 3.214305), 2) - 0.0357145)/pow(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475, 2) - (0.00520829166675*pow(tanh(-0.083333*Xvalp1[0] - 3.3749865), 2) - 0.00520829166675)*(-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)*tanh(-0.083333*Xvalp1[0] - 3.3749865)/pow(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475, 2) + (-0.06249975*pow(tanh(-0.083333*Xvalp1[0] - 3.3749865), 2) + 0.06249975)*(-0.031249875*pow(tanh(-0.083333*Xvalp1[0] - 3.3749865), 2) + 0.031249875)*(-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475, 3)) + 0.25*hstep*(0.5*Xval[4] - Xval2[4] + 0.5*Xvalp1[4] + 0.125*hstep*((Xvalp1[4] - 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) - 0.5)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475)))*((-0.005102102041*pow(tanh(-0.071429*Xvalp1[0] - 3.214305), 2) + 0.005102102041)*tanh(-0.071429*Xvalp1[0] - 3.214305)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) - (-0.06249975*pow(tanh(-0.083333*Xvalp1[0] - 3.3749865), 2) + 0.06249975)*(0.0357145*pow(tanh(-0.071429*Xvalp1[0] - 3.214305), 2) - 0.0357145)/pow(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475, 2) + (0.00520829166675*pow(tanh(-0.083333*Xvalp1[0] - 3.3749865), 2) - 0.00520829166675)*(-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)*tanh(-0.083333*Xvalp1[0] - 3.3749865)/pow(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475, 2) - (-0.06249975*pow(tanh(-0.083333*Xvalp1[0] - 3.3749865), 2) + 0.06249975)*(-0.031249875*pow(tanh(-0.083333*Xvalp1[0] - 3.3749865), 2) + 0.031249875)*(-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475, 3)))/(2*Time+1);

   values[10*(Time+1)+16*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + 2)*((0.0357145*pow(tanh(-0.071429*Xval[0] - 3.214305), 2) - 0.0357145)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + (-0.031249875*pow(tanh(-0.083333*Xval[0] - 3.3749865), 2) + 0.031249875)*(-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475, 2)) + 0.125*hstep*(-0.25*hstep/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + 1.0)*((0.0357145*pow(tanh(-0.071429*Xval[0] - 3.214305), 2) - 0.0357145)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + (-0.031249875*pow(tanh(-0.083333*Xval[0] - 3.3749865), 2) + 0.031249875)*(-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475, 2)) - 0.333333333333333*hstep*(-0.031249875*pow(tanh(-0.083333*Xval[0] - 3.3749865), 2) + 0.031249875)*(Xval[4] - Xvalp1[4] + 0.166666666666667*hstep*((-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-4*Xval2[4] + 2.0*tanh(-0.071429*Xval2[0] - 3.214305) + 2.0)/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) + (-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475)))/pow(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475, 2) - 0.25*hstep*(-0.031249875*pow(tanh(-0.083333*Xval[0] - 3.3749865), 2) + 0.031249875)*(0.5*Xval[4] - Xval2[4] + 0.5*Xvalp1[4] + 0.125*hstep*((Xvalp1[4] - 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) - 0.5)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475)))/pow(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475, 2))/(2*Time+1);

   values[11*(Time+1)+16*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + 2)*((0.0357145*pow(tanh(-0.071429*Xvalp1[0] - 3.214305), 2) - 0.0357145)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-0.031249875*pow(tanh(-0.083333*Xvalp1[0] - 3.3749865), 2) + 0.031249875)*(-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475, 2)) + 0.125*hstep*(-0.25*hstep/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + 1.0)*((-0.0357145*pow(tanh(-0.071429*Xvalp1[0] - 3.214305), 2) + 0.0357145)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) - (-0.031249875*pow(tanh(-0.083333*Xvalp1[0] - 3.3749865), 2) + 0.031249875)*(-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475, 2)))/(2*Time+1);

   values[14*(Time+1)+20*(Time)+0+jt] += obj_factor*bounds[4][3]*((-0.333333333333333*hstep/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + 2)*(-0.166666666666667*hstep/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + 1) + (-0.25*hstep/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + 1.0)*(-0.125*hstep/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + 0.5))/(2*Time+1);

   values[15*(Time+1)+20*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) - 2)*((0.0357145*pow(tanh(-0.071429*Xval[0] - 3.214305), 2) - 0.0357145)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + (-0.031249875*pow(tanh(-0.083333*Xval[0] - 3.3749865), 2) + 0.031249875)*(-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475, 2)) + 0.125*hstep*(0.25*hstep/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + 1.0)*((0.0357145*pow(tanh(-0.071429*Xval[0] - 3.214305), 2) - 0.0357145)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + (-0.031249875*pow(tanh(-0.083333*Xval[0] - 3.3749865), 2) + 0.031249875)*(-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475, 2)))/(2*Time+1);

   values[10*(Time+1)+16*(Time)+0+ 1+jt] += obj_factor*bounds[4][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) - 2)*((0.0357145*pow(tanh(-0.071429*Xvalp1[0] - 3.214305), 2) - 0.0357145)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-0.031249875*pow(tanh(-0.083333*Xvalp1[0] - 3.3749865), 2) + 0.031249875)*(-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475, 2)) + 0.125*hstep*(0.25*hstep/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + 1.0)*((-0.0357145*pow(tanh(-0.071429*Xvalp1[0] - 3.214305), 2) + 0.0357145)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) - (-0.031249875*pow(tanh(-0.083333*Xvalp1[0] - 3.3749865), 2) + 0.031249875)*(-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475, 2)) - 0.333333333333333*hstep*(-0.031249875*pow(tanh(-0.083333*Xvalp1[0] - 3.3749865), 2) + 0.031249875)*(Xval[4] - Xvalp1[4] + 0.166666666666667*hstep*((-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-4*Xval2[4] + 2.0*tanh(-0.071429*Xval2[0] - 3.214305) + 2.0)/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) + (-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475)))/pow(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475, 2) + 0.25*hstep*(-0.031249875*pow(tanh(-0.083333*Xvalp1[0] - 3.3749865), 2) + 0.031249875)*(0.5*Xval[4] - Xval2[4] + 0.5*Xvalp1[4] + 0.125*hstep*((Xvalp1[4] - 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) - 0.5)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475)))/pow(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475, 2))/(2*Time+1);

   values[15*(Time+1)+24*(Time)+0+jt] += obj_factor*bounds[4][3]*((-0.166666666666667*hstep/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + 1)*(-0.333333333333333*hstep/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) - 2) + (-0.125*hstep/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + 0.5)*(0.25*hstep/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + 1.0))/(2*Time+1);

   values[14*(Time+1)+20*(Time)+0+ 1+jt] += obj_factor*bounds[4][3]*((-0.333333333333333*hstep/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) - 2)*(-0.166666666666667*hstep/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) - 1) + (0.125*hstep/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + 0.5)*(0.25*hstep/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + 1.0))/(2*Time+1);

   values[22*(Time+1)+37*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.0555555555555556*pow(hstep, 2)*((0.0357145*pow(tanh(-0.071429*Xval[0] - 3.214305), 2) - 0.0357145)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + (-0.031249875*pow(tanh(-0.083333*Xval[0] - 3.3749865), 2) + 0.031249875)*(-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475, 2))*((0.142858*pow(tanh(-0.071429*Xval2[0] - 3.214305), 2) - 0.142858)/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) + (-0.1249995*pow(tanh(-0.083333*Xval2[0] - 3.3749865), 2) + 0.1249995)*(-Xval2[4] + 0.5*tanh(-0.071429*Xval2[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475, 2)))/(2*Time+1);

   values[22*(Time+1)+38*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.0555555555555556*pow(hstep, 2)*((0.142858*pow(tanh(-0.071429*Xval2[0] - 3.214305), 2) - 0.142858)/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) + (-0.1249995*pow(tanh(-0.083333*Xval2[0] - 3.3749865), 2) + 0.1249995)*(-Xval2[4] + 0.5*tanh(-0.071429*Xval2[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475, 2))*((0.0357145*pow(tanh(-0.071429*Xvalp1[0] - 3.214305), 2) - 0.0357145)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-0.031249875*pow(tanh(-0.083333*Xvalp1[0] - 3.3749865), 2) + 0.031249875)*(-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475, 2)))/(2*Time+1);

   values[22*(Time+1)+45*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.333333333333333*hstep*(-0.166666666666667*hstep/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + 1)*((0.142858*pow(tanh(-0.071429*Xval2[0] - 3.214305), 2) - 0.142858)/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) + (-0.1249995*pow(tanh(-0.083333*Xval2[0] - 3.3749865), 2) + 0.1249995)*(-Xval2[4] + 0.5*tanh(-0.071429*Xval2[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475, 2)))/(2*Time+1);

   values[22*(Time+1)+46*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.333333333333333*hstep*(-0.166666666666667*hstep/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) - 1)*((0.142858*pow(tanh(-0.071429*Xval2[0] - 3.214305), 2) - 0.142858)/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) + (-0.1249995*pow(tanh(-0.083333*Xval2[0] - 3.3749865), 2) + 0.1249995)*(-Xval2[4] + 0.5*tanh(-0.071429*Xval2[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475, 2)))/(2*Time+1);

   values[22*(Time+1)+51*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.0555555555555556*pow(hstep, 2)*pow((0.142858*pow(tanh(-0.071429*Xval2[0] - 3.214305), 2) - 0.142858)/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) + (-0.1249995*pow(tanh(-0.083333*Xval2[0] - 3.3749865), 2) + 0.1249995)*(-Xval2[4] + 0.5*tanh(-0.071429*Xval2[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475, 2), 2) + 0.333333333333333*hstep*(Xval[4] - Xvalp1[4] + 0.166666666666667*hstep*((-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-4*Xval2[4] + 2.0*tanh(-0.071429*Xval2[0] - 3.214305) + 2.0)/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) + (-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475)))*((0.020408408164*pow(tanh(-0.071429*Xval2[0] - 3.214305), 2) - 0.020408408164)*tanh(-0.071429*Xval2[0] - 3.214305)/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) + (-0.249999*pow(tanh(-0.083333*Xval2[0] - 3.3749865), 2) + 0.249999)*(0.0357145*pow(tanh(-0.071429*Xval2[0] - 3.214305), 2) - 0.0357145)/pow(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475, 2) - (0.020833166667*pow(tanh(-0.083333*Xval2[0] - 3.3749865), 2) - 0.020833166667)*(-Xval2[4] + 0.5*tanh(-0.071429*Xval2[0] - 3.214305) + 0.5)*tanh(-0.083333*Xval2[0] - 3.3749865)/pow(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475, 2) + (-0.249999*pow(tanh(-0.083333*Xval2[0] - 3.3749865), 2) + 0.249999)*(-0.031249875*pow(tanh(-0.083333*Xval2[0] - 3.3749865), 2) + 0.031249875)*(-Xval2[4] + 0.5*tanh(-0.071429*Xval2[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475, 3)))/(2*Time+1);

   values[22*(Time+1)+95*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.222222222222222*pow(hstep, 2)*((0.0357145*pow(tanh(-0.071429*Xval[0] - 3.214305), 2) - 0.0357145)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + (-0.031249875*pow(tanh(-0.083333*Xval[0] - 3.3749865), 2) + 0.031249875)*(-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475, 2))/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) - 0.25*hstep*((0.0357145*pow(tanh(-0.071429*Xval[0] - 3.214305), 2) - 0.0357145)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + (-0.031249875*pow(tanh(-0.083333*Xval[0] - 3.3749865), 2) + 0.031249875)*(-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475, 2)))/(2*Time+1);

   values[22*(Time+1)+96*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.222222222222222*pow(hstep, 2)*((0.0357145*pow(tanh(-0.071429*Xvalp1[0] - 3.214305), 2) - 0.0357145)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-0.031249875*pow(tanh(-0.083333*Xvalp1[0] - 3.3749865), 2) + 0.031249875)*(-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475, 2))/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) - 0.25*hstep*((-0.0357145*pow(tanh(-0.071429*Xvalp1[0] - 3.214305), 2) + 0.0357145)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) - (-0.031249875*pow(tanh(-0.083333*Xvalp1[0] - 3.3749865), 2) + 0.031249875)*(-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475, 2)))/(2*Time+1);

   values[22*(Time+1)+103*(Time)+0+jt] += obj_factor*bounds[4][3]*(-1.33333333333333*hstep*(-0.166666666666667*hstep/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) + 1)/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) + 0.25*hstep/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475) - 1.0)/(2*Time+1);

   values[22*(Time+1)+104*(Time)+0+jt] += obj_factor*bounds[4][3]*(-1.33333333333333*hstep*(-0.166666666666667*hstep/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) - 1)/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) - 0.25*hstep/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) - 1.0)/(2*Time+1);

   values[22*(Time+1)+105*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.222222222222222*pow(hstep, 2)*((0.142858*pow(tanh(-0.071429*Xval2[0] - 3.214305), 2) - 0.142858)/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) + (-0.1249995*pow(tanh(-0.083333*Xval2[0] - 3.3749865), 2) + 0.1249995)*(-Xval2[4] + 0.5*tanh(-0.071429*Xval2[0] - 3.214305) + 0.5)/pow(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475, 2))/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) - 1.33333333333333*hstep*(-0.031249875*pow(tanh(-0.083333*Xval2[0] - 3.3749865), 2) + 0.031249875)*(Xval[4] - Xvalp1[4] + 0.166666666666667*hstep*((-Xvalp1[4] + 0.5*tanh(-0.071429*Xvalp1[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xvalp1[0] - 3.3749865) + 0.475) + (-4*Xval2[4] + 2.0*tanh(-0.071429*Xval2[0] - 3.214305) + 2.0)/(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475) + (-Xval[4] + 0.5*tanh(-0.071429*Xval[0] - 3.214305) + 0.5)/(0.375*tanh(-0.083333*Xval[0] - 3.3749865) + 0.475)))/pow(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475, 2))/(2*Time+1);

   values[22*(Time+1)+109*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.888888888888889*pow(hstep, 2)/pow(0.375*tanh(-0.083333*Xval2[0] - 3.3749865) + 0.475, 2) + 2)/(2*Time+1);

   values[2*(Time+1)+2*(Time)+0+jt] += obj_factor*bounds[5][3]*(0.0868055555555556*pow(hstep, 2)*pow(-0.025*pow(tanh(0.05*Xval[1] + 2.0), 2) + 0.025, 2) - 0.00833333333333333*hstep*(-0.1*pow(tanh(0.05*Xval[1] + 2.0), 2) + 0.1)*(0.166666666666667*hstep*(-Xval[5] - 4*Xval2[5] - Xvalp1[5] + 0.5*tanh(0.05*Xval[1] + 2.0) + 2.0*tanh(0.05*Xval2[1] + 2.0) + 0.5*tanh(0.05*Xvalp1[1] + 2.0) + 3.0) + Xval[5] - Xvalp1[5])*tanh(0.05*Xval[1] + 2.0) - 0.00625*hstep*(-0.1*pow(tanh(0.05*Xval[1] + 2.0), 2) + 0.1)*(0.125*hstep*(-Xval[5] + Xvalp1[5] + 0.5*tanh(0.05*Xval[1] + 2.0) - 0.5*tanh(0.05*Xvalp1[1] + 2.0)) + 0.5*Xval[5] - Xval2[5] + 0.5*Xvalp1[5])*tanh(0.05*Xval[1] + 2.0))/(2*Time+1);

   values[3*(Time+1)+3*(Time)+0+jt] += obj_factor*bounds[5][3]*(0.0555555555555556*pow(hstep, 2)*(-0.025*pow(tanh(0.05*Xval[1] + 2.0), 2) + 0.025)*(-0.025*pow(tanh(0.05*Xvalp1[1] + 2.0), 2) + 0.025) + 0.03125*pow(hstep, 2)*(-0.025*pow(tanh(0.05*Xval[1] + 2.0), 2) + 0.025)*(0.025*pow(tanh(0.05*Xvalp1[1] + 2.0), 2) - 0.025))/(2*Time+1);

   values[2*(Time+1)+2*(Time)+0+ 1+jt] += obj_factor*bounds[5][3]*(0.0555555555555556*pow(hstep, 2)*pow(-0.025*pow(tanh(0.05*Xvalp1[1] + 2.0), 2) + 0.025, 2) + 0.03125*pow(hstep, 2)*pow(0.025*pow(tanh(0.05*Xvalp1[1] + 2.0), 2) - 0.025, 2) - 0.00833333333333333*hstep*(-0.1*pow(tanh(0.05*Xvalp1[1] + 2.0), 2) + 0.1)*(0.166666666666667*hstep*(-Xval[5] - 4*Xval2[5] - Xvalp1[5] + 0.5*tanh(0.05*Xval[1] + 2.0) + 2.0*tanh(0.05*Xval2[1] + 2.0) + 0.5*tanh(0.05*Xvalp1[1] + 2.0) + 3.0) + Xval[5] - Xvalp1[5])*tanh(0.05*Xvalp1[1] + 2.0) + 0.00625*hstep*(-0.1*pow(tanh(0.05*Xvalp1[1] + 2.0), 2) + 0.1)*(0.125*hstep*(-Xval[5] + Xvalp1[5] + 0.5*tanh(0.05*Xval[1] + 2.0) - 0.5*tanh(0.05*Xvalp1[1] + 2.0)) + 0.5*Xval[5] - Xval2[5] + 0.5*Xvalp1[5])*tanh(0.05*Xvalp1[1] + 2.0))/(2*Time+1);

   values[16*(Time+1)+26*(Time)+0+jt] += obj_factor*bounds[5][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep + 2)*(-0.025*pow(tanh(0.05*Xval[1] + 2.0), 2) + 0.025) + 0.125*hstep*(-0.25*hstep + 1.0)*(-0.025*pow(tanh(0.05*Xval[1] + 2.0), 2) + 0.025))/(2*Time+1);

   values[17*(Time+1)+26*(Time)+0+jt] += obj_factor*bounds[5][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep + 2)*(-0.025*pow(tanh(0.05*Xvalp1[1] + 2.0), 2) + 0.025) + 0.125*hstep*(-0.25*hstep + 1.0)*(0.025*pow(tanh(0.05*Xvalp1[1] + 2.0), 2) - 0.025))/(2*Time+1);

   values[17*(Time+1)+27*(Time)+0+jt] += obj_factor*bounds[5][3]*((-0.333333333333333*hstep + 2)*(-0.166666666666667*hstep + 1) + (-0.25*hstep + 1.0)*(-0.125*hstep + 0.5))/(2*Time+1);

   values[18*(Time+1)+28*(Time)+0+jt] += obj_factor*bounds[5][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep - 2)*(-0.025*pow(tanh(0.05*Xval[1] + 2.0), 2) + 0.025) + 0.125*hstep*(0.25*hstep + 1.0)*(-0.025*pow(tanh(0.05*Xval[1] + 2.0), 2) + 0.025))/(2*Time+1);

   values[16*(Time+1)+26*(Time)+0+ 1+jt] += obj_factor*bounds[5][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep - 2)*(-0.025*pow(tanh(0.05*Xvalp1[1] + 2.0), 2) + 0.025) + 0.125*hstep*(0.25*hstep + 1.0)*(0.025*pow(tanh(0.05*Xvalp1[1] + 2.0), 2) - 0.025))/(2*Time+1);

   values[18*(Time+1)+29*(Time)+0+jt] += obj_factor*bounds[5][3]*((-0.333333333333333*hstep - 2)*(-0.166666666666667*hstep + 1) + (-0.125*hstep + 0.5)*(0.25*hstep + 1.0))/(2*Time+1);

   values[17*(Time+1)+27*(Time)+0+ 1+jt] += obj_factor*bounds[5][3]*((-0.333333333333333*hstep - 2)*(-0.166666666666667*hstep - 1) + (0.125*hstep + 0.5)*(0.25*hstep + 1.0))/(2*Time+1);

   values[22*(Time+1)+54*(Time)+0+jt] += obj_factor*bounds[5][3]*(0.0555555555555556*pow(hstep, 2)*(-0.025*pow(tanh(0.05*Xval[1] + 2.0), 2) + 0.025)*(-0.1*pow(tanh(0.05*Xval2[1] + 2.0), 2) + 0.1))/(2*Time+1);

   values[22*(Time+1)+55*(Time)+0+jt] += obj_factor*bounds[5][3]*(0.0555555555555556*pow(hstep, 2)*(-0.1*pow(tanh(0.05*Xval2[1] + 2.0), 2) + 0.1)*(-0.025*pow(tanh(0.05*Xvalp1[1] + 2.0), 2) + 0.025))/(2*Time+1);

   values[22*(Time+1)+62*(Time)+0+jt] += obj_factor*bounds[5][3]*(0.333333333333333*hstep*(-0.166666666666667*hstep + 1)*(-0.1*pow(tanh(0.05*Xval2[1] + 2.0), 2) + 0.1))/(2*Time+1);

   values[22*(Time+1)+63*(Time)+0+jt] += obj_factor*bounds[5][3]*(0.333333333333333*hstep*(-0.166666666666667*hstep - 1)*(-0.1*pow(tanh(0.05*Xval2[1] + 2.0), 2) + 0.1))/(2*Time+1);

   values[22*(Time+1)+67*(Time)+0+jt] += obj_factor*bounds[5][3]*(0.0555555555555556*pow(hstep, 2)*pow(-0.1*pow(tanh(0.05*Xval2[1] + 2.0), 2) + 0.1, 2) - 0.0333333333333333*hstep*(-0.1*pow(tanh(0.05*Xval2[1] + 2.0), 2) + 0.1)*(0.166666666666667*hstep*(-Xval[5] - 4*Xval2[5] - Xvalp1[5] + 0.5*tanh(0.05*Xval[1] + 2.0) + 2.0*tanh(0.05*Xval2[1] + 2.0) + 0.5*tanh(0.05*Xvalp1[1] + 2.0) + 3.0) + Xval[5] - Xvalp1[5])*tanh(0.05*Xval2[1] + 2.0))/(2*Time+1);

   values[22*(Time+1)+112*(Time)+0+jt] += obj_factor*bounds[5][3]*(-0.222222222222222*pow(hstep, 2)*(-0.025*pow(tanh(0.05*Xval[1] + 2.0), 2) + 0.025) - 0.25*hstep*(-0.025*pow(tanh(0.05*Xval[1] + 2.0), 2) + 0.025))/(2*Time+1);

   values[22*(Time+1)+113*(Time)+0+jt] += obj_factor*bounds[5][3]*(-0.222222222222222*pow(hstep, 2)*(-0.025*pow(tanh(0.05*Xvalp1[1] + 2.0), 2) + 0.025) - 0.25*hstep*(0.025*pow(tanh(0.05*Xvalp1[1] + 2.0), 2) - 0.025))/(2*Time+1);

   values[22*(Time+1)+114*(Time)+0+jt] += obj_factor*bounds[5][3]*(-1.33333333333333*hstep*(-0.166666666666667*hstep + 1) + 0.25*hstep - 1.0)/(2*Time+1);

   values[22*(Time+1)+115*(Time)+0+jt] += obj_factor*bounds[5][3]*(-1.33333333333333*hstep*(-0.166666666666667*hstep - 1) - 0.25*hstep - 1.0)/(2*Time+1);

   values[22*(Time+1)+119*(Time)+0+jt] += obj_factor*bounds[5][3]*(-0.222222222222222*pow(hstep, 2)*(-0.1*pow(tanh(0.05*Xval2[1] + 2.0), 2) + 0.1))/(2*Time+1);

   values[22*(Time+1)+120*(Time)+0+jt] += obj_factor*bounds[5][3]*(0.888888888888889*pow(hstep, 2) + 2)/(2*Time+1);

   values[2*(Time+1)+2*(Time)+0+jt] += obj_factor*bounds[6][3]*(// Not C:
// ghkjac(Xval[1],1)
// ghkhes(Xval[1],1,0)
0.0868055555555556*pow(hstep, 2)*pow(-1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.000222222222222222*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]), 2) + 0.333333333333333*hstep*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1))*(-1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkhes(Xval[1],1,0) - 1.64609053497942e-5*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) + 0.000444444444444444*pow(Xval[5], 2)*exp(-0.0740740740740741*Xval[1])*ghkjac(Xval[1],1)) + 0.25*hstep*(-1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkhes(Xval[1],1,0) - 1.64609053497942e-5*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) + 0.000444444444444444*pow(Xval[5], 2)*exp(-0.0740740740740741*Xval[1])*ghkjac(Xval[1],1))*(0.5*Xval[6] - Xval2[6] + 0.5*Xvalp1[6] + 0.125*hstep*(-0.0303030303030303*Xval[6] + 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]))))/(2*Time+1);

   values[3*(Time+1)+3*(Time)+0+jt] += obj_factor*bounds[6][3]*(// Not C:
// ghkjac(Xval[1],1)
// ghkjac(Xvalp1[1],1)
0.0555555555555556*pow(hstep, 2)*(-1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.000222222222222222*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]))*(-1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) + 0.000222222222222222*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1])) + 0.03125*pow(hstep, 2)*(-1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.000222222222222222*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]))*(1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) - 0.000222222222222222*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1])))/(2*Time+1);

   values[2*(Time+1)+2*(Time)+0+ 1+jt] += obj_factor*bounds[6][3]*(// Not C:
// ghkjac(Xvalp1[1],1)
// ghkhes(Xvalp1[1],1,0)
0.0555555555555556*pow(hstep, 2)*pow(-1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) + 0.000222222222222222*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]), 2) + 0.03125*pow(hstep, 2)*pow(1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) - 0.000222222222222222*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]), 2) + 0.333333333333333*hstep*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1))*(-1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkhes(Xvalp1[1],1,0) - 1.64609053497942e-5*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) + 0.000444444444444444*pow(Xvalp1[5], 2)*exp(-0.0740740740740741*Xvalp1[1])*ghkjac(Xvalp1[1],1)) + 0.25*hstep*(1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkhes(Xvalp1[1],1,0) + 1.64609053497942e-5*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) - 0.000444444444444444*pow(Xvalp1[5], 2)*exp(-0.0740740740740741*Xvalp1[1])*ghkjac(Xvalp1[1],1))*(0.5*Xval[6] - Xval2[6] + 0.5*Xvalp1[6] + 0.125*hstep*(-0.0303030303030303*Xval[6] + 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]))))/(2*Time+1);

   values[16*(Time+1)+26*(Time)+0+jt] += obj_factor*bounds[6][3]*(// Not C:
// ghkjac(Xval[1],1)
-2.08333333333333e-7*pow(hstep, 2)*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(-1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.000222222222222222*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 8.0e-7*hstep*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1))*ghkjac(Xval[1],1) - 6.0e-7*hstep*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(0.5*Xval[6] - Xval2[6] + 0.5*Xvalp1[6] + 0.125*hstep*(-0.0303030303030303*Xval[6] + 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])))*ghkjac(Xval[1],1) + 0.000148148148148148*hstep*Xval[5]*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1))*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]) + 0.000111111111111111*hstep*Xval[5]*(0.5*Xval[6] - Xval2[6] + 0.5*Xvalp1[6] + 0.125*hstep*(-0.0303030303030303*Xval[6] + 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])))*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]))/(2*Time+1);

   values[17*(Time+1)+26*(Time)+0+jt] += obj_factor*bounds[6][3]*(// Not C:
// ghkjac(Xvalp1[1],1)
-1.33333333333333e-7*pow(hstep, 2)*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(-1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) + 0.000222222222222222*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xval[1]) - 7.5e-8*pow(hstep, 2)*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) - 0.000222222222222222*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xval[1]))/(2*Time+1);

   values[17*(Time+1)+27*(Time)+0+jt] += obj_factor*bounds[6][3]*(5.0e-13*pow(hstep, 2)*pow(Xval[5], 2)*pow(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]), 2)*pow(ghk(Xval[1]), 2) - 8.0e-7*hstep*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1))*ghk(Xval[1]) - 6.0e-7*hstep*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(0.5*Xval[6] - Xval2[6] + 0.5*Xvalp1[6] + 0.125*hstep*(-0.0303030303030303*Xval[6] + 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])))*ghk(Xval[1]))/(2*Time+1);

   values[18*(Time+1)+28*(Time)+0+jt] += obj_factor*bounds[6][3]*(// Not C:
// ghkjac(Xval[1],1)
-5.83333333333333e-8*pow(hstep, 2)*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(-1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.000222222222222222*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]))*ghk(Xvalp1[1]))/(2*Time+1);

   values[16*(Time+1)+26*(Time)+0+ 1+jt] += obj_factor*bounds[6][3]*(// Not C:
// ghkjac(Xvalp1[1],1)
-1.33333333333333e-7*pow(hstep, 2)*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(-1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) + 0.000222222222222222*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 7.5e-8*pow(hstep, 2)*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) - 0.000222222222222222*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) - 8.0e-7*hstep*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1))*ghkjac(Xvalp1[1],1) + 6.0e-7*hstep*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(0.5*Xval[6] - Xval2[6] + 0.5*Xvalp1[6] + 0.125*hstep*(-0.0303030303030303*Xval[6] + 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])))*ghkjac(Xvalp1[1],1) + 0.000148148148148148*hstep*Xvalp1[5]*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1))*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]) - 0.000111111111111111*hstep*Xvalp1[5]*(0.5*Xval[6] - Xval2[6] + 0.5*Xvalp1[6] + 0.125*hstep*(-0.0303030303030303*Xval[6] + 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])))*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]))/(2*Time+1);

   values[18*(Time+1)+29*(Time)+0+jt] += obj_factor*bounds[6][3]*(1.4e-13*pow(hstep, 2)*Xval[5]*Xvalp1[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xval[1])*ghk(Xvalp1[1]))/(2*Time+1);

   values[17*(Time+1)+27*(Time)+0+ 1+jt] += obj_factor*bounds[6][3]*(5.0e-13*pow(hstep, 2)*pow(Xvalp1[5], 2)*pow(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]), 2)*pow(ghk(Xvalp1[1]), 2) - 8.0e-7*hstep*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1))*ghk(Xvalp1[1]) + 6.0e-7*hstep*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(0.5*Xval[6] - Xval2[6] + 0.5*Xvalp1[6] + 0.125*hstep*(-0.0303030303030303*Xval[6] + 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])))*ghk(Xvalp1[1]))/(2*Time+1);

   values[19*(Time+1)+31*(Time)+0+jt] += obj_factor*bounds[6][3]*(// Not C:
// ghkjac(Xval[1],1)
4.0e-7*hstep*pow(Xval[5], 2)*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1))*ghkjac(Xval[1],1) + 3.0e-7*hstep*pow(Xval[5], 2)*(0.5*Xval[6] - Xval2[6] + 0.5*Xvalp1[6] + 0.125*hstep*(-0.0303030303030303*Xval[6] + 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])))*ghkjac(Xval[1],1) + 0.125*hstep*(0.25*hstep*(1.2e-6*pow(Xval[5], 2)*ghk(Xval[1]) - 0.0303030303030303) + 1.0)*(-1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.000222222222222222*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1])) + 0.166666666666667*hstep*(0.333333333333333*hstep*(1.2e-6*pow(Xval[5], 2)*ghk(Xval[1]) - 0.0303030303030303) + 2)*(-1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.000222222222222222*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1])))/(2*Time+1);

   values[20*(Time+1)+31*(Time)+0+jt] += obj_factor*bounds[6][3]*(// Not C:
// ghkjac(Xvalp1[1],1)
0.125*hstep*(0.25*hstep*(1.2e-6*pow(Xval[5], 2)*ghk(Xval[1]) - 0.0303030303030303) + 1.0)*(1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) - 0.000222222222222222*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1])) + 0.166666666666667*hstep*(0.333333333333333*hstep*(1.2e-6*pow(Xval[5], 2)*ghk(Xval[1]) - 0.0303030303030303) + 2)*(-1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) + 0.000222222222222222*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1])))/(2*Time+1);

   values[20*(Time+1)+32*(Time)+0+jt] += obj_factor*bounds[6][3]*(-3.0e-7*hstep*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(0.25*hstep*(1.2e-6*pow(Xval[5], 2)*ghk(Xval[1]) - 0.0303030303030303) + 1.0)*ghk(Xval[1]) - 4.0e-7*hstep*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(0.333333333333333*hstep*(1.2e-6*pow(Xval[5], 2)*ghk(Xval[1]) - 0.0303030303030303) + 2)*ghk(Xval[1]) + 8.0e-7*hstep*Xval[5]*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1))*ghk(Xval[1]) + 6.0e-7*hstep*Xval[5]*(0.5*Xval[6] - Xval2[6] + 0.5*Xvalp1[6] + 0.125*hstep*(-0.0303030303030303*Xval[6] + 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])))*ghk(Xval[1]))/(2*Time+1);

   values[21*(Time+1)+32*(Time)+0+jt] += obj_factor*bounds[6][3]*(3.0e-7*hstep*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(0.25*hstep*(1.2e-6*pow(Xval[5], 2)*ghk(Xval[1]) - 0.0303030303030303) + 1.0)*ghk(Xvalp1[1]) - 4.0e-7*hstep*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(0.333333333333333*hstep*(1.2e-6*pow(Xval[5], 2)*ghk(Xval[1]) - 0.0303030303030303) + 2)*ghk(Xvalp1[1]))/(2*Time+1);

   values[21*(Time+1)+33*(Time)+0+jt] += obj_factor*bounds[6][3]*((0.125*hstep*(1.2e-6*pow(Xval[5], 2)*ghk(Xval[1]) - 0.0303030303030303) + 0.5)*(0.25*hstep*(1.2e-6*pow(Xval[5], 2)*ghk(Xval[1]) - 0.0303030303030303) + 1.0) + (0.166666666666667*hstep*(1.2e-6*pow(Xval[5], 2)*ghk(Xval[1]) - 0.0303030303030303) + 1)*(0.333333333333333*hstep*(1.2e-6*pow(Xval[5], 2)*ghk(Xval[1]) - 0.0303030303030303) + 2))/(2*Time+1);

   values[22*(Time+1)+34*(Time)+0+jt] += obj_factor*bounds[6][3]*(// Not C:
// ghkjac(Xval[1],1)
0.125*hstep*(0.25*hstep*(-1.2e-6*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]) + 0.0303030303030303) + 1.0)*(-1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.000222222222222222*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1])) + 0.166666666666667*hstep*(0.333333333333333*hstep*(1.2e-6*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]) - 0.0303030303030303) - 2)*(-1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.000222222222222222*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1])))/(2*Time+1);

   values[19*(Time+1)+31*(Time)+0+ 1+jt] += obj_factor*bounds[6][3]*(// Not C:
// ghkjac(Xvalp1[1],1)
4.0e-7*hstep*pow(Xvalp1[5], 2)*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1))*ghkjac(Xvalp1[1],1) - 3.0e-7*hstep*pow(Xvalp1[5], 2)*(0.5*Xval[6] - Xval2[6] + 0.5*Xvalp1[6] + 0.125*hstep*(-0.0303030303030303*Xval[6] + 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])))*ghkjac(Xvalp1[1],1) + 0.125*hstep*(0.25*hstep*(-1.2e-6*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]) + 0.0303030303030303) + 1.0)*(1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) - 0.000222222222222222*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1])) + 0.166666666666667*hstep*(0.333333333333333*hstep*(1.2e-6*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]) - 0.0303030303030303) - 2)*(-1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) + 0.000222222222222222*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1])))/(2*Time+1);

   values[22*(Time+1)+35*(Time)+0+jt] += obj_factor*bounds[6][3]*(-3.0e-7*hstep*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(0.25*hstep*(-1.2e-6*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]) + 0.0303030303030303) + 1.0)*ghk(Xval[1]) - 4.0e-7*hstep*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(0.333333333333333*hstep*(1.2e-6*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]) - 0.0303030303030303) - 2)*ghk(Xval[1]))/(2*Time+1);

   values[20*(Time+1)+32*(Time)+0+ 1+jt] += obj_factor*bounds[6][3]*(3.0e-7*hstep*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(0.25*hstep*(-1.2e-6*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]) + 0.0303030303030303) + 1.0)*ghk(Xvalp1[1]) - 4.0e-7*hstep*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(0.333333333333333*hstep*(1.2e-6*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]) - 0.0303030303030303) - 2)*ghk(Xvalp1[1]) + 8.0e-7*hstep*Xvalp1[5]*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1))*ghk(Xvalp1[1]) - 6.0e-7*hstep*Xvalp1[5]*(0.5*Xval[6] - Xval2[6] + 0.5*Xvalp1[6] + 0.125*hstep*(-0.0303030303030303*Xval[6] + 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) + 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1])))*ghk(Xvalp1[1]))/(2*Time+1);

   values[22*(Time+1)+36*(Time)+0+jt] += obj_factor*bounds[6][3]*((0.125*hstep*(1.2e-6*pow(Xval[5], 2)*ghk(Xval[1]) - 0.0303030303030303) + 0.5)*(0.25*hstep*(-1.2e-6*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]) + 0.0303030303030303) + 1.0) + (0.166666666666667*hstep*(1.2e-6*pow(Xval[5], 2)*ghk(Xval[1]) - 0.0303030303030303) + 1)*(0.333333333333333*hstep*(1.2e-6*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]) - 0.0303030303030303) - 2))/(2*Time+1);

   values[21*(Time+1)+33*(Time)+0+ 1+jt] += obj_factor*bounds[6][3]*((0.125*hstep*(-1.2e-6*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]) + 0.0303030303030303) + 0.5)*(0.25*hstep*(-1.2e-6*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]) + 0.0303030303030303) + 1.0) + (0.166666666666667*hstep*(1.2e-6*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]) - 0.0303030303030303) - 1)*(0.333333333333333*hstep*(1.2e-6*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]) - 0.0303030303030303) - 2))/(2*Time+1);

   values[22*(Time+1)+54*(Time)+0+jt] += obj_factor*bounds[6][3]*(// Not C:
// ghkjac(Xval[1],1)
// ghkjac(Xval2[1],1)
0.0555555555555556*pow(hstep, 2)*(-1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.000222222222222222*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]))*(-4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkjac(Xval2[1],1) + 0.000888888888888889*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1])))/(2*Time+1);

   values[22*(Time+1)+55*(Time)+0+jt] += obj_factor*bounds[6][3]*(// Not C:
// ghkjac(Xval2[1],1)
// ghkjac(Xvalp1[1],1)
0.0555555555555556*pow(hstep, 2)*(-4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkjac(Xval2[1],1) + 0.000888888888888889*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1]))*(-1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) + 0.000222222222222222*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1])))/(2*Time+1);

   values[22*(Time+1)+62*(Time)+0+jt] += obj_factor*bounds[6][3]*(// Not C:
// ghkjac(Xval2[1],1)
-1.33333333333333e-7*pow(hstep, 2)*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(-4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkjac(Xval2[1],1) + 0.000888888888888889*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval[1]))/(2*Time+1);

   values[22*(Time+1)+63*(Time)+0+jt] += obj_factor*bounds[6][3]*(// Not C:
// ghkjac(Xval2[1],1)
-1.33333333333333e-7*pow(hstep, 2)*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(-4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkjac(Xval2[1],1) + 0.000888888888888889*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1]))*ghk(Xvalp1[1]))/(2*Time+1);

   values[22*(Time+1)+64*(Time)+0+jt] += obj_factor*bounds[6][3]*(// Not C:
// ghkjac(Xval2[1],1)
0.333333333333333*hstep*(0.166666666666667*hstep*(1.2e-6*pow(Xval[5], 2)*ghk(Xval[1]) - 0.0303030303030303) + 1)*(-4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkjac(Xval2[1],1) + 0.000888888888888889*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1])))/(2*Time+1);

   values[22*(Time+1)+65*(Time)+0+jt] += obj_factor*bounds[6][3]*(// Not C:
// ghkjac(Xval2[1],1)
0.333333333333333*hstep*(0.166666666666667*hstep*(1.2e-6*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]) - 0.0303030303030303) - 1)*(-4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkjac(Xval2[1],1) + 0.000888888888888889*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1])))/(2*Time+1);

   values[22*(Time+1)+67*(Time)+0+jt] += obj_factor*bounds[6][3]*(// Not C:
// ghkjac(Xval2[1],1)
// ghkhes(Xval2[1],1,0)
0.0555555555555556*pow(hstep, 2)*pow(-4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkjac(Xval2[1],1) + 0.000888888888888889*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1]), 2) + 0.333333333333333*hstep*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1))*(-4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkhes(Xval2[1],1,0) - 6.5843621399177e-5*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1]) + 0.00177777777777778*pow(Xval2[5], 2)*exp(-0.0740740740740741*Xval2[1])*ghkjac(Xval2[1],1)))/(2*Time+1);

   values[22*(Time+1)+112*(Time)+0+jt] += obj_factor*bounds[6][3]*(// Not C:
// ghkjac(Xval[1],1)
-5.33333333333333e-7*pow(hstep, 2)*Xval2[5]*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*(-1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.000222222222222222*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1]))*ghk(Xval2[1]))/(2*Time+1);

   values[22*(Time+1)+113*(Time)+0+jt] += obj_factor*bounds[6][3]*(// Not C:
// ghkjac(Xvalp1[1],1)
-5.33333333333333e-7*pow(hstep, 2)*Xval2[5]*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*(-1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) + 0.000222222222222222*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xval2[1]))/(2*Time+1);

   values[22*(Time+1)+114*(Time)+0+jt] += obj_factor*bounds[6][3]*(1.28e-12*pow(hstep, 2)*Xval[5]*Xval2[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval[1])*ghk(Xval2[1]))/(2*Time+1);

   values[22*(Time+1)+115*(Time)+0+jt] += obj_factor*bounds[6][3]*(1.28e-12*pow(hstep, 2)*Xval2[5]*Xvalp1[5]*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xval2[1])*ghk(Xvalp1[1]))/(2*Time+1);

   values[22*(Time+1)+116*(Time)+0+jt] += obj_factor*bounds[6][3]*(-3.2e-6*hstep*Xval2[5]*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*(0.166666666666667*hstep*(1.2e-6*pow(Xval[5], 2)*ghk(Xval[1]) - 0.0303030303030303) + 1)*ghk(Xval2[1]))/(2*Time+1);

   values[22*(Time+1)+117*(Time)+0+jt] += obj_factor*bounds[6][3]*(-3.2e-6*hstep*Xval2[5]*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*(0.166666666666667*hstep*(1.2e-6*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]) - 0.0303030303030303) - 1)*ghk(Xval2[1]))/(2*Time+1);

   values[22*(Time+1)+119*(Time)+0+jt] += obj_factor*bounds[6][3]*(// Not C:
// ghkjac(Xval2[1],1)
-5.33333333333333e-7*pow(hstep, 2)*Xval2[5]*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*(-4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkjac(Xval2[1],1) + 0.000888888888888889*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 3.2e-6*hstep*Xval2[5]*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1))*ghkjac(Xval2[1],1) + 0.000592592592592593*hstep*Xval2[5]*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1))*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1]))/(2*Time+1);

   values[22*(Time+1)+120*(Time)+0+jt] += obj_factor*bounds[6][3]*(5.12e-12*pow(hstep, 2)*pow(Xval2[5], 2)*pow(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]), 2)*pow(ghk(Xval2[1]), 2) - 3.2e-6*hstep*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1))*ghk(Xval2[1]))/(2*Time+1);

   values[22*(Time+1)+123*(Time)+0+jt] += obj_factor*bounds[6][3]*(// Not C:
// ghkjac(Xval[1],1)
0.0555555555555556*pow(hstep, 2)*(4.8e-6*pow(Xval2[5], 2)*ghk(Xval2[1]) - 0.121212121212121)*(-1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.000222222222222222*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1])) - 0.25*hstep*(-1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghkjac(Xval[1],1) + 0.000222222222222222*pow(Xval[5], 2)*ghk(Xval[1])*exp(-0.0740740740740741*Xval[1])))/(2*Time+1);

   values[22*(Time+1)+124*(Time)+0+jt] += obj_factor*bounds[6][3]*(// Not C:
// ghkjac(Xvalp1[1],1)
0.0555555555555556*pow(hstep, 2)*(4.8e-6*pow(Xval2[5], 2)*ghk(Xval2[1]) - 0.121212121212121)*(-1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) + 0.000222222222222222*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1])) - 0.25*hstep*(1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghkjac(Xvalp1[1],1) - 0.000222222222222222*pow(Xvalp1[5], 2)*ghk(Xvalp1[1])*exp(-0.0740740740740741*Xvalp1[1])))/(2*Time+1);

   values[22*(Time+1)+125*(Time)+0+jt] += obj_factor*bounds[6][3]*(-1.33333333333333e-7*pow(hstep, 2)*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*(4.8e-6*pow(Xval2[5], 2)*ghk(Xval2[1]) - 0.121212121212121)*ghk(Xval[1]) + 6.0e-7*hstep*Xval[5]*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]))/(2*Time+1);

   values[22*(Time+1)+126*(Time)+0+jt] += obj_factor*bounds[6][3]*(-1.33333333333333e-7*pow(hstep, 2)*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*(4.8e-6*pow(Xval2[5], 2)*ghk(Xval2[1]) - 0.121212121212121)*ghk(Xvalp1[1]) - 6.0e-7*hstep*Xvalp1[5]*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]))/(2*Time+1);

   values[22*(Time+1)+127*(Time)+0+jt] += obj_factor*bounds[6][3]*(0.333333333333333*hstep*(0.166666666666667*hstep*(1.2e-6*pow(Xval[5], 2)*ghk(Xval[1]) - 0.0303030303030303) + 1)*(4.8e-6*pow(Xval2[5], 2)*ghk(Xval2[1]) - 0.121212121212121) - 0.25*hstep*(1.2e-6*pow(Xval[5], 2)*ghk(Xval[1]) - 0.0303030303030303) - 1.0)/(2*Time+1);

   values[22*(Time+1)+128*(Time)+0+jt] += obj_factor*bounds[6][3]*(0.333333333333333*hstep*(0.166666666666667*hstep*(1.2e-6*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]) - 0.0303030303030303) - 1)*(4.8e-6*pow(Xval2[5], 2)*ghk(Xval2[1]) - 0.121212121212121) - 0.25*hstep*(-1.2e-6*pow(Xvalp1[5], 2)*ghk(Xvalp1[1]) + 0.0303030303030303) - 1.0)/(2*Time+1);

   values[22*(Time+1)+130*(Time)+0+jt] += obj_factor*bounds[6][3]*(// Not C:
// ghkjac(Xval2[1],1)
0.0555555555555556*pow(hstep, 2)*(4.8e-6*pow(Xval2[5], 2)*ghk(Xval2[1]) - 0.121212121212121)*(-4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghkjac(Xval2[1],1) + 0.000888888888888889*pow(Xval2[5], 2)*ghk(Xval2[1])*exp(-0.0740740740740741*Xval2[1])) + 1.6e-6*hstep*pow(Xval2[5], 2)*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1))*ghkjac(Xval2[1],1))/(2*Time+1);

   values[22*(Time+1)+131*(Time)+0+jt] += obj_factor*bounds[6][3]*(-5.33333333333333e-7*pow(hstep, 2)*Xval2[5]*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*(4.8e-6*pow(Xval2[5], 2)*ghk(Xval2[1]) - 0.121212121212121)*ghk(Xval2[1]) + 3.2e-6*hstep*Xval2[5]*(Xval[6] - Xvalp1[6] + 0.166666666666667*hstep*(-0.0303030303030303*Xval[6] - 0.121212121212121*Xval2[6] - 0.0303030303030303*Xvalp1[6] - 1.2e-6*pow(Xval[5], 2)*(-Xval[6] + 2500*exp(-0.0740740740740741*Xval[1]))*ghk(Xval[1]) - 4.8e-6*pow(Xval2[5], 2)*(-Xval2[6] + 2500*exp(-0.0740740740740741*Xval2[1]))*ghk(Xval2[1]) - 1.2e-6*pow(Xvalp1[5], 2)*(-Xvalp1[6] + 2500*exp(-0.0740740740740741*Xvalp1[1]))*ghk(Xvalp1[1]) + 0.1))*ghk(Xval2[1]))/(2*Time+1);

   values[22*(Time+1)+132*(Time)+0+jt] += obj_factor*bounds[6][3]*(0.0555555555555556*pow(hstep, 2)*pow(4.8e-6*pow(Xval2[5], 2)*ghk(Xval2[1]) - 0.121212121212121, 2) + 2)/(2*Time+1);

   } // end for loop 

// Add last element
     for(Index i=0;i<nY;i++) {
        Xval[i] = x[Time + i*(Time+1)];
        Xvalp1[i] = 0;
        Xval2[i] = 0;
     } //end for loop

     for(Index i=0;i<nU;i++) {
        K11val[i] = x[Time + nY*(Time+1) + 2*i*(Time+1)];
        K11valp1[i] = 0;
        K11val2[i] = 0;
        dK11val[i] = x[Time + (nY+2*i+1)*(Time+1)];
        dK11valp1[i] = 0;
        dK11val2[i] = 0;
     } //end for loop

     Xdval[0] = VDATA0[2*Time];
     Xdval2[0] = 0;
     Xdvalp1[0] = 0;
     Ival[0] = Iinj[2*Time];
     Ival2[0] = 0;
     Ivalp1[0] = 0;

     for(Index i=0;i<nP;i++) {
        Pval[i] = x[(2*Time+1)*(nY+2*nU)+i];
     } //end for loop


   values[0*(Time+1)+0*(Time)+0+ Time] += obj_factor*(2)/(2*Time+1);

  } // end else 

   return true;
}


void SIMPLE_NAKL_NLP::finalize_solution(SolverReturn status,
                        Index n, const Number* x, const Number* z_L, const Number* z_U,
                        Index m, const Number* g, const Number* lambda,
                        Number obj_value,
                        const IpoptData* ip_data,
                        IpoptCalculatedQuantities* ip_cq)
{
  // here is where the solution is written to file


    FILE *OUTPUT1;
    char filename[20];
	sprintf(filename,"D%d_M%d_PATH%d.dat", nY,nM,taskid);
  	if(beta==0){
		OUTPUT1 = fopen (filename,"w");
	}
	else
		OUTPUT1 = fopen (filename,"a");

	fprintf(OUTPUT1, "%d %d %e ",beta, (status == SUCCESS), obj_value);
	for(Index i=0;i<Ntotal;i++){
		solution[i] = x[i];
	}
	for (Index i=0;i<Time;i++) {
     	for (Index j=0;j<nY;j++) {
        	fprintf(OUTPUT1,"%e ", x[j*(Time+1)+i]);
        }
      	for (Index j=0;j<nY;j++) {
        	fprintf(OUTPUT1,"%e ", x[(nY+2*nU)*(Time+1)+j*Time+i]);
		}
    }
  	for (Index j=0;j<nY;j++) {
     	fprintf(OUTPUT1,"%e ", x[j*(Time+1)+Time]);
    }
  	for (Index j=0;j<nP;j++) {
     	fprintf(OUTPUT1,"%e ", x[(2*Time+1)*(nY+2*nU)+j]);
    }
	fprintf(OUTPUT1,"\n");
	fclose (OUTPUT1);
	

}
