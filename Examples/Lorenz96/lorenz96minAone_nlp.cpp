// lorenz96.cpp
// Nonlinear Ipopt program

// Author: Bryan A. Toth
// btoth@physics.ucsd.edu

#include "lorenz96minAone_nlp.hpp"
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

// constructor
LORENZ96_NLP::LORENZ96_NLP(int id)
{
  nU=0;
  nP=1;
  nY=5;
  nM=1;
  nI=0;

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
    data1 = new double[2*Time+1];
    data1dummy = new double[skip];

    Ntotal = (2*Time+1)*nY+nP;
    solution = new double[Ntotal];
    FILE *pFile0;
    filename = specs[3];
    pFile0 = fopen(filename.c_str(),"r");

    for(Index jt=0;jt<skip;jt++)
	{
	ret = fscanf (pFile0, "%lf", &data1dummy[jt]);
	if (ret == EOF) break;
	}
    for(Index jt=0;jt<2*Time+1;jt++)
	{
	ret = fscanf (pFile0, "%lf", &data1[jt]);
	if (ret == EOF) break;
	}
    fclose (pFile0);
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
LORENZ96_NLP::~LORENZ96_NLP()
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
    delete [] data1;
    delete [] data1dummy;
  int rows = nY+2*nU+nP;
  for (Index i=0;i<rows;i++) delete [] bounds[i];
  delete [] bounds;

}


bool LORENZ96_NLP::changeRf(){
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
bool LORENZ96_NLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
				Index& nnz_h_lag, IndexStyleEnum& index_style)

{
  // Number of variables
  n = 10*Time+6;

  // Number of equality constraints
  m = 0;

  // Number of Jacobian nonzero entries
  nnz_jac_g = 0;

  // Number of Hessian nonzero entries
  nnz_h_lag = 20*(Time+1)+95*Time+1;

  // use the C style indexing (0-based)
  index_style = TNLP::C_STYLE;

  return true;
}


// returns the variable bounds
bool LORENZ96_NLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{
  // Here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.
  assert(n == 10*Time+6);
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
bool LORENZ96_NLP::get_starting_point(Index n, bool init_x, Number* x,
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


if(beta==0){

      int skipROWS = skip;
      int ROWS = 2*Time+1;
      int COLS = nY;
    
      double **skipinit = new double* [skipROWS];
      double **init = new double* [ROWS];
      for(Index i=0;i<skipROWS;i++) skipinit[i] = new double[COLS];
      for(Index i=0;i<ROWS;i++) init[i] = new double[COLS];
      
      string filename;
      filename = specs[4+nM+nI];

      if (specs[3+nM+nI] =="1")
      {
      FILE *initFILE;
      int ret;
      initFILE = fopen(filename.c_str(),"r");
    
      for(Index jt=0;jt<skip;jt++)
          {
	  ret = fscanf (initFILE,"%lf %lf %lf %lf %lf ",&skipinit[jt][0],&skipinit[jt][1],&skipinit[jt][2],&skipinit[jt][3],&skipinit[jt][4]);
	  if (ret == EOF) break;
          }
      for(Index jt=0;jt<2*Time+1;jt++)
          {
	  ret = fscanf (initFILE,"%lf %lf %lf %lf %lf ",&init[jt][0],&init[jt][1],&init[jt][2],&init[jt][3],&init[jt][4]);
	  if (ret == EOF) break;
          }
    fclose (initFILE);
    }

  for(Index jt=0;jt<Time+1;jt++) {
     for(Index var=0;var<nY;var++) {
       // Initial conditions for x
       if (specs[3+nM+nI] == "1")
                  {
		  x[(Time+1)*var+jt] = init[2*jt][var];
		  }
	        else
	          {
		  for(int i=0; i<taskid;i++) x[(Time+1)*var+jt] = rand()*1.0/RAND_MAX*(bounds[var][1]-bounds[var][0])+bounds[var][0];
		  }
       // Initial conditions for midpoints
       if(jt<Time) {
                  if (specs[3+nM+nI] == "1")
		    {
		    x[(Time+1)*(nY+2*nU)+Time*var+jt] = init[2*jt+1][var];
		    }
		  else
		    {
		    for(int i=0; i<taskid;i++) x[(Time+1)*(nY+2*nU)+Time*var+jt] = rand()*1.0/RAND_MAX*(bounds[var][1]-bounds[var][0])+bounds[var][0];
		    }
		  }
		}
     for(Index cup=0;cup<2*nU;cup++) {
       // Initial conditions for k
       x[(Time+1)*(cup+nY)+jt]=bounds[cup+nY][2];
       // Initial conditions for midpoints
       if(jt<Time) {
                  x[(Time+1)*(nY+2*nU)+Time*(cup+nY)+jt]=bounds[cup+nY][2];
	          }
	      }
  } // End for loop

     for(Index par=0;par<nP;par++) {
     // Initial conditions for p5
     if (specs[3+nM+nI] == "1"){
					x[2*Time*(nY+2*nU)+nY+2*nU+par]=bounds[nY+2*nU+par][2];
			  }
			  else{
			  		 for(int i=0; i<taskid;i++) x[2*Time*(nY+2*nU)+nY+2*nU+par]=rand()*1.0/RAND_MAX*(bounds[nY+2*nU+par][1]-bounds[nY+2*nU+par][0])+bounds[nY+2*nU+par][0];
			  }
              }

  for(Index i=0;i<ROWS;i++) delete [] init[i];
   delete [] init;
}
else{
	for(Index i=0;i<Ntotal;i++) x[i] = solution[i];
	}
  return true;
}


// returns the value of the objective function
bool LORENZ96_NLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  assert(n == 10*Time+6);
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

     Xdval[0] = data1[2*jt];
     Xdval2[0] = data1[2*jt+1];
     Xdvalp1[0] = data1[2*jt+2];

     for(Index i=0;i<nP;i++) {
        Pval[i] = x[(2*Time+1)*(nY+2*nU)+i];
     } //end for loop


    obj_value += (Xdval[0] - Xval[0])*(4*Xdval[0] - 4*Xval[0]) + (Xdval2[0] - Xval2[0])*(4*Xdval2[0] - 4*Xval2[0]); 

    obj_value += bounds[0][3]*(pow(0.166666666666667*hstep*(6*Pval[0] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) + 4*Xval2[4]*(Xval2[1] - Xval2[3]) + Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + Xval[0] - Xvalp1[0], 2) + pow(0.125*hstep*(-Xval[0] + Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) - Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0], 2)); 

    obj_value += bounds[1][3]*(pow(0.166666666666667*hstep*(6*Pval[0] + Xval[0]*(Xval[2] - Xval[4]) + 4*Xval2[0]*(Xval2[2] - Xval2[4]) + Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] - 4*Xval2[1] - Xvalp1[1]) + Xval[1] - Xvalp1[1], 2) + pow(0.125*hstep*(Xval[0]*(Xval[2] - Xval[4]) - Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] + Xvalp1[1]) + 0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1], 2)); 

    obj_value += bounds[2][3]*(pow(0.166666666666667*hstep*(6*Pval[0] + Xval[1]*(-Xval[0] + Xval[3]) + 4*Xval2[1]*(-Xval2[0] + Xval2[3]) + Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] - 4*Xval2[2] - Xvalp1[2]) + Xval[2] - Xvalp1[2], 2) + pow(0.125*hstep*(Xval[1]*(-Xval[0] + Xval[3]) - Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] + Xvalp1[2]) + 0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2], 2)); 

    obj_value += bounds[3][3]*(pow(0.166666666666667*hstep*(6*Pval[0] + Xval[2]*(-Xval[1] + Xval[4]) + 4*Xval2[2]*(-Xval2[1] + Xval2[4]) + Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] - 4*Xval2[3] - Xvalp1[3]) + Xval[3] - Xvalp1[3], 2) + pow(0.125*hstep*(Xval[2]*(-Xval[1] + Xval[4]) - Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] + Xvalp1[3]) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3], 2)); 

    obj_value += bounds[4][3]*(pow(0.166666666666667*hstep*(6*Pval[0] + Xval[3]*(Xval[0] - Xval[2]) + 4*Xval2[3]*(Xval2[0] - Xval2[2]) + Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] - 4*Xval2[4] - Xvalp1[4]) + Xval[4] - Xvalp1[4], 2) + pow(0.125*hstep*(Xval[3]*(Xval[0] - Xval[2]) - Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] + Xvalp1[4]) + 0.5*Xval[4] - Xval2[4] + 0.5*Xvalp1[4], 2)); 

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

     Xdval[0] = data1[2*Time];
     Xdval2[0] = 0;
     Xdvalp1[0] = 0;

     for(Index i=0;i<nP;i++) {
        Pval[i] = x[(2*Time+1)*(nY+2*nU)+i];
     } //end for loop


  obj_value += (Xdval[0] - Xval[0])*(4*Xdval[0] - 4*Xval[0]) + (Xdval2[0] - Xval2[0])*(4*Xdval2[0] - 4*Xval2[0]);

  obj_value = obj_value/(2*Time+1);

  return true;
}


// return the gradient of the objective function grad_{x} f(x)
bool LORENZ96_NLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  assert(n == 10*Time+6);

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

     Xdval[0] = data1[2*jt];
     Xdval2[0] = data1[2*jt+1];
     Xdvalp1[0] = data1[2*jt+2];

     for(Index i=0;i<nP;i++) {
        Pval[i] = x[(2*Time+1)*(nY+2*nU)+i];
     } //end for loop

    grad_f[jt+0*(Time+1)] += (-8*Xdval[0] + 8*Xval[0])/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 0*Time + jt] += (-8*Xdval2[0] + 8*Xval2[0])/(2*Time+1);

    grad_f[jt+0*(Time+1)] += bounds[0][3]*((-0.333333333333333*hstep + 2)*(0.166666666666667*hstep*(6*Pval[0] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) + 4*Xval2[4]*(Xval2[1] - Xval2[3]) + Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + Xval[0] - Xvalp1[0]) + (-0.25*hstep + 1.0)*(0.125*hstep*(-Xval[0] + Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) - Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
    grad_f[jt+1*(Time+1)] += bounds[0][3]*(0.333333333333333*hstep*Xval[4]*(0.166666666666667*hstep*(6*Pval[0] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) + 4*Xval2[4]*(Xval2[1] - Xval2[3]) + Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + Xval[0] - Xvalp1[0]) + 0.25*hstep*Xval[4]*(0.125*hstep*(-Xval[0] + Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) - Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
    grad_f[jt+3*(Time+1)] += bounds[0][3]*(-0.333333333333333*hstep*Xval[4]*(0.166666666666667*hstep*(6*Pval[0] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) + 4*Xval2[4]*(Xval2[1] - Xval2[3]) + Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + Xval[0] - Xvalp1[0]) - 0.25*hstep*Xval[4]*(0.125*hstep*(-Xval[0] + Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) - Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
    grad_f[jt+4*(Time+1)] += bounds[0][3]*(0.333333333333333*hstep*(Xval[1] - Xval[3])*(0.166666666666667*hstep*(6*Pval[0] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) + 4*Xval2[4]*(Xval2[1] - Xval2[3]) + Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + Xval[0] - Xvalp1[0]) + 0.25*hstep*(Xval[1] - Xval[3])*(0.125*hstep*(-Xval[0] + Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) - Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
    grad_f[jt+1+0*(Time+1)] += bounds[0][3]*((-0.333333333333333*hstep - 2)*(0.166666666666667*hstep*(6*Pval[0] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) + 4*Xval2[4]*(Xval2[1] - Xval2[3]) + Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + Xval[0] - Xvalp1[0]) + (0.25*hstep + 1.0)*(0.125*hstep*(-Xval[0] + Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) - Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
    grad_f[jt+1+1*(Time+1)] += bounds[0][3]*(0.333333333333333*hstep*Xvalp1[4]*(0.166666666666667*hstep*(6*Pval[0] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) + 4*Xval2[4]*(Xval2[1] - Xval2[3]) + Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + Xval[0] - Xvalp1[0]) - 0.25*hstep*Xvalp1[4]*(0.125*hstep*(-Xval[0] + Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) - Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
    grad_f[jt+1+3*(Time+1)] += bounds[0][3]*(-0.333333333333333*hstep*Xvalp1[4]*(0.166666666666667*hstep*(6*Pval[0] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) + 4*Xval2[4]*(Xval2[1] - Xval2[3]) + Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + Xval[0] - Xvalp1[0]) + 0.25*hstep*Xvalp1[4]*(0.125*hstep*(-Xval[0] + Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) - Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
    grad_f[jt+1+4*(Time+1)] += bounds[0][3]*(0.25*hstep*(-Xvalp1[1] + Xvalp1[3])*(0.125*hstep*(-Xval[0] + Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) - Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]) + 0.333333333333333*hstep*(Xvalp1[1] - Xvalp1[3])*(0.166666666666667*hstep*(6*Pval[0] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) + 4*Xval2[4]*(Xval2[1] - Xval2[3]) + Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 0*Time + jt] += bounds[0][3]*(-1.33333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) + 4*Xval2[4]*(Xval2[1] - Xval2[3]) + Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + Xval[0] - Xvalp1[0]) - 0.25*hstep*(-Xval[0] + Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) - Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) - 1.0*Xval[0] + 2*Xval2[0] - 1.0*Xvalp1[0])/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 1*Time + jt] += bounds[0][3]*(1.33333333333333*hstep*Xval2[4]*(0.166666666666667*hstep*(6*Pval[0] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) + 4*Xval2[4]*(Xval2[1] - Xval2[3]) + Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 3*Time + jt] += bounds[0][3]*(-1.33333333333333*hstep*Xval2[4]*(0.166666666666667*hstep*(6*Pval[0] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) + 4*Xval2[4]*(Xval2[1] - Xval2[3]) + Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 4*Time + jt] += bounds[0][3]*(0.333333333333333*hstep*(4*Xval2[1] - 4*Xval2[3])*(0.166666666666667*hstep*(6*Pval[0] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) + 4*Xval2[4]*(Xval2[1] - Xval2[3]) + Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
     grad_f[(2*Time+1)*(nY+2*nU)+0] += bounds[0][3]*(2.0*hstep*(0.166666666666667*hstep*(6*Pval[0] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) + 4*Xval2[4]*(Xval2[1] - Xval2[3]) + Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
    grad_f[jt+0*(Time+1)] += bounds[1][3]*(0.333333333333333*hstep*(Xval[2] - Xval[4])*(0.166666666666667*hstep*(6*Pval[0] + Xval[0]*(Xval[2] - Xval[4]) + 4*Xval2[0]*(Xval2[2] - Xval2[4]) + Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] - 4*Xval2[1] - Xvalp1[1]) + Xval[1] - Xvalp1[1]) + 0.25*hstep*(Xval[2] - Xval[4])*(0.125*hstep*(Xval[0]*(Xval[2] - Xval[4]) - Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] + Xvalp1[1]) + 0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1]))/(2*Time+1);
    grad_f[jt+1*(Time+1)] += bounds[1][3]*((-0.333333333333333*hstep + 2)*(0.166666666666667*hstep*(6*Pval[0] + Xval[0]*(Xval[2] - Xval[4]) + 4*Xval2[0]*(Xval2[2] - Xval2[4]) + Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] - 4*Xval2[1] - Xvalp1[1]) + Xval[1] - Xvalp1[1]) + (-0.25*hstep + 1.0)*(0.125*hstep*(Xval[0]*(Xval[2] - Xval[4]) - Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] + Xvalp1[1]) + 0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1]))/(2*Time+1);
    grad_f[jt+2*(Time+1)] += bounds[1][3]*(0.333333333333333*hstep*Xval[0]*(0.166666666666667*hstep*(6*Pval[0] + Xval[0]*(Xval[2] - Xval[4]) + 4*Xval2[0]*(Xval2[2] - Xval2[4]) + Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] - 4*Xval2[1] - Xvalp1[1]) + Xval[1] - Xvalp1[1]) + 0.25*hstep*Xval[0]*(0.125*hstep*(Xval[0]*(Xval[2] - Xval[4]) - Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] + Xvalp1[1]) + 0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1]))/(2*Time+1);
    grad_f[jt+4*(Time+1)] += bounds[1][3]*(-0.333333333333333*hstep*Xval[0]*(0.166666666666667*hstep*(6*Pval[0] + Xval[0]*(Xval[2] - Xval[4]) + 4*Xval2[0]*(Xval2[2] - Xval2[4]) + Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] - 4*Xval2[1] - Xvalp1[1]) + Xval[1] - Xvalp1[1]) - 0.25*hstep*Xval[0]*(0.125*hstep*(Xval[0]*(Xval[2] - Xval[4]) - Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] + Xvalp1[1]) + 0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1]))/(2*Time+1);
    grad_f[jt+1+0*(Time+1)] += bounds[1][3]*(0.25*hstep*(-Xvalp1[2] + Xvalp1[4])*(0.125*hstep*(Xval[0]*(Xval[2] - Xval[4]) - Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] + Xvalp1[1]) + 0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1]) + 0.333333333333333*hstep*(Xvalp1[2] - Xvalp1[4])*(0.166666666666667*hstep*(6*Pval[0] + Xval[0]*(Xval[2] - Xval[4]) + 4*Xval2[0]*(Xval2[2] - Xval2[4]) + Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] - 4*Xval2[1] - Xvalp1[1]) + Xval[1] - Xvalp1[1]))/(2*Time+1);
    grad_f[jt+1+1*(Time+1)] += bounds[1][3]*((-0.333333333333333*hstep - 2)*(0.166666666666667*hstep*(6*Pval[0] + Xval[0]*(Xval[2] - Xval[4]) + 4*Xval2[0]*(Xval2[2] - Xval2[4]) + Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] - 4*Xval2[1] - Xvalp1[1]) + Xval[1] - Xvalp1[1]) + (0.25*hstep + 1.0)*(0.125*hstep*(Xval[0]*(Xval[2] - Xval[4]) - Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] + Xvalp1[1]) + 0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1]))/(2*Time+1);
    grad_f[jt+1+2*(Time+1)] += bounds[1][3]*(0.333333333333333*hstep*Xvalp1[0]*(0.166666666666667*hstep*(6*Pval[0] + Xval[0]*(Xval[2] - Xval[4]) + 4*Xval2[0]*(Xval2[2] - Xval2[4]) + Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] - 4*Xval2[1] - Xvalp1[1]) + Xval[1] - Xvalp1[1]) - 0.25*hstep*Xvalp1[0]*(0.125*hstep*(Xval[0]*(Xval[2] - Xval[4]) - Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] + Xvalp1[1]) + 0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1]))/(2*Time+1);
    grad_f[jt+1+4*(Time+1)] += bounds[1][3]*(-0.333333333333333*hstep*Xvalp1[0]*(0.166666666666667*hstep*(6*Pval[0] + Xval[0]*(Xval[2] - Xval[4]) + 4*Xval2[0]*(Xval2[2] - Xval2[4]) + Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] - 4*Xval2[1] - Xvalp1[1]) + Xval[1] - Xvalp1[1]) + 0.25*hstep*Xvalp1[0]*(0.125*hstep*(Xval[0]*(Xval[2] - Xval[4]) - Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] + Xvalp1[1]) + 0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1]))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 0*Time + jt] += bounds[1][3]*(0.333333333333333*hstep*(4*Xval2[2] - 4*Xval2[4])*(0.166666666666667*hstep*(6*Pval[0] + Xval[0]*(Xval[2] - Xval[4]) + 4*Xval2[0]*(Xval2[2] - Xval2[4]) + Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] - 4*Xval2[1] - Xvalp1[1]) + Xval[1] - Xvalp1[1]))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 1*Time + jt] += bounds[1][3]*(-1.33333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[0]*(Xval[2] - Xval[4]) + 4*Xval2[0]*(Xval2[2] - Xval2[4]) + Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] - 4*Xval2[1] - Xvalp1[1]) + Xval[1] - Xvalp1[1]) - 0.25*hstep*(Xval[0]*(Xval[2] - Xval[4]) - Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] + Xvalp1[1]) - 1.0*Xval[1] + 2*Xval2[1] - 1.0*Xvalp1[1])/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 2*Time + jt] += bounds[1][3]*(1.33333333333333*hstep*Xval2[0]*(0.166666666666667*hstep*(6*Pval[0] + Xval[0]*(Xval[2] - Xval[4]) + 4*Xval2[0]*(Xval2[2] - Xval2[4]) + Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] - 4*Xval2[1] - Xvalp1[1]) + Xval[1] - Xvalp1[1]))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 4*Time + jt] += bounds[1][3]*(-1.33333333333333*hstep*Xval2[0]*(0.166666666666667*hstep*(6*Pval[0] + Xval[0]*(Xval[2] - Xval[4]) + 4*Xval2[0]*(Xval2[2] - Xval2[4]) + Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] - 4*Xval2[1] - Xvalp1[1]) + Xval[1] - Xvalp1[1]))/(2*Time+1);
     grad_f[(2*Time+1)*(nY+2*nU)+0] += bounds[1][3]*(2.0*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[0]*(Xval[2] - Xval[4]) + 4*Xval2[0]*(Xval2[2] - Xval2[4]) + Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] - 4*Xval2[1] - Xvalp1[1]) + Xval[1] - Xvalp1[1]))/(2*Time+1);
    grad_f[jt+0*(Time+1)] += bounds[2][3]*(-0.333333333333333*hstep*Xval[1]*(0.166666666666667*hstep*(6*Pval[0] + Xval[1]*(-Xval[0] + Xval[3]) + 4*Xval2[1]*(-Xval2[0] + Xval2[3]) + Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] - 4*Xval2[2] - Xvalp1[2]) + Xval[2] - Xvalp1[2]) - 0.25*hstep*Xval[1]*(0.125*hstep*(Xval[1]*(-Xval[0] + Xval[3]) - Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] + Xvalp1[2]) + 0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2]))/(2*Time+1);
    grad_f[jt+1*(Time+1)] += bounds[2][3]*(0.333333333333333*hstep*(-Xval[0] + Xval[3])*(0.166666666666667*hstep*(6*Pval[0] + Xval[1]*(-Xval[0] + Xval[3]) + 4*Xval2[1]*(-Xval2[0] + Xval2[3]) + Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] - 4*Xval2[2] - Xvalp1[2]) + Xval[2] - Xvalp1[2]) + 0.25*hstep*(-Xval[0] + Xval[3])*(0.125*hstep*(Xval[1]*(-Xval[0] + Xval[3]) - Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] + Xvalp1[2]) + 0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2]))/(2*Time+1);
    grad_f[jt+2*(Time+1)] += bounds[2][3]*((-0.333333333333333*hstep + 2)*(0.166666666666667*hstep*(6*Pval[0] + Xval[1]*(-Xval[0] + Xval[3]) + 4*Xval2[1]*(-Xval2[0] + Xval2[3]) + Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] - 4*Xval2[2] - Xvalp1[2]) + Xval[2] - Xvalp1[2]) + (-0.25*hstep + 1.0)*(0.125*hstep*(Xval[1]*(-Xval[0] + Xval[3]) - Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] + Xvalp1[2]) + 0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2]))/(2*Time+1);
    grad_f[jt+3*(Time+1)] += bounds[2][3]*(0.333333333333333*hstep*Xval[1]*(0.166666666666667*hstep*(6*Pval[0] + Xval[1]*(-Xval[0] + Xval[3]) + 4*Xval2[1]*(-Xval2[0] + Xval2[3]) + Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] - 4*Xval2[2] - Xvalp1[2]) + Xval[2] - Xvalp1[2]) + 0.25*hstep*Xval[1]*(0.125*hstep*(Xval[1]*(-Xval[0] + Xval[3]) - Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] + Xvalp1[2]) + 0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2]))/(2*Time+1);
    grad_f[jt+1+0*(Time+1)] += bounds[2][3]*(-0.333333333333333*hstep*Xvalp1[1]*(0.166666666666667*hstep*(6*Pval[0] + Xval[1]*(-Xval[0] + Xval[3]) + 4*Xval2[1]*(-Xval2[0] + Xval2[3]) + Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] - 4*Xval2[2] - Xvalp1[2]) + Xval[2] - Xvalp1[2]) + 0.25*hstep*Xvalp1[1]*(0.125*hstep*(Xval[1]*(-Xval[0] + Xval[3]) - Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] + Xvalp1[2]) + 0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2]))/(2*Time+1);
    grad_f[jt+1+1*(Time+1)] += bounds[2][3]*(0.333333333333333*hstep*(-Xvalp1[0] + Xvalp1[3])*(0.166666666666667*hstep*(6*Pval[0] + Xval[1]*(-Xval[0] + Xval[3]) + 4*Xval2[1]*(-Xval2[0] + Xval2[3]) + Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] - 4*Xval2[2] - Xvalp1[2]) + Xval[2] - Xvalp1[2]) + 0.25*hstep*(Xvalp1[0] - Xvalp1[3])*(0.125*hstep*(Xval[1]*(-Xval[0] + Xval[3]) - Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] + Xvalp1[2]) + 0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2]))/(2*Time+1);
    grad_f[jt+1+2*(Time+1)] += bounds[2][3]*((-0.333333333333333*hstep - 2)*(0.166666666666667*hstep*(6*Pval[0] + Xval[1]*(-Xval[0] + Xval[3]) + 4*Xval2[1]*(-Xval2[0] + Xval2[3]) + Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] - 4*Xval2[2] - Xvalp1[2]) + Xval[2] - Xvalp1[2]) + (0.25*hstep + 1.0)*(0.125*hstep*(Xval[1]*(-Xval[0] + Xval[3]) - Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] + Xvalp1[2]) + 0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2]))/(2*Time+1);
    grad_f[jt+1+3*(Time+1)] += bounds[2][3]*(0.333333333333333*hstep*Xvalp1[1]*(0.166666666666667*hstep*(6*Pval[0] + Xval[1]*(-Xval[0] + Xval[3]) + 4*Xval2[1]*(-Xval2[0] + Xval2[3]) + Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] - 4*Xval2[2] - Xvalp1[2]) + Xval[2] - Xvalp1[2]) - 0.25*hstep*Xvalp1[1]*(0.125*hstep*(Xval[1]*(-Xval[0] + Xval[3]) - Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] + Xvalp1[2]) + 0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2]))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 0*Time + jt] += bounds[2][3]*(-1.33333333333333*hstep*Xval2[1]*(0.166666666666667*hstep*(6*Pval[0] + Xval[1]*(-Xval[0] + Xval[3]) + 4*Xval2[1]*(-Xval2[0] + Xval2[3]) + Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] - 4*Xval2[2] - Xvalp1[2]) + Xval[2] - Xvalp1[2]))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 1*Time + jt] += bounds[2][3]*(0.333333333333333*hstep*(-4*Xval2[0] + 4*Xval2[3])*(0.166666666666667*hstep*(6*Pval[0] + Xval[1]*(-Xval[0] + Xval[3]) + 4*Xval2[1]*(-Xval2[0] + Xval2[3]) + Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] - 4*Xval2[2] - Xvalp1[2]) + Xval[2] - Xvalp1[2]))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 2*Time + jt] += bounds[2][3]*(-1.33333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[1]*(-Xval[0] + Xval[3]) + 4*Xval2[1]*(-Xval2[0] + Xval2[3]) + Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] - 4*Xval2[2] - Xvalp1[2]) + Xval[2] - Xvalp1[2]) - 0.25*hstep*(Xval[1]*(-Xval[0] + Xval[3]) - Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] + Xvalp1[2]) - 1.0*Xval[2] + 2*Xval2[2] - 1.0*Xvalp1[2])/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 3*Time + jt] += bounds[2][3]*(1.33333333333333*hstep*Xval2[1]*(0.166666666666667*hstep*(6*Pval[0] + Xval[1]*(-Xval[0] + Xval[3]) + 4*Xval2[1]*(-Xval2[0] + Xval2[3]) + Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] - 4*Xval2[2] - Xvalp1[2]) + Xval[2] - Xvalp1[2]))/(2*Time+1);
     grad_f[(2*Time+1)*(nY+2*nU)+0] += bounds[2][3]*(2.0*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[1]*(-Xval[0] + Xval[3]) + 4*Xval2[1]*(-Xval2[0] + Xval2[3]) + Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] - 4*Xval2[2] - Xvalp1[2]) + Xval[2] - Xvalp1[2]))/(2*Time+1);
    grad_f[jt+1*(Time+1)] += bounds[3][3]*(-0.333333333333333*hstep*Xval[2]*(0.166666666666667*hstep*(6*Pval[0] + Xval[2]*(-Xval[1] + Xval[4]) + 4*Xval2[2]*(-Xval2[1] + Xval2[4]) + Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] - 4*Xval2[3] - Xvalp1[3]) + Xval[3] - Xvalp1[3]) - 0.25*hstep*Xval[2]*(0.125*hstep*(Xval[2]*(-Xval[1] + Xval[4]) - Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] + Xvalp1[3]) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3]))/(2*Time+1);
    grad_f[jt+2*(Time+1)] += bounds[3][3]*(0.333333333333333*hstep*(-Xval[1] + Xval[4])*(0.166666666666667*hstep*(6*Pval[0] + Xval[2]*(-Xval[1] + Xval[4]) + 4*Xval2[2]*(-Xval2[1] + Xval2[4]) + Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] - 4*Xval2[3] - Xvalp1[3]) + Xval[3] - Xvalp1[3]) + 0.25*hstep*(-Xval[1] + Xval[4])*(0.125*hstep*(Xval[2]*(-Xval[1] + Xval[4]) - Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] + Xvalp1[3]) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3]))/(2*Time+1);
    grad_f[jt+3*(Time+1)] += bounds[3][3]*((-0.333333333333333*hstep + 2)*(0.166666666666667*hstep*(6*Pval[0] + Xval[2]*(-Xval[1] + Xval[4]) + 4*Xval2[2]*(-Xval2[1] + Xval2[4]) + Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] - 4*Xval2[3] - Xvalp1[3]) + Xval[3] - Xvalp1[3]) + (-0.25*hstep + 1.0)*(0.125*hstep*(Xval[2]*(-Xval[1] + Xval[4]) - Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] + Xvalp1[3]) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3]))/(2*Time+1);
    grad_f[jt+4*(Time+1)] += bounds[3][3]*(0.333333333333333*hstep*Xval[2]*(0.166666666666667*hstep*(6*Pval[0] + Xval[2]*(-Xval[1] + Xval[4]) + 4*Xval2[2]*(-Xval2[1] + Xval2[4]) + Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] - 4*Xval2[3] - Xvalp1[3]) + Xval[3] - Xvalp1[3]) + 0.25*hstep*Xval[2]*(0.125*hstep*(Xval[2]*(-Xval[1] + Xval[4]) - Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] + Xvalp1[3]) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3]))/(2*Time+1);
    grad_f[jt+1+1*(Time+1)] += bounds[3][3]*(-0.333333333333333*hstep*Xvalp1[2]*(0.166666666666667*hstep*(6*Pval[0] + Xval[2]*(-Xval[1] + Xval[4]) + 4*Xval2[2]*(-Xval2[1] + Xval2[4]) + Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] - 4*Xval2[3] - Xvalp1[3]) + Xval[3] - Xvalp1[3]) + 0.25*hstep*Xvalp1[2]*(0.125*hstep*(Xval[2]*(-Xval[1] + Xval[4]) - Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] + Xvalp1[3]) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3]))/(2*Time+1);
    grad_f[jt+1+2*(Time+1)] += bounds[3][3]*(0.333333333333333*hstep*(-Xvalp1[1] + Xvalp1[4])*(0.166666666666667*hstep*(6*Pval[0] + Xval[2]*(-Xval[1] + Xval[4]) + 4*Xval2[2]*(-Xval2[1] + Xval2[4]) + Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] - 4*Xval2[3] - Xvalp1[3]) + Xval[3] - Xvalp1[3]) + 0.25*hstep*(Xvalp1[1] - Xvalp1[4])*(0.125*hstep*(Xval[2]*(-Xval[1] + Xval[4]) - Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] + Xvalp1[3]) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3]))/(2*Time+1);
    grad_f[jt+1+3*(Time+1)] += bounds[3][3]*((-0.333333333333333*hstep - 2)*(0.166666666666667*hstep*(6*Pval[0] + Xval[2]*(-Xval[1] + Xval[4]) + 4*Xval2[2]*(-Xval2[1] + Xval2[4]) + Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] - 4*Xval2[3] - Xvalp1[3]) + Xval[3] - Xvalp1[3]) + (0.25*hstep + 1.0)*(0.125*hstep*(Xval[2]*(-Xval[1] + Xval[4]) - Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] + Xvalp1[3]) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3]))/(2*Time+1);
    grad_f[jt+1+4*(Time+1)] += bounds[3][3]*(0.333333333333333*hstep*Xvalp1[2]*(0.166666666666667*hstep*(6*Pval[0] + Xval[2]*(-Xval[1] + Xval[4]) + 4*Xval2[2]*(-Xval2[1] + Xval2[4]) + Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] - 4*Xval2[3] - Xvalp1[3]) + Xval[3] - Xvalp1[3]) - 0.25*hstep*Xvalp1[2]*(0.125*hstep*(Xval[2]*(-Xval[1] + Xval[4]) - Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] + Xvalp1[3]) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3]))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 1*Time + jt] += bounds[3][3]*(-1.33333333333333*hstep*Xval2[2]*(0.166666666666667*hstep*(6*Pval[0] + Xval[2]*(-Xval[1] + Xval[4]) + 4*Xval2[2]*(-Xval2[1] + Xval2[4]) + Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] - 4*Xval2[3] - Xvalp1[3]) + Xval[3] - Xvalp1[3]))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 2*Time + jt] += bounds[3][3]*(0.333333333333333*hstep*(-4*Xval2[1] + 4*Xval2[4])*(0.166666666666667*hstep*(6*Pval[0] + Xval[2]*(-Xval[1] + Xval[4]) + 4*Xval2[2]*(-Xval2[1] + Xval2[4]) + Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] - 4*Xval2[3] - Xvalp1[3]) + Xval[3] - Xvalp1[3]))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 3*Time + jt] += bounds[3][3]*(-1.33333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[2]*(-Xval[1] + Xval[4]) + 4*Xval2[2]*(-Xval2[1] + Xval2[4]) + Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] - 4*Xval2[3] - Xvalp1[3]) + Xval[3] - Xvalp1[3]) - 0.25*hstep*(Xval[2]*(-Xval[1] + Xval[4]) - Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] + Xvalp1[3]) - 1.0*Xval[3] + 2*Xval2[3] - 1.0*Xvalp1[3])/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 4*Time + jt] += bounds[3][3]*(1.33333333333333*hstep*Xval2[2]*(0.166666666666667*hstep*(6*Pval[0] + Xval[2]*(-Xval[1] + Xval[4]) + 4*Xval2[2]*(-Xval2[1] + Xval2[4]) + Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] - 4*Xval2[3] - Xvalp1[3]) + Xval[3] - Xvalp1[3]))/(2*Time+1);
     grad_f[(2*Time+1)*(nY+2*nU)+0] += bounds[3][3]*(2.0*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[2]*(-Xval[1] + Xval[4]) + 4*Xval2[2]*(-Xval2[1] + Xval2[4]) + Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] - 4*Xval2[3] - Xvalp1[3]) + Xval[3] - Xvalp1[3]))/(2*Time+1);
    grad_f[jt+0*(Time+1)] += bounds[4][3]*(0.333333333333333*hstep*Xval[3]*(0.166666666666667*hstep*(6*Pval[0] + Xval[3]*(Xval[0] - Xval[2]) + 4*Xval2[3]*(Xval2[0] - Xval2[2]) + Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] - 4*Xval2[4] - Xvalp1[4]) + Xval[4] - Xvalp1[4]) + 0.25*hstep*Xval[3]*(0.125*hstep*(Xval[3]*(Xval[0] - Xval[2]) - Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] + Xvalp1[4]) + 0.5*Xval[4] - Xval2[4] + 0.5*Xvalp1[4]))/(2*Time+1);
    grad_f[jt+2*(Time+1)] += bounds[4][3]*(-0.333333333333333*hstep*Xval[3]*(0.166666666666667*hstep*(6*Pval[0] + Xval[3]*(Xval[0] - Xval[2]) + 4*Xval2[3]*(Xval2[0] - Xval2[2]) + Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] - 4*Xval2[4] - Xvalp1[4]) + Xval[4] - Xvalp1[4]) - 0.25*hstep*Xval[3]*(0.125*hstep*(Xval[3]*(Xval[0] - Xval[2]) - Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] + Xvalp1[4]) + 0.5*Xval[4] - Xval2[4] + 0.5*Xvalp1[4]))/(2*Time+1);
    grad_f[jt+3*(Time+1)] += bounds[4][3]*(0.333333333333333*hstep*(Xval[0] - Xval[2])*(0.166666666666667*hstep*(6*Pval[0] + Xval[3]*(Xval[0] - Xval[2]) + 4*Xval2[3]*(Xval2[0] - Xval2[2]) + Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] - 4*Xval2[4] - Xvalp1[4]) + Xval[4] - Xvalp1[4]) + 0.25*hstep*(Xval[0] - Xval[2])*(0.125*hstep*(Xval[3]*(Xval[0] - Xval[2]) - Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] + Xvalp1[4]) + 0.5*Xval[4] - Xval2[4] + 0.5*Xvalp1[4]))/(2*Time+1);
    grad_f[jt+4*(Time+1)] += bounds[4][3]*((-0.333333333333333*hstep + 2)*(0.166666666666667*hstep*(6*Pval[0] + Xval[3]*(Xval[0] - Xval[2]) + 4*Xval2[3]*(Xval2[0] - Xval2[2]) + Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] - 4*Xval2[4] - Xvalp1[4]) + Xval[4] - Xvalp1[4]) + (-0.25*hstep + 1.0)*(0.125*hstep*(Xval[3]*(Xval[0] - Xval[2]) - Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] + Xvalp1[4]) + 0.5*Xval[4] - Xval2[4] + 0.5*Xvalp1[4]))/(2*Time+1);
    grad_f[jt+1+0*(Time+1)] += bounds[4][3]*(0.333333333333333*hstep*Xvalp1[3]*(0.166666666666667*hstep*(6*Pval[0] + Xval[3]*(Xval[0] - Xval[2]) + 4*Xval2[3]*(Xval2[0] - Xval2[2]) + Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] - 4*Xval2[4] - Xvalp1[4]) + Xval[4] - Xvalp1[4]) - 0.25*hstep*Xvalp1[3]*(0.125*hstep*(Xval[3]*(Xval[0] - Xval[2]) - Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] + Xvalp1[4]) + 0.5*Xval[4] - Xval2[4] + 0.5*Xvalp1[4]))/(2*Time+1);
    grad_f[jt+1+2*(Time+1)] += bounds[4][3]*(-0.333333333333333*hstep*Xvalp1[3]*(0.166666666666667*hstep*(6*Pval[0] + Xval[3]*(Xval[0] - Xval[2]) + 4*Xval2[3]*(Xval2[0] - Xval2[2]) + Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] - 4*Xval2[4] - Xvalp1[4]) + Xval[4] - Xvalp1[4]) + 0.25*hstep*Xvalp1[3]*(0.125*hstep*(Xval[3]*(Xval[0] - Xval[2]) - Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] + Xvalp1[4]) + 0.5*Xval[4] - Xval2[4] + 0.5*Xvalp1[4]))/(2*Time+1);
    grad_f[jt+1+3*(Time+1)] += bounds[4][3]*(0.25*hstep*(-Xvalp1[0] + Xvalp1[2])*(0.125*hstep*(Xval[3]*(Xval[0] - Xval[2]) - Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] + Xvalp1[4]) + 0.5*Xval[4] - Xval2[4] + 0.5*Xvalp1[4]) + 0.333333333333333*hstep*(Xvalp1[0] - Xvalp1[2])*(0.166666666666667*hstep*(6*Pval[0] + Xval[3]*(Xval[0] - Xval[2]) + 4*Xval2[3]*(Xval2[0] - Xval2[2]) + Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] - 4*Xval2[4] - Xvalp1[4]) + Xval[4] - Xvalp1[4]))/(2*Time+1);
    grad_f[jt+1+4*(Time+1)] += bounds[4][3]*((-0.333333333333333*hstep - 2)*(0.166666666666667*hstep*(6*Pval[0] + Xval[3]*(Xval[0] - Xval[2]) + 4*Xval2[3]*(Xval2[0] - Xval2[2]) + Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] - 4*Xval2[4] - Xvalp1[4]) + Xval[4] - Xvalp1[4]) + (0.25*hstep + 1.0)*(0.125*hstep*(Xval[3]*(Xval[0] - Xval[2]) - Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] + Xvalp1[4]) + 0.5*Xval[4] - Xval2[4] + 0.5*Xvalp1[4]))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 0*Time + jt] += bounds[4][3]*(1.33333333333333*hstep*Xval2[3]*(0.166666666666667*hstep*(6*Pval[0] + Xval[3]*(Xval[0] - Xval[2]) + 4*Xval2[3]*(Xval2[0] - Xval2[2]) + Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] - 4*Xval2[4] - Xvalp1[4]) + Xval[4] - Xvalp1[4]))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 2*Time + jt] += bounds[4][3]*(-1.33333333333333*hstep*Xval2[3]*(0.166666666666667*hstep*(6*Pval[0] + Xval[3]*(Xval[0] - Xval[2]) + 4*Xval2[3]*(Xval2[0] - Xval2[2]) + Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] - 4*Xval2[4] - Xvalp1[4]) + Xval[4] - Xvalp1[4]))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 3*Time + jt] += bounds[4][3]*(0.333333333333333*hstep*(4*Xval2[0] - 4*Xval2[2])*(0.166666666666667*hstep*(6*Pval[0] + Xval[3]*(Xval[0] - Xval[2]) + 4*Xval2[3]*(Xval2[0] - Xval2[2]) + Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] - 4*Xval2[4] - Xvalp1[4]) + Xval[4] - Xvalp1[4]))/(2*Time+1);
    grad_f[(Time+1)*(2*nU+nY) + 4*Time + jt] += bounds[4][3]*(-1.33333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[3]*(Xval[0] - Xval[2]) + 4*Xval2[3]*(Xval2[0] - Xval2[2]) + Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] - 4*Xval2[4] - Xvalp1[4]) + Xval[4] - Xvalp1[4]) - 0.25*hstep*(Xval[3]*(Xval[0] - Xval[2]) - Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] + Xvalp1[4]) - 1.0*Xval[4] + 2*Xval2[4] - 1.0*Xvalp1[4])/(2*Time+1);
     grad_f[(2*Time+1)*(nY+2*nU)+0] += bounds[4][3]*(2.0*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[3]*(Xval[0] - Xval[2]) + 4*Xval2[3]*(Xval2[0] - Xval2[2]) + Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] - 4*Xval2[4] - Xvalp1[4]) + Xval[4] - Xvalp1[4]))/(2*Time+1);

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

     Xdval[0] = data1[2*Time];
     Xdval2[0] = 0;
     Xdvalp1[0] = 0;

     for(Index i=0;i<nP;i++) {
        Pval[i] = x[(2*Time+1)*(nY+2*nU)+i];
     } //end for loop


    grad_f[Time+0*(Time+1)] += (-8*Xdval[0] + 8*Xval[0])/(2*Time+1);

  return true;
}


// return the value of the constraints: g(x)
bool LORENZ96_NLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  assert(n == 10*Time+6);
  assert(m == 0);

  return true;
}


// return the structure or values of the jacobian
bool LORENZ96_NLP::eval_jac_g(Index n, const Number* x, bool new_x,
                            Index m, Index nele_jac, Index* iRow, Index* jCol,
                            Number* values)
{

return true;
}


// return the structure or values of the hessian
bool LORENZ96_NLP::eval_h(Index n, const Number* x, bool new_x,
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

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+25*(Time)+0+jt] = (Time+1)*5+Time*0+jt;
     jCol[15*(Time+1)+25*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+26*(Time)+0+jt] = (Time+1)*5+Time*0+jt;
     jCol[15*(Time+1)+26*(Time)+0+jt] = (Time+1)*0+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+27*(Time)+0+jt] = (Time+1)*5+Time*0+jt;
     jCol[15*(Time+1)+27*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+28*(Time)+0+jt] = (Time+1)*5+Time*0+jt;
     jCol[15*(Time+1)+28*(Time)+0+jt] = (Time+1)*1+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+29*(Time)+0+jt] = (Time+1)*5+Time*0+jt;
     jCol[15*(Time+1)+29*(Time)+0+jt] = (Time+1)*2+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+30*(Time)+0+jt] = (Time+1)*5+Time*0+jt;
     jCol[15*(Time+1)+30*(Time)+0+jt] = (Time+1)*2+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+31*(Time)+0+jt] = (Time+1)*5+Time*0+jt;
     jCol[15*(Time+1)+31*(Time)+0+jt] = (Time+1)*3+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+32*(Time)+0+jt] = (Time+1)*5+Time*0+jt;
     jCol[15*(Time+1)+32*(Time)+0+jt] = (Time+1)*3+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+33*(Time)+0+jt] = (Time+1)*5+Time*0+jt;
     jCol[15*(Time+1)+33*(Time)+0+jt] = (Time+1)*4+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+34*(Time)+0+jt] = (Time+1)*5+Time*0+jt;
     jCol[15*(Time+1)+34*(Time)+0+jt] = (Time+1)*4+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+35*(Time)+0+jt] = (Time+1)*5+Time*0+jt;
     jCol[15*(Time+1)+35*(Time)+0+jt] = (Time+1)*5+Time*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+36*(Time)+0+jt] = (Time+1)*5+Time*1+jt;
     jCol[15*(Time+1)+36*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+37*(Time)+0+jt] = (Time+1)*5+Time*1+jt;
     jCol[15*(Time+1)+37*(Time)+0+jt] = (Time+1)*0+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+38*(Time)+0+jt] = (Time+1)*5+Time*1+jt;
     jCol[15*(Time+1)+38*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+39*(Time)+0+jt] = (Time+1)*5+Time*1+jt;
     jCol[15*(Time+1)+39*(Time)+0+jt] = (Time+1)*1+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+40*(Time)+0+jt] = (Time+1)*5+Time*1+jt;
     jCol[15*(Time+1)+40*(Time)+0+jt] = (Time+1)*2+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+41*(Time)+0+jt] = (Time+1)*5+Time*1+jt;
     jCol[15*(Time+1)+41*(Time)+0+jt] = (Time+1)*2+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+42*(Time)+0+jt] = (Time+1)*5+Time*1+jt;
     jCol[15*(Time+1)+42*(Time)+0+jt] = (Time+1)*3+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+43*(Time)+0+jt] = (Time+1)*5+Time*1+jt;
     jCol[15*(Time+1)+43*(Time)+0+jt] = (Time+1)*3+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+44*(Time)+0+jt] = (Time+1)*5+Time*1+jt;
     jCol[15*(Time+1)+44*(Time)+0+jt] = (Time+1)*4+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+45*(Time)+0+jt] = (Time+1)*5+Time*1+jt;
     jCol[15*(Time+1)+45*(Time)+0+jt] = (Time+1)*4+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+46*(Time)+0+jt] = (Time+1)*5+Time*1+jt;
     jCol[15*(Time+1)+46*(Time)+0+jt] = (Time+1)*5+Time*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+47*(Time)+0+jt] = (Time+1)*5+Time*1+jt;
     jCol[15*(Time+1)+47*(Time)+0+jt] = (Time+1)*5+Time*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+48*(Time)+0+jt] = (Time+1)*5+Time*2+jt;
     jCol[15*(Time+1)+48*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+49*(Time)+0+jt] = (Time+1)*5+Time*2+jt;
     jCol[15*(Time+1)+49*(Time)+0+jt] = (Time+1)*0+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+50*(Time)+0+jt] = (Time+1)*5+Time*2+jt;
     jCol[15*(Time+1)+50*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+51*(Time)+0+jt] = (Time+1)*5+Time*2+jt;
     jCol[15*(Time+1)+51*(Time)+0+jt] = (Time+1)*1+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+52*(Time)+0+jt] = (Time+1)*5+Time*2+jt;
     jCol[15*(Time+1)+52*(Time)+0+jt] = (Time+1)*2+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+53*(Time)+0+jt] = (Time+1)*5+Time*2+jt;
     jCol[15*(Time+1)+53*(Time)+0+jt] = (Time+1)*2+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+54*(Time)+0+jt] = (Time+1)*5+Time*2+jt;
     jCol[15*(Time+1)+54*(Time)+0+jt] = (Time+1)*3+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+55*(Time)+0+jt] = (Time+1)*5+Time*2+jt;
     jCol[15*(Time+1)+55*(Time)+0+jt] = (Time+1)*3+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+56*(Time)+0+jt] = (Time+1)*5+Time*2+jt;
     jCol[15*(Time+1)+56*(Time)+0+jt] = (Time+1)*4+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+57*(Time)+0+jt] = (Time+1)*5+Time*2+jt;
     jCol[15*(Time+1)+57*(Time)+0+jt] = (Time+1)*4+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+58*(Time)+0+jt] = (Time+1)*5+Time*2+jt;
     jCol[15*(Time+1)+58*(Time)+0+jt] = (Time+1)*5+Time*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+59*(Time)+0+jt] = (Time+1)*5+Time*2+jt;
     jCol[15*(Time+1)+59*(Time)+0+jt] = (Time+1)*5+Time*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+60*(Time)+0+jt] = (Time+1)*5+Time*2+jt;
     jCol[15*(Time+1)+60*(Time)+0+jt] = (Time+1)*5+Time*2+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+61*(Time)+0+jt] = (Time+1)*5+Time*3+jt;
     jCol[15*(Time+1)+61*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+62*(Time)+0+jt] = (Time+1)*5+Time*3+jt;
     jCol[15*(Time+1)+62*(Time)+0+jt] = (Time+1)*0+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+63*(Time)+0+jt] = (Time+1)*5+Time*3+jt;
     jCol[15*(Time+1)+63*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+64*(Time)+0+jt] = (Time+1)*5+Time*3+jt;
     jCol[15*(Time+1)+64*(Time)+0+jt] = (Time+1)*1+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+65*(Time)+0+jt] = (Time+1)*5+Time*3+jt;
     jCol[15*(Time+1)+65*(Time)+0+jt] = (Time+1)*2+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+66*(Time)+0+jt] = (Time+1)*5+Time*3+jt;
     jCol[15*(Time+1)+66*(Time)+0+jt] = (Time+1)*2+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+67*(Time)+0+jt] = (Time+1)*5+Time*3+jt;
     jCol[15*(Time+1)+67*(Time)+0+jt] = (Time+1)*3+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+68*(Time)+0+jt] = (Time+1)*5+Time*3+jt;
     jCol[15*(Time+1)+68*(Time)+0+jt] = (Time+1)*3+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+69*(Time)+0+jt] = (Time+1)*5+Time*3+jt;
     jCol[15*(Time+1)+69*(Time)+0+jt] = (Time+1)*4+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+70*(Time)+0+jt] = (Time+1)*5+Time*3+jt;
     jCol[15*(Time+1)+70*(Time)+0+jt] = (Time+1)*4+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+71*(Time)+0+jt] = (Time+1)*5+Time*3+jt;
     jCol[15*(Time+1)+71*(Time)+0+jt] = (Time+1)*5+Time*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+72*(Time)+0+jt] = (Time+1)*5+Time*3+jt;
     jCol[15*(Time+1)+72*(Time)+0+jt] = (Time+1)*5+Time*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+73*(Time)+0+jt] = (Time+1)*5+Time*3+jt;
     jCol[15*(Time+1)+73*(Time)+0+jt] = (Time+1)*5+Time*2+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+74*(Time)+0+jt] = (Time+1)*5+Time*3+jt;
     jCol[15*(Time+1)+74*(Time)+0+jt] = (Time+1)*5+Time*3+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+75*(Time)+0+jt] = (Time+1)*5+Time*4+jt;
     jCol[15*(Time+1)+75*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+76*(Time)+0+jt] = (Time+1)*5+Time*4+jt;
     jCol[15*(Time+1)+76*(Time)+0+jt] = (Time+1)*0+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+77*(Time)+0+jt] = (Time+1)*5+Time*4+jt;
     jCol[15*(Time+1)+77*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+78*(Time)+0+jt] = (Time+1)*5+Time*4+jt;
     jCol[15*(Time+1)+78*(Time)+0+jt] = (Time+1)*1+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+79*(Time)+0+jt] = (Time+1)*5+Time*4+jt;
     jCol[15*(Time+1)+79*(Time)+0+jt] = (Time+1)*2+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+80*(Time)+0+jt] = (Time+1)*5+Time*4+jt;
     jCol[15*(Time+1)+80*(Time)+0+jt] = (Time+1)*2+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+81*(Time)+0+jt] = (Time+1)*5+Time*4+jt;
     jCol[15*(Time+1)+81*(Time)+0+jt] = (Time+1)*3+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+82*(Time)+0+jt] = (Time+1)*5+Time*4+jt;
     jCol[15*(Time+1)+82*(Time)+0+jt] = (Time+1)*3+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+83*(Time)+0+jt] = (Time+1)*5+Time*4+jt;
     jCol[15*(Time+1)+83*(Time)+0+jt] = (Time+1)*4+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+84*(Time)+0+jt] = (Time+1)*5+Time*4+jt;
     jCol[15*(Time+1)+84*(Time)+0+jt] = (Time+1)*4+jt+1;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+85*(Time)+0+jt] = (Time+1)*5+Time*4+jt;
     jCol[15*(Time+1)+85*(Time)+0+jt] = (Time+1)*5+Time*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+86*(Time)+0+jt] = (Time+1)*5+Time*4+jt;
     jCol[15*(Time+1)+86*(Time)+0+jt] = (Time+1)*5+Time*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+87*(Time)+0+jt] = (Time+1)*5+Time*4+jt;
     jCol[15*(Time+1)+87*(Time)+0+jt] = (Time+1)*5+Time*2+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+88*(Time)+0+jt] = (Time+1)*5+Time*4+jt;
     jCol[15*(Time+1)+88*(Time)+0+jt] = (Time+1)*5+Time*3+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[15*(Time+1)+89*(Time)+0+jt] = (Time+1)*5+Time*4+jt;
     jCol[15*(Time+1)+89*(Time)+0+jt] = (Time+1)*5+Time*4+jt;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[15*(Time+1)+90*(Time)+0+jt] = 2*Time*5+5;
     jCol[15*(Time+1)+90*(Time)+0+jt] = (Time+1)*0+jt;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[16*(Time+1)+90*(Time)+0+jt] = 2*Time*5+5;
     jCol[16*(Time+1)+90*(Time)+0+jt] = (Time+1)*1+jt;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[17*(Time+1)+90*(Time)+0+jt] = 2*Time*5+5;
     jCol[17*(Time+1)+90*(Time)+0+jt] = (Time+1)*2+jt;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[18*(Time+1)+90*(Time)+0+jt] = 2*Time*5+5;
     jCol[18*(Time+1)+90*(Time)+0+jt] = (Time+1)*3+jt;
   }

   for(Index jt=0;jt<Time+1;jt++) {
     iRow[19*(Time+1)+90*(Time)+0+jt] = 2*Time*5+5;
     jCol[19*(Time+1)+90*(Time)+0+jt] = (Time+1)*4+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[20*(Time+1)+90*(Time)+0+jt] = 2*Time*5+5;
     jCol[20*(Time+1)+90*(Time)+0+jt] = (Time+1)*5+Time*0+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[20*(Time+1)+91*(Time)+0+jt] = 2*Time*5+5;
     jCol[20*(Time+1)+91*(Time)+0+jt] = (Time+1)*5+Time*1+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[20*(Time+1)+92*(Time)+0+jt] = 2*Time*5+5;
     jCol[20*(Time+1)+92*(Time)+0+jt] = (Time+1)*5+Time*2+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[20*(Time+1)+93*(Time)+0+jt] = 2*Time*5+5;
     jCol[20*(Time+1)+93*(Time)+0+jt] = (Time+1)*5+Time*3+jt;
   }

   for(Index jt=0;jt<Time;jt++) {
     iRow[20*(Time+1)+94*(Time)+0+jt] = 2*Time*5+5;
     jCol[20*(Time+1)+94*(Time)+0+jt] = (Time+1)*5+Time*4+jt;
   }

   for(Index jt=0;jt<1;jt++) {
     iRow[20*(Time+1)+95*(Time)+0+jt] = 2*Time*5+5;
     jCol[20*(Time+1)+95*(Time)+0+jt] = 2*Time*5+5;
   }
}

else {
  // return the values.  This is a symmetric matrix, fill the lower left
  // triangle only
  // initialize the values array
  // Point to the initial starting spot for the Hessian elements

  for(Index jt=0;jt<20*(Time+1)+95*Time+1;jt++) values[jt] = 0.; // Initialize matrix

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

     Xdval[0] = data1[2*jt];
     Xdval2[0] = data1[2*jt+1];
     Xdvalp1[0] = data1[2*jt+2];

     for(Index i=0;i<nP;i++) {
        Pval[i] = x[(2*Time+1)*(nY+2*nU)+i];
     } //end for loop

    values[20*(Time+1)+95*(Time)+0] += obj_factor*bounds[0][3]*(2.0*pow(hstep, 2))/(2*Time+1);

    values[20*(Time+1)+95*(Time)+0] += obj_factor*bounds[1][3]*(2.0*pow(hstep, 2))/(2*Time+1);

    values[20*(Time+1)+95*(Time)+0] += obj_factor*bounds[2][3]*(2.0*pow(hstep, 2))/(2*Time+1);

    values[20*(Time+1)+95*(Time)+0] += obj_factor*bounds[3][3]*(2.0*pow(hstep, 2))/(2*Time+1);

    values[20*(Time+1)+95*(Time)+0] += obj_factor*bounds[4][3]*(2.0*pow(hstep, 2))/(2*Time+1);

   values[0*(Time+1)+0*(Time)+0+jt] += obj_factor*1*(8)/(2*Time+1);

   values[15*(Time+1)+35*(Time)+0+jt] += obj_factor*1*(8)/(2*Time+1);

   values[0*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[0][3]*((-0.333333333333333*hstep + 2)*(-0.166666666666667*hstep + 1) + (-0.25*hstep + 1.0)*(-0.125*hstep + 0.5))/(2*Time+1);

   values[1*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[0][3]*((-0.333333333333333*hstep - 2)*(-0.166666666666667*hstep + 1) + (-0.125*hstep + 0.5)*(0.25*hstep + 1.0))/(2*Time+1);

   values[0*(Time+1)+0*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*((-0.333333333333333*hstep - 2)*(-0.166666666666667*hstep - 1) + (0.125*hstep + 0.5)*(0.25*hstep + 1.0))/(2*Time+1);

   values[1*(Time+1)+1*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.333333333333333*hstep*Xval[4]*(-0.166666666666667*hstep + 1) + 0.25*hstep*Xval[4]*(-0.125*hstep + 0.5))/(2*Time+1);

   values[2*(Time+1)+1*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.333333333333333*hstep*Xval[4]*(-0.166666666666667*hstep - 1) + 0.25*hstep*Xval[4]*(0.125*hstep + 0.5))/(2*Time+1);

   values[2*(Time+1)+2*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0868055555555556*pow(hstep, 2)*pow(Xval[4], 2))/(2*Time+1);

   values[3*(Time+1)+2*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.333333333333333*hstep*Xvalp1[4]*(-0.166666666666667*hstep + 1) - 0.25*hstep*Xvalp1[4]*(-0.125*hstep + 0.5))/(2*Time+1);

   values[1*(Time+1)+1*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(0.333333333333333*hstep*Xvalp1[4]*(-0.166666666666667*hstep - 1) - 0.25*hstep*Xvalp1[4]*(0.125*hstep + 0.5))/(2*Time+1);

   values[3*(Time+1)+3*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0243055555555556*pow(hstep, 2)*Xval[4]*Xvalp1[4])/(2*Time+1);

   values[2*(Time+1)+2*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(0.0868055555555556*pow(hstep, 2)*pow(Xvalp1[4], 2))/(2*Time+1);

   values[6*(Time+1)+9*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.333333333333333*hstep*Xval[4]*(-0.166666666666667*hstep + 1) - 0.25*hstep*Xval[4]*(-0.125*hstep + 0.5))/(2*Time+1);

   values[7*(Time+1)+9*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.333333333333333*hstep*Xval[4]*(-0.166666666666667*hstep - 1) - 0.25*hstep*Xval[4]*(0.125*hstep + 0.5))/(2*Time+1);

   values[7*(Time+1)+10*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.0868055555555556*pow(hstep, 2)*pow(Xval[4], 2))/(2*Time+1);

   values[8*(Time+1)+10*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.0243055555555556*pow(hstep, 2)*Xval[4]*Xvalp1[4])/(2*Time+1);

   values[9*(Time+1)+12*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0868055555555556*pow(hstep, 2)*pow(Xval[4], 2))/(2*Time+1);

   values[10*(Time+1)+12*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.333333333333333*hstep*Xvalp1[4]*(-0.166666666666667*hstep + 1) + 0.25*hstep*Xvalp1[4]*(-0.125*hstep + 0.5))/(2*Time+1);

   values[6*(Time+1)+9*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(-0.333333333333333*hstep*Xvalp1[4]*(-0.166666666666667*hstep - 1) + 0.25*hstep*Xvalp1[4]*(0.125*hstep + 0.5))/(2*Time+1);

   values[10*(Time+1)+13*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.0243055555555556*pow(hstep, 2)*Xval[4]*Xvalp1[4])/(2*Time+1);

   values[7*(Time+1)+10*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(-0.0868055555555556*pow(hstep, 2)*pow(Xvalp1[4], 2))/(2*Time+1);

   values[10*(Time+1)+15*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0243055555555556*pow(hstep, 2)*Xval[4]*Xvalp1[4])/(2*Time+1);

   values[9*(Time+1)+12*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(0.0868055555555556*pow(hstep, 2)*pow(Xvalp1[4], 2))/(2*Time+1);

   values[10*(Time+1)+16*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.333333333333333*hstep*(-0.166666666666667*hstep + 1)*(Xval[1] - Xval[3]) + 0.25*hstep*(-0.125*hstep + 0.5)*(Xval[1] - Xval[3]))/(2*Time+1);

   values[11*(Time+1)+16*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.333333333333333*hstep*(-0.166666666666667*hstep - 1)*(Xval[1] - Xval[3]) + 0.25*hstep*(0.125*hstep + 0.5)*(Xval[1] - Xval[3]))/(2*Time+1);

   values[11*(Time+1)+17*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0868055555555556*pow(hstep, 2)*Xval[4]*(Xval[1] - Xval[3]) + 0.333333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) + 4*Xval2[4]*(Xval2[1] - Xval2[3]) + Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + Xval[0] - Xvalp1[0]) + 0.25*hstep*(0.125*hstep*(-Xval[0] + Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) - Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);

   values[12*(Time+1)+17*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0243055555555556*pow(hstep, 2)*Xvalp1[4]*(Xval[1] - Xval[3]))/(2*Time+1);

   values[13*(Time+1)+19*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.0868055555555556*pow(hstep, 2)*Xval[4]*(Xval[1] - Xval[3]) - 0.333333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) + 4*Xval2[4]*(Xval2[1] - Xval2[3]) + Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + Xval[0] - Xvalp1[0]) - 0.25*hstep*(0.125*hstep*(-Xval[0] + Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) - Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);

   values[14*(Time+1)+19*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.0243055555555556*pow(hstep, 2)*Xvalp1[4]*(Xval[1] - Xval[3]))/(2*Time+1);

   values[14*(Time+1)+20*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0868055555555556*pow(hstep, 2)*pow(Xval[1] - Xval[3], 2))/(2*Time+1);

   values[15*(Time+1)+20*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.333333333333333*hstep*(-0.166666666666667*hstep + 1)*(Xvalp1[1] - Xvalp1[3]) + 0.25*hstep*(-0.125*hstep + 0.5)*(-Xvalp1[1] + Xvalp1[3]))/(2*Time+1);

   values[10*(Time+1)+16*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(0.333333333333333*hstep*(-0.166666666666667*hstep - 1)*(Xvalp1[1] - Xvalp1[3]) + 0.25*hstep*(0.125*hstep + 0.5)*(-Xvalp1[1] + Xvalp1[3]))/(2*Time+1);

   values[15*(Time+1)+21*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.03125*pow(hstep, 2)*Xval[4]*(-Xvalp1[1] + Xvalp1[3]) + 0.0555555555555556*pow(hstep, 2)*Xval[4]*(Xvalp1[1] - Xvalp1[3]))/(2*Time+1);

   values[11*(Time+1)+17*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(-0.03125*pow(hstep, 2)*Xvalp1[4]*(-Xvalp1[1] + Xvalp1[3]) + 0.0555555555555556*pow(hstep, 2)*Xvalp1[4]*(Xvalp1[1] - Xvalp1[3]) + 0.333333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) + 4*Xval2[4]*(Xval2[1] - Xval2[3]) + Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + Xval[0] - Xvalp1[0]) - 0.25*hstep*(0.125*hstep*(-Xval[0] + Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) - Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);

   values[15*(Time+1)+23*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.03125*pow(hstep, 2)*Xval[4]*(-Xvalp1[1] + Xvalp1[3]) - 0.0555555555555556*pow(hstep, 2)*Xval[4]*(Xvalp1[1] - Xvalp1[3]))/(2*Time+1);

   values[13*(Time+1)+19*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(0.03125*pow(hstep, 2)*Xvalp1[4]*(-Xvalp1[1] + Xvalp1[3]) - 0.0555555555555556*pow(hstep, 2)*Xvalp1[4]*(Xvalp1[1] - Xvalp1[3]) - 0.333333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) + 4*Xval2[4]*(Xval2[1] - Xval2[3]) + Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + Xval[0] - Xvalp1[0]) + 0.25*hstep*(0.125*hstep*(-Xval[0] + Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) - Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);

   values[15*(Time+1)+24*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.03125*pow(hstep, 2)*(Xval[1] - Xval[3])*(-Xvalp1[1] + Xvalp1[3]) + 0.0555555555555556*pow(hstep, 2)*(Xval[1] - Xval[3])*(Xvalp1[1] - Xvalp1[3]))/(2*Time+1);

   values[14*(Time+1)+20*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(0.03125*pow(hstep, 2)*pow(-Xvalp1[1] + Xvalp1[3], 2) + 0.0555555555555556*pow(hstep, 2)*pow(Xvalp1[1] - Xvalp1[3], 2))/(2*Time+1);

   values[15*(Time+1)+25*(Time)+0+jt] += obj_factor*bounds[0][3]*(-1.33333333333333*hstep*(-0.166666666666667*hstep + 1) + 0.25*hstep - 1.0)/(2*Time+1);

   values[15*(Time+1)+26*(Time)+0+jt] += obj_factor*bounds[0][3]*(-1.33333333333333*hstep*(-0.166666666666667*hstep - 1) - 0.25*hstep - 1.0)/(2*Time+1);

   values[15*(Time+1)+27*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.222222222222222*pow(hstep, 2)*Xval[4] - 0.25*hstep*Xval[4])/(2*Time+1);

   values[15*(Time+1)+28*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.222222222222222*pow(hstep, 2)*Xvalp1[4] + 0.25*hstep*Xvalp1[4])/(2*Time+1);

   values[15*(Time+1)+31*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.222222222222222*pow(hstep, 2)*Xval[4] + 0.25*hstep*Xval[4])/(2*Time+1);

   values[15*(Time+1)+32*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.222222222222222*pow(hstep, 2)*Xvalp1[4] - 0.25*hstep*Xvalp1[4])/(2*Time+1);

   values[15*(Time+1)+33*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.222222222222222*pow(hstep, 2)*(Xval[1] - Xval[3]) - 0.25*hstep*(Xval[1] - Xval[3]))/(2*Time+1);

   values[15*(Time+1)+34*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.222222222222222*pow(hstep, 2)*(Xvalp1[1] - Xvalp1[3]) - 0.25*hstep*(-Xvalp1[1] + Xvalp1[3]))/(2*Time+1);

   values[15*(Time+1)+35*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.888888888888889*pow(hstep, 2) + 2)/(2*Time+1);

   values[15*(Time+1)+36*(Time)+0+jt] += obj_factor*bounds[0][3]*(1.33333333333333*hstep*Xval2[4]*(-0.166666666666667*hstep + 1))/(2*Time+1);

   values[15*(Time+1)+37*(Time)+0+jt] += obj_factor*bounds[0][3]*(1.33333333333333*hstep*Xval2[4]*(-0.166666666666667*hstep - 1))/(2*Time+1);

   values[15*(Time+1)+38*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.222222222222222*pow(hstep, 2)*Xval[4]*Xval2[4])/(2*Time+1);

   values[15*(Time+1)+39*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.222222222222222*pow(hstep, 2)*Xval2[4]*Xvalp1[4])/(2*Time+1);

   values[15*(Time+1)+42*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.222222222222222*pow(hstep, 2)*Xval[4]*Xval2[4])/(2*Time+1);

   values[15*(Time+1)+43*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[4]*Xvalp1[4])/(2*Time+1);

   values[15*(Time+1)+44*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.222222222222222*pow(hstep, 2)*Xval2[4]*(Xval[1] - Xval[3]))/(2*Time+1);

   values[15*(Time+1)+45*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.222222222222222*pow(hstep, 2)*Xval2[4]*(Xvalp1[1] - Xvalp1[3]))/(2*Time+1);

   values[15*(Time+1)+46*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.888888888888889*pow(hstep, 2)*Xval2[4])/(2*Time+1);

   values[15*(Time+1)+47*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.888888888888889*pow(hstep, 2)*pow(Xval2[4], 2))/(2*Time+1);

   values[15*(Time+1)+61*(Time)+0+jt] += obj_factor*bounds[0][3]*(-1.33333333333333*hstep*Xval2[4]*(-0.166666666666667*hstep + 1))/(2*Time+1);

   values[15*(Time+1)+62*(Time)+0+jt] += obj_factor*bounds[0][3]*(-1.33333333333333*hstep*Xval2[4]*(-0.166666666666667*hstep - 1))/(2*Time+1);

   values[15*(Time+1)+63*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.222222222222222*pow(hstep, 2)*Xval[4]*Xval2[4])/(2*Time+1);

   values[15*(Time+1)+64*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[4]*Xvalp1[4])/(2*Time+1);

   values[15*(Time+1)+67*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.222222222222222*pow(hstep, 2)*Xval[4]*Xval2[4])/(2*Time+1);

   values[15*(Time+1)+68*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.222222222222222*pow(hstep, 2)*Xval2[4]*Xvalp1[4])/(2*Time+1);

   values[15*(Time+1)+69*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[4]*(Xval[1] - Xval[3]))/(2*Time+1);

   values[15*(Time+1)+70*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[4]*(Xvalp1[1] - Xvalp1[3]))/(2*Time+1);

   values[15*(Time+1)+71*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.888888888888889*pow(hstep, 2)*Xval2[4])/(2*Time+1);

   values[15*(Time+1)+72*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.888888888888889*pow(hstep, 2)*pow(Xval2[4], 2))/(2*Time+1);

   values[15*(Time+1)+74*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.888888888888889*pow(hstep, 2)*pow(Xval2[4], 2))/(2*Time+1);

   values[15*(Time+1)+75*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.333333333333333*hstep*(-0.166666666666667*hstep + 1)*(4*Xval2[1] - 4*Xval2[3]))/(2*Time+1);

   values[15*(Time+1)+76*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.333333333333333*hstep*(-0.166666666666667*hstep - 1)*(4*Xval2[1] - 4*Xval2[3]))/(2*Time+1);

   values[15*(Time+1)+77*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0555555555555556*pow(hstep, 2)*Xval[4]*(4*Xval2[1] - 4*Xval2[3]))/(2*Time+1);

   values[15*(Time+1)+78*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0555555555555556*pow(hstep, 2)*Xvalp1[4]*(4*Xval2[1] - 4*Xval2[3]))/(2*Time+1);

   values[15*(Time+1)+81*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.0555555555555556*pow(hstep, 2)*Xval[4]*(4*Xval2[1] - 4*Xval2[3]))/(2*Time+1);

   values[15*(Time+1)+82*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.0555555555555556*pow(hstep, 2)*Xvalp1[4]*(4*Xval2[1] - 4*Xval2[3]))/(2*Time+1);

   values[15*(Time+1)+83*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0555555555555556*pow(hstep, 2)*(Xval[1] - Xval[3])*(4*Xval2[1] - 4*Xval2[3]))/(2*Time+1);

   values[15*(Time+1)+84*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0555555555555556*pow(hstep, 2)*(4*Xval2[1] - 4*Xval2[3])*(Xvalp1[1] - Xvalp1[3]))/(2*Time+1);

   values[15*(Time+1)+85*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.222222222222222*pow(hstep, 2)*(4*Xval2[1] - 4*Xval2[3]))/(2*Time+1);

   values[15*(Time+1)+86*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.222222222222222*pow(hstep, 2)*Xval2[4]*(4*Xval2[1] - 4*Xval2[3]) + 1.33333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) + 4*Xval2[4]*(Xval2[1] - Xval2[3]) + Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + Xval[0] - Xvalp1[0]))/(2*Time+1);

   values[15*(Time+1)+88*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[4]*(4*Xval2[1] - 4*Xval2[3]) - 1.33333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] - Xval[0] - 4*Xval2[0] - Xvalp1[0] + Xval[4]*(Xval[1] - Xval[3]) + 4*Xval2[4]*(Xval2[1] - Xval2[3]) + Xvalp1[4]*(Xvalp1[1] - Xvalp1[3])) + Xval[0] - Xvalp1[0]))/(2*Time+1);

   values[15*(Time+1)+89*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0555555555555556*pow(hstep, 2)*pow(4*Xval2[1] - 4*Xval2[3], 2))/(2*Time+1);

   values[15*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[0][3]*(2.0*hstep*(-0.166666666666667*hstep + 1))/(2*Time+1);

   values[15*(Time+1)+90*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(2.0*hstep*(-0.166666666666667*hstep - 1))/(2*Time+1);

   values[16*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.333333333333333*pow(hstep, 2)*Xval[4])/(2*Time+1);

   values[16*(Time+1)+90*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(0.333333333333333*pow(hstep, 2)*Xvalp1[4])/(2*Time+1);

   values[18*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.333333333333333*pow(hstep, 2)*Xval[4])/(2*Time+1);

   values[18*(Time+1)+90*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(-0.333333333333333*pow(hstep, 2)*Xvalp1[4])/(2*Time+1);

   values[19*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.333333333333333*pow(hstep, 2)*(Xval[1] - Xval[3]))/(2*Time+1);

   values[19*(Time+1)+90*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(0.333333333333333*pow(hstep, 2)*(Xvalp1[1] - Xvalp1[3]))/(2*Time+1);

   values[20*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[0][3]*(-1.33333333333333*pow(hstep, 2))/(2*Time+1);

   values[20*(Time+1)+91*(Time)+0+jt] += obj_factor*bounds[0][3]*(1.33333333333333*pow(hstep, 2)*Xval2[4])/(2*Time+1);

   values[20*(Time+1)+93*(Time)+0+jt] += obj_factor*bounds[0][3]*(-1.33333333333333*pow(hstep, 2)*Xval2[4])/(2*Time+1);

   values[20*(Time+1)+94*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.333333333333333*pow(hstep, 2)*(4*Xval2[1] - 4*Xval2[3]))/(2*Time+1);

   values[0*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.0868055555555556*pow(hstep, 2)*pow(Xval[2] - Xval[4], 2))/(2*Time+1);

   values[1*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.03125*pow(hstep, 2)*(Xval[2] - Xval[4])*(-Xvalp1[2] + Xvalp1[4]) + 0.0555555555555556*pow(hstep, 2)*(Xval[2] - Xval[4])*(Xvalp1[2] - Xvalp1[4]))/(2*Time+1);

   values[0*(Time+1)+0*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(0.03125*pow(hstep, 2)*pow(-Xvalp1[2] + Xvalp1[4], 2) + 0.0555555555555556*pow(hstep, 2)*pow(Xvalp1[2] - Xvalp1[4], 2))/(2*Time+1);

   values[1*(Time+1)+1*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep + 2)*(Xval[2] - Xval[4]) + 0.125*hstep*(-0.25*hstep + 1.0)*(Xval[2] - Xval[4]))/(2*Time+1);

   values[2*(Time+1)+1*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep + 2)*(Xvalp1[2] - Xvalp1[4]) + 0.125*hstep*(-0.25*hstep + 1.0)*(-Xvalp1[2] + Xvalp1[4]))/(2*Time+1);

   values[2*(Time+1)+2*(Time)+0+jt] += obj_factor*bounds[1][3]*((-0.333333333333333*hstep + 2)*(-0.166666666666667*hstep + 1) + (-0.25*hstep + 1.0)*(-0.125*hstep + 0.5))/(2*Time+1);

   values[3*(Time+1)+2*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep - 2)*(Xval[2] - Xval[4]) + 0.125*hstep*(0.25*hstep + 1.0)*(Xval[2] - Xval[4]))/(2*Time+1);

   values[1*(Time+1)+1*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep - 2)*(Xvalp1[2] - Xvalp1[4]) + 0.125*hstep*(0.25*hstep + 1.0)*(-Xvalp1[2] + Xvalp1[4]))/(2*Time+1);

   values[3*(Time+1)+3*(Time)+0+jt] += obj_factor*bounds[1][3]*((-0.333333333333333*hstep - 2)*(-0.166666666666667*hstep + 1) + (-0.125*hstep + 0.5)*(0.25*hstep + 1.0))/(2*Time+1);

   values[2*(Time+1)+2*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*((-0.333333333333333*hstep - 2)*(-0.166666666666667*hstep - 1) + (0.125*hstep + 0.5)*(0.25*hstep + 1.0))/(2*Time+1);

   values[3*(Time+1)+4*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.0868055555555556*pow(hstep, 2)*Xval[0]*(Xval[2] - Xval[4]) + 0.333333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[0]*(Xval[2] - Xval[4]) + 4*Xval2[0]*(Xval2[2] - Xval2[4]) + Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] - 4*Xval2[1] - Xvalp1[1]) + Xval[1] - Xvalp1[1]) + 0.25*hstep*(0.125*hstep*(Xval[0]*(Xval[2] - Xval[4]) - Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] + Xvalp1[1]) + 0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1]))/(2*Time+1);

   values[4*(Time+1)+4*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.03125*pow(hstep, 2)*Xval[0]*(-Xvalp1[2] + Xvalp1[4]) + 0.0555555555555556*pow(hstep, 2)*Xval[0]*(Xvalp1[2] - Xvalp1[4]))/(2*Time+1);

   values[4*(Time+1)+5*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.333333333333333*hstep*Xval[0]*(-0.166666666666667*hstep + 1) + 0.25*hstep*Xval[0]*(-0.125*hstep + 0.5))/(2*Time+1);

   values[5*(Time+1)+5*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.333333333333333*hstep*Xval[0]*(-0.166666666666667*hstep - 1) + 0.25*hstep*Xval[0]*(0.125*hstep + 0.5))/(2*Time+1);

   values[5*(Time+1)+6*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.0868055555555556*pow(hstep, 2)*pow(Xval[0], 2))/(2*Time+1);

   values[6*(Time+1)+6*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.0243055555555556*pow(hstep, 2)*Xvalp1[0]*(Xval[2] - Xval[4]))/(2*Time+1);

   values[3*(Time+1)+4*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(-0.03125*pow(hstep, 2)*Xvalp1[0]*(-Xvalp1[2] + Xvalp1[4]) + 0.0555555555555556*pow(hstep, 2)*Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) + 0.333333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[0]*(Xval[2] - Xval[4]) + 4*Xval2[0]*(Xval2[2] - Xval2[4]) + Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] - 4*Xval2[1] - Xvalp1[1]) + Xval[1] - Xvalp1[1]) - 0.25*hstep*(0.125*hstep*(Xval[0]*(Xval[2] - Xval[4]) - Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] + Xvalp1[1]) + 0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1]))/(2*Time+1);

   values[6*(Time+1)+7*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.333333333333333*hstep*Xvalp1[0]*(-0.166666666666667*hstep + 1) - 0.25*hstep*Xvalp1[0]*(-0.125*hstep + 0.5))/(2*Time+1);

   values[4*(Time+1)+5*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(0.333333333333333*hstep*Xvalp1[0]*(-0.166666666666667*hstep - 1) - 0.25*hstep*Xvalp1[0]*(0.125*hstep + 0.5))/(2*Time+1);

   values[6*(Time+1)+8*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.0243055555555556*pow(hstep, 2)*Xval[0]*Xvalp1[0])/(2*Time+1);

   values[5*(Time+1)+6*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(0.0868055555555556*pow(hstep, 2)*pow(Xvalp1[0], 2))/(2*Time+1);

   values[10*(Time+1)+16*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.0868055555555556*pow(hstep, 2)*Xval[0]*(Xval[2] - Xval[4]) - 0.333333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[0]*(Xval[2] - Xval[4]) + 4*Xval2[0]*(Xval2[2] - Xval2[4]) + Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] - 4*Xval2[1] - Xvalp1[1]) + Xval[1] - Xvalp1[1]) - 0.25*hstep*(0.125*hstep*(Xval[0]*(Xval[2] - Xval[4]) - Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] + Xvalp1[1]) + 0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1]))/(2*Time+1);

   values[11*(Time+1)+16*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.03125*pow(hstep, 2)*Xval[0]*(-Xvalp1[2] + Xvalp1[4]) - 0.0555555555555556*pow(hstep, 2)*Xval[0]*(Xvalp1[2] - Xvalp1[4]))/(2*Time+1);

   values[11*(Time+1)+17*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.333333333333333*hstep*Xval[0]*(-0.166666666666667*hstep + 1) - 0.25*hstep*Xval[0]*(-0.125*hstep + 0.5))/(2*Time+1);

   values[12*(Time+1)+17*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.333333333333333*hstep*Xval[0]*(-0.166666666666667*hstep - 1) - 0.25*hstep*Xval[0]*(0.125*hstep + 0.5))/(2*Time+1);

   values[12*(Time+1)+18*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.0868055555555556*pow(hstep, 2)*pow(Xval[0], 2))/(2*Time+1);

   values[13*(Time+1)+18*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.0243055555555556*pow(hstep, 2)*Xval[0]*Xvalp1[0])/(2*Time+1);

   values[14*(Time+1)+20*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.0868055555555556*pow(hstep, 2)*pow(Xval[0], 2))/(2*Time+1);

   values[15*(Time+1)+20*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.0243055555555556*pow(hstep, 2)*Xvalp1[0]*(Xval[2] - Xval[4]))/(2*Time+1);

   values[10*(Time+1)+16*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(0.03125*pow(hstep, 2)*Xvalp1[0]*(-Xvalp1[2] + Xvalp1[4]) - 0.0555555555555556*pow(hstep, 2)*Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - 0.333333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[0]*(Xval[2] - Xval[4]) + 4*Xval2[0]*(Xval2[2] - Xval2[4]) + Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] - 4*Xval2[1] - Xvalp1[1]) + Xval[1] - Xvalp1[1]) + 0.25*hstep*(0.125*hstep*(Xval[0]*(Xval[2] - Xval[4]) - Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] + Xvalp1[1]) + 0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1]))/(2*Time+1);

   values[15*(Time+1)+21*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.333333333333333*hstep*Xvalp1[0]*(-0.166666666666667*hstep + 1) + 0.25*hstep*Xvalp1[0]*(-0.125*hstep + 0.5))/(2*Time+1);

   values[11*(Time+1)+17*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(-0.333333333333333*hstep*Xvalp1[0]*(-0.166666666666667*hstep - 1) + 0.25*hstep*Xvalp1[0]*(0.125*hstep + 0.5))/(2*Time+1);

   values[15*(Time+1)+22*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.0243055555555556*pow(hstep, 2)*Xval[0]*Xvalp1[0])/(2*Time+1);

   values[12*(Time+1)+18*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(-0.0868055555555556*pow(hstep, 2)*pow(Xvalp1[0], 2))/(2*Time+1);

   values[15*(Time+1)+24*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.0243055555555556*pow(hstep, 2)*Xval[0]*Xvalp1[0])/(2*Time+1);

   values[14*(Time+1)+20*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(0.0868055555555556*pow(hstep, 2)*pow(Xvalp1[0], 2))/(2*Time+1);

   values[15*(Time+1)+25*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.0555555555555556*pow(hstep, 2)*(Xval[2] - Xval[4])*(4*Xval2[2] - 4*Xval2[4]))/(2*Time+1);

   values[15*(Time+1)+26*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.0555555555555556*pow(hstep, 2)*(4*Xval2[2] - 4*Xval2[4])*(Xvalp1[2] - Xvalp1[4]))/(2*Time+1);

   values[15*(Time+1)+27*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.333333333333333*hstep*(-0.166666666666667*hstep + 1)*(4*Xval2[2] - 4*Xval2[4]))/(2*Time+1);

   values[15*(Time+1)+28*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.333333333333333*hstep*(-0.166666666666667*hstep - 1)*(4*Xval2[2] - 4*Xval2[4]))/(2*Time+1);

   values[15*(Time+1)+29*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.0555555555555556*pow(hstep, 2)*Xval[0]*(4*Xval2[2] - 4*Xval2[4]))/(2*Time+1);

   values[15*(Time+1)+30*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.0555555555555556*pow(hstep, 2)*Xvalp1[0]*(4*Xval2[2] - 4*Xval2[4]))/(2*Time+1);

   values[15*(Time+1)+33*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.0555555555555556*pow(hstep, 2)*Xval[0]*(4*Xval2[2] - 4*Xval2[4]))/(2*Time+1);

   values[15*(Time+1)+34*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.0555555555555556*pow(hstep, 2)*Xvalp1[0]*(4*Xval2[2] - 4*Xval2[4]))/(2*Time+1);

   values[15*(Time+1)+35*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.0555555555555556*pow(hstep, 2)*pow(4*Xval2[2] - 4*Xval2[4], 2))/(2*Time+1);

   values[15*(Time+1)+36*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.222222222222222*pow(hstep, 2)*(Xval[2] - Xval[4]) - 0.25*hstep*(Xval[2] - Xval[4]))/(2*Time+1);

   values[15*(Time+1)+37*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.222222222222222*pow(hstep, 2)*(Xvalp1[2] - Xvalp1[4]) - 0.25*hstep*(-Xvalp1[2] + Xvalp1[4]))/(2*Time+1);

   values[15*(Time+1)+38*(Time)+0+jt] += obj_factor*bounds[1][3]*(-1.33333333333333*hstep*(-0.166666666666667*hstep + 1) + 0.25*hstep - 1.0)/(2*Time+1);

   values[15*(Time+1)+39*(Time)+0+jt] += obj_factor*bounds[1][3]*(-1.33333333333333*hstep*(-0.166666666666667*hstep - 1) - 0.25*hstep - 1.0)/(2*Time+1);

   values[15*(Time+1)+40*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.222222222222222*pow(hstep, 2)*Xval[0] - 0.25*hstep*Xval[0])/(2*Time+1);

   values[15*(Time+1)+41*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.222222222222222*pow(hstep, 2)*Xvalp1[0] + 0.25*hstep*Xvalp1[0])/(2*Time+1);

   values[15*(Time+1)+44*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.222222222222222*pow(hstep, 2)*Xval[0] + 0.25*hstep*Xval[0])/(2*Time+1);

   values[15*(Time+1)+45*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.222222222222222*pow(hstep, 2)*Xvalp1[0] - 0.25*hstep*Xvalp1[0])/(2*Time+1);

   values[15*(Time+1)+46*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.222222222222222*pow(hstep, 2)*(4*Xval2[2] - 4*Xval2[4]))/(2*Time+1);

   values[15*(Time+1)+47*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.888888888888889*pow(hstep, 2) + 2)/(2*Time+1);

   values[15*(Time+1)+48*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.222222222222222*pow(hstep, 2)*Xval2[0]*(Xval[2] - Xval[4]))/(2*Time+1);

   values[15*(Time+1)+49*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.222222222222222*pow(hstep, 2)*Xval2[0]*(Xvalp1[2] - Xvalp1[4]))/(2*Time+1);

   values[15*(Time+1)+50*(Time)+0+jt] += obj_factor*bounds[1][3]*(1.33333333333333*hstep*Xval2[0]*(-0.166666666666667*hstep + 1))/(2*Time+1);

   values[15*(Time+1)+51*(Time)+0+jt] += obj_factor*bounds[1][3]*(1.33333333333333*hstep*Xval2[0]*(-0.166666666666667*hstep - 1))/(2*Time+1);

   values[15*(Time+1)+52*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.222222222222222*pow(hstep, 2)*Xval[0]*Xval2[0])/(2*Time+1);

   values[15*(Time+1)+53*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.222222222222222*pow(hstep, 2)*Xval2[0]*Xvalp1[0])/(2*Time+1);

   values[15*(Time+1)+56*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.222222222222222*pow(hstep, 2)*Xval[0]*Xval2[0])/(2*Time+1);

   values[15*(Time+1)+57*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[0]*Xvalp1[0])/(2*Time+1);

   values[15*(Time+1)+58*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.222222222222222*pow(hstep, 2)*Xval2[0]*(4*Xval2[2] - 4*Xval2[4]) + 1.33333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[0]*(Xval[2] - Xval[4]) + 4*Xval2[0]*(Xval2[2] - Xval2[4]) + Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] - 4*Xval2[1] - Xvalp1[1]) + Xval[1] - Xvalp1[1]))/(2*Time+1);

   values[15*(Time+1)+59*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.888888888888889*pow(hstep, 2)*Xval2[0])/(2*Time+1);

   values[15*(Time+1)+60*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.888888888888889*pow(hstep, 2)*pow(Xval2[0], 2))/(2*Time+1);

   values[15*(Time+1)+75*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[0]*(Xval[2] - Xval[4]))/(2*Time+1);

   values[15*(Time+1)+76*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[0]*(Xvalp1[2] - Xvalp1[4]))/(2*Time+1);

   values[15*(Time+1)+77*(Time)+0+jt] += obj_factor*bounds[1][3]*(-1.33333333333333*hstep*Xval2[0]*(-0.166666666666667*hstep + 1))/(2*Time+1);

   values[15*(Time+1)+78*(Time)+0+jt] += obj_factor*bounds[1][3]*(-1.33333333333333*hstep*Xval2[0]*(-0.166666666666667*hstep - 1))/(2*Time+1);

   values[15*(Time+1)+79*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.222222222222222*pow(hstep, 2)*Xval[0]*Xval2[0])/(2*Time+1);

   values[15*(Time+1)+80*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[0]*Xvalp1[0])/(2*Time+1);

   values[15*(Time+1)+83*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.222222222222222*pow(hstep, 2)*Xval[0]*Xval2[0])/(2*Time+1);

   values[15*(Time+1)+84*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.222222222222222*pow(hstep, 2)*Xval2[0]*Xvalp1[0])/(2*Time+1);

   values[15*(Time+1)+85*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[0]*(4*Xval2[2] - 4*Xval2[4]) - 1.33333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[0]*(Xval[2] - Xval[4]) + 4*Xval2[0]*(Xval2[2] - Xval2[4]) + Xvalp1[0]*(Xvalp1[2] - Xvalp1[4]) - Xval[1] - 4*Xval2[1] - Xvalp1[1]) + Xval[1] - Xvalp1[1]))/(2*Time+1);

   values[15*(Time+1)+86*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.888888888888889*pow(hstep, 2)*Xval2[0])/(2*Time+1);

   values[15*(Time+1)+87*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.888888888888889*pow(hstep, 2)*pow(Xval2[0], 2))/(2*Time+1);

   values[15*(Time+1)+89*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.888888888888889*pow(hstep, 2)*pow(Xval2[0], 2))/(2*Time+1);

   values[15*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.333333333333333*pow(hstep, 2)*(Xval[2] - Xval[4]))/(2*Time+1);

   values[15*(Time+1)+90*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(0.333333333333333*pow(hstep, 2)*(Xvalp1[2] - Xvalp1[4]))/(2*Time+1);

   values[16*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[1][3]*(2.0*hstep*(-0.166666666666667*hstep + 1))/(2*Time+1);

   values[16*(Time+1)+90*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(2.0*hstep*(-0.166666666666667*hstep - 1))/(2*Time+1);

   values[17*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.333333333333333*pow(hstep, 2)*Xval[0])/(2*Time+1);

   values[17*(Time+1)+90*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(0.333333333333333*pow(hstep, 2)*Xvalp1[0])/(2*Time+1);

   values[19*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[1][3]*(-0.333333333333333*pow(hstep, 2)*Xval[0])/(2*Time+1);

   values[19*(Time+1)+90*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(-0.333333333333333*pow(hstep, 2)*Xvalp1[0])/(2*Time+1);

   values[20*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[1][3]*(0.333333333333333*pow(hstep, 2)*(4*Xval2[2] - 4*Xval2[4]))/(2*Time+1);

   values[20*(Time+1)+91*(Time)+0+jt] += obj_factor*bounds[1][3]*(-1.33333333333333*pow(hstep, 2))/(2*Time+1);

   values[20*(Time+1)+92*(Time)+0+jt] += obj_factor*bounds[1][3]*(1.33333333333333*pow(hstep, 2)*Xval2[0])/(2*Time+1);

   values[20*(Time+1)+94*(Time)+0+jt] += obj_factor*bounds[1][3]*(-1.33333333333333*pow(hstep, 2)*Xval2[0])/(2*Time+1);

   values[0*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.0868055555555556*pow(hstep, 2)*pow(Xval[1], 2))/(2*Time+1);

   values[1*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.0243055555555556*pow(hstep, 2)*Xval[1]*Xvalp1[1])/(2*Time+1);

   values[0*(Time+1)+0*(Time)+0+ 1+jt] += obj_factor*bounds[2][3]*(0.0868055555555556*pow(hstep, 2)*pow(Xvalp1[1], 2))/(2*Time+1);

   values[1*(Time+1)+1*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.0868055555555556*pow(hstep, 2)*Xval[1]*(-Xval[0] + Xval[3]) - 0.333333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[1]*(-Xval[0] + Xval[3]) + 4*Xval2[1]*(-Xval2[0] + Xval2[3]) + Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] - 4*Xval2[2] - Xvalp1[2]) + Xval[2] - Xvalp1[2]) - 0.25*hstep*(0.125*hstep*(Xval[1]*(-Xval[0] + Xval[3]) - Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] + Xvalp1[2]) + 0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2]))/(2*Time+1);

   values[2*(Time+1)+1*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.0243055555555556*pow(hstep, 2)*Xvalp1[1]*(-Xval[0] + Xval[3]))/(2*Time+1);

   values[2*(Time+1)+2*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.0868055555555556*pow(hstep, 2)*pow(-Xval[0] + Xval[3], 2))/(2*Time+1);

   values[3*(Time+1)+2*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.0555555555555556*pow(hstep, 2)*Xval[1]*(-Xvalp1[0] + Xvalp1[3]) - 0.03125*pow(hstep, 2)*Xval[1]*(Xvalp1[0] - Xvalp1[3]))/(2*Time+1);

   values[1*(Time+1)+1*(Time)+0+ 1+jt] += obj_factor*bounds[2][3]*(-0.0555555555555556*pow(hstep, 2)*Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) + 0.03125*pow(hstep, 2)*Xvalp1[1]*(Xvalp1[0] - Xvalp1[3]) - 0.333333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[1]*(-Xval[0] + Xval[3]) + 4*Xval2[1]*(-Xval2[0] + Xval2[3]) + Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] - 4*Xval2[2] - Xvalp1[2]) + Xval[2] - Xvalp1[2]) + 0.25*hstep*(0.125*hstep*(Xval[1]*(-Xval[0] + Xval[3]) - Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] + Xvalp1[2]) + 0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2]))/(2*Time+1);

   values[3*(Time+1)+3*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.0555555555555556*pow(hstep, 2)*(-Xval[0] + Xval[3])*(-Xvalp1[0] + Xvalp1[3]) + 0.03125*pow(hstep, 2)*(-Xval[0] + Xval[3])*(Xvalp1[0] - Xvalp1[3]))/(2*Time+1);

   values[2*(Time+1)+2*(Time)+0+ 1+jt] += obj_factor*bounds[2][3]*(0.0555555555555556*pow(hstep, 2)*pow(-Xvalp1[0] + Xvalp1[3], 2) + 0.03125*pow(hstep, 2)*pow(Xvalp1[0] - Xvalp1[3], 2))/(2*Time+1);

   values[3*(Time+1)+4*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.166666666666667*hstep*Xval[1]*(-0.333333333333333*hstep + 2) - 0.125*hstep*Xval[1]*(-0.25*hstep + 1.0))/(2*Time+1);

   values[4*(Time+1)+4*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.166666666666667*hstep*Xvalp1[1]*(-0.333333333333333*hstep + 2) + 0.125*hstep*Xvalp1[1]*(-0.25*hstep + 1.0))/(2*Time+1);

   values[4*(Time+1)+5*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep + 2)*(-Xval[0] + Xval[3]) + 0.125*hstep*(-0.25*hstep + 1.0)*(-Xval[0] + Xval[3]))/(2*Time+1);

   values[5*(Time+1)+5*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep + 2)*(-Xvalp1[0] + Xvalp1[3]) + 0.125*hstep*(-0.25*hstep + 1.0)*(Xvalp1[0] - Xvalp1[3]))/(2*Time+1);

   values[5*(Time+1)+6*(Time)+0+jt] += obj_factor*bounds[2][3]*((-0.333333333333333*hstep + 2)*(-0.166666666666667*hstep + 1) + (-0.25*hstep + 1.0)*(-0.125*hstep + 0.5))/(2*Time+1);

   values[6*(Time+1)+6*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.166666666666667*hstep*Xval[1]*(-0.333333333333333*hstep - 2) - 0.125*hstep*Xval[1]*(0.25*hstep + 1.0))/(2*Time+1);

   values[3*(Time+1)+4*(Time)+0+ 1+jt] += obj_factor*bounds[2][3]*(-0.166666666666667*hstep*Xvalp1[1]*(-0.333333333333333*hstep - 2) + 0.125*hstep*Xvalp1[1]*(0.25*hstep + 1.0))/(2*Time+1);

   values[6*(Time+1)+7*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep - 2)*(-Xval[0] + Xval[3]) + 0.125*hstep*(0.25*hstep + 1.0)*(-Xval[0] + Xval[3]))/(2*Time+1);

   values[4*(Time+1)+5*(Time)+0+ 1+jt] += obj_factor*bounds[2][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep - 2)*(-Xvalp1[0] + Xvalp1[3]) + 0.125*hstep*(0.25*hstep + 1.0)*(Xvalp1[0] - Xvalp1[3]))/(2*Time+1);

   values[6*(Time+1)+8*(Time)+0+jt] += obj_factor*bounds[2][3]*((-0.333333333333333*hstep - 2)*(-0.166666666666667*hstep + 1) + (-0.125*hstep + 0.5)*(0.25*hstep + 1.0))/(2*Time+1);

   values[5*(Time+1)+6*(Time)+0+ 1+jt] += obj_factor*bounds[2][3]*((-0.333333333333333*hstep - 2)*(-0.166666666666667*hstep - 1) + (0.125*hstep + 0.5)*(0.25*hstep + 1.0))/(2*Time+1);

   values[6*(Time+1)+9*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.0868055555555556*pow(hstep, 2)*pow(Xval[1], 2))/(2*Time+1);

   values[7*(Time+1)+9*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.0243055555555556*pow(hstep, 2)*Xval[1]*Xvalp1[1])/(2*Time+1);

   values[7*(Time+1)+10*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.0868055555555556*pow(hstep, 2)*Xval[1]*(-Xval[0] + Xval[3]) + 0.333333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[1]*(-Xval[0] + Xval[3]) + 4*Xval2[1]*(-Xval2[0] + Xval2[3]) + Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] - 4*Xval2[2] - Xvalp1[2]) + Xval[2] - Xvalp1[2]) + 0.25*hstep*(0.125*hstep*(Xval[1]*(-Xval[0] + Xval[3]) - Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] + Xvalp1[2]) + 0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2]))/(2*Time+1);

   values[8*(Time+1)+10*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.0555555555555556*pow(hstep, 2)*Xval[1]*(-Xvalp1[0] + Xvalp1[3]) + 0.03125*pow(hstep, 2)*Xval[1]*(Xvalp1[0] - Xvalp1[3]))/(2*Time+1);

   values[8*(Time+1)+11*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.333333333333333*hstep*Xval[1]*(-0.166666666666667*hstep + 1) + 0.25*hstep*Xval[1]*(-0.125*hstep + 0.5))/(2*Time+1);

   values[9*(Time+1)+11*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.333333333333333*hstep*Xval[1]*(-0.166666666666667*hstep - 1) + 0.25*hstep*Xval[1]*(0.125*hstep + 0.5))/(2*Time+1);

   values[9*(Time+1)+12*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.0868055555555556*pow(hstep, 2)*pow(Xval[1], 2))/(2*Time+1);

   values[10*(Time+1)+12*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.0243055555555556*pow(hstep, 2)*Xval[1]*Xvalp1[1])/(2*Time+1);

   values[6*(Time+1)+9*(Time)+0+ 1+jt] += obj_factor*bounds[2][3]*(-0.0868055555555556*pow(hstep, 2)*pow(Xvalp1[1], 2))/(2*Time+1);

   values[10*(Time+1)+13*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.0243055555555556*pow(hstep, 2)*Xvalp1[1]*(-Xval[0] + Xval[3]))/(2*Time+1);

   values[7*(Time+1)+10*(Time)+0+ 1+jt] += obj_factor*bounds[2][3]*(0.0555555555555556*pow(hstep, 2)*Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - 0.03125*pow(hstep, 2)*Xvalp1[1]*(Xvalp1[0] - Xvalp1[3]) + 0.333333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[1]*(-Xval[0] + Xval[3]) + 4*Xval2[1]*(-Xval2[0] + Xval2[3]) + Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] - 4*Xval2[2] - Xvalp1[2]) + Xval[2] - Xvalp1[2]) - 0.25*hstep*(0.125*hstep*(Xval[1]*(-Xval[0] + Xval[3]) - Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] + Xvalp1[2]) + 0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2]))/(2*Time+1);

   values[10*(Time+1)+14*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.333333333333333*hstep*Xvalp1[1]*(-0.166666666666667*hstep + 1) - 0.25*hstep*Xvalp1[1]*(-0.125*hstep + 0.5))/(2*Time+1);

   values[8*(Time+1)+11*(Time)+0+ 1+jt] += obj_factor*bounds[2][3]*(0.333333333333333*hstep*Xvalp1[1]*(-0.166666666666667*hstep - 1) - 0.25*hstep*Xvalp1[1]*(0.125*hstep + 0.5))/(2*Time+1);

   values[10*(Time+1)+15*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.0243055555555556*pow(hstep, 2)*Xval[1]*Xvalp1[1])/(2*Time+1);

   values[9*(Time+1)+12*(Time)+0+ 1+jt] += obj_factor*bounds[2][3]*(0.0868055555555556*pow(hstep, 2)*pow(Xvalp1[1], 2))/(2*Time+1);

   values[15*(Time+1)+25*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.222222222222222*pow(hstep, 2)*Xval[1]*Xval2[1])/(2*Time+1);

   values[15*(Time+1)+26*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.222222222222222*pow(hstep, 2)*Xval2[1]*Xvalp1[1])/(2*Time+1);

   values[15*(Time+1)+27*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[1]*(-Xval[0] + Xval[3]))/(2*Time+1);

   values[15*(Time+1)+28*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[1]*(-Xvalp1[0] + Xvalp1[3]))/(2*Time+1);

   values[15*(Time+1)+29*(Time)+0+jt] += obj_factor*bounds[2][3]*(-1.33333333333333*hstep*Xval2[1]*(-0.166666666666667*hstep + 1))/(2*Time+1);

   values[15*(Time+1)+30*(Time)+0+jt] += obj_factor*bounds[2][3]*(-1.33333333333333*hstep*Xval2[1]*(-0.166666666666667*hstep - 1))/(2*Time+1);

   values[15*(Time+1)+31*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.222222222222222*pow(hstep, 2)*Xval[1]*Xval2[1])/(2*Time+1);

   values[15*(Time+1)+32*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[1]*Xvalp1[1])/(2*Time+1);

   values[15*(Time+1)+35*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.888888888888889*pow(hstep, 2)*pow(Xval2[1], 2))/(2*Time+1);

   values[15*(Time+1)+36*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.0555555555555556*pow(hstep, 2)*Xval[1]*(-4*Xval2[0] + 4*Xval2[3]))/(2*Time+1);

   values[15*(Time+1)+37*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.0555555555555556*pow(hstep, 2)*Xvalp1[1]*(-4*Xval2[0] + 4*Xval2[3]))/(2*Time+1);

   values[15*(Time+1)+38*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.0555555555555556*pow(hstep, 2)*(-Xval[0] + Xval[3])*(-4*Xval2[0] + 4*Xval2[3]))/(2*Time+1);

   values[15*(Time+1)+39*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.0555555555555556*pow(hstep, 2)*(-4*Xval2[0] + 4*Xval2[3])*(-Xvalp1[0] + Xvalp1[3]))/(2*Time+1);

   values[15*(Time+1)+40*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.333333333333333*hstep*(-0.166666666666667*hstep + 1)*(-4*Xval2[0] + 4*Xval2[3]))/(2*Time+1);

   values[15*(Time+1)+41*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.333333333333333*hstep*(-0.166666666666667*hstep - 1)*(-4*Xval2[0] + 4*Xval2[3]))/(2*Time+1);

   values[15*(Time+1)+42*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.0555555555555556*pow(hstep, 2)*Xval[1]*(-4*Xval2[0] + 4*Xval2[3]))/(2*Time+1);

   values[15*(Time+1)+43*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.0555555555555556*pow(hstep, 2)*Xvalp1[1]*(-4*Xval2[0] + 4*Xval2[3]))/(2*Time+1);

   values[15*(Time+1)+46*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[1]*(-4*Xval2[0] + 4*Xval2[3]) - 1.33333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[1]*(-Xval[0] + Xval[3]) + 4*Xval2[1]*(-Xval2[0] + Xval2[3]) + Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] - 4*Xval2[2] - Xvalp1[2]) + Xval[2] - Xvalp1[2]))/(2*Time+1);

   values[15*(Time+1)+47*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.0555555555555556*pow(hstep, 2)*pow(-4*Xval2[0] + 4*Xval2[3], 2))/(2*Time+1);

   values[15*(Time+1)+48*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.222222222222222*pow(hstep, 2)*Xval[1] + 0.25*hstep*Xval[1])/(2*Time+1);

   values[15*(Time+1)+49*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.222222222222222*pow(hstep, 2)*Xvalp1[1] - 0.25*hstep*Xvalp1[1])/(2*Time+1);

   values[15*(Time+1)+50*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.222222222222222*pow(hstep, 2)*(-Xval[0] + Xval[3]) - 0.25*hstep*(-Xval[0] + Xval[3]))/(2*Time+1);

   values[15*(Time+1)+51*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.222222222222222*pow(hstep, 2)*(-Xvalp1[0] + Xvalp1[3]) - 0.25*hstep*(Xvalp1[0] - Xvalp1[3]))/(2*Time+1);

   values[15*(Time+1)+52*(Time)+0+jt] += obj_factor*bounds[2][3]*(-1.33333333333333*hstep*(-0.166666666666667*hstep + 1) + 0.25*hstep - 1.0)/(2*Time+1);

   values[15*(Time+1)+53*(Time)+0+jt] += obj_factor*bounds[2][3]*(-1.33333333333333*hstep*(-0.166666666666667*hstep - 1) - 0.25*hstep - 1.0)/(2*Time+1);

   values[15*(Time+1)+54*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.222222222222222*pow(hstep, 2)*Xval[1] - 0.25*hstep*Xval[1])/(2*Time+1);

   values[15*(Time+1)+55*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.222222222222222*pow(hstep, 2)*Xvalp1[1] + 0.25*hstep*Xvalp1[1])/(2*Time+1);

   values[15*(Time+1)+58*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.888888888888889*pow(hstep, 2)*Xval2[1])/(2*Time+1);

   values[15*(Time+1)+59*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.222222222222222*pow(hstep, 2)*(-4*Xval2[0] + 4*Xval2[3]))/(2*Time+1);

   values[15*(Time+1)+60*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.888888888888889*pow(hstep, 2) + 2)/(2*Time+1);

   values[15*(Time+1)+61*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.222222222222222*pow(hstep, 2)*Xval[1]*Xval2[1])/(2*Time+1);

   values[15*(Time+1)+62*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[1]*Xvalp1[1])/(2*Time+1);

   values[15*(Time+1)+63*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.222222222222222*pow(hstep, 2)*Xval2[1]*(-Xval[0] + Xval[3]))/(2*Time+1);

   values[15*(Time+1)+64*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.222222222222222*pow(hstep, 2)*Xval2[1]*(-Xvalp1[0] + Xvalp1[3]))/(2*Time+1);

   values[15*(Time+1)+65*(Time)+0+jt] += obj_factor*bounds[2][3]*(1.33333333333333*hstep*Xval2[1]*(-0.166666666666667*hstep + 1))/(2*Time+1);

   values[15*(Time+1)+66*(Time)+0+jt] += obj_factor*bounds[2][3]*(1.33333333333333*hstep*Xval2[1]*(-0.166666666666667*hstep - 1))/(2*Time+1);

   values[15*(Time+1)+67*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.222222222222222*pow(hstep, 2)*Xval[1]*Xval2[1])/(2*Time+1);

   values[15*(Time+1)+68*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.222222222222222*pow(hstep, 2)*Xval2[1]*Xvalp1[1])/(2*Time+1);

   values[15*(Time+1)+71*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.888888888888889*pow(hstep, 2)*pow(Xval2[1], 2))/(2*Time+1);

   values[15*(Time+1)+72*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.222222222222222*pow(hstep, 2)*Xval2[1]*(-4*Xval2[0] + 4*Xval2[3]) + 1.33333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[1]*(-Xval[0] + Xval[3]) + 4*Xval2[1]*(-Xval2[0] + Xval2[3]) + Xvalp1[1]*(-Xvalp1[0] + Xvalp1[3]) - Xval[2] - 4*Xval2[2] - Xvalp1[2]) + Xval[2] - Xvalp1[2]))/(2*Time+1);

   values[15*(Time+1)+73*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.888888888888889*pow(hstep, 2)*Xval2[1])/(2*Time+1);

   values[15*(Time+1)+74*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.888888888888889*pow(hstep, 2)*pow(Xval2[1], 2))/(2*Time+1);

   values[15*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[2][3]*(-0.333333333333333*pow(hstep, 2)*Xval[1])/(2*Time+1);

   values[15*(Time+1)+90*(Time)+0+ 1+jt] += obj_factor*bounds[2][3]*(-0.333333333333333*pow(hstep, 2)*Xvalp1[1])/(2*Time+1);

   values[16*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.333333333333333*pow(hstep, 2)*(-Xval[0] + Xval[3]))/(2*Time+1);

   values[16*(Time+1)+90*(Time)+0+ 1+jt] += obj_factor*bounds[2][3]*(0.333333333333333*pow(hstep, 2)*(-Xvalp1[0] + Xvalp1[3]))/(2*Time+1);

   values[17*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[2][3]*(2.0*hstep*(-0.166666666666667*hstep + 1))/(2*Time+1);

   values[17*(Time+1)+90*(Time)+0+ 1+jt] += obj_factor*bounds[2][3]*(2.0*hstep*(-0.166666666666667*hstep - 1))/(2*Time+1);

   values[18*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.333333333333333*pow(hstep, 2)*Xval[1])/(2*Time+1);

   values[18*(Time+1)+90*(Time)+0+ 1+jt] += obj_factor*bounds[2][3]*(0.333333333333333*pow(hstep, 2)*Xvalp1[1])/(2*Time+1);

   values[20*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[2][3]*(-1.33333333333333*pow(hstep, 2)*Xval2[1])/(2*Time+1);

   values[20*(Time+1)+91*(Time)+0+jt] += obj_factor*bounds[2][3]*(0.333333333333333*pow(hstep, 2)*(-4*Xval2[0] + 4*Xval2[3]))/(2*Time+1);

   values[20*(Time+1)+92*(Time)+0+jt] += obj_factor*bounds[2][3]*(-1.33333333333333*pow(hstep, 2))/(2*Time+1);

   values[20*(Time+1)+93*(Time)+0+jt] += obj_factor*bounds[2][3]*(1.33333333333333*pow(hstep, 2)*Xval2[1])/(2*Time+1);

   values[2*(Time+1)+2*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.0868055555555556*pow(hstep, 2)*pow(Xval[2], 2))/(2*Time+1);

   values[3*(Time+1)+3*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.0243055555555556*pow(hstep, 2)*Xval[2]*Xvalp1[2])/(2*Time+1);

   values[2*(Time+1)+2*(Time)+0+ 1+jt] += obj_factor*bounds[3][3]*(0.0868055555555556*pow(hstep, 2)*pow(Xvalp1[2], 2))/(2*Time+1);

   values[4*(Time+1)+5*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.0868055555555556*pow(hstep, 2)*Xval[2]*(-Xval[1] + Xval[4]) - 0.333333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[2]*(-Xval[1] + Xval[4]) + 4*Xval2[2]*(-Xval2[1] + Xval2[4]) + Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] - 4*Xval2[3] - Xvalp1[3]) + Xval[3] - Xvalp1[3]) - 0.25*hstep*(0.125*hstep*(Xval[2]*(-Xval[1] + Xval[4]) - Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] + Xvalp1[3]) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3]))/(2*Time+1);

   values[5*(Time+1)+5*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.0243055555555556*pow(hstep, 2)*Xvalp1[2]*(-Xval[1] + Xval[4]))/(2*Time+1);

   values[5*(Time+1)+6*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.0868055555555556*pow(hstep, 2)*pow(-Xval[1] + Xval[4], 2))/(2*Time+1);

   values[6*(Time+1)+7*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.0555555555555556*pow(hstep, 2)*Xval[2]*(-Xvalp1[1] + Xvalp1[4]) - 0.03125*pow(hstep, 2)*Xval[2]*(Xvalp1[1] - Xvalp1[4]))/(2*Time+1);

   values[4*(Time+1)+5*(Time)+0+ 1+jt] += obj_factor*bounds[3][3]*(-0.0555555555555556*pow(hstep, 2)*Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) + 0.03125*pow(hstep, 2)*Xvalp1[2]*(Xvalp1[1] - Xvalp1[4]) - 0.333333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[2]*(-Xval[1] + Xval[4]) + 4*Xval2[2]*(-Xval2[1] + Xval2[4]) + Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] - 4*Xval2[3] - Xvalp1[3]) + Xval[3] - Xvalp1[3]) + 0.25*hstep*(0.125*hstep*(Xval[2]*(-Xval[1] + Xval[4]) - Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] + Xvalp1[3]) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3]))/(2*Time+1);

   values[6*(Time+1)+8*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.0555555555555556*pow(hstep, 2)*(-Xval[1] + Xval[4])*(-Xvalp1[1] + Xvalp1[4]) + 0.03125*pow(hstep, 2)*(-Xval[1] + Xval[4])*(Xvalp1[1] - Xvalp1[4]))/(2*Time+1);

   values[5*(Time+1)+6*(Time)+0+ 1+jt] += obj_factor*bounds[3][3]*(0.0555555555555556*pow(hstep, 2)*pow(-Xvalp1[1] + Xvalp1[4], 2) + 0.03125*pow(hstep, 2)*pow(Xvalp1[1] - Xvalp1[4], 2))/(2*Time+1);

   values[7*(Time+1)+10*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.166666666666667*hstep*Xval[2]*(-0.333333333333333*hstep + 2) - 0.125*hstep*Xval[2]*(-0.25*hstep + 1.0))/(2*Time+1);

   values[8*(Time+1)+10*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.166666666666667*hstep*Xvalp1[2]*(-0.333333333333333*hstep + 2) + 0.125*hstep*Xvalp1[2]*(-0.25*hstep + 1.0))/(2*Time+1);

   values[8*(Time+1)+11*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep + 2)*(-Xval[1] + Xval[4]) + 0.125*hstep*(-0.25*hstep + 1.0)*(-Xval[1] + Xval[4]))/(2*Time+1);

   values[9*(Time+1)+11*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep + 2)*(-Xvalp1[1] + Xvalp1[4]) + 0.125*hstep*(-0.25*hstep + 1.0)*(Xvalp1[1] - Xvalp1[4]))/(2*Time+1);

   values[9*(Time+1)+12*(Time)+0+jt] += obj_factor*bounds[3][3]*((-0.333333333333333*hstep + 2)*(-0.166666666666667*hstep + 1) + (-0.25*hstep + 1.0)*(-0.125*hstep + 0.5))/(2*Time+1);

   values[10*(Time+1)+13*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.166666666666667*hstep*Xval[2]*(-0.333333333333333*hstep - 2) - 0.125*hstep*Xval[2]*(0.25*hstep + 1.0))/(2*Time+1);

   values[7*(Time+1)+10*(Time)+0+ 1+jt] += obj_factor*bounds[3][3]*(-0.166666666666667*hstep*Xvalp1[2]*(-0.333333333333333*hstep - 2) + 0.125*hstep*Xvalp1[2]*(0.25*hstep + 1.0))/(2*Time+1);

   values[10*(Time+1)+14*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep - 2)*(-Xval[1] + Xval[4]) + 0.125*hstep*(0.25*hstep + 1.0)*(-Xval[1] + Xval[4]))/(2*Time+1);

   values[8*(Time+1)+11*(Time)+0+ 1+jt] += obj_factor*bounds[3][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep - 2)*(-Xvalp1[1] + Xvalp1[4]) + 0.125*hstep*(0.25*hstep + 1.0)*(Xvalp1[1] - Xvalp1[4]))/(2*Time+1);

   values[10*(Time+1)+15*(Time)+0+jt] += obj_factor*bounds[3][3]*((-0.333333333333333*hstep - 2)*(-0.166666666666667*hstep + 1) + (-0.125*hstep + 0.5)*(0.25*hstep + 1.0))/(2*Time+1);

   values[9*(Time+1)+12*(Time)+0+ 1+jt] += obj_factor*bounds[3][3]*((-0.333333333333333*hstep - 2)*(-0.166666666666667*hstep - 1) + (0.125*hstep + 0.5)*(0.25*hstep + 1.0))/(2*Time+1);

   values[11*(Time+1)+17*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.0868055555555556*pow(hstep, 2)*pow(Xval[2], 2))/(2*Time+1);

   values[12*(Time+1)+17*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.0243055555555556*pow(hstep, 2)*Xval[2]*Xvalp1[2])/(2*Time+1);

   values[12*(Time+1)+18*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.0868055555555556*pow(hstep, 2)*Xval[2]*(-Xval[1] + Xval[4]) + 0.333333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[2]*(-Xval[1] + Xval[4]) + 4*Xval2[2]*(-Xval2[1] + Xval2[4]) + Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] - 4*Xval2[3] - Xvalp1[3]) + Xval[3] - Xvalp1[3]) + 0.25*hstep*(0.125*hstep*(Xval[2]*(-Xval[1] + Xval[4]) - Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] + Xvalp1[3]) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3]))/(2*Time+1);

   values[13*(Time+1)+18*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.0555555555555556*pow(hstep, 2)*Xval[2]*(-Xvalp1[1] + Xvalp1[4]) + 0.03125*pow(hstep, 2)*Xval[2]*(Xvalp1[1] - Xvalp1[4]))/(2*Time+1);

   values[13*(Time+1)+19*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.333333333333333*hstep*Xval[2]*(-0.166666666666667*hstep + 1) + 0.25*hstep*Xval[2]*(-0.125*hstep + 0.5))/(2*Time+1);

   values[14*(Time+1)+19*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.333333333333333*hstep*Xval[2]*(-0.166666666666667*hstep - 1) + 0.25*hstep*Xval[2]*(0.125*hstep + 0.5))/(2*Time+1);

   values[14*(Time+1)+20*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.0868055555555556*pow(hstep, 2)*pow(Xval[2], 2))/(2*Time+1);

   values[15*(Time+1)+21*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.0243055555555556*pow(hstep, 2)*Xval[2]*Xvalp1[2])/(2*Time+1);

   values[11*(Time+1)+17*(Time)+0+ 1+jt] += obj_factor*bounds[3][3]*(-0.0868055555555556*pow(hstep, 2)*pow(Xvalp1[2], 2))/(2*Time+1);

   values[15*(Time+1)+22*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.0243055555555556*pow(hstep, 2)*Xvalp1[2]*(-Xval[1] + Xval[4]))/(2*Time+1);

   values[12*(Time+1)+18*(Time)+0+ 1+jt] += obj_factor*bounds[3][3]*(0.0555555555555556*pow(hstep, 2)*Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - 0.03125*pow(hstep, 2)*Xvalp1[2]*(Xvalp1[1] - Xvalp1[4]) + 0.333333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[2]*(-Xval[1] + Xval[4]) + 4*Xval2[2]*(-Xval2[1] + Xval2[4]) + Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] - 4*Xval2[3] - Xvalp1[3]) + Xval[3] - Xvalp1[3]) - 0.25*hstep*(0.125*hstep*(Xval[2]*(-Xval[1] + Xval[4]) - Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] + Xvalp1[3]) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3]))/(2*Time+1);

   values[15*(Time+1)+23*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.333333333333333*hstep*Xvalp1[2]*(-0.166666666666667*hstep + 1) - 0.25*hstep*Xvalp1[2]*(-0.125*hstep + 0.5))/(2*Time+1);

   values[13*(Time+1)+19*(Time)+0+ 1+jt] += obj_factor*bounds[3][3]*(0.333333333333333*hstep*Xvalp1[2]*(-0.166666666666667*hstep - 1) - 0.25*hstep*Xvalp1[2]*(0.125*hstep + 0.5))/(2*Time+1);

   values[15*(Time+1)+24*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.0243055555555556*pow(hstep, 2)*Xval[2]*Xvalp1[2])/(2*Time+1);

   values[14*(Time+1)+20*(Time)+0+ 1+jt] += obj_factor*bounds[3][3]*(0.0868055555555556*pow(hstep, 2)*pow(Xvalp1[2], 2))/(2*Time+1);

   values[15*(Time+1)+38*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.222222222222222*pow(hstep, 2)*Xval[2]*Xval2[2])/(2*Time+1);

   values[15*(Time+1)+39*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.222222222222222*pow(hstep, 2)*Xval2[2]*Xvalp1[2])/(2*Time+1);

   values[15*(Time+1)+40*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[2]*(-Xval[1] + Xval[4]))/(2*Time+1);

   values[15*(Time+1)+41*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[2]*(-Xvalp1[1] + Xvalp1[4]))/(2*Time+1);

   values[15*(Time+1)+42*(Time)+0+jt] += obj_factor*bounds[3][3]*(-1.33333333333333*hstep*Xval2[2]*(-0.166666666666667*hstep + 1))/(2*Time+1);

   values[15*(Time+1)+43*(Time)+0+jt] += obj_factor*bounds[3][3]*(-1.33333333333333*hstep*Xval2[2]*(-0.166666666666667*hstep - 1))/(2*Time+1);

   values[15*(Time+1)+44*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.222222222222222*pow(hstep, 2)*Xval[2]*Xval2[2])/(2*Time+1);

   values[15*(Time+1)+45*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[2]*Xvalp1[2])/(2*Time+1);

   values[15*(Time+1)+47*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.888888888888889*pow(hstep, 2)*pow(Xval2[2], 2))/(2*Time+1);

   values[15*(Time+1)+50*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.0555555555555556*pow(hstep, 2)*Xval[2]*(-4*Xval2[1] + 4*Xval2[4]))/(2*Time+1);

   values[15*(Time+1)+51*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.0555555555555556*pow(hstep, 2)*Xvalp1[2]*(-4*Xval2[1] + 4*Xval2[4]))/(2*Time+1);

   values[15*(Time+1)+52*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.0555555555555556*pow(hstep, 2)*(-Xval[1] + Xval[4])*(-4*Xval2[1] + 4*Xval2[4]))/(2*Time+1);

   values[15*(Time+1)+53*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.0555555555555556*pow(hstep, 2)*(-4*Xval2[1] + 4*Xval2[4])*(-Xvalp1[1] + Xvalp1[4]))/(2*Time+1);

   values[15*(Time+1)+54*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.333333333333333*hstep*(-0.166666666666667*hstep + 1)*(-4*Xval2[1] + 4*Xval2[4]))/(2*Time+1);

   values[15*(Time+1)+55*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.333333333333333*hstep*(-0.166666666666667*hstep - 1)*(-4*Xval2[1] + 4*Xval2[4]))/(2*Time+1);

   values[15*(Time+1)+56*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.0555555555555556*pow(hstep, 2)*Xval[2]*(-4*Xval2[1] + 4*Xval2[4]))/(2*Time+1);

   values[15*(Time+1)+57*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.0555555555555556*pow(hstep, 2)*Xvalp1[2]*(-4*Xval2[1] + 4*Xval2[4]))/(2*Time+1);

   values[15*(Time+1)+59*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[2]*(-4*Xval2[1] + 4*Xval2[4]) - 1.33333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[2]*(-Xval[1] + Xval[4]) + 4*Xval2[2]*(-Xval2[1] + Xval2[4]) + Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] - 4*Xval2[3] - Xvalp1[3]) + Xval[3] - Xvalp1[3]))/(2*Time+1);

   values[15*(Time+1)+60*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.0555555555555556*pow(hstep, 2)*pow(-4*Xval2[1] + 4*Xval2[4], 2))/(2*Time+1);

   values[15*(Time+1)+63*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.222222222222222*pow(hstep, 2)*Xval[2] + 0.25*hstep*Xval[2])/(2*Time+1);

   values[15*(Time+1)+64*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.222222222222222*pow(hstep, 2)*Xvalp1[2] - 0.25*hstep*Xvalp1[2])/(2*Time+1);

   values[15*(Time+1)+65*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.222222222222222*pow(hstep, 2)*(-Xval[1] + Xval[4]) - 0.25*hstep*(-Xval[1] + Xval[4]))/(2*Time+1);

   values[15*(Time+1)+66*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.222222222222222*pow(hstep, 2)*(-Xvalp1[1] + Xvalp1[4]) - 0.25*hstep*(Xvalp1[1] - Xvalp1[4]))/(2*Time+1);

   values[15*(Time+1)+67*(Time)+0+jt] += obj_factor*bounds[3][3]*(-1.33333333333333*hstep*(-0.166666666666667*hstep + 1) + 0.25*hstep - 1.0)/(2*Time+1);

   values[15*(Time+1)+68*(Time)+0+jt] += obj_factor*bounds[3][3]*(-1.33333333333333*hstep*(-0.166666666666667*hstep - 1) - 0.25*hstep - 1.0)/(2*Time+1);

   values[15*(Time+1)+69*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.222222222222222*pow(hstep, 2)*Xval[2] - 0.25*hstep*Xval[2])/(2*Time+1);

   values[15*(Time+1)+70*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.222222222222222*pow(hstep, 2)*Xvalp1[2] + 0.25*hstep*Xvalp1[2])/(2*Time+1);

   values[15*(Time+1)+72*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.888888888888889*pow(hstep, 2)*Xval2[2])/(2*Time+1);

   values[15*(Time+1)+73*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.222222222222222*pow(hstep, 2)*(-4*Xval2[1] + 4*Xval2[4]))/(2*Time+1);

   values[15*(Time+1)+74*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.888888888888889*pow(hstep, 2) + 2)/(2*Time+1);

   values[15*(Time+1)+77*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.222222222222222*pow(hstep, 2)*Xval[2]*Xval2[2])/(2*Time+1);

   values[15*(Time+1)+78*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[2]*Xvalp1[2])/(2*Time+1);

   values[15*(Time+1)+79*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.222222222222222*pow(hstep, 2)*Xval2[2]*(-Xval[1] + Xval[4]))/(2*Time+1);

   values[15*(Time+1)+80*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.222222222222222*pow(hstep, 2)*Xval2[2]*(-Xvalp1[1] + Xvalp1[4]))/(2*Time+1);

   values[15*(Time+1)+81*(Time)+0+jt] += obj_factor*bounds[3][3]*(1.33333333333333*hstep*Xval2[2]*(-0.166666666666667*hstep + 1))/(2*Time+1);

   values[15*(Time+1)+82*(Time)+0+jt] += obj_factor*bounds[3][3]*(1.33333333333333*hstep*Xval2[2]*(-0.166666666666667*hstep - 1))/(2*Time+1);

   values[15*(Time+1)+83*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.222222222222222*pow(hstep, 2)*Xval[2]*Xval2[2])/(2*Time+1);

   values[15*(Time+1)+84*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.222222222222222*pow(hstep, 2)*Xval2[2]*Xvalp1[2])/(2*Time+1);

   values[15*(Time+1)+86*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.888888888888889*pow(hstep, 2)*pow(Xval2[2], 2))/(2*Time+1);

   values[15*(Time+1)+87*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.222222222222222*pow(hstep, 2)*Xval2[2]*(-4*Xval2[1] + 4*Xval2[4]) + 1.33333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[2]*(-Xval[1] + Xval[4]) + 4*Xval2[2]*(-Xval2[1] + Xval2[4]) + Xvalp1[2]*(-Xvalp1[1] + Xvalp1[4]) - Xval[3] - 4*Xval2[3] - Xvalp1[3]) + Xval[3] - Xvalp1[3]))/(2*Time+1);

   values[15*(Time+1)+88*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.888888888888889*pow(hstep, 2)*Xval2[2])/(2*Time+1);

   values[15*(Time+1)+89*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.888888888888889*pow(hstep, 2)*pow(Xval2[2], 2))/(2*Time+1);

   values[16*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[3][3]*(-0.333333333333333*pow(hstep, 2)*Xval[2])/(2*Time+1);

   values[16*(Time+1)+90*(Time)+0+ 1+jt] += obj_factor*bounds[3][3]*(-0.333333333333333*pow(hstep, 2)*Xvalp1[2])/(2*Time+1);

   values[17*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.333333333333333*pow(hstep, 2)*(-Xval[1] + Xval[4]))/(2*Time+1);

   values[17*(Time+1)+90*(Time)+0+ 1+jt] += obj_factor*bounds[3][3]*(0.333333333333333*pow(hstep, 2)*(-Xvalp1[1] + Xvalp1[4]))/(2*Time+1);

   values[18*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[3][3]*(2.0*hstep*(-0.166666666666667*hstep + 1))/(2*Time+1);

   values[18*(Time+1)+90*(Time)+0+ 1+jt] += obj_factor*bounds[3][3]*(2.0*hstep*(-0.166666666666667*hstep - 1))/(2*Time+1);

   values[19*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.333333333333333*pow(hstep, 2)*Xval[2])/(2*Time+1);

   values[19*(Time+1)+90*(Time)+0+ 1+jt] += obj_factor*bounds[3][3]*(0.333333333333333*pow(hstep, 2)*Xvalp1[2])/(2*Time+1);

   values[20*(Time+1)+91*(Time)+0+jt] += obj_factor*bounds[3][3]*(-1.33333333333333*pow(hstep, 2)*Xval2[2])/(2*Time+1);

   values[20*(Time+1)+92*(Time)+0+jt] += obj_factor*bounds[3][3]*(0.333333333333333*pow(hstep, 2)*(-4*Xval2[1] + 4*Xval2[4]))/(2*Time+1);

   values[20*(Time+1)+93*(Time)+0+jt] += obj_factor*bounds[3][3]*(-1.33333333333333*pow(hstep, 2))/(2*Time+1);

   values[20*(Time+1)+94*(Time)+0+jt] += obj_factor*bounds[3][3]*(1.33333333333333*pow(hstep, 2)*Xval2[2])/(2*Time+1);

   values[0*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.0868055555555556*pow(hstep, 2)*pow(Xval[3], 2))/(2*Time+1);

   values[1*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.0243055555555556*pow(hstep, 2)*Xval[3]*Xvalp1[3])/(2*Time+1);

   values[0*(Time+1)+0*(Time)+0+ 1+jt] += obj_factor*bounds[4][3]*(0.0868055555555556*pow(hstep, 2)*pow(Xvalp1[3], 2))/(2*Time+1);

   values[3*(Time+1)+4*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.0868055555555556*pow(hstep, 2)*pow(Xval[3], 2))/(2*Time+1);

   values[4*(Time+1)+4*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.0243055555555556*pow(hstep, 2)*Xval[3]*Xvalp1[3])/(2*Time+1);

   values[5*(Time+1)+6*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.0868055555555556*pow(hstep, 2)*pow(Xval[3], 2))/(2*Time+1);

   values[6*(Time+1)+6*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.0243055555555556*pow(hstep, 2)*Xval[3]*Xvalp1[3])/(2*Time+1);

   values[3*(Time+1)+4*(Time)+0+ 1+jt] += obj_factor*bounds[4][3]*(-0.0868055555555556*pow(hstep, 2)*pow(Xvalp1[3], 2))/(2*Time+1);

   values[6*(Time+1)+8*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.0243055555555556*pow(hstep, 2)*Xval[3]*Xvalp1[3])/(2*Time+1);

   values[5*(Time+1)+6*(Time)+0+ 1+jt] += obj_factor*bounds[4][3]*(0.0868055555555556*pow(hstep, 2)*pow(Xvalp1[3], 2))/(2*Time+1);

   values[6*(Time+1)+9*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.0868055555555556*pow(hstep, 2)*Xval[3]*(Xval[0] - Xval[2]) + 0.333333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[3]*(Xval[0] - Xval[2]) + 4*Xval2[3]*(Xval2[0] - Xval2[2]) + Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] - 4*Xval2[4] - Xvalp1[4]) + Xval[4] - Xvalp1[4]) + 0.25*hstep*(0.125*hstep*(Xval[3]*(Xval[0] - Xval[2]) - Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] + Xvalp1[4]) + 0.5*Xval[4] - Xval2[4] + 0.5*Xvalp1[4]))/(2*Time+1);

   values[7*(Time+1)+9*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.0243055555555556*pow(hstep, 2)*Xvalp1[3]*(Xval[0] - Xval[2]))/(2*Time+1);

   values[8*(Time+1)+11*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.0868055555555556*pow(hstep, 2)*Xval[3]*(Xval[0] - Xval[2]) - 0.333333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[3]*(Xval[0] - Xval[2]) + 4*Xval2[3]*(Xval2[0] - Xval2[2]) + Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] - 4*Xval2[4] - Xvalp1[4]) + Xval[4] - Xvalp1[4]) - 0.25*hstep*(0.125*hstep*(Xval[3]*(Xval[0] - Xval[2]) - Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] + Xvalp1[4]) + 0.5*Xval[4] - Xval2[4] + 0.5*Xvalp1[4]))/(2*Time+1);

   values[9*(Time+1)+11*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.0243055555555556*pow(hstep, 2)*Xvalp1[3]*(Xval[0] - Xval[2]))/(2*Time+1);

   values[9*(Time+1)+12*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.0868055555555556*pow(hstep, 2)*pow(Xval[0] - Xval[2], 2))/(2*Time+1);

   values[10*(Time+1)+12*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.03125*pow(hstep, 2)*Xval[3]*(-Xvalp1[0] + Xvalp1[2]) + 0.0555555555555556*pow(hstep, 2)*Xval[3]*(Xvalp1[0] - Xvalp1[2]))/(2*Time+1);

   values[6*(Time+1)+9*(Time)+0+ 1+jt] += obj_factor*bounds[4][3]*(-0.03125*pow(hstep, 2)*Xvalp1[3]*(-Xvalp1[0] + Xvalp1[2]) + 0.0555555555555556*pow(hstep, 2)*Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) + 0.333333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[3]*(Xval[0] - Xval[2]) + 4*Xval2[3]*(Xval2[0] - Xval2[2]) + Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] - 4*Xval2[4] - Xvalp1[4]) + Xval[4] - Xvalp1[4]) - 0.25*hstep*(0.125*hstep*(Xval[3]*(Xval[0] - Xval[2]) - Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] + Xvalp1[4]) + 0.5*Xval[4] - Xval2[4] + 0.5*Xvalp1[4]))/(2*Time+1);

   values[10*(Time+1)+14*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.03125*pow(hstep, 2)*Xval[3]*(-Xvalp1[0] + Xvalp1[2]) - 0.0555555555555556*pow(hstep, 2)*Xval[3]*(Xvalp1[0] - Xvalp1[2]))/(2*Time+1);

   values[8*(Time+1)+11*(Time)+0+ 1+jt] += obj_factor*bounds[4][3]*(0.03125*pow(hstep, 2)*Xvalp1[3]*(-Xvalp1[0] + Xvalp1[2]) - 0.0555555555555556*pow(hstep, 2)*Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - 0.333333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[3]*(Xval[0] - Xval[2]) + 4*Xval2[3]*(Xval2[0] - Xval2[2]) + Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] - 4*Xval2[4] - Xvalp1[4]) + Xval[4] - Xvalp1[4]) + 0.25*hstep*(0.125*hstep*(Xval[3]*(Xval[0] - Xval[2]) - Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] + Xvalp1[4]) + 0.5*Xval[4] - Xval2[4] + 0.5*Xvalp1[4]))/(2*Time+1);

   values[10*(Time+1)+15*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.03125*pow(hstep, 2)*(Xval[0] - Xval[2])*(-Xvalp1[0] + Xvalp1[2]) + 0.0555555555555556*pow(hstep, 2)*(Xval[0] - Xval[2])*(Xvalp1[0] - Xvalp1[2]))/(2*Time+1);

   values[9*(Time+1)+12*(Time)+0+ 1+jt] += obj_factor*bounds[4][3]*(0.03125*pow(hstep, 2)*pow(-Xvalp1[0] + Xvalp1[2], 2) + 0.0555555555555556*pow(hstep, 2)*pow(Xvalp1[0] - Xvalp1[2], 2))/(2*Time+1);

   values[10*(Time+1)+16*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.166666666666667*hstep*Xval[3]*(-0.333333333333333*hstep + 2) + 0.125*hstep*Xval[3]*(-0.25*hstep + 1.0))/(2*Time+1);

   values[11*(Time+1)+16*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.166666666666667*hstep*Xvalp1[3]*(-0.333333333333333*hstep + 2) - 0.125*hstep*Xvalp1[3]*(-0.25*hstep + 1.0))/(2*Time+1);

   values[12*(Time+1)+18*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.166666666666667*hstep*Xval[3]*(-0.333333333333333*hstep + 2) - 0.125*hstep*Xval[3]*(-0.25*hstep + 1.0))/(2*Time+1);

   values[13*(Time+1)+18*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.166666666666667*hstep*Xvalp1[3]*(-0.333333333333333*hstep + 2) + 0.125*hstep*Xvalp1[3]*(-0.25*hstep + 1.0))/(2*Time+1);

   values[13*(Time+1)+19*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep + 2)*(Xval[0] - Xval[2]) + 0.125*hstep*(-0.25*hstep + 1.0)*(Xval[0] - Xval[2]))/(2*Time+1);

   values[14*(Time+1)+19*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep + 2)*(Xvalp1[0] - Xvalp1[2]) + 0.125*hstep*(-0.25*hstep + 1.0)*(-Xvalp1[0] + Xvalp1[2]))/(2*Time+1);

   values[14*(Time+1)+20*(Time)+0+jt] += obj_factor*bounds[4][3]*((-0.333333333333333*hstep + 2)*(-0.166666666666667*hstep + 1) + (-0.25*hstep + 1.0)*(-0.125*hstep + 0.5))/(2*Time+1);

   values[15*(Time+1)+20*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.166666666666667*hstep*Xval[3]*(-0.333333333333333*hstep - 2) + 0.125*hstep*Xval[3]*(0.25*hstep + 1.0))/(2*Time+1);

   values[10*(Time+1)+16*(Time)+0+ 1+jt] += obj_factor*bounds[4][3]*(0.166666666666667*hstep*Xvalp1[3]*(-0.333333333333333*hstep - 2) - 0.125*hstep*Xvalp1[3]*(0.25*hstep + 1.0))/(2*Time+1);

   values[15*(Time+1)+22*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.166666666666667*hstep*Xval[3]*(-0.333333333333333*hstep - 2) - 0.125*hstep*Xval[3]*(0.25*hstep + 1.0))/(2*Time+1);

   values[12*(Time+1)+18*(Time)+0+ 1+jt] += obj_factor*bounds[4][3]*(-0.166666666666667*hstep*Xvalp1[3]*(-0.333333333333333*hstep - 2) + 0.125*hstep*Xvalp1[3]*(0.25*hstep + 1.0))/(2*Time+1);

   values[15*(Time+1)+23*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep - 2)*(Xval[0] - Xval[2]) + 0.125*hstep*(0.25*hstep + 1.0)*(Xval[0] - Xval[2]))/(2*Time+1);

   values[13*(Time+1)+19*(Time)+0+ 1+jt] += obj_factor*bounds[4][3]*(0.166666666666667*hstep*(-0.333333333333333*hstep - 2)*(Xvalp1[0] - Xvalp1[2]) + 0.125*hstep*(0.25*hstep + 1.0)*(-Xvalp1[0] + Xvalp1[2]))/(2*Time+1);

   values[15*(Time+1)+24*(Time)+0+jt] += obj_factor*bounds[4][3]*((-0.333333333333333*hstep - 2)*(-0.166666666666667*hstep + 1) + (-0.125*hstep + 0.5)*(0.25*hstep + 1.0))/(2*Time+1);

   values[14*(Time+1)+20*(Time)+0+ 1+jt] += obj_factor*bounds[4][3]*((-0.333333333333333*hstep - 2)*(-0.166666666666667*hstep - 1) + (0.125*hstep + 0.5)*(0.25*hstep + 1.0))/(2*Time+1);

   values[15*(Time+1)+25*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.222222222222222*pow(hstep, 2)*Xval[3]*Xval2[3])/(2*Time+1);

   values[15*(Time+1)+26*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.222222222222222*pow(hstep, 2)*Xval2[3]*Xvalp1[3])/(2*Time+1);

   values[15*(Time+1)+29*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.222222222222222*pow(hstep, 2)*Xval[3]*Xval2[3])/(2*Time+1);

   values[15*(Time+1)+30*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[3]*Xvalp1[3])/(2*Time+1);

   values[15*(Time+1)+31*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.222222222222222*pow(hstep, 2)*Xval2[3]*(Xval[0] - Xval[2]))/(2*Time+1);

   values[15*(Time+1)+32*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.222222222222222*pow(hstep, 2)*Xval2[3]*(Xvalp1[0] - Xvalp1[2]))/(2*Time+1);

   values[15*(Time+1)+33*(Time)+0+jt] += obj_factor*bounds[4][3]*(1.33333333333333*hstep*Xval2[3]*(-0.166666666666667*hstep + 1))/(2*Time+1);

   values[15*(Time+1)+34*(Time)+0+jt] += obj_factor*bounds[4][3]*(1.33333333333333*hstep*Xval2[3]*(-0.166666666666667*hstep - 1))/(2*Time+1);

   values[15*(Time+1)+35*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.888888888888889*pow(hstep, 2)*pow(Xval2[3], 2))/(2*Time+1);

   values[15*(Time+1)+48*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.222222222222222*pow(hstep, 2)*Xval[3]*Xval2[3])/(2*Time+1);

   values[15*(Time+1)+49*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[3]*Xvalp1[3])/(2*Time+1);

   values[15*(Time+1)+52*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.222222222222222*pow(hstep, 2)*Xval[3]*Xval2[3])/(2*Time+1);

   values[15*(Time+1)+53*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.222222222222222*pow(hstep, 2)*Xval2[3]*Xvalp1[3])/(2*Time+1);

   values[15*(Time+1)+54*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[3]*(Xval[0] - Xval[2]))/(2*Time+1);

   values[15*(Time+1)+55*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[3]*(Xvalp1[0] - Xvalp1[2]))/(2*Time+1);

   values[15*(Time+1)+56*(Time)+0+jt] += obj_factor*bounds[4][3]*(-1.33333333333333*hstep*Xval2[3]*(-0.166666666666667*hstep + 1))/(2*Time+1);

   values[15*(Time+1)+57*(Time)+0+jt] += obj_factor*bounds[4][3]*(-1.33333333333333*hstep*Xval2[3]*(-0.166666666666667*hstep - 1))/(2*Time+1);

   values[15*(Time+1)+58*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.888888888888889*pow(hstep, 2)*pow(Xval2[3], 2))/(2*Time+1);

   values[15*(Time+1)+60*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.888888888888889*pow(hstep, 2)*pow(Xval2[3], 2))/(2*Time+1);

   values[15*(Time+1)+61*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.0555555555555556*pow(hstep, 2)*Xval[3]*(4*Xval2[0] - 4*Xval2[2]))/(2*Time+1);

   values[15*(Time+1)+62*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.0555555555555556*pow(hstep, 2)*Xvalp1[3]*(4*Xval2[0] - 4*Xval2[2]))/(2*Time+1);

   values[15*(Time+1)+65*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.0555555555555556*pow(hstep, 2)*Xval[3]*(4*Xval2[0] - 4*Xval2[2]))/(2*Time+1);

   values[15*(Time+1)+66*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.0555555555555556*pow(hstep, 2)*Xvalp1[3]*(4*Xval2[0] - 4*Xval2[2]))/(2*Time+1);

   values[15*(Time+1)+67*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.0555555555555556*pow(hstep, 2)*(Xval[0] - Xval[2])*(4*Xval2[0] - 4*Xval2[2]))/(2*Time+1);

   values[15*(Time+1)+68*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.0555555555555556*pow(hstep, 2)*(4*Xval2[0] - 4*Xval2[2])*(Xvalp1[0] - Xvalp1[2]))/(2*Time+1);

   values[15*(Time+1)+69*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.333333333333333*hstep*(-0.166666666666667*hstep + 1)*(4*Xval2[0] - 4*Xval2[2]))/(2*Time+1);

   values[15*(Time+1)+70*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.333333333333333*hstep*(-0.166666666666667*hstep - 1)*(4*Xval2[0] - 4*Xval2[2]))/(2*Time+1);

   values[15*(Time+1)+71*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.222222222222222*pow(hstep, 2)*Xval2[3]*(4*Xval2[0] - 4*Xval2[2]) + 1.33333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[3]*(Xval[0] - Xval[2]) + 4*Xval2[3]*(Xval2[0] - Xval2[2]) + Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] - 4*Xval2[4] - Xvalp1[4]) + Xval[4] - Xvalp1[4]))/(2*Time+1);

   values[15*(Time+1)+73*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.222222222222222*pow(hstep, 2)*Xval2[3]*(4*Xval2[0] - 4*Xval2[2]) - 1.33333333333333*hstep*(0.166666666666667*hstep*(6*Pval[0] + Xval[3]*(Xval[0] - Xval[2]) + 4*Xval2[3]*(Xval2[0] - Xval2[2]) + Xvalp1[3]*(Xvalp1[0] - Xvalp1[2]) - Xval[4] - 4*Xval2[4] - Xvalp1[4]) + Xval[4] - Xvalp1[4]))/(2*Time+1);

   values[15*(Time+1)+74*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.0555555555555556*pow(hstep, 2)*pow(4*Xval2[0] - 4*Xval2[2], 2))/(2*Time+1);

   values[15*(Time+1)+75*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.222222222222222*pow(hstep, 2)*Xval[3] - 0.25*hstep*Xval[3])/(2*Time+1);

   values[15*(Time+1)+76*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.222222222222222*pow(hstep, 2)*Xvalp1[3] + 0.25*hstep*Xvalp1[3])/(2*Time+1);

   values[15*(Time+1)+79*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.222222222222222*pow(hstep, 2)*Xval[3] + 0.25*hstep*Xval[3])/(2*Time+1);

   values[15*(Time+1)+80*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.222222222222222*pow(hstep, 2)*Xvalp1[3] - 0.25*hstep*Xvalp1[3])/(2*Time+1);

   values[15*(Time+1)+81*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.222222222222222*pow(hstep, 2)*(Xval[0] - Xval[2]) - 0.25*hstep*(Xval[0] - Xval[2]))/(2*Time+1);

   values[15*(Time+1)+82*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.222222222222222*pow(hstep, 2)*(Xvalp1[0] - Xvalp1[2]) - 0.25*hstep*(-Xvalp1[0] + Xvalp1[2]))/(2*Time+1);

   values[15*(Time+1)+83*(Time)+0+jt] += obj_factor*bounds[4][3]*(-1.33333333333333*hstep*(-0.166666666666667*hstep + 1) + 0.25*hstep - 1.0)/(2*Time+1);

   values[15*(Time+1)+84*(Time)+0+jt] += obj_factor*bounds[4][3]*(-1.33333333333333*hstep*(-0.166666666666667*hstep - 1) - 0.25*hstep - 1.0)/(2*Time+1);

   values[15*(Time+1)+85*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.888888888888889*pow(hstep, 2)*Xval2[3])/(2*Time+1);

   values[15*(Time+1)+87*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.888888888888889*pow(hstep, 2)*Xval2[3])/(2*Time+1);

   values[15*(Time+1)+88*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.222222222222222*pow(hstep, 2)*(4*Xval2[0] - 4*Xval2[2]))/(2*Time+1);

   values[15*(Time+1)+89*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.888888888888889*pow(hstep, 2) + 2)/(2*Time+1);

   values[15*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.333333333333333*pow(hstep, 2)*Xval[3])/(2*Time+1);

   values[15*(Time+1)+90*(Time)+0+ 1+jt] += obj_factor*bounds[4][3]*(0.333333333333333*pow(hstep, 2)*Xvalp1[3])/(2*Time+1);

   values[17*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[4][3]*(-0.333333333333333*pow(hstep, 2)*Xval[3])/(2*Time+1);

   values[17*(Time+1)+90*(Time)+0+ 1+jt] += obj_factor*bounds[4][3]*(-0.333333333333333*pow(hstep, 2)*Xvalp1[3])/(2*Time+1);

   values[18*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.333333333333333*pow(hstep, 2)*(Xval[0] - Xval[2]))/(2*Time+1);

   values[18*(Time+1)+90*(Time)+0+ 1+jt] += obj_factor*bounds[4][3]*(0.333333333333333*pow(hstep, 2)*(Xvalp1[0] - Xvalp1[2]))/(2*Time+1);

   values[19*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[4][3]*(2.0*hstep*(-0.166666666666667*hstep + 1))/(2*Time+1);

   values[19*(Time+1)+90*(Time)+0+ 1+jt] += obj_factor*bounds[4][3]*(2.0*hstep*(-0.166666666666667*hstep - 1))/(2*Time+1);

   values[20*(Time+1)+90*(Time)+0+jt] += obj_factor*bounds[4][3]*(1.33333333333333*pow(hstep, 2)*Xval2[3])/(2*Time+1);

   values[20*(Time+1)+92*(Time)+0+jt] += obj_factor*bounds[4][3]*(-1.33333333333333*pow(hstep, 2)*Xval2[3])/(2*Time+1);

   values[20*(Time+1)+93*(Time)+0+jt] += obj_factor*bounds[4][3]*(0.333333333333333*pow(hstep, 2)*(4*Xval2[0] - 4*Xval2[2]))/(2*Time+1);

   values[20*(Time+1)+94*(Time)+0+jt] += obj_factor*bounds[4][3]*(-1.33333333333333*pow(hstep, 2))/(2*Time+1);

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

     Xdval[0] = data1[2*Time];
     Xdval2[0] = 0;
     Xdvalp1[0] = 0;

     for(Index i=0;i<nP;i++) {
        Pval[i] = x[(2*Time+1)*(nY+2*nU)+i];
     } //end for loop


   values[0*(Time+1)+0*(Time)+0+ Time] += obj_factor*(8)/(2*Time+1);

  } // end else 

   return true;
}


void LORENZ96_NLP::finalize_solution(SolverReturn status,
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
