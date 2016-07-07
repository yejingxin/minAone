#!/usr/bin/python

#####################################################
#
#  10 January 2011
#  Orginally written by Bryan A. Toth
#  University of California, San Diego
#  btoth@physics.ucsd.edu
# 
#
#  30 September 2014
#  Modified by Jingxin Ye 
#  University of California, San Diego
#  j9ye@physics.ucsd.edu
#
#
#  5 June 2016
#  Modified by Nirag Kadakia 
#  University of California, San Diego
#  nkadakia@physics.ucsd.edu
#  Added functionality for different types of observation and 
#    stimulus files so many jobs can be run in parallel with 
#    different data sets. Also added difrent output file 
#    formats. See user guide updates for further information.
#  
#####################################################
#
#  This script builds C++ code to minimize a path integral
#  formulation of a parameter estimation problem with
#  the optimization software IPOPT.
#
#  Specifically, given a vector-field (model) of the form:
#
#           dx_1(t) = G_1(x_1(t),x_p(t),q)
#           dx_p(t) = G_p(x_1(t),x_p(t),q)
#
#          where x_p denotes 1 or more equations in the model,
#
#  this code takes the discAzerod vector field and objective 
#  function (discAzerod in companion script discAzero.py),
#  and builds the requisite IPOPT functions to solve the 
#  resulting optimization problem.
#
#  This script has been developed as part of a suite of 
#  python scripts to define a dynamic parameter estimation
#  problem using the optimization software IPOPT, but could
#  easily be modified for use with other optimization software.
#
######################################################


#  For ease of use, all necessary C++ files are built with one keyboard
#   command.  The following four scripts each write a necessary C++ file.
#   See the individual scripts for more information.


import discAone
import makecppAone
import makehppAone
import makemakeAone
import makeoptAone
# DiscAone.py reads the file equations.txt and sets up the given vector field
#  in the correct format.
print("Writing the problem files\n")
# Import the problem name and change to upper and lower case
prob = discAone.Problem
probu = prob.upper()
probl = prob.lower()

FILE = probl + 'minAone_nlp.cpp'

nU = discAone.nU
nP = discAone.nP
nY = discAone.nY
nI = discAone.nI
nF = discAone.nF
nM = discAone.nM
# The name of the IPOPT file to be written to
f = open(FILE,'w')

# The following write commands are writing C++ code.

# Front matter
f.write('// %s.cpp\n' % probl)
f.write('// Nonlinear Ipopt program\n')
f.write('\n// Author: Bryan A. Toth\n// btoth@physics.ucsd.edu\n\n')
f.write('#include \"%sminAone_nlp.hpp\"\n' % probl)
f.write('#include <cmath>\n\
#include <cstdio>\n\
#include <iostream>\n\
#include <fstream>\n\
#include <string>\n\
#include <stdlib.h>\n\
#include <cstring>\n')

f.write('#ifndef HAVE_CSTDIO\n\
#define HAVE_CSTDIO\n\
# include <cstdio>\n\
# include <iostream>\n')
if nF>0:
    f.write('# include "myfunctions.hpp"\n')
f.write('# include <fstream>\n\
# include <string>\n\
# include <stdlib.h>\n\
#else\n\
# ifndef HAVE_STDIO_H\n\
#define HAVE_STDIO_H\n\
#  include <stdio.h>\n\
# else\n\
#  error \"don\'t have header file for stdio\"\n\
# endif\n\
#endif\n\n')

f.write('using namespace Ipopt;\n\n')
f.write('using namespace std;\n\n')
for i in range(nF):
  args = discAone.Funcarg[i]
  f.write('double %s(' % discAone.Sf[i])
  for j in range(args-1):
     f.write('double, ')
  f.write('double);\n')
  f.write('double %sjac(' % discAone.Sf[i])
  for j in range(args-1):
     f.write('double, ')
  f.write('double, int);\n')
  f.write('double %shes(' % discAone.Sf[i])
  for j in range(args-1):
     f.write('double, ')
  f.write('double, int, int);\n')
f.write('// constructor\n\
%s_NLP::%s_NLP(int id)\n\
{\n' % (probu, probu))

# Define problem parameters
f.write('  nU=%d;\n' % nU)
f.write('  nP=%d;\n' % nP)
f.write('  nY=%d;\n' % nY)
f.write('  nM=%d;\n' % nM)
f.write('  nI=%d;\n\n' % nI)

# Define variable names 
f.write('\
  K11val = new double[nU];\n\
  K11val2 = new double[nU];\n\
  K11valp1 = new double[nU];\n\
  dK11val = new double[nU];\n\
  dK11val2 = new double[nU];\n\
  dK11valp1 = new double[nU];\n\
  Xdval = new double[nM];\n\
  Xdval2 = new double[nM];\n\
  Xdvalp1 = new double[nM];\n')
f.write('\
  Xval = new double[nY];\n\
  Xval2 = new double[nY];\n\
  Xvalp1 = new double[nY];\n')
f.write('\
  Pval = new double[nP];\n')
f.write('\
  Ival = new double[nI];\n\
  Ival2 = new double[nI];\n\
  Ivalp1 = new double[nI];\n\
  Rf0 = new double[nY];\n\n')

# Commands to read specs.txt file
# The string length for specs is an overshoot; different for different file_fmt's
f.write('\
  string buffer;\n\
  specs = new string[9+nP+nY+nI+nU+nM+1+1];\n\n\
  int count;\n\
  count = 0;\n\
  \n\
  ifstream fin ("specs.txt");\n')
f.write("\
  if (fin.is_open()){\n\
    while (! fin.eof()){\n\
      getline (fin,buffer);\n\
      if (buffer[0] !='#'){\n\
        specs[count] = buffer;\n\
        count++;\n\
      }\n\
    }\n\
    fin.close();\n\
  }\n\
  else cout << \"Unable to open file\";\n\n")

# Write Time to file
#  Time is a misnomer - this is a measure of the the number
#  of time steps used in the problem
f.write('  Time = atoi(specs[0].c_str());\n')
# Write skip to file
#  Skip is a dummy variable to allow the use of various parts
#  of a given data file
f.write('  skip = atoi(specs[1].c_str());\n')
# Write hstep to file
#  Hstep is the time-step of the discretization
f.write('  hstep = atof(specs[2].c_str());\n\n')


# Declare variables for input files and states

f.write('\
  string filename;\n\
  int ret;\n\
  double tempdata;\n\
  Ntotal = (2*Time+1)*nY +nP;\n\
  solution = new double[Ntotal];\n\
  int toggle = 0;\n\
  file_fmt = atoi(specs[3].c_str());')


# Read in Observation files; depending on file structure
# The type of input is tagged in the 4th line of specs
# 0: No initial data file; stimulus and observations in separate files
# 1: Initial data file; stimulus and observations in separate files
# 2: No initial data file; stimulus and observations in single files; 
#    observation files must have dummy colums for unmeasured variables too;
#    stimulus files must have nI columns. This is to allow for ease of including
#    more measurements without having to write separate measurement files.
# 3: Initial data file; stimulus and observations in single files; 
#    observation files must have dummy colums for unmeasured variables too;
#    stimulus files must have nI columns

if nM != 0:
  f.write('\n\n  // Read in observation files\n\n')
  f.write('\
    if (file_fmt == 1) toggle = 1;\n\
    if (file_fmt == 2) toggle = 2;\n\
    if (file_fmt == 3) toggle = 3;\n\n')

  # Read multiple files if tag = 0 or 1
  f.write('\
    if (file_fmt < 2) {\n\
      // Read in files one by one if file_fmt = 0,1\n\n\
      taskid = id;\n\n')
  for i in range(nM):
    f.write('\
      %s = new double[2*Time+1];\n\
      %s = new double[skip];\n\
      FILE *pFile%d;\n\
      filename = specs[%d+toggle];\n\
      pFile%d = fopen(filename.c_str(),"r");\n' % (discAone.Ldata[i],discAone.Ldata[i]+'dummy',i,4+i,i))
    temp1 = "%lf"
    f.write('\
      for(Index jt=0;jt<skip;jt++)\n\
      {\n\
        ret = fscanf (pFile%d, "%s", &%s[jt]);\n\
        if (ret == EOF) break;\n\
      }\n\
      for(Index jt=0;jt<2*Time+1;jt++)\n\
      {\n\
        ret = fscanf (pFile%d, "%s", &%s[jt]);\n\
        if (ret == EOF) break;\n\
      }\n\
      fclose (pFile%d);\n\n' % (i,temp1, discAone.Ldata[i]+'dummy',i,temp1,discAone.Ldata[i],i))
  f.write('\
    }\n')

  # Read single indexed file if tag is 2 or 3
  f.write('\
    else if (file_fmt >= 2){\n\
      // Read in files from single file if file_fmt = 2,3\n\n') 

  temp1 = "%d"
  temp2 = "%"
  temp3 = "%s"
  for i in range(nM):
    f.write('\
      %s = new double[2*Time+1];\n\
      %s = new double[skip];\n' % (discAone.Ldata[i],discAone.Ldata[i]+'dummy'))
  f.write('\n')
  f.write('\
      FILE *pFile0;\n\
      filename = specs[4+toggle];\n\
      char init_idx[21]; \n\
      modid = atoi(specs[4].c_str());\n\
      if (modid != 0){\n\
        pathid = id / modid;\n\
        taskid = id %s modid;\n\
        sprintf(init_idx, "%s.%s", pathid, specs[5].c_str());\n\
      }\n\
      else{\n\
        pathid = 0;\n\
        taskid = id;\n\
        sprintf(init_idx, ".%s", specs[5].c_str());\n\
      }\n\
      filename += init_idx;\n\
      pFile0 = fopen(filename.c_str(),"r");\n\n' % (temp2, temp1, temp3, temp3))
  f.write('\
      for(Index rows=0;rows<skip;rows++){\n\
        for (Index cols=0;cols<nY;cols++){\n')
  for i in range(nM):
    temp1 = "%lf"
    if i == 0:
      f.write('\
          if (cols == %d){\n\
            ret = fscanf (pFile0, "%s", &%sdummy[rows]);\n\
          }\n' % (i,temp1,discAone.Ldata[i]))
    else: 
      f.write('\
          else if (cols == %d){\n\
            ret = fscanf (pFile0, "%s", &%sdummy[rows]);\n\
          }\n' % (i,temp1,discAone.Ldata[i]))
  f.write('\
          else {\n\
            ret = fscanf (pFile0, "%s", &tempdata);\n\
          }          \n\
           if (ret == EOF) break;\n\
        }\n\
      }\n' % temp1)

  f.write('\
      for(Index rows=0;rows<2*Time+1;rows++){\n\
        for (Index cols=0;cols<nY;cols++){\n')
  for i in range(nM):
    temp1 = "%lf"
    if i == 0:  
      f.write('\
          if (cols == %d){\n\
            ret = fscanf (pFile0, "%s", &%s[rows]);\n\
          }\n' % (i,temp1,discAone.Ldata[i]))
    else:
      f.write('\
          else if (cols == %d){\n\
            ret = fscanf (pFile0, "%s", &%s[rows]);\n\
          }\n' % (i,temp1,discAone.Ldata[i]))
  f.write('\
          else {\n\
            ret = fscanf (pFile0, "%s", &tempdata);\n\
          }          \n\
           if (ret == EOF) break;\n\
        }\n\
      }\n' % temp1)
  f.write('\
      fclose (pFile0);\n\
    }\n')


# Read in Stimulus files; depending on file structure
if nI != 0:
  f.write('\n  // Read in stimulus files\n\n')

  # Read multiple files if tag = 0 or 1
  f.write('\
    if (file_fmt < 2) {\n\
      // Read in files one by one if file_fmt = 0,1\n\n') 
  for i in range(nI):
    f.write('\
      %s = new double[2*Time+1];\n\
      %s = new double[skip];\n\
      FILE *qFile%d;\n\
      filename = specs[%d+toggle];\n\
      qFile%d = fopen(filename.c_str(),"r");\n' % (discAone.Lstimuli[i],discAone.Lstimuli[i]+'dummy',i,4+nM+i,i))
    temp1 = "%lf"
    f.write('\
      for(Index jt=0;jt<skip;jt++)\n\
      {\n\
        ret = fscanf (qFile%d, "%s", &%s[jt]);\n\
        if (ret == EOF) break;\n\
      }\n\
      for(Index jt=0;jt<2*Time+1;jt++)\n\
      {\n\
        ret = fscanf (qFile%d, "%s", &%s[jt]);\n\
        if (ret == EOF) break;\n\
      }\n\
      fclose (qFile%d);\n\n' % (i,temp1, discAone.Lstimuli[i]+'dummy',i,temp1,discAone.Lstimuli[i],i))
  f.write('\
    }\n\n')

  # Read single indexed file if tag is 2 or 3
  f.write('\
    else if (file_fmt >= 2){\n\
      // Read in data from a single file if file_fmt = 2,3\n\n') 
  temp1 = "%d"
  temp2 = "%"
  temp3 = "%s"

  for i in range(nI):
    f.write('\
      %s = new double[2*Time+1];\n\
      %s = new double[skip];\n' % (discAone.Lstimuli[i],discAone.Lstimuli[i]+'dummy'))
  f.write('\n')
  f.write('\
      FILE *qFile0;\n\
      filename = specs[5+toggle];\n\
      char init_idx[21]; \n\
      modid = atoi(specs[4].c_str());\n\
      if (modid != 0){\n\
        pathid = id / modid;\n\
        taskid = id %s modid;\n\
        sprintf(init_idx, "%s.%s", pathid, specs[5].c_str());\n\
      }\n\
      else{\n\
        pathid = 0;\n\
        taskid = id;\n\
        sprintf(init_idx, ".%s", specs[5].c_str());\n\
      }\n\
      filename += init_idx;\n\
      qFile0 = fopen(filename.c_str(),"r");\n\n' % (temp2, temp1, temp3, temp3))
  f.write('\
      for(Index rows=0;rows<skip;rows++){\n\
        for (Index cols=0;cols<nI;cols++){\n')
  for i in range(nI):
    temp1 = "%lf"
    if i == 0:
      f.write('\
          if (cols == %d){\n\
            ret = fscanf (qFile0, "%s", &%sdummy[rows]);\n\
          }\n' % (i,temp1,discAone.Lstimuli[i]))
    else: 
      f.write('\
          else if (cols == %d){\n\
            ret = fscanf (qFile0, "%s", &%sdummy[rows]);\n\
          }\n' % (i,temp1,discAone.Lstimuli[i]))
  f.write('\
        }\n\
      }\n')

  f.write('\
      for(Index rows=0;rows<2*Time+1;rows++){\n\
        for (Index cols=0;cols<nI;cols++){\n')
  for i in range(nI):
    temp1 = "%lf"
    if i == 0:  
      f.write('\
          if (cols == %d){\n\
            ret = fscanf (qFile0, "%s", &%s[rows]);\n\
          }\n' % (i,temp1,discAone.Lstimuli[i]))
    else:
      f.write('\
          else if (cols == %d){\n\
            ret = fscanf (qFile0, "%s", &%s[rows]);\n\
          }\n' % (i,temp1,discAone.Lstimuli[i]))
  f.write('\
        }\n\
      }\n')

  f.write('\
      fclose (qFile0);\n\
    }\n')


# Read in the initial and boundary conditions for all variables into arrays
f.write('\n  //Read in boundary conditions\n\n')

f.write('\
  if (file_fmt == 0) toggle = nM + nI;\n\
  if (file_fmt == 1) toggle = nM + nI + 1;\n\
  if (file_fmt == 2) toggle = (nM > 0) + (nI > 0) + 2;\n\
  if (file_fmt == 3) toggle = (nM > 0) + (nI > 0) + 3;\n\n\
  int rows = nY+nU+nP+1;\n\
  bounds = new double*[rows];\n\
  for (Index i=0;i<rows;i++) bounds[i] = new double[4];\n\
  int counter;\n\
  for(Index k=0;k<rows;k++){\n\
    counter=0;\n\
    char* tmp = new char[specs[4+toggle+k].size()+1];\n\
    strcpy( tmp, specs[4+toggle+k].c_str() );\n\
    char *ptr = strtok(tmp,",");\n\
    bounds[k][3] = 0.0;\n\
    while(ptr != 0){\n\
      if(counter<3){\n\
        bounds[k][counter] = atof(ptr);\n\
      }\n\
      if(counter==3) {\n\
        bounds[k][counter] = atof(ptr);\n\
      }\n\
      ptr = strtok(0,",");\n\
      counter++;\n\
    }\n\
  }\n\
  for (Index i=0;i<nY;i++){\n\
    Rf0[i]=bounds[i][2];\n\
    bounds[i][3]=bounds[i][2];\n\
  }\n\n\
  beta=0;\n\
  alpha = bounds[nY+nU+nP][0];\n\
  delta_beta=(int) bounds[nY+nU+nP][1];\n\
  max_beta=(int) bounds[nY+nU+nP][2];\n\n\
  //Read in output_fmt (zero is default)\n\
  output_fmt = atoi(specs[4+toggle+rows].c_str());\n')

f.write('\
}\n\n')

### Write the destructor
f.write('// destructor\n\
%s_NLP::~%s_NLP()\n\
{\n\
  delete [] K11val;\n\
  delete [] K11val2;\n\
  delete [] K11valp1;\n\
  delete [] dK11val;\n\
  delete [] dK11val2;\n\
  delete [] dK11valp1;\n\
  delete [] Xdval;\n\
  delete [] Xdval2;\n\
  delete [] Xdvalp1;\n\
  delete [] Xval;\n\
  delete [] Xval2;\n\
  delete [] Xvalp1;\n\
  delete [] Pval;\n\
  delete [] Ival;\n\
  delete [] Ival2;\n\
  delete [] Ivalp1;\n\
  delete [] specs;\n' % (probu, probu))
for i in range(nM):
  f.write('\
  delete [] %s;\n\
  delete [] %s;\n' % (discAone.Ldata[i],discAone.Ldata[i]+'dummy'))
for i in range(nI):
  f.write('\
  delete [] %s;\n\
  delete [] %s;\n' % (discAone.Lstimuli[i],discAone.Lstimuli[i]+'dummy'))
f.write('\
  int rows = nY+nU+nP;\n\
  for (Index i=0;i<rows;i++) delete [] bounds[i];\n\
  delete [] bounds;\n\
}\n')


### Below are individual functions needed including ipopt functions

# change RF
f.write('\n\
bool %s_NLP::changeRf(){\n\
  if((beta+delta_beta)>(max_beta-1))\n\
    return false;\n\
  else\n\
    beta = beta + delta_beta;\n\
  for (Index i=0;i<nY;i++) {\n\
    bounds[i][3]=pow(alpha,beta)*Rf0[i];\n\
  }\n\
  printf("\\n\\n\\n\\n<------       Beta=%%d       ------->\\n\\n\\n\\n",beta);\n\
  return true;\n\
}\n\n'% probu)


# GET_NLP_INFO 
f.write('// returns the size of the problem\n\
bool %s_NLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,\n\
        Index& nnz_h_lag, IndexStyleEnum& index_style)\n\n' % probu)

# Number of variables 
alpha = 2*(nY+nU)
beta = nY+nU+nP
f.write('{\n\
  // Number of variables\n\
  n = %d*Time+%d;\n' % (alpha,beta))

# Number of equality constraints


f.write('\n\
  // Number of equality constraints\n\
  m = 0;\n')

# Number of Jacobian nonone entries
f.write('\n\
  // Number of Jacobian nonzero entries\n\
  nnz_jac_g = 0;\n')

# Number of Hessian non-zeros in lower left of diagonal
# Determine number of large, medium, small elements in Hessian
Hessize = discAone.Hessize
Hesindex = discAone.Hesindex
cat = len(discAone.Smod)
sma = 0
med = 0
lar = 0
for i in range(cat):
  for j in range(i+1):
      if Hesindex[i][j] == -1:
         temp = bytes.decode(Hessize[i][j])
         if temp == 's':
            sma = sma + 1
         elif temp == 'm':
            med = med + 1
         elif temp == 'l':
            lar = lar + 1 

f.write('\n\
  // Number of Hessian nonzero entries\n\
  nnz_h_lag = %d*(Time+1)+%d*Time+%d;\n' % (lar, med, sma))

f.write('\n\
  // use the C style indexing (0-based)\n\
  index_style = TNLP::C_STYLE;\n\n\
  return true;\n\
}\n\n\n')



# GET_BOUNDS_INFO
f.write('// returns the variable bounds\n\
bool %s_NLP::get_bounds_info(Index n, Number* x_l, Number* x_u,\n\
                                Index m, Number* g_l, Number* g_u)\n\
{\n\
  // Here, the n and m we gave IPOPT in get_nlp_info are passed back to us.\n' % probu)
f.write('  assert(n == %d*Time+%d);\n' % (alpha, beta))
f.write('  assert(m == 0);\n\n')

# Set bounds on state variables
f.write('\
  for(Index jt=0;jt<Time+1;jt++) {\n\
    for(Index var=0;var<nY;var++) {\n\
      // Bounds for x\n\
      x_l[(Time+1)*var+jt]=bounds[var][0];\n\
      x_u[(Time+1)*var+jt]=bounds[var][1];\n\
      // Bounds for midpoints\n\
      if(jt<Time){\n\
        x_l[(Time+1)*(nY+nU)+Time*var+jt]=bounds[var][0];\n\
        x_u[(Time+1)*(nY+nU)+Time*var+jt]=bounds[var][1];\n\
      }\n\
    }\n')

# Set bounds on controls 
f.write('\
    for(Index cup=0;cup<nU;cup++) {\n\
      // Bounds for k\n\
      x_l[(Time+1)*(nY+cup)+jt]=bounds[nY+cup][0];\n\
      x_u[(Time+1)*(nY+cup)+jt]=bounds[nY+cup][1];\n\
      // Bounds for midpoints\n\
      if(jt<Time) {\n\
        x_l[(Time+1)*(nY+nU)+Time*(nY+cup)+jt]=bounds[nY+cup][0];\n\
        x_u[(Time+1)*(nY+nU)+Time*(nY+cup)+jt]=bounds[nY+cup][1];\n\
      }\n\
    }\n\
  } \n\n')

# Set bounds on parameters
f.write('\
  for(Index par=0;par<nP;par++) {\n\
    // Bounds for parameters\n\
    x_l[2*Time*(nY+nU)+nY+nU+par]=bounds[nY+nU+par][0];\n\
    x_u[2*Time*(nY+nU)+nY+nU+par]=bounds[nY+nU+par][1];\n\
  }\n\n\
  return true;\n\
}\n')


# GET_STARTING_POINT
f.write('// returns the initial point for the problem\n\
bool %s_NLP::get_starting_point(Index n, bool init_x, Number* x,\n\
                                   bool init_z, Number* z_L, Number* z_U,\n\
                                   Index m, bool init_lambda,\n\
                                   Number* lambda)\n\
{\n\
  assert(init_x == true);\n\
  assert(init_z == false);\n\
  assert(init_lambda == false);\n\n\
  for (Index i=0; i<n; i++) {\n\
    x[i] = 0.0;\n\
  }\n\n' % probu)

temp1 = "%lf"
temp2 = "%d"
temp3 = "%s"

f.write('\
  if(beta==0){\n\
    double **skipinit = new double* [skip];\n\
    double **init = new double* [2*Time + 1];\n\
    double *param_init = new double [nP];\n\n\
    for(Index i=0;i<skip;i++) skipinit[i] = new double[nY];\n\
    for(Index i=0;i<2*Time+1;i++) init[i] = new double[nY];\n\n\
    string filename;\n\
    int ret;\n\n\
    // Read in initial data files for file_fmt = 1,3\n\n\
    if (file_fmt == 1){\n\
      filename = specs[4];\n\
      FILE *initFILE;\n\
      initFILE = fopen(filename.c_str(),"r");\n\
      for(Index jt=0;jt<skip;jt++){\n\
        for (Index jy=0;jy<nY;jy++){\n\
          ret = fscanf (initFILE,"%s",&skipinit[jt][jy]);\n\
          if (ret == EOF) break;\n\
        }\n\
      }\n\
      for(Index jt=0;jt<(2*Time+1);jt++){\n\
        for (Index jy=0;jy<nY;jy++){\n\
          ret = fscanf (initFILE,"%s",&init[jt][jy]);\n\
          if (ret == EOF) break;\n\
        }\n\
      }\n\
      for(Index i=0;i<nP;i++){\n\
        ret = fscanf (initFILE,"%s",&param_init[i]);\n\
        if (ret == EOF) break;\n\
      }\n\
      fclose (initFILE);\n\
    }\n' % (temp1, temp1, temp1))

f.write('\
    else if (file_fmt == 3){\n\
      filename = specs[6];\n\
      FILE *initFILE;\n\
      char init_idx[21];\n\
      if (modid != 0){\n\
        sprintf(init_idx, "%s.%s", taskid, specs[5].c_str());\n\
      }\n\
      else{\n\
        sprintf(init_idx, ".%s", specs[5].c_str());\n\
      }\n\
      filename += init_idx;\n\
      initFILE = fopen(filename.c_str(),"r");\n\n' % (temp2, temp3, temp3))

f.write('\
      for(Index jt=0;jt<skip;jt++){\n\
        for (Index jy=0;jy<nY;jy++){\n\
          ret = fscanf (initFILE,"%s",&skipinit[jt][jy]);\n\
          if (ret == EOF) break;\n\
        }\n\
      }\n\
      for(Index jt=0;jt<(2*Time+1);jt++){\n\
        for (Index jy=0;jy<nY;jy++){\n\
          ret = fscanf (initFILE,"%s",&init[jt][jy]);\n\
          if (ret == EOF) break;\n\
        }\n\
      }\n\
      for(Index i=0;i<nP;i++){\n\
        ret = fscanf (initFILE,"%s",&param_init[i]);\n\
        if (ret == EOF) break;\n\
      }\n\
      fclose (initFILE);\n\
    }\n\n' % (temp1, temp1, temp1))

f.write('\
    // Save init_data to variables for file_fmt = 0,1,2,3\n\n\
    if ((file_fmt == 1) || (file_fmt == 3)){\n\
      for (Index jt=0;jt<Time+1;jt++){\n\
        for (Index var=0;var<nY;var++){\n\n\
        // Initial conditions for x and midpoints\n\
          x[(Time+1)*var+jt] = init[2*jt][var];\n\
          if (jt<Time) x[(Time+1)*(nY+nU)+Time*var+jt] = init[2*jt+1][var];\n\
        }\n\n\
        // Initial conditions for controls and midpoints\n\
        for (Index cup=0;cup<nU;cup++){\n\
          x[(Time+1)*(cup+nY)+jt]=bounds[cup+nY][2];\n\
          if (jt<Time) x[(Time+1)*(nY+nU)+Time*(cup+nY)+jt]=bounds[cup+nY][2];\n\
        }\n\
      }\n\n\
      // Initial conditions for parameters\n\
      for (Index par=0;par<nP;par++){\n\
        x[(2*Time+1)*(nY+nU)+par] = param_init[par];\n\
      }\n\
    }\n\
    else if ((file_fmt == 0) || (file_fmt == 2)){\n\
      srand(taskid);\n\
      for (Index jt=0;jt<Time+1;jt++){\n\
        for(Index var=0;var<nY;var++){\n\
          x[(Time+1)*var+jt] = rand()*1.0/RAND_MAX*(bounds[var][1]-bounds[var][0])+bounds[var][0];\n\
          if (jt<Time) x[(Time+1)*(nY+nU)+Time*var+jt] = rand()*1.0/RAND_MAX*(bounds[var][1]-bounds[var][0])+bounds[var][0];\n\
        }\n\
        for (Index cup=0;cup<nU;cup++){\n\
          x[(Time+1)*(cup+nY)+jt]=bounds[cup+nY][2];\n\
          if (jt<Time) x[(Time+1)*(nY+nU)+Time*(cup+nY)+jt]=bounds[cup+nY][2];\n\
        }\n\
        for (Index par=0;par<nP;par++){\n\
          x[(2*Time+1)*(nY+nU)+par] = rand()*1.0/RAND_MAX*(bounds[nY+nU+par][1]-bounds[nY+nU+par][0])+bounds[nY+nU+par][0];\n\
        }\n\
      }\n\
    }\n\
    for(Index i=0;i<2*Time+1;i++) delete [] init[i];\n\
    delete [] init;\n\
  }\n\
  else{\n\
    for(Index i=0;i<Ntotal-nP;i++){\n\
      x[i] = solution[i];\n\
    }\n\
    for (Index jt=0;jt<Time+1;jt++){\n\
      for (Index cup=0;cup<nU;cup++){\n\
        x[(Time+1)*(cup+nY)+jt]=bounds[cup+nY][2];\n\
        if (jt<Time) x[(Time+1)*(nY+nU)+Time*(cup+nY)+jt]=bounds[cup+nY][2];\n\
      }\n\
    }\n\
    for(Index i=0;i<nP;i++){\n\
      x[(2*Time+1)*(nU+nY)+i] = solution[Ntotal-nP+i];\n\
    }\n\
  }\n\
  return true;\n\
}\n\n')


# EVAL_F
# Subroutine to calculate the objective value
# Here, and in the following subroutines, strings from the result
#  of symbolic discretization and differentiation in discAzero.py

#  are inserted to the code.

f.write('// returns the value of the objective function\n\
bool %s_NLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)\n\
{\n' % probu)
f.write('  assert(n == %d*Time+%d);\n' % (alpha, beta))
f.write('  obj_value = 0;\n\n')
f.write('  for(Index jt=0;jt<Time;jt++) {\n\n')
f.write('     for(Index i=0;i<nY;i++) {\n')
f.write('        Xval[i] = x[i*(Time+1) + jt];\n')
f.write('        Xvalp1[i] = x[i*(Time+1) + jt + 1];\n')
f.write('        Xval2[i] = x[(Time+1)*(nY+nU) + i*(Time) + jt];\n')
f.write('     } //end for loop\n\n')

f.write('\n')
f.write('     for(Index cup=0;cup<nU;cup++) {\n')
f.write('        K11val[cup] = x[nY*(Time+1) + cup*(Time+1) + jt];\n')
f.write('        K11valp1[cup] = x[nY*(Time+1) + cup*(Time+1) + jt + 1];\n')
f.write('        K11val2[cup] = x[(Time+1)*(nY+nU) + Time*nY + cup*(Time) + jt];\n')

f.write('     } //end for loop\n\n')

for i in range(nM):
    f.write('     Xdval[%d] = %s[2*jt];\n' % (i,discAone.Ldata[i]))
    f.write('     Xdval2[%d] = %s[2*jt+1];\n' % (i,discAone.Ldata[i]))
    f.write('     Xdvalp1[%d] = %s[2*jt+2];\n' % (i,discAone.Ldata[i]))

for i in range(nI):
    f.write('     Ival[%d] = %s[2*jt];\n' % (i,discAone.Lstimuli[i]))
    f.write('     Ival2[%d] = %s[2*jt+1];\n' % (i,discAone.Lstimuli[i]))
    f.write('     Ivalp1[%d] = %s[2*jt+2];\n' % (i,discAone.Lstimuli[i]))

f.write('\n')
f.write('     for(Index i=0;i<nP;i++) {\n')
f.write('        Pval[i] = x[(2*Time+1)*(nY+nU)+i];\n')
f.write('     } //end for loop\n\n')
f.write('\n')
f.write('    obj_value += %s; \n\n' % (discAone.strObj[0]))
for i in range(nY):
   f.write('    obj_value += bounds[%s][3]*(%s); \n\n' % (i,discAone.strObj[1][i]))
  
f.write('  } //end for loop\n\n')

# Add code for last element

f.write('// Add last element\n')
f.write('     for(Index i=0;i<nY;i++) {\n')
f.write('        Xval[i] = x[Time + i*(Time+1)];\n')
f.write('        Xvalp1[i] = 0;\n')
f.write('        Xval2[i] = 0;\n')
f.write('     } //end for loop\n\n')

f.write('\n')
f.write('     for(Index cup=0;cup<nU;cup++) {\n')
f.write('        K11val[cup] = x[nY*(Time+1) + Time + cup*(Time+1)];\n')
f.write('        K11valp1[cup] = 0;\n')
f.write('        K11val2[cup] = 0;\n')

f.write('     } //end for loop\n\n')

for i in range(nM):
    f.write('     Xdval[%d] = %s[2*Time];\n' % (i,discAone.Ldata[i]))
    f.write('     Xdval2[%d] = 0;\n' % i)
    f.write('     Xdvalp1[%d] = 0;\n' % i)

for i in range(nI):
    f.write('     Ival[%d] = %s[2*Time];\n' % (i,discAone.Lstimuli[i]))
    f.write('     Ival2[%d] = 0;\n' % i)
    f.write('     Ivalp1[%d] = 0;\n' % i)

f.write('\n')
f.write('     for(Index i=0;i<nP;i++) {\n')
f.write('        Pval[i] = x[(2*Time+1)*(nY+nU)+i];\n')
f.write('     } //end for loop\n\n')
f.write('\n')

f.write('  obj_value += %s;\n\n' % (discAone.strObj[0]))

# Adding a line to divide the overall objective function by 2T + 1
# This normalizes the objective function

f.write('  obj_value = obj_value/(2*Time+1);\n\n')

f.write('  return true;\n\
}\n\n\n')



# EVAL_GRAD_F

f.write('// return the gradient of the objective function grad_{x} f(x)\n\
bool %s_NLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)\n\
{\n' % probu)
f.write('  assert(n == %d*Time+%d);\n\n' % (alpha, beta))

f.write('  for(Index i=0;i<n;i++) {\n')
f.write('     grad_f[i] = 0;\n')
f.write('  }\n\n')

f.write('  for(Index jt=0;jt<Time;jt++) {\n\n')
f.write('     for(Index i=0;i<nY;i++) {\n')
f.write('        Xval[i] = x[jt + i*(Time+1)];\n')
f.write('        Xvalp1[i] = x[jt + i*(Time+1) + 1];\n')
f.write('        Xval2[i] = x[(Time+1)*(nY+nU) + i*(Time)+ jt];\n')
f.write('     } //end for loop\n\n')

f.write('     for(Index cup=0;cup<nU;cup++) {\n')
f.write('        K11val[cup] = x[nY*(Time+1) + cup*(Time+1) + jt];\n')
f.write('        K11valp1[cup] = x[nY*(Time+1) + cup*(Time+1) + jt + 1];\n')
f.write('        K11val2[cup] = x[(Time+1)*(nY+nU) + (nY)*Time + cup*(Time) +  jt];\n')

f.write('     } //end for loop\n\n')

for i in range(nM):
    f.write('     Xdval[%d] = %s[2*jt];\n' % (i,discAone.Ldata[i]))
    f.write('     Xdval2[%d] = %s[2*jt+1];\n' % (i,discAone.Ldata[i]))
    f.write('     Xdvalp1[%d] = %s[2*jt+2];\n' % (i,discAone.Ldata[i]))

for i in range(nI):
    f.write('     Ival[%d] = %s[2*jt];\n' % (i,discAone.Lstimuli[i]))
    f.write('     Ival2[%d] = %s[2*jt+1];\n' % (i,discAone.Lstimuli[i]))
    f.write('     Ivalp1[%d] = %s[2*jt+2];\n' % (i,discAone.Lstimuli[i]))

f.write('\n')
f.write('     for(Index i=0;i<nP;i++) {\n')
f.write('        Pval[i] = x[(2*Time+1)*(nY+nU)+i];\n')
f.write('     } //end for loop\n')
f.write('\n')

VObj = discAone.VObj
data = len(VObj[0])
model = len(VObj[1])
# Do the gradient elements for the data aspect of the objective function
for i in range(data):
    if VObj[0][i][2] < nY+nU:
      f.write('    grad_f[jt+%d*(Time+1)] += (%s)/(2*Time+1);\n' % (VObj[0][i][2], VObj[0][i][0]))
    elif VObj[0][i][2] < 2*(nY+nU):
      f.write('    grad_f[jt+1+%d*(Time+1)] += (%s)/(2*Time+1);\n' % (VObj[0][i][2] - (nY+nU), VObj[0][i][0]))
    elif VObj[0][i][2] < 3*(nY+nU):
      f.write('    grad_f[(Time+1)*(nU+nY) + %d*Time + jt] += (%s)/(2*Time+1);\n' % (VObj[0][i][2] - 2*(nY+nU), VObj[0][i][0]))
    else:
      f.write('     grad_f[(2*Time+1)*(nY+nU)+%d] += (%s)/(2*Time+1);\n' % (VObj[0][i][2]-3*(nY+nU), VObj[0][i][0]))
f.write('\n')
#  Now do the model aspect gradient
for i in range(model):
    eqnno = VObj[1][i][2]
    if eqnno < nY+nU:
      f.write('    grad_f[jt+%d*(Time+1)] += bounds[%s][3]*(%s)/(2*Time+1);\n' % (VObj[1][i][2], VObj[1][i][1]-1, VObj[1][i][0]))
    elif eqnno < 2*(nY+nU):
      f.write('    grad_f[jt+1+%d*(Time+1)] += bounds[%s][3]*(%s)/(2*Time+1);\n' % (VObj[1][i][2]-(nY+nU), VObj[1][i][1]-1, VObj[1][i][0]))
    elif eqnno < 3*(nY+nU):
      f.write('    grad_f[(Time+1)*(nU+nY) + %d*Time + jt] += bounds[%s][3]*(%s)/(2*Time+1);\n' % (VObj[1][i][2]-2*(nY+nU), VObj[1][i][1]-1, VObj[1][i][0]))
    else:
      f.write('     grad_f[(2*Time+1)*(nY+nU)+%d] += bounds[%s][3]*(%s)/(2*Time+1);\n' % (VObj[1][i][2]-3*(nY+nU), VObj[1][i][1]-1,VObj[1][i][0]))
f.write('\n')
f.write('  } //end for loop\n\n')

# Add code for last gradient element

f.write('// Add last element\n')
f.write('     for(Index i=0;i<nY;i++) {\n')
f.write('        Xval[i] = x[i*(Time+1) + Time];\n')
f.write('        Xvalp1[i] = 0;\n')
f.write('        Xval2[i] = 0;\n')
f.write('     } //end for loop\n\n')

f.write('     for(Index cup=0;cup<nU;cup++) {\n')
f.write('        K11val[cup] = x[nY*(Time+1) + cup*(Time+1) + Time];\n')
f.write('        K11valp1[cup] = 0;\n')
f.write('        K11val2[cup] = 0;\n')

f.write('     } //end for loop\n\n')

for i in range(nM):
    f.write('     Xdval[%d] = %s[2*Time];\n' % (i,discAone.Ldata[i]))
    f.write('     Xdval2[%d] = 0;\n' % i)
    f.write('     Xdvalp1[%d] = 0;\n' % i)

for i in range(nI):
    f.write('     Ival[%d] = %s[2*Time];\n' % (i,discAone.Lstimuli[i]))
    f.write('     Ival2[%d] = 0;\n' % i)
    f.write('     Ivalp1[%d] = 0;\n' % i)

f.write('\n')
f.write('     for(Index i=0;i<nP;i++) {\n')
f.write('        Pval[i] = x[(2*Time+1)*(nY+nU)+i];\n')
f.write('     } //end for loop\n\n')
f.write('\n')

for i in range(data):
    if VObj[0][i][2] < nY+nU:
      f.write('    grad_f[Time+%d*(Time+1)] += (%s)/(2*Time+1);\n' % (VObj[0][i][2], VObj[0][i][0]))

f.write('\n')
f.write('  return true;\n\
}\n\n\n')



# EVAL_G

f.write('// return the value of the constraints: g(x)\n\
bool %s_NLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)\n\
{\n' % probu)
f.write('  assert(n == %d*Time+%d);\n' % (alpha, beta))
f.write('  assert(m == 0);\n\n')


f.write('  return true;\n\
}\n\n\n')



# EVAL_JAC_G

f.write('// return the structure or values of the jacobian\n\
bool %s_NLP::eval_jac_g(Index n, const Number* x, bool new_x,\n\
                            Index m, Index nele_jac, Index* iRow, Index* jCol,\n\
                            Number* values)\n\
{\n\n\
return true;\n\
}\n\n\n' % probu)



# EVAL_H

f.write('// return the structure or values of the hessian\n\
bool %s_NLP::eval_h(Index n, const Number* x, bool new_x,\n\
                       Number obj_factor, Index m, const Number* lambda,\n\
                       bool new_lambda, Index nele_hess, Index* iRow,\n\
                       Index* jCol, Number* values)\n\
{\n\n\
if (values == NULL) {\n\
   // return the structure.  This is a symmetric matrix, fill in the lower left\n\
   // triangle only.\n\n' % probu)

# Set up dictionaries to mark start point and length of each Hessian entry
# dictlength will be a string of either 1, Time, or Time +1 depending on how
# many entries in the Hessian there are for an element
# dictstart determines where each of the parts start based on what has come before

f.write('   // Each non-one Hessian element has its own explicit loop\n\
   // since each element needs a different number of matrix elements\n\n')

dictstart = {}
dictlength = {}

VHes = discAone.VHes

# sma, med, lar track how many 1, T, T+1 there are
sma = 0
med = 0
lar = 0
# count keeps track of the unique index for each entry
# Start at 10 to be sure it is unique
count = 1
for i in range(cat):
  for j in range(i+1):

    # Only take the unique entries: no xp1/xp1 for instance
    if Hesindex[i][j] == -1:
      row = i
      col = j
      temp = bytes.decode(Hessize[i][j])
      start = '%d*(Time+1)+%d*(Time)+%d' % (lar,med,sma)
      if temp == 's':
         length = '1'
         sma = sma + 1
      elif temp == 'm':
         length = 'Time'
         med = med + 1
      elif temp == 'l':
         length = 'Time+1'
         lar = lar + 1 
      d1 = {count:start}  
      d2 = {count:length}

      # Update the Hesindex matrix to reflect count
      Hesindex[i][j] = count
      count = count + 1
      dictstart.update(d1)
      dictlength.update(d2)
      f.write('     for(Index jt=0;jt<%s;jt++) {\n' % length)
      f.write('     iRow[%s+jt] = ' % start)
      
      # Now set out which row and column the derivatives are in respect to
      if row < 2*(nY+nU):
        if row % 2 == 0:
           f.write('(Time+1)*%d+jt;\n' % (row/2))
        else:
           f.write('(Time+1)*%d+jt+1;\n' % (row/2))
      elif row < 3*(nY+nU):
        f.write('(Time+1)*%d+Time*%d+jt;\n' % (nY+nU,row-2*(nY+nU)))
      else:
        f.write('2*Time*%d+%d;\n' % (nY+nU, row-2*(nY+nU)))

      f.write('     jCol[%s+jt] = ' % start)
      if col < 2*(nY+nU):
        if col % 2 == 0:
          f.write('(Time+1)*%d+jt;\n' % (col/2))
        else:
          f.write('(Time+1)*%d+jt+1;\n' % (col/2))
      elif col < 3*(nY+nU):
        f.write('(Time+1)*%d+Time*%d+jt;\n' % (nY+nU,col-2*(nY+nU)))
      else:
        f.write('2*Time*%d+%d;\n' % (nY+nU, col-2*(nY+nU)))
      f.write('   }\n')

    else:

      # Change Hesindex to reflect the correct index
      row = i
      col = j
      if i < 2*(nY+nU):
        Hesindex[i][j] = Hesindex[i-1][j-1]
      else:
        Hesindex[i][j] = Hesindex[i][j-1]

f.write('}\n\
else {\n\
  // return the values.  This is a symmetric matrix, fill the lower left\n\
  // triangle only\n\
  // initialize the values array\n\
  // Point to the initial starting spot for the Hessian elements\n\n\
  for(Index jt=0;jt<%d*(Time+1)+%d*Time+%d;jt++) values[jt] = 0.; // Initialize matrix' % (lar, med ,sma))
f.write('\n\n   // fill the objective portion\n\n')

# Loop over all other elements

f.write('  for(Index jt=0;jt<Time;jt++) {\n\n')
f.write('     for(Index i=0;i<nY;i++) {\n')
f.write('        Xval[i] = x[jt + i*(Time+1)];\n')
f.write('        Xvalp1[i] = x[jt + i*(Time+1) + 1];\n')
f.write('        Xval2[i] = x[(Time+1)*(nY+nU) + jt + i*(Time)];\n')
f.write('     } //end for loop\n\n')

f.write('     for(Index cup=0;cup<nU;cup++) {\n')
f.write('        K11val[cup] = x[jt + nY*(Time+1) + cup*(Time+1)];\n')
f.write('        K11valp1[cup] = x[jt + nY*(Time+1) + cup*(Time+1) + 1];\n')
f.write('        K11val2[cup] = x[(Time+1)*(nY+nU) + (nY+cup)*Time + jt];\n')

f.write('     } //end for loop\n\n')

for i in range(nM):
    f.write('     Xdval[%d] = %s[2*jt];\n' % (i,discAone.Ldata[i]))
    f.write('     Xdval2[%d] = %s[2*jt+1];\n' % (i,discAone.Ldata[i]))
    f.write('     Xdvalp1[%d] = %s[2*jt+2];\n' % (i,discAone.Ldata[i]))

for i in range(nI):
    f.write('     Ival[%d] = %s[2*jt];\n' % (i,discAone.Lstimuli[i]))
    f.write('     Ival2[%d] = %s[2*jt+1];\n' % (i,discAone.Lstimuli[i]))
    f.write('     Ivalp1[%d] = %s[2*jt+2];\n' % (i,discAone.Lstimuli[i]))

f.write('\n')
f.write('     for(Index i=0;i<nP;i++) {\n')
f.write('        Pval[i] = x[(2*Time+1)*(nY+nU)+i];\n')
f.write('     } //end for loop\n')
f.write('\n')
# Do singletons first
for i in range(len(VHes)):
   Hesrow = VHes[i][1]
   Hescol = VHes[i][2]
   if bytes.decode(Hessize[Hesrow][Hescol]) == 's':
      constraint = VHes[i][0]
      count = Hesindex[Hesrow][Hescol]
      string = VHes[i][4]
      start = dictstart[count]
      if constraint == 0:
          R = '1'
      else:
          R = 'bounds[%s][3]' % (constraint - 1) 
      f.write('    values[%s] += ' % start)
      f.write('obj_factor*%s*(%s)/(2*Time+1);\n' % (R,string))
#  Check if this works correctly!
for i in range(len(VHes)):
    Hesrow = VHes[i][1]
    Hescol = VHes[i][2]
    if bytes.decode(Hessize[Hesrow][Hescol]) != 's':
       constraint = VHes[i][0]
       count = Hesindex[Hesrow][Hescol]
       string = VHes[i][4]
       start = dictstart[count]
       if constraint == 0:
          R = '1'
       else:
          R = 'bounds[%s][3]' % (constraint - 1) 
       toggle = VHes[i][3]
       if toggle == 0:      
          f.write('   values[%s+jt] += ' % start)
          f.write('obj_factor*%s*(%s)/(2*Time+1);\n' % (R, string))
       elif toggle == 1:
          newstart = start + '+ 1'
          f.write('   values[%s+jt] += ' % newstart)
          f.write('obj_factor*%s*(%s)/(2*Time+1);\n' % (R, string))

f.write('   } // end for loop \n\n')

#  Add in code for last element from data part
f.write('// Add last element\n')
f.write('     for(Index i=0;i<nY;i++) {\n')
f.write('        Xval[i] = x[Time + i*(Time+1)];\n')
f.write('        Xvalp1[i] = 0;\n')
f.write('        Xval2[i] = 0;\n')
f.write('     } //end for loop\n\n')

f.write('     for(Index cup=0;cup<nU;cup++) {\n')
f.write('        K11val[cup] = x[Time + nY*(Time+1) + cup*(Time+1)];\n')
f.write('        K11valp1[cup] = 0;\n')
f.write('        K11val2[cup] = 0;\n')

f.write('     } //end for loop\n\n')

for i in range(nM):
    f.write('     Xdval[%d] = %s[2*Time];\n' % (i,discAone.Ldata[i]))
    f.write('     Xdval2[%d] = 0;\n' % i)
    f.write('     Xdvalp1[%d] = 0;\n' % i)

for i in range(nI):
    f.write('     Ival[%d] = %s[2*Time];\n' % (i,discAone.Lstimuli[i]))
    f.write('     Ival2[%d] = 0;\n' % i)
    f.write('     Ivalp1[%d] = 0;\n' % i)

f.write('\n')
f.write('     for(Index i=0;i<nP;i++) {\n')
f.write('        Pval[i] = x[(2*Time+1)*(nY+nU)+i];\n')
f.write('     } //end for loop\n\n')
f.write('\n')

for i in range(len(VHes)):
    constraint = VHes[i][0]
    if constraint == 0:
        Hesrow = VHes[i][1]
        Hescol = VHes[i][2]
        if bytes.decode(Hessize[Hesrow][Hescol]) == 'l':
           string = VHes[i][4]
           count = Hesindex[Hesrow][Hescol]
           start = dictstart[count]
           newstart = start + '+ Time'
           f.write('   values[%s] += ' % newstart)
           f.write('obj_factor*(%s)/(2*Time+1);\n' % string)

f.write('  } // end else \n\n')



f.write('   return true;\n\
}\n\n\n')


# FINALIZE_SOLUTION


f.write('\
void %s_NLP::finalize_solution(SolverReturn status,\n\
                        Index n, const Number* x, const Number* z_L, const Number* z_U,\n\
                        Index m, const Number* g, const Number* lambda,\n\
                        Number obj_value,\n\
                        const IpoptData* ip_data,\n\
                        IpoptCalculatedQuantities* ip_cq)\n\
{\n' % probu)
f.write('\
  FILE *OUTPUT1;\n\
  char filename[20];\n\
  if ((file_fmt < 2) || (modid == 0)){\n\
    sprintf(filename,"D%d_M%d_IC%d.dat", nY,nM,taskid);\n\
  }\n\
  else if ((file_fmt >= 2) && (modid != 0)){\n\
    sprintf(filename,"D%d_M%d_PATH%d_IC%d.dat", nY,nM,pathid,taskid);\n\
  }\n\
  if(beta==0){\n\
    OUTPUT1 = fopen (filename,"w");\n\
  }\n\
  else{\n\
    if ((output_fmt == 1) || (output_fmt < 0)){\n\
      OUTPUT1 = fopen (filename,"w");\n\
    }\n\
    else \n\
      OUTPUT1 = fopen (filename,"a");\n\
  }\n')
f.write('\
  // Write solution for annealing\n\
  for(Index i=0;i<Ntotal-nP;i++){\n\
    solution[i] = x[i];\n\
  }\n\
  for(Index i=0;i<nP;i++){\n\
    solution[Ntotal-nP+i] = x[(2*Time+1)*(nU+nY)+i];\n\
  }\n\n\
  // Write to file\n\
  if (output_fmt != 1){\n\
    fprintf(OUTPUT1, "%d %d %e ",beta, (status == SUCCESS), obj_value);\n\
  }\n\
  for (Index i=0;i<Time;i++) {\n\
    for (Index j=0;j<nY;j++) {\n\
      fprintf(OUTPUT1,"%e ", x[j*(Time+1)+i]);\n\
    }\n\
    if (abs(output_fmt) == 2){\n\
      for (Index cup=0;cup<nU;cup++) {\n\
        fprintf(OUTPUT1,"%e ", x[nY*(Time+1)+ cup*(Time+1)+ i]);\n\
      }\n\
    }\n\
    for (Index j=0;j<nY;j++) {\n\
      fprintf(OUTPUT1,"%e ", x[(nY+nU)*(Time+1) + j*Time + i]);\n\
    }\n\
    if (abs(output_fmt) == 2){\n\
      for (Index cup=0;cup<nU;cup++) {\n\
        fprintf(OUTPUT1,"%e ", x[(nY+nU)*(Time+1) + nY*Time + cup*Time + i]);\n\
      }\n\
    }\n\
  }\n\
  for (Index j=0;j<nY;j++) {\n\
     fprintf(OUTPUT1,"%e ", x[j*(Time+1) + Time]);\n\
  }\n\
  if (abs(output_fmt) == 2){\n\
    for (Index cup=0;cup<nU;cup++) {\n\
       fprintf(OUTPUT1,"%e ", x[nY*(Time+1) + cup*(Time+1) + Time]);\n\
    }\n\
  }\n\
  for (Index j=0;j<nP;j++) {\n\
     fprintf(OUTPUT1,"%e ", x[(2*Time+1)*(nY+nU)+j]);\n\
  }\n\
  fprintf(OUTPUT1,"\\n");\n\
  fclose (OUTPUT1);\n\
  \n\n')
f.write('}\n')       



f.close( )
