#####################################################
#
#  20 October 2009
#  Bryan A. Toth
#  University of California, San Diego
#  btoth@physics.ucsd.edu
# 
#  This script writes the class header file for a C++ IPOPT 
#  program defined by the vector field in the file
#  equations.txt and written by the script makecode.py.
#  The file written by this script defines a class that is 
#  further described in the file written by makecode.py.
#
#  This script has been developed as part of a suite of 
#  python scripts to define a dynamic parameter estimation
#  problem using the optimization software IPOPT, but is 
#  generally applicable to any application needing
#  discAzerod derivatives of a vector field.
#
######################################################

import discAone

prob = discAone.Problem
probu = prob.upper()
probl = prob.lower()

FILE = probl + 'minAone_nlp.hpp'

# The name of the IPOPT header file
f = open(FILE,'w')

f.write('// %s.hpp\n\
// Header file for %s.cpp\n\
// For use with IPOPT\n' % (probl, probl))

f.write('\n\
\n\
#ifndef __%s_NLP_HPP__\n\
#define __%s_NLP_HPP__\n\
\n\
#include "IpTNLP.hpp"\n\
#include <iostream>\n\
#include <fstream>\n\
#include <string>\n\
#include <stdlib.h>\n\
#include <cstring>\n\
#include "assert.h"\n\
\n\
using namespace std;\n\
using namespace Ipopt;\n' % (probu, probu))


f.write('\n\n\
class %s_NLP : public TNLP\n\
{\n\
public:\n\
  /** default constructor */\n\
  %s_NLP(int id);\n\
\n\
  /** default destructor */\n\
  virtual ~%s_NLP();\n\
\n\
  /**@name Overloaded from TNLP */\n\
  //@{\n\
  /** Method to return some info about the nlp */\n\
  virtual bool changeRf();\n\
  virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,\n\
                            Index& nnz_h_lag, IndexStyleEnum& index_style);\n\
\n\
  /** Method to return the bounds for my problem */\n\
  virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,\n\
                               Index m, Number* g_l, Number* g_u);\n\
\n\
  /** Method to return the starting point for the algorithm */\n\
  virtual bool get_starting_point(Index n, bool init_x, Number* x,\n\
                                  bool init_z, Number* z_L, Number* z_U,\n\
                                  Index m, bool init_lambda,\n\
                                  Number* lambda);\n\
\n\
  /** Method to return the objective value */\n\
  virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);\n\
\n\
  /** Method to return the gradient of the objective */\n\
  virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);\n\
\n\
  /** Method to return the constraint residuals */\n\
  virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);\n\
\n\
  /** Method to return:\n\
   *   1) The structure of the jacobian (if "values" is NULL)\n\
   *   2) The values of the jacobian (if "values" is not NULL)\n\
   */\n\
  virtual bool eval_jac_g(Index n, const Number* x, bool new_x,\n\
                          Index m, Index nele_jac, Index* iRow, Index *jCol,\n\
                          Number* values);\n\
\n\
  /** Method to return:\n\
   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)\n\
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)\n\
   */\n\
  virtual bool eval_h(Index n, const Number* x, bool new_x,\n\
                      Number obj_factor, Index m, const Number* lambda,\n\
                      bool new_lambda, Index nele_hess, Index* iRow,\n\
                      Index* jCol, Number* values);\n\
\n\
  //@}\n\
\n\
  /** @name Solution Methods */\n\
  //@{\n\
  /** This method is called when the algorithm is complete so the TNLP can store/write the solution */\n\
  virtual void finalize_solution(SolverReturn status,\n\
                                 Index n, const Number* x, const Number* z_L, const Number* z_U,\n\
                                 Index m, const Number* g, const Number* lambda,\n\
                                 Number obj_value,\n\
				 const IpoptData* ip_data,\n\
				 IpoptCalculatedQuantities* ip_cq);\n\
  //@}\n\
\n\
private:\n\
  /**@name Methods to block default compiler methods.\n\
   * The compiler automatically generates the following three methods.\n\
   *  Since the default compiler implementation is generally not what\n\
   *  you want (for all but the most simple classes), we usually \n\
   *  put the declarations of these methods in the private section\n\
   *  and never implement them. This prevents the compiler from\n\
   *  implementing an incorrect "default" behavior without us\n\
   *  knowing. (See Scott Meyers book, "Effective C++")\n\
   *  \n\
   */\n\
  //@{\n\
  //  %s_NLP();\n' % (probu, probu, probu, probu))

for i in range(discAone.nM):
  f.write('\
  double* %s;\n\
  double* %s;\n' % (discAone.Ldata[i], discAone.Ldata[i]+'dummy'))

for i in range(discAone.nI):
  f.write('\
  double* %s;\n\
  double* %s;\n' % (discAone.Lstimuli[i], discAone.Lstimuli[i]+'dummy'))

f.write('\
  double tmpdata;')

#  int %s;\n' % (discAzero.Ldata[i], discAzero.Ldata[i]+'dummy',discAzero.Ldata[i]+'skip'))
  
f.write('\
\n\
  int nU;\n\
  int nP;\n\
  int nY;\n\
  int nI;\n\
  int nF;\n\
  int nM;\n\
  int skip;\n\
  int Ntotal;\n\
  double* K11val;\n\
  double* K11val2;\n\
  double* K11valp1;\n\
  double* dK11val;\n\
  double* dK11val2;\n\
  double* dK11valp1;\n\
  double* Xdval;\n\
  double* Xdval2;\n\
  double* Xdvalp1;\n\
  double* Xval;\n\
  double* Xval2;\n\
  double* Xvalp1;\n\
  double* Pval;\n\
  double* Ival;\n\
  double* Ival2;\n\
  double* Ivalp1;\n\
  int Time;\n\
  double hstep;\n\
  string buffer;\n\
  string* specs;\n\
  double** bounds;\n\
  double* Rf0;\n\
  double* solution;\n\
  double alpha;\n\
  int beta, delta_beta, max_beta;\n\
  int taskid;\n\
  int pathid;\n\
  int modid;\n\
  int file_fmt;\n\
  int output_fmt;\n\
  %s_NLP(const %s_NLP&);\n\
  %s_NLP& operator=(const %s_NLP&);\n\
  //@}\n\
};\n\
\n\
\n\
#endif\n' % (probu, probu, probu, probu))

f.close()
