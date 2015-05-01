#####################################################
#
#  20 October 2009
#  Bryan A. Toth
#  University of California, San Diego
#  btoth@physics.ucsd.edu
# 
#
#  30 September 2014
#  Modified by Jingxin Ye 
#  University of California, San Diego
#  j9ye@physics.ucsd.edu
#
#  This script writes the options file for a C++ IPOPT 
#  program defined by the vector field in the file
#  equations.txt and written by the script makecode.py.
#  Thie file produces a "bare bones" options file -
#  explanation of the options and additional options
#  can be found using the IPOPT documentation, available
#  at http://www.coin-or.org/Ipopt/documentation/
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
probl = prob.lower()
# Name of opt file
FILE = probl + '.opt'

f = open(FILE,'w')

f.write('\
# Set the max number of iterations\n\
max_iter 10000\n\
# set derivative test\n\
#derivative_test second-order\n\
# set termination criteria\n\
tol 1.0e-12\n\
#dual_inf_tol 0.001\n\
#compl_inf_tol 1.0e-12\n\
#constr_viol_tol 1.0e-8\n\
#acceptable_tol 1.0e-10\n\
#acceptable_iter \n\
#turn off the NLP scaling\n\
#nlp_scaling_method none\n\
#mehrotra_algorithm yes\n\
mu_strategy adaptive\n\
adaptive_mu_globalization never-monotone-mode\n\
linear_solver ma97\n\
#linear_system_scaling none\n\
bound_relax_factor 0\n\
#ma27_pivtol 1.0e-6\n\
\n')

f.close()
