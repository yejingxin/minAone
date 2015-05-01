// simple_nakl_main.cpp
// Main file for use with IPOPT

#include "IpIpoptApplication.hpp"
#include "simple_naklminAone_nlp.hpp"

// for printf
#ifndef HAVE_CSTDIO
#define HAVE_CSTDIO
# include <cstdio>
#else
# ifndef HAVE_STDIO_H
#  define HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif

using namespace Ipopt;

int main(int argv, char* argc[])
{

  // Create a new instance of your nlp
  //  (use a SmartPtr, not raw)
  SmartPtr<TNLP> mynlp = new SIMPLE_NAKL_NLP((argv<2)?0:atoi(argc[1]));

  // Create a new instance of IpoptApplication
  //  (use a SmartPtr, not raw)
  SmartPtr<IpoptApplication> app = new IpoptApplication();

  // Change some options
  // Note: The following choices are only examples, they might not be
  //       suitable for your optimization problem.
  app->Options()->SetNumericValue("tol", 1e-12);
  app->Options()->SetStringValue("mu_strategy", "adaptive");
  app->Options()->SetStringValue("output_file", "ipopt.out");
  // The following overwrites the default name (ipopt.opt) of the
  // options file
  app->Options()->SetStringValue("option_file_name", "simple_nakl.opt");

  // Intialize the IpoptApplication and process the options
  ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Solve_Succeeded) {
    printf("\n\n *** Error during initialization!\n");
    return (int) status;
  }

  // Ask Ipopt to solve the problem
    for (Index i=0;i<100;i++) { 
	status = app->OptimizeTNLP(mynlp);
	if(!(*dynamic_cast<SIMPLE_NAKL_NLP*>(GetRawPtr(mynlp))).changeRf()) break; 
  }


  if (status == Solve_Succeeded) {
    printf("\n\n*** The problem solved!\n");
  }
  else {
    printf("\n\n*** The problem FAILED!\n");
  }

  // As the SmartPtrs go out of scope, the reference count
  // will be decremented and the objects will automatically
  // be deleted.

  return (int) status;
}
