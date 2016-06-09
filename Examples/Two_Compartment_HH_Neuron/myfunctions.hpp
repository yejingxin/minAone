#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

double tol = 5.0;

double ghk(double V)
{
  double Theta = 13.5;
  double value=0.0;
  if (abs(V) < tol){
     value = -Theta - 0.5*V - V*V/(12.0*Theta)+ V*V*V*V/(720.0*Theta*Theta*Theta) - V*V*V*V*V*V/(30240.0*Theta*Theta*Theta*Theta*Theta) + V*V*V*V*V*V*V*V/(1209600.0*Theta*Theta*Theta*Theta*Theta*Theta*Theta) - V*V*V*V*V*V*V*V*V*V/(47900160.0*Theta*Theta*Theta*Theta*Theta*Theta*Theta*Theta*Theta);
     }
  else {
     value = V/(exp(-V/Theta) - 1.0);
     }
  return value;
}

double ghkjac(double V, int n)
{
  double Theta = 13.5;
  double value=0.0;
        if (abs(V) < tol){
           value = - 0.5- V/(6.0*Theta)+ V*V*V/(180.0*Theta*Theta*Theta)- V*V*V*V*V/(5040.0*Theta*Theta*Theta*Theta*Theta) + V*V*V*V*V*V*V/(151200.0*Theta*Theta*Theta*Theta*Theta*Theta*Theta) - V*V*V*V*V*V*V*V*V/(4790016.0*Theta*Theta*Theta*Theta*Theta*Theta*Theta*Theta*Theta);
           }
        else {
           value = V*exp(-V/Theta)/(Theta*(-1.0 + exp(-V/Theta))*(-1.0 + exp(-V/Theta))) + 1.0/(-1.0 + exp(-V/Theta));
           }
  return value;
}

double ghkhes(double V, int n, int m)
  {
  double Theta = 13.5;
  double value=0.0;
               if (abs(V) < tol){
                  value = -1.0/(6.0*Theta)+ V*V/(60.0*Theta*Theta*Theta)-V*V*V*V/(1008.0*Theta*Theta*Theta*Theta*Theta) + V*V*V*V*V*V/(21600.0*Theta*Theta*Theta*Theta*Theta*Theta*Theta) - V*V*V*V*V*V*V*V/(532224.0*Theta*Theta*Theta*Theta*Theta*Theta*Theta*Theta*Theta);
  		  }
               else {
                    value = -V*exp(-V/Theta)/(Theta*Theta*(-1.0 + exp(-V/Theta))*(-1.0 + exp(-V/Theta))) + 2.0*V*exp(-2.0*V/Theta)/(Theta*Theta*(-1.0 + exp(-V/Theta))*(-1.0 + exp(-V/Theta))*(-1.0 + exp(-V/Theta))) + 2.0*exp(-V/Theta)/(Theta*(-1.0 + exp(-V/Theta))*(-1.0 + exp(-V/Theta)));
                   }
  return value;
}

