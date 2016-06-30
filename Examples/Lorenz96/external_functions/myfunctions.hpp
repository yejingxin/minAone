#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

double Lorenzvectorfield(double x, double x_p1, double x_m1, double x_m2, double forcing)
{
  // Lorenz model equations; x = x_n, x_m1 = x_(n-1), x_m2 = x_(n-2), x_p1 = x_(n+1)
  double value=0.0;
  value = x_m1*(x_p1-x_m2) - x + forcing;
  return value;
}

double Lorenzvectorfieldjac(double x, double x_p1, double x_m1, double x_m2, double forcing, int n)
{
  double value=0.0;
  switch (n) {
    case 1:
    // derivative with respect to x
      value = -1; break;
    case 2:
    // derivative with respect to x_p1
      value = x_m1; break;
    case 3:
    // derivative with respect to x_m1
      value = x_p1 - x_m2; break;
    case 4:
    // derivative with respect to x_m2
      value = -x_m1; break;
    case 5:
    // derivative with respect to forcing
      value = 1; break;
    default:
      cout << "Error in user-defined function jacobian definition";break;
  } // end switch
  return value;
}

double Lorenzvectorfieldhes(double x, double x_p1, double x_m1, double x_m2, double forcing, int n, int m)
{
  double value=0.0;
  switch (n) {
    case 1:
    // derivative with respect to x
      value = 0; break;
    case 2:
    // derivative with respect to x_p1
      if (m == 3){
         value = 1; break;
      }
      else{
         value = 0; break;
      }      
    case 3:
    // derivative with respect to x_m1
        if (m == 2){
           value = 1; break;
        }
        else if (m == 4){
           value = -1; break;
        }      
        else{
           value = 0; break;
        }      
    case 4:
    // derivative with respect to x_m2
        if (m == 3){
           value = -1; break;
        }
        else{
           value = 0; break;
        }      
    case 5:
    // derivative with respect to forcing
      value = 0; break;
    default:
      cout << "Error in user-defined function Hessian";break;
  } 
  return value;
}

