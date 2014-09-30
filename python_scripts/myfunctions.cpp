#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;
double tol = .001;

double dingo(double a, double b, double c, double d)
{
  double value=0.0;
  value = 1.0/(a + exp(b*(c+d)));
  return value;
}

double dingojac(double a, double b, double c, double d, int n)
{
  double value=0.0;
  switch (n) {
    case 1:
    // derivative with respect to a
      value = -1.0*pow(dingo(a,b,c,d),2); break;
    case 2:
    // derivative with respect to b
      value = -(c+d)/(a*a*exp(-b*(c+d))+2*a+exp(b*(c+d))); break;
    case 3:
    // derivative with respect to c
      value = -b/(a*a*exp(-b*(c+d))+2*a+exp(b*(c+d))); break;
    case 4:
    // derivative with respect to d
      value = -b/(a*a*exp(-b*(c+d))+2*a+exp(b*(c+d))); break;
    default:
      cout << "Error";break;
  } // end switch
  return value;
}

double dingohes(double a, double b, double c, double d, int n, int m)
{
  double value=0.0;
  switch (n) {
    case 1:
    // derivative with respect to a
        switch (m) {
              case 1:
              // derivative with respect to a
               value = 2.0*pow(dingo(a,b,c,d),3); break;
              case 2:
              // derivative with respect to b
               value = 2.0*(c+d)/(a*a*a*exp(-b*(c+d))+exp(2*b*(c+d))+3*a*a+3*a*exp(b*(c+d))); break;
              case 3:
              // derivative with respect to c
               value = 2.0*b/(a*a*a*exp(-b*(c+d))+exp(2*b*(c+d))+3*a*a+3*a*exp(b*(c+d))); break;
              case 4:
              // derivative with respect to d
               value = 2.0*b/(a*a*a*exp(-b*(c+d))+exp(2*b*(c+d))+3*a*a+3*a*exp(b*(c+d))); break;
              default:
               cout << "Error";break;
           } // end switch
           break;
    case 2:
    // derivative with respect to b
        switch (m) {
              case 1:
              // derivative with respect to a
               value = 2.0*(c+d)/(a*a*a*exp(-b*(c+d))+exp(2*b*(c+d))+3*a*a+3*a*exp(b*(c+d))); break;
              case 2:
              // derivative with respect to b
               value = pow(c+d,2)*(1.0-a*exp(-b*(c+d)))/(a*a*a*exp(-2*b*(c+d))+exp(b*(c+d))+3*a*a*exp(-b*(c+d))+3*a); break;
              case 3: 
              // derivative with respect to c
               value = (-1+b*(c+d)-a*exp(-b*(c+d))*(1+b*(c+d)))/(a*a*a*exp(-2*b*(c+d))+exp(b*(c+d))+3*a*a*exp(-b*(c+d))+3*a); break;
              case 4:
              // derivative with respect to d
               value = (-1+b*(c+d)-a*exp(-b*(c+d))*(1+b*(c+d)))/(a*a*a*exp(-2*b*(c+d))+exp(b*(c+d))+3*a*a*exp(-b*(c+d))+3*a); break;
              default:
               cout << "Error";break;
          } // end switch
          break;
    case 3:
    // derivative with respect to c
        switch (m) {
              case 1:
              // derivative with respect to a
               value = 2.0*b/(a*a*a*exp(-b*(c+d))+exp(2*b*(c+d))+3*a*a+3*a*exp(b*(c+d))); break;
              case 2:
              // derivative with respect to b
               value = (-1+b*(c+d)-a*exp(-b*(c+d))*(1+b*(c+d)))/(a*a*a*exp(-2*b*(c+d))+exp(b*(c+d))+3*a*a*exp(-b*(c+d))+3*a); break;
              case 3:
              // derivative with respect to c
               value = (b*b*(1.0-a*exp(-b*(c+d))))/(a*a*a*exp(-2*b*(c+d))+exp(b*(c+d))+3*a*a*exp(-b*(c+d))+3*a); break;
              case 4:
              // derivative with respect to d
               value = (b*b*(1.0-a*exp(-b*(c+d))))/(a*a*a*exp(-2*b*(c+d))+exp(b*(c+d))+3*a*a*exp(-b*(c+d))+3*a); break;
              default:
               cout << "Error";break;
          } // end switch
          break;
    case 4:
    // derivative with respect to d
        switch (m) {
              case 1:
              // derivative with respect to a
               value = 2.0*b/(a*a*a*exp(-b*(c+d))+exp(2*b*(c+d))+3*a*a+3*a*exp(b*(c+d))); break;
              case 2:
              // derivative with respect to b
               value = (-1+b*(c+d)-a*exp(-b*(c+d))*(1+b*(c+d)))/(a*a*a*exp(-2*b*(c+d))+exp(b*(c+d))+3*a*a*exp(-b*(c+d))+3*a); break;
              case 3:
              // derivative with respect to c
               value = (b*b*(1.0-a*exp(-b*(c+d))))/(a*a*a*exp(-2*b*(c+d))+exp(b*(c+d))+3*a*a*exp(-b*(c+d))+3*a); break;
              case 4:
              // derivative with respect to d
               value = (b*b*(1.0-a*exp(-b*(c+d))))/(a*a*a*exp(-2*b*(c+d))+exp(b*(c+d))+3*a*a*exp(-b*(c+d))+3*a); break;
              default:
               cout << "Error";break;
         } // end switch
         break;
    default:
      cout << "Error";break;
  } // end switch
  return value;
}

double efunc(double a, double b, double c, double d)
{
   double value=0.0;
   if (abs(b+d) < tol) {
     value = a/c - 0.5*a*(b+d)+1.0/12*a*c*pow(b+d,2)-1.0/720*pow(b+d,4);
     }
   else  value = a*(b+d)/(exp((b+d)*c)-1.0);
   return value;
}

double efuncjac(double a, double b, double c, double d, int n)
{
   double value=0.0;
   switch (n) {
     case 1:
     // derivative with respect to a
       value = efunc(a,b,c,d)/a; break;
     case 2:
     // derivative with respect to b
       if (abs(b+d) < tol){
         value = -a/2+1.0/6*a*c*(b+d)-1.0/180*a*pow(c,3)*pow((b+d),3); break;
         }
       else {
         value = -(1+(-1+c*(b+d))*exp((b+d)*c))/(a*pow((b+d),2))*pow(efunc(a,b,c,d),2); break;
	 }
     case 3:
     // derivative with respect to c
       value = -pow(efunc(a,b,c,d),2)*exp((b+d)*c)/a; break; 
     case 4:
     // derivative with respect to d 
       if (abs(b+d) < tol){
         value = -a/2+1.0/6*a*c*(b+d)-1.0/180*a*pow(c,3)*pow((b+d),3); break;
         }
       else {
         value = -(1+(-1+c*(b+d))*exp((b+d)*c))/(a*pow((b+d),2))*pow(efunc(a,b,c,d),2); break;
         }	 

     default:
       cout << "I don't know what to do jac";break;

   } // end switch
   return value;
}

double efunches(double a, double b, double c, double d, int n, int m)
{
   double value=0.0;
   switch (n) {
     case 1:
     // derivative with respect to a
       switch (m) {
         case 1:
         // derivative with respect to a
	 value = 0; break;
         case 2:
         // derivative with respect to b
           if (abs(b+d) < tol){
             value = -1.0/2+1.0/6*c*(b+d)-1.0/180*pow(c,3)*pow((b+d),3); break;
             }
           else {
             value = -(1+(-1+c*(b+d))*exp((b+d)*c))/(pow(a,2)*pow((b+d),2))*pow(efunc(a,b,c,d),2); break;
	     }
         case 3:
         // derivative with respect to c
           if (abs(b+d) < tol){
             value = -1.0/c/c+1.0/12*pow((b+d),2)-1.0/240*c*c*pow((b+d),4); break;
	     } 
	   else {
	     value = -(exp((b+d)*c))/a/a*pow(efunc(a,b,c,d),2); break; 
             }
	 case 4:
         // derivative with respect to d 
           if (abs(b+d) < tol){
             value = -1.0/2+1.0/6*c*(b+d)-1.0/180*pow(c,3)*pow((b+d),3); break;
             }
           else {
             value = -(1+(-1+c*(b+d))*exp((b+d)*c))/(pow(a,2)*pow((b+d),2))*pow(efunc(a,b,c,d),2); break;
             }	 

         default:
           cout << "I don't know what to do";break;
         } // end switch
         break;
     case 2:
     // derivative with respect to b
       switch (m) {
         case 1:
         // derivative with respect to a
           if (abs(b+d) < tol){
             value = -1.0/2+1.0/6*c*(b+d)-1.0/180*pow(c,3)*pow((b+d),3); break;
             }
           else {
             value = -(1+(-1+c*(b+d))*exp((b+d)*c))/(pow(a,2)*pow((b+d),2))*pow(efunc(a,b,c,d),2); break;
	     }
         case 2:
         // derivative with respect to b
           if (abs(b+d) < tol){
             value = a*c/6-1.0/60*a*pow(c,3)*pow((b+d),2); break;
             }
           else {
             value = c*exp((b+d)*c)*(2.0-2.0*exp((b+d)*c)+b*c*(1+exp((b+d)*c))+c*d*(1+exp((b+d)*c)))/a/a/pow((b+d),3)*pow(efunc(a,b,c,d),3); break;
	     }
         case 3:
         // derivative with respect to c
           if (abs(b+d) < tol){
             value = 1.0/6*a*(b+d)-1.0/60*a*c*c*pow((b+d),3); break;      
	     } 
	   else {
             value = exp((b+d)*c)*(2.0-2.0*exp((b+d)*c)+b*c*(1+exp((b+d)*c))+c*d*(1+exp((b+d)*c)))/a/a/pow((b+d),2)*pow(efunc(a,b,c,d),3); break;
             }
	 case 4:
         // derivative with respect to d 
           if (abs(b+d) < tol){
             value = a*c/6-1.0/60*a*pow(c,3)*pow((b+d),2); break;
             }
           else {
             value = c*exp((b+d)*c)*(2.0-2.0*exp((b+d)*c)+b*c*(1+exp((b+d)*c))+c*d*(1+exp((b+d)*c)))/a/a/pow((b+d),3)*pow(efunc(a,b,c,d),3); break;
	     }

         default:
           cout << "I don't know what to do";break;
         } // end switch
         break;
     case 3:
     // derivative with respect to c
       switch (m) {
         case 1:
         // derivative with respect to a
           if (abs(b+d) < tol){
             value = -1.0/c/c+1.0/12*pow((b+d),2)-1.0/240*c*c*pow((b+d),4); break;
	     } 
	   else {
	     value = -(exp((b+d)*c))/a/a*pow(efunc(a,b,c,d),2); break; 
             }
         case 2:
         // derivative with respect to b
           if (abs(b+d) < tol){
             value = 1.0/6*a*(b+d)-1.0/60*a*c*c*pow((b+d),3); break;      
	     } 
	   else {
             value = exp((b+d)*c)*(2.0-2.0*exp((b+d)*c)+b*c*(1+exp((b+d)*c))+c*d*(1+exp((b+d)*c)))/a/a/pow((b+d),2)*pow(efunc(a,b,c,d),3); break;
             }
         case 3:
         // derivative with respect to c
           value = exp((b+d)*c)*(1+exp((b+d)*c))/pow(a,2)*pow(efunc(a,b,c,d),3); break; 
         case 4:
         // derivative with respect to d 
           if (abs(b+d) < tol){
             value = 1.0/6*a*(b+d)-1.0/60*a*c*c*pow((b+d),3); break;      
	     } 
	   else {
             value = exp((b+d)*c)*(2.0-2.0*exp((b+d)*c)+b*c*(1+exp((b+d)*c))+c*d*(1+exp((b+d)*c)))/a/a/pow((b+d),2)*pow(efunc(a,b,c,d),3); break;
             }
         
	 default:
           cout << "I don't know what to do";break;
         } // end switch
         break;
     case 4:
     // derivative with respect to d 
       switch (m) {
         case 1:
           if (abs(b+d) < tol){
             value = -1.0/2+1.0/6*c*(b+d)-1.0/180*pow(c,3)*pow((b+d),3); break;
             }
           else {
             value = -(1+(-1+c*(b+d))*exp((b+d)*c))/(pow(a,2)*pow((b+d),2))*pow(efunc(a,b,c,d),2); break;
             }	 
         // derivative with respect to a
         case 2:
         // derivative with respect to b
           if (abs(b+d) < tol){
             value = a*c/6-1.0/60*a*pow(c,3)*pow((b+d),2); break;
             }
           else {
             value = c*exp((b+d)*c)*(2.0-2.0*exp((b+d)*c)+b*c*(1+exp((b+d)*c))+c*d*(1+exp((b+d)*c)))/a/a/pow((b+d),3)*pow(efunc(a,b,c,d),3); break;
	     }
         case 3:
         // derivative with respect to c
           if (abs(b+d) < tol){
             value = 1.0/6*a*(b+d)-1.0/60*a*c*c*pow((b+d),3); break;      
	     } 
	   else {
             value = exp((b+d)*c)*(2.0-2.0*exp((b+d)*c)+b*c*(1+exp((b+d)*c))+c*d*(1+exp((b+d)*c)))/a/a/pow((b+d),2)*pow(efunc(a,b,c,d),3); break;
             }
         case 4:
         // derivative with respect to d 
           if (abs(b+d) < tol){
             value = a*c/6-1.0/60*a*pow(c,3)*pow((b+d),2); break;
             }
           else {
             value = c*exp((b+d)*c)*(2.0-2.0*exp((b+d)*c)+b*c*(1+exp((b+d)*c))+c*d*(1+exp((b+d)*c)))/a/a/pow((b+d),3)*pow(efunc(a,b,c,d),3); break;
	     }

         default:
           cout << "I don't know what to do";break;
         } // end switch
         break;

     default:
       cout << "I don't know what to do hes"; break;

   } // end switch
   return value;
}
