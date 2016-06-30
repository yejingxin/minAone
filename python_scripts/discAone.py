#####################################################
#
#  5 January 2011
#  Bryan A. Toth
#  University of California, San Diego
#  btoth@physics.ucsd.edu
# 
#  This script performs symbolic Hermite-Simpson 
#  integration on a vector field given in the text file 
#  equations.txt, takes the Jacobian and Hessian of the
#  vector field, and stores the results in arrays that 
#  are used by other python scripts in this directory.
#
#  This script has been developed as part of a suite of 
#  python scripts to define a dynamic parameter estimation
#  problem using the optimization software IPOPT, but is 
#  generally applicable to any application needing
#  discretized derivatives of a vector field.
#
#  
#  20 May 2016
#  Daniel Breen
#  University of California, San Diego
#  dlbreen@physics.ucsd.edu
#  Fixed a bug in subfunc, which creates the Jacobian and 
#  Hessian terms. Now properly substitutes the values of 
#  n and m as defined in the examples in myfunctions.cpp.
#  
#
######################################################

import sympy as sym 
from sympy import *
import re

#  Opening and reading the text file with the vector field information
print("Importing equations.txt\n")
file = open('equations.txt','r')
temp=[]  # Array to hold equations.txt information
for line in file:
  if line.startswith('#'): # Pound used as comment in text file
     continue
  elif line.startswith('\\'): # In case file has UTF-8 markers
     continue
  else:
     temp.append(line)
file.close()

h=[]  # Array to hold unformatted equations.txt information
for i in range(len(temp)):
  temp1=temp[i].rstrip( )
  h.append(temp1)

# Initialize problem variables
nY=0
nP=0
nU=0
nI=0
nF=0
nM=0

# Problem name
Problem = h[0]

# Problem variables
a=h[1].split(',')  
nY=int(a[0])
nP=int(a[1])
nU=int(a[2])
nI=int(a[3])
if len(a) > 4:
  nF=int(a[4])
  nM=int(a[5])
else:
  nM=nU

# Import equations as strings
Feqnstr = []
for k in range(nY):
   Feqnstr.append(h[k+2])

# Import objective function as string
Fobjstr = []
Fobjstr.append(h[nY+2])

# Import variable, parameter, control, data, and stimuli names as strings
Lvars = []
for k in range(nY):
   Lvars.append(h[k+3+nY])

Lcouple = []
for k in range(nU):
    Lcouple.append(h[k+3+nY+nY])

Lparams = []
for k in range(nP):
   Lparams.append(h[k+3+nY+nY+nU])

Ldata = []
for k in range(nM):
    Ldata.append(h[k+3+nY+nY+nU+nP])

Lstimuli = []
for k in range(nI):
    Lstimuli.append(h[k+3+2*nY+nU+nP+nM])

# Import function names as strings
Funcstr = []
Funcarg = []
for k in range(nF):
    temp = h[k+3+2*nY+nP+nU+nM+nI].split(',')
    Funcstr.append(temp[0])
    Funcarg.append(int(temp[1]))
#Lvars.reverse()
#Lcouple.reverse()
Fdim = len(Feqnstr)
Pdim = len(Lparams)
print("Making symbols\n")
#  Make symbols using sympy module
#  Make symbols for next time and mid time
#    for all but parameters
Sv = []
Sp = []
Sk = []
Sd = []
Si = []
Svp1 = []
Sp = []
Skp1 = []
Sdp1 = []
Sip1 = []
Sv2 = []
Sp = []
Sk2 = []
Sd2 = []
Si2 = []
Smod = []  # Smod gets the symbols in the correct order for the Hessian
for i in range(len(Lvars)):
  Sv.append(sym.Symbol(Lvars[i]))
  Svp1.append(sym.Symbol(Lvars[i]+'p1'))
  Sv2.append(sym.Symbol(Lvars[i]+'mid'))
for i in range(len(Lparams)):
  Sp.append(sym.Symbol(Lparams[i]))
for i in range(nM):
  Sd.append(sym.Symbol(Ldata[i]))
  Sdp1.append(sym.Symbol(Ldata[i]+'p1'))
  Sd2.append(sym.Symbol(Ldata[i]+'mid'))
for i in range(nU):
  Sk.append(sym.Symbol(Lcouple[i]))
  Skp1.append(sym.Symbol(Lcouple[i]+'p1'))
  Sk2.append(sym.Symbol(Lcouple[i]+'mid'))
for i in range(nI):
  Si.append(sym.Symbol(Lstimuli[i]))
  Sip1.append(sym.Symbol(Lstimuli[i]+'p1'))
  Si2.append(sym.Symbol(Lstimuli[i]+'mid'))

# Assemble Smod order x,xp1,y,yp1, ... , xnp1, k1, k1p1, k2 ..., dkap1, p1, p2, ...
for i in range(len(Lvars)):
  Smod.append(sym.Symbol(Lvars[i]))
  Smod.append(sym.Symbol(Lvars[i]+'p1'))
for i in range(nU):
  Smod.append(sym.Symbol(Lcouple[i]))
  Smod.append(sym.Symbol(Lcouple[i]+'p1'))
for i in range(len(Lvars)):
  Smod.append(sym.Symbol(Lvars[i]+'mid'))
for i in range(nU):
  Smod.append(sym.Symbol(Lcouple[i]+'mid'))
for i in range(len(Lparams)):
  Smod.append(sym.Symbol(Lparams[i]))
Sall = Sv + Sk + Svp1 + Skp1 + Sv2 +  Sk2 + Sp

# Make symbols for functions
Sf = []
for i in range(nF):
  Sf.append(sym.Function(Funcstr[i]))

hstep = sym.Symbol("hstep")

print("Defining symbolic objective function\n")
# Define symbolic objective function
Fobj = []
sTemp1 = Fobjstr[0]
sTemp1a = Fobjstr[0]
for i in range(len(Lvars)):
  sTemp2 = "Sv[%d]" % i
  sTemp1 = sTemp1.replace(Lvars[i],sTemp2)
  sTemp2a = "Sv2[%d]" % i
  sTemp1a = sTemp1a.replace(Lvars[i],sTemp2a)
for i in range(len(Lparams)):
  sTemp2 = "Sp[%d]" % i
  sTemp1 = sTemp1.replace(Lparams[i],sTemp2)
  sTemp1a = sTemp1a.replace(Lparams[i],sTemp2)
for i in range(len(Lcouple)):
  sTemp2 = "Sk[%d]" % i
  sTemp1 = sTemp1.replace(Lcouple[i],sTemp2)
  sTemp2a = "Sk2[%d]" % i
  sTemp1a = sTemp1a.replace(Lcouple[i],sTemp2a)
for i in range(len(Lstimuli)):
  sTemp2 = "Si[%d]" % i
  sTemp1 = sTemp1.replace(Lstimuli[i],sTemp2)
  sTemp2a = "Si2[%d]" % i
  sTemp1a = sTemp1a.replace(Lstimuli[i],sTemp2a)
for i in range(nM):
  sTemp2 = "Sd[%d]" % i
  sTemp1 = sTemp1.replace(Ldata[i],sTemp2)
  sTemp2a = "Sd2[%d]" % i
  sTemp1a = sTemp1a.replace(Ldata[i],sTemp2a)
sTemp2 = "Fobj.append("
sTemp2 = sTemp2 + sTemp1 + " + " + sTemp1a + ")"
exec(sTemp2)

# Define symbolic vector field
# This vector field is for the continuous case
# Add it as a second element of objective function
# i.e. the other summation ...
Ftemp = []
for k in range(Fdim):
  sTemp1 = Feqnstr[k]
  sTemp1a = Feqnstr[k]
  sTemp1b = Feqnstr[k]
  for i in range(len(Lvars)):
    sTemp2 = "Sv[%d]" % i
    sTemp1 = sTemp1.replace(Lvars[i],sTemp2)
    sTemp2a = "Svp1[%d]" % i
    sTemp1a = sTemp1a.replace(Lvars[i],sTemp2a)
    sTemp2b = "Sv2[%d]" % i
    sTemp1b = sTemp1b.replace(Lvars[i],sTemp2b)
  for i in range(len(Lparams)):
    sTemp2 = "Sp[%d]" % i
    sTemp1 = sTemp1.replace(Lparams[i],sTemp2)
    sTemp1a = sTemp1a.replace(Lparams[i],sTemp2)
    sTemp1b = sTemp1b.replace(Lparams[i],sTemp2)
  for i in range(len(Lcouple)):
    sTemp2 = "Sk[%d]" % i
    sTemp1 = sTemp1.replace(Lcouple[i],sTemp2)
    sTemp2a = "Skp1[%d]" % i
    sTemp1a = sTemp1a.replace(Lcouple[i],sTemp2a)
    sTemp2b = "Sk2[%d]" % i
    sTemp1b = sTemp1b.replace(Lcouple[i],sTemp2b)
  for i in range(nM):
    sTemp2 = "Sd[%d]" % i
    sTemp1 = sTemp1.replace(Ldata[i],sTemp2)
    sTemp2a = "Sdp1[%d]" % i
    sTemp1a = sTemp1a.replace(Ldata[i],sTemp2a)
    sTemp2b = "Sd2[%d]" % i
    sTemp1b = sTemp1b.replace(Ldata[i],sTemp2b)
  for i in range(nI):
    sTemp2 = "Si[%d]" % i
    sTemp1 = sTemp1.replace(Lstimuli[i],sTemp2)
    sTemp2a = "Sip1[%d]" % i
    sTemp1a = sTemp1a.replace(Lstimuli[i],sTemp2a)
    sTemp2b = "Si2[%d]" % i
    sTemp1b = sTemp1b.replace(Lstimuli[i],sTemp2b)
  for i in range(nF):
    sTemp2 = "Sf[%d]" % i
    sTemp1 = sTemp1.replace(Funcstr[i],sTemp2)
    sTemp1a = sTemp1a.replace(Funcstr[i],sTemp2)
    sTemp1b = sTemp1b.replace(Funcstr[i],sTemp2)
  sTemp2 = "Ftemp.append((("
  sTemp2 = sTemp2 + sTemp1 + " + " + sTemp1a + " + 4*(" + sTemp1b + "))*hstep/6.0 + Sv[%d]- Svp1[%d])**2)" % (k,k)
  exec(sTemp2)
  sTemp2 = "Ftemp.append((("
  sTemp2 = sTemp2 + sTemp1 + " - (" + sTemp1a + "))*hstep/8.0 + 0.5*Sv[%d] + 0.5*Svp1[%d] - Sv2[%d])**2)" % (k,k,k)
  exec(sTemp2)

# At this point, Feqns is a list of each discretized equation, with Simpson
#  followed by Hermite constraint, for each state variable
# Need to add in the constrains for the coupling midpoint?  That seems unnecessary.
# Combine the terms of Feqns, since it is a summation
# DO NOT COMBINE ... these terms are a summation, but will have different R-factors in front.
#   Combine in generated code.
Feqns = []
for i in range(nY):
  temp1 = Ftemp[2*i] + Ftemp[2*i+1]
  Feqns.append(temp1) 

dict1 = {0:"Xval",1:"Xvalp1",2:"Xval2"}
dict2 = {0:"K11val",1:"K11valp1",2:"K11val2"}
dict3 = {0:"Xdval",1:"Xdvalp1",2:"Xdval2"}
dict4 = {0:"Ival",1:"Ivalp1",2:"Ival2"}

# SUBVARS function turns the output into c++ code.  Modified from discretize.py
#  to account for plus-one and midpoints

def subvars(mystr,myi):
   mytemp = mystr
#  The following two lines are very important
#  Sympy converts all inputs into a simplified form
#  for calculations.  Specifically, this involves
#  any exponentials (a^b) put in the form a**b.
#  Whereas this form is acceptable for fortran outputs,
#  this needs to change to pow(a,b) for C++ outputs.
#  Sympify and ccode combine to make this transformation.
#  This transformation is done at this point in the code,
#  since sympify will not operate on an
#  expression that includes brackets - which are added
#  in this function.
   mytemp = sym.sympify(mytemp)
   mytemp = sym.ccode(mytemp)
    
   n = myi
   while (n > -1):
     for j in range(len(Sv)):
       if (n == 0):
         Srep = dict1[n] + "[%d]" % j
         Sfind = Lvars[j]
         mytemp = mytemp.replace(Sfind,Srep)
       elif (n == 1):
         Srep = dict1[n] + "[%d]" % j
         Sfind = Lvars[j] + "p1"
         mytemp = mytemp.replace(Sfind,Srep)
       elif (n == 2):
         Srep = dict1[n] + "[%d]" % j
         Sfind = Lvars[j] + "mid"
         mytemp = mytemp.replace(Sfind,Srep)

     for j in range(len(Sp)):
       Srep = "Pval[%d]" % j
       Sfind = Lparams[j]
       mytemp = mytemp.replace(Sfind,Srep)

     for j in range(len(Sk)):
         if (n == 0):
           Srep = dict2[n] + "[%d]" % (j)
           Sfind = Lcouple[j]
           mytemp = mytemp.replace(Sfind,Srep)
         elif (n == 1):
           Srep = dict2[n] + "[%d]" % (j)
           Sfind = Lcouple[j] + "p1"
           mytemp = mytemp.replace(Sfind,Srep)
         elif (n == 2):
           Srep = dict2[n] + "[%d]" % (j)
           Sfind = Lcouple[j] + "mid"
           mytemp = mytemp.replace(Sfind,Srep)
 
     for j in range(len(Sd)):
       if (n == 0):
         Srep = dict3[n] + "[%d]" % j
         Sfind = Ldata[j]
         mytemp = mytemp.replace(Sfind,Srep)
       elif (n == 1):
         Srep = dict3[n] + "[%d]" % j
         Sfind = Ldata[j] + "p1"
         mytemp = mytemp.replace(Sfind,Srep)
       elif (n == 2):
         Srep = dict3[n] + "[%d]" % j
         Sfind = Ldata[j] + "mid"
         mytemp = mytemp.replace(Sfind,Srep)

     for j in range(len(Si)):
       if (n == 0):
         Srep = dict4[n] + "[%d]" % j
         Sfind = Lstimuli[j]
         mytemp = mytemp.replace(Sfind,Srep)
       elif (n == 1):
         Srep = dict4[n] + "[%d]" % j
         Sfind = Lstimuli[j] + "p1"
         mytemp = mytemp.replace(Sfind,Srep)
       elif (n == 2):
         Srep = dict4[n] + "[%d]" % j
         Sfind = Lstimuli[j] + "mid"
         mytemp = mytemp.replace(Sfind,Srep)
     n = n - 1
   return mytemp
# END subvars

def subfunc(mystr,myi):
   mytemp = mystr
   Dsearch = re.findall('Derivative\(', mytemp)
   for num in range(len(Dsearch)):
     jacsearch = re.search('Derivative\(([A-Za-z]+)\(([A-Za-z0-9]+\[[0-9]+\])+(, [A-Za-z0-9]+\[[0-9]+\])*\), (-?[0-9]*(\.[0-9]+)?|[A-Za-z0-9]+\[[0-9]+\])\)', mytemp)
     hessearch = re.search('Derivative\(([A-Za-z]+)\(([A-Za-z0-9]+\[[0-9]+\])+(, [A-Za-z0-9]+\[[0-9]+\])*\), (-?[0-9]*(\.[0-9]+)?|[A-Za-z0-9]+\[[0-9]+\]), (-?[0-9]*(\.[0-9]+)?|[A-Za-z0-9]+\[[0-9]+\])\)', mytemp)
     if jacsearch:
          
       dervar = jacsearch.group(4)
       var = re.findall('[A-za-z0-9]+\[[0-9]+\]|-?[0-9]*\.?[0-9]+', jacsearch.group(0))
       for i in range(len(var)-1):
         if var[i] == dervar:
           
           jacnum = i+1
       rep = jacsearch.group(1) + 'jac('
       for i in range(len(var)-1):
         temp = var[i] + ','
         rep += temp
       rep += str(jacnum) +')'
       mytemp = re.sub('Derivative\(([A-Za-z]+)\(([A-Za-z0-9]+\[[0-9]+\])+(, [A-Za-z0-9]+\[[0-9]+\])*\), (-?[0-9]+(\.[0-9]+)?|[A-Za-z0-9]+\[[0-9]+\])\)', rep,mytemp,1)
     if hessearch:
       dervar1 = hessearch.group(4)
       dervar2 = hessearch.group(6)
       hesnum1 = 0
       hesnum2 = 0
       var = re.findall('[A-za-z0-9]+\[[0-9]+\]|-?[0-9]*\.?[0-9]+', hessearch.group(0))
       for i in range(len(var)-2):
         if var[i] == dervar1:
           hesnum1 = i+1
         if var[i] == dervar2:
           hesnum2 = i+1
       rep = hessearch.group(1) + 'hes('
       
       for i in range(len(var)-2):
         temp = var[i] + ','
         rep += temp
       rep += str(hesnum1) +','+ str(hesnum2)+ ')'
       mytemp = re.sub('Derivative\(([A-Za-z]+)\(([A-Za-z0-9]+\[[0-9]+\])+(, [A-Za-z0-9]+\[[0-9]+\])*\), (-?[0-9]+(\.[0-9]+)?|[A-Za-z0-9]+\[[0-9]+\]), (-?[0-9]+(\.[0-9]+)?|[A-Za-z0-9]+\[[0-9]+\])\)', rep, mytemp, 1)
   return mytemp
#end subfunc
print("Building objective function strings\n")
# Build objective function strings
#  In format (data,model) where model is an array
strObj = []
Stemp = str(Fobj[0])
Stemp = subvars(Stemp,2)
Stemp = subfunc(Stemp,2)
strObj.append(Stemp) 
temp2 = []
for j in range(len(Feqns)):
  Stemp = str(Feqns[j])
  Stemp = subvars(Stemp,2)
  Stemp = subfunc(Stemp,2)
  temp2.append(Stemp)
strObj.append(temp2)

print("Building objective gradient strings\n")
# Build objective gradient strings
strGrad = []
temp1 = []
for jvar in range(len(Sall)):
   Stemp = str(sym.diff(Fobj[0],Sall[jvar]))
   Stemp = subvars(Stemp,2)
   Stemp = subfunc(Stemp,2)
   temp1.append(Stemp)
strGrad.append(temp1)

for icon in range(len(Feqns)):
   temp1 = []
   for jvar in range(len(Sall)):
      Stemp = str(sym.diff(Feqns[icon],Sall[jvar]))
      Stemp = subvars(Stemp,2)
      Stemp = subfunc(Stemp,2)
      temp1.append(Stemp)
   strGrad.append(temp1)

#def sub_user_defined_functions(mystr, func, var):
     

print("Building Hessian strings\n")
# Build Hessian strings
strHes = []
temp1=[]
for jvar in range(len(Smod)):
    temp2 = []
    for kvar in range(jvar+1):
    #for kvar in range(len(Smod)):
        Stemp = str(sym.diff(sym.diff(Fobj[0],Smod[jvar]),Smod[kvar]))
        Stemp = subvars(Stemp,2)
        Stemp = subfunc(Stemp,2)
        temp2.append(Stemp)
    temp1.append(temp2)
strHes.append(temp1)

# Build Hessian strings of constraint
for icon in range(len(Feqns)):
   temp1 = []
   for jvar in range(len(Smod)):
      temp2 = []
      for kvar in range(jvar+1):
      #for kvar in range(len(Smod)):
         Stemp = str(sym.diff(sym.diff(Feqns[icon],Smod[jvar]),Smod[kvar]))
         Stemp = subvars(Stemp,2)
         Stemp = subfunc(Stemp,2)
         temp2.append(Stemp)
      temp1.append(temp2)
   strHes.append(temp1)
# Derivatives taken above.
# Now, remove all non-zero entries
# Fill out objective gradient
# Do data part first
print("Filling out objective gradient non-zero elements\n")
VObj = []
temp2 = []
for j in range(len(Sall)):
   F = strGrad[0][j]
   if F != '0' :
       temp1 = []
       temp1.append(F)
       temp1.append(0)       
       temp1.append(j)
       temp2.append(temp1)
VObj.append(temp2)
# Now do model part
temp2 = []
for i in range(len(strGrad)-1):
 for j in range(len(Sall)):
   F = strGrad[i+1][j]
   if F != '0' :
       temp1 = []
       temp1.append(F)
       temp1.append(i+1)
       temp1.append(j)
       temp2.append(temp1)
VObj.append(temp2)
# Nothing needed for VJac, since it does not exist
# Fill out Hessian vector, keeping track of the index counter
#  for specific row/column combinations
print("Filling out Hessian non-zero elements\n")
from numpy import *
# Set up indexing for VHes
# Hesindex will keep track of which elements are used
# Hessize defines the number of elements in each entry
#  's' = 1, 'm' = N, 'l' = N+1
Hesindex = zeros([len(Sall),len(Sall)], integer)
Hessize = zeros([len(Sall),len(Sall)], character)
for i in range(len(Sall)):
   for j in range(i+1):
       if i > (3*(nY+nU)-1):
          if j > (3*(nY+nU)-1):
            # Parameter/parameter derivatives only have 1 element 
             Hessize[i][j] = 's'
       # Two instances where derivative have N+1 elements
       # First instance is x/par derivatives
          elif j < (2*(nY+nU)):
             if j % 2 == 0:
                Hessize[i][j] = 'l'
             else:
                Hessize[i][j] = 'm'
          else: 
             Hessize[i][j] = 'm'
       # Second instance is x/x, x/y derivatives since these
       #  include xp1/xp1 and xp1/yp1 etc ... terms
       elif i < (2*(nY+nU)):
          if i % 2 == 0:
             if j % 2 == 0:
                Hessize[i][j] = 'l'
             else:
                Hessize[i][j] = 'm'
          else: 
             Hessize[i][j] = 'm'
       # All other instances have N elements
       else:
         Hessize[i][j] = 'm'
VHes = []

# Symmetrical matrix - filling out lower half diagonal only.
# Fill out constraint for data part first
for i in range(nY+1):
  for j in range(len(Smod)):
   for k in range(j+1):
     H = strHes[i][j][k]
     if j < (3*(nY+nU)):
         n = 0
         if H != '0':
             temp1 = []
             Hesindex[j][k] = -1
             if j < (2*(nY+nU)):
               if j%2 != 0:
                if k%2 != 0: 
                  Hesindex[j][k]=0
                  n = 1
             temp1.append(i)
             temp1.append(j)
             temp1.append(k)
             temp1.append(n)
             temp1.append(H)
             VHes.append(temp1)
     else:  # These are the parameter/parameter derivatives
       n = 0
       if H != '0':
         temp2 = []
         Hesindex[j][k] = -1
         if k < (2*(nY+nU)):
           if k%2 != 0:
              Hesindex[j][k]=0
              n = 1
         temp2.append(i)
         temp2.append(j)
         temp2.append(k)
         temp2.append(n)
         temp2.append(H)
         VHes.append(temp2)
