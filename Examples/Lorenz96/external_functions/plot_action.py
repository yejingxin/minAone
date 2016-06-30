import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

num_states = 10
num_meas = 5
init_conds = range(1,100)
min_val = 1e8
min_path = 1
num_steps = 30

for iy in init_conds:
  print (iy, end = " ", flush=True)
  tmp = np.loadtxt('results/D%s_M%s_IC%d.dat' % (num_states, num_meas, iy))
  plt.scatter(np.arange(num_steps), np.log(tmp[:,2])/np.log(10))
  if tmp[-1,2] < min_val:
    min_val = tmp[-1,2]
    min_path = iy
    min_param = tmp[-1,-1]

print ("\n\nMinimum path is IC = %d; parameter value = %f" %(min_path, min_param))
plt.xlabel(r'$\beta = \log_{\alpha}[R_f/R_{f0}]$',fontsize=15)
plt.ylabel(r'$\log[A]$',fontsize=20)
plt.savefig('results/action.png')
plt.show()

