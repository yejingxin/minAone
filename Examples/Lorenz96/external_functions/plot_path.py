import numpy as np
import matplotlib.pyplot as plt

IC = 1
num_states = 10
num_meas = 5
num_con = 0
skip = 100
dt = 0.01

tmp = np.loadtxt('results/D%s_M%s_IC%d.dat' % (num_states, num_meas, IC))
data = np.reshape(tmp[-1,3:-1], ((-1, num_states + num_con)), order = 'c')
true = np.loadtxt('observations/true_path.txt')

N  = len(data[:,0])
time = np.arange(0, N*dt, dt)
fig = plt.figure()
fig.set_size_inches(12,12)

for idx in range(num_states + num_con):
  plt.subplot(np.floor((num_states + num_con)/2)+1, 2, idx+1)
  plt.plot(time, data[:,idx],label = "%s, est" % idx)
  if idx < num_states: plt.plot(time, true[skip:len(data[:,idx])+skip,idx], label = "%s, true" % idx)
  plt.xlim(0,N*dt)
  plt.legend(fontsize=10)


plt.savefig('results/estimation.png',bbox_inches='tight')
plt.show()
