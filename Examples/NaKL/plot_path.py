import numpy as np
import matplotlib.pyplot as plt

IC = 1
num_states = 4
num_meas = 1
skip = 0
dt = 0.02
names = ["V", "n", "h", "m"]

tmp = np.loadtxt('results/D%s_M%s_IC%d.dat' % (num_states, num_meas, IC))
data = np.reshape(tmp[-1,3:-19], ((-1, num_states)), order = 'c')
true = np.loadtxt('input_data/allstates.dat')

N  = len(data[:,0])
time = np.arange(0, N*dt, dt)
fig = plt.figure()
fig.set_size_inches(8,12)

for idx in range(num_states):
  plt.subplot(num_states, 1, idx+1)
  plt.plot(time, data[:,idx],label = "%s, est" % names[idx])
  plt.plot(time, true[skip:len(data[:,idx])+skip,idx], label = "%s, true" % names[idx])
  plt.xlim(0,N*dt)
  plt.legend(fontsize=10)

plt.savefig('results/estimation.png',bbox_inches='tight')
plt.show()
