import numpy as np
import matplotlib.pyplot as plt

data_set = 1
IC = 0
num_states = 10
num_meas = 4
skip = 100
dt = 0.01

tmp = np.loadtxt('results/D%s_M%s_PATH%s_IC%d.dat' % (num_states, num_meas, data_set, IC))
data = np.reshape(tmp[-1,3:-1], ((-1, num_states)), order = 'c')
true = np.loadtxt('observations/true_path_%s.txt' % data_set)

N  = len(data[:,0])
time = np.arange(0, N*dt, dt)
fig = plt.figure()
fig.set_size_inches(12,12)

for idx in range(num_states):
  plt.subplot(np.floor(num_states/2)+1, 2, idx+1)
  plt.plot(time, data[:,idx],label = "%s, est" % idx)
  plt.plot(time, true[skip:len(data[:,idx])+skip,idx], label = "%s, true" % idx)
  plt.xlim(0,N*dt)
  plt.legend(fontsize=10)

plt.savefig('results/estimation_%s.png' % data_set,bbox_inches='tight')
plt.show()
