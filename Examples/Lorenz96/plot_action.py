import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use("Agg")
font = {'size'   : 18}
matplotlib.rc('font', **font)
D = 5
M = 1
PATH=100
B=30
minima = np.zeros((PATH*B,2))
for p in range(PATH):
    data = np.loadtxt('D%d_M%d_PATH%d.dat'%(D,M,p+1))
    minima[B*p:B*(p+1),:] = data[:,[0,2]]
    print p,
np.savetxt('minima0.dat', minima, fmt='%e')
fig = plt.figure(figsize=(12,9))
ax1 = fig.add_subplot(111)
ax1.plot(minima[:,0],minima[:,1],'bo')
plt.axhline(y=1, xmin=-5, xmax=35, linewidth=2, color = 'r',linestyle='--')
ax1.set_title('Lorenz96 D=5 L=1 NT=161*0.025')
ax1.set_ylabel('Minimized Action')
ax1.set_xlabel('beta')
ax1.set_yscale('log')
#ax1.set_ylim(())
plt.savefig('Action.png',bbox_inches='tight')
