import numpy as np
data = np.loadtxt('./measured.dat')
NT = np.shape(data)[0]
print NT
data = data + np.random.randn(NT)
np.savetxt('noise_measured.dat', data, fmt='%e')
