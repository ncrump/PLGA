import numpy as np
import matplotlib.pyplot as plt

# read data file
t,PE,KE,E = np.loadtxt('energy_nvt0_h2o_bigbox3x.xvg',skiprows=21,unpack=True)


# make time axis
indx = len(t)
t = t/1000    # time in ns

# make plots
plt.figure()
plt.subplot(3,1,1)
plt.plot(t,PE,'r')
plt.ylabel('<PE> (kJ/mol)')
plt.subplot(3,1,2)
plt.plot(t,KE,'g')
plt.ylabel('<KE> (kJ/mol)')
plt.subplot(3,1,3)
plt.plot(t,E,'b')
plt.ylabel('<E> (kJ/mol)')
plt.xlabel('Time (ns)')