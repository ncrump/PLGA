import numpy as np
import matplotlib.pyplot as plt

# read data file
t,msd = np.loadtxt('msd_gromacs.xvg',skiprows=15,unpack=True)

# make time axis
indx = len(t)
t = np.linspace(0,0.1*(indx-1),len(t))

# plot mean square deviation (MSD) vs time
plt.figure()
fontsize = 16
plt.plot(t,msd,'r-',label='Gromacs System MSD')
plt.xlabel('time (ns)',fontsize=fontsize)
plt.ylabel('MSD (nm$^2$)',fontsize=fontsize)
plt.legend(loc=2)