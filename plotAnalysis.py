import numpy as np
import matplotlib.pyplot as plt

# read data file
t,rmsd = np.loadtxt('rmsd_T900_LGA.xvg',skiprows=13,unpack=True)
t,msd = np.loadtxt('msd_T900_LGA.xvg',skiprows=15,unpack=True)
r,rdf = np.loadtxt('rdf_T900_LGA.xvg',skiprows=14,unpack=True)
t,rg,rgx,rgy,rgz = np.loadtxt('gyrate_T900_LGA.xvg',skiprows=22,unpack=True)

# make time axis
indx = len(t)
t = np.linspace(0,0.005*(indx-1),len(t))

# plot root mean square deviation (RMSD) vs time
plt.figure()
plt.plot(t,rmsd,label='System RMSD')
plt.xlabel('Time (ns)')
plt.ylabel('RMSD (nm)')
plt.legend(loc=2)

# plot mean square deviation (MSD) vs time
plt.figure()
plt.plot(t,msd,label='System MSD')
plt.xlabel('Time (ns)')
plt.ylabel('MSD (nm$^2$)')
plt.legend(loc=2)

# plot radial distribution function (Gr) vs distance
plt.figure()
plt.plot(r,rdf,label='G(r) Mol-COM')
plt.plot([0,max(r)],[1,1],'k--')
plt.xlabel('r (nm)')
plt.ylabel('G(r)')
plt.legend(loc=2)

# plot radius of gyration vs time
plt.figure()
plt.plot(t,rg,'k',label='Rg')
plt.plot(t,rgx,'r',label='RgX')
plt.plot(t,rgy,'g',label='RgY')
plt.plot(t,rgz,'b',label='RgZ')
plt.ylabel('Radius of Gyration (nm)')
plt.xlabel('Time (ns)')
plt.legend(loc=2)