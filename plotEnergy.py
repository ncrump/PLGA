import numpy as np
import matplotlib.pyplot as plt

# read data file
t,GA,GP,GI,LJ14,C14,LJSR,DisCor,CSR,CR,P = np.loadtxt('energy_EA.xvg',skiprows=28,unpack=True)


# make time axis
indx = len(t)
t = t/1000    # time in ns

# plot bonded and non-bonded energy terms
fig = plt.figure()
ax = plt.subplot(111)
plt.plot(t,GA,'-.',label='G96Angle')
plt.plot(t,GP,'-.',label='Proper Dih')
plt.plot(t,GI,'-.',label='Improper Dih')
plt.plot(t,LJ14,label='LJ-14')
plt.plot(t,C14,label='Coulomb-14')
plt.plot(t,LJSR,label='LJ (SR)')
plt.plot(t,DisCor,label='Disper Corr')
plt.plot(t,CSR,label='Coulomb (SR)')
plt.plot(t,CR,label='Coulomb Recip')
plt.plot(t,P,'--',label='Potential')
plt.xlabel('Time (ns)')
plt.ylabel('Potential Energy (kJ/mol)')
plt.ticklabel_format(style='sci',useOffset=False,axis='y',scilimits=(-3,3))
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
plt.legend(loc=6,bbox_to_anchor=(1,0.5),fontsize=11)

# plot potential energy
plt.figure()
plt.plot(t,P,label='Potential')
plt.xlabel('Time (ns)')
plt.ylabel('Potential (kJ/mol)')
plt.ticklabel_format(style='sci',useOffset=False,axis='y',scilimits=(-3,3))
plt.legend(loc=2)