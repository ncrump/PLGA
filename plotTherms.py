import numpy as np
import matplotlib.pyplot as plt

# read data file
t,T,P,V,D = np.loadtxt('therms0_EA_testnpt.xvg',skiprows=22,unpack=True)

# plot running averages (0=no, 1=yes)
avg = 1

# make time axis
indx = len(t)
t = t/1000    # time in ns

# plot temperature vs time
plt.figure()
plt.subplot(2,1,1)
plt.plot(t,T,label='Temperature')
plt.ylabel('Temperature (K)')
plt.ticklabel_format(style='plain',useOffset=False,axis='both')
plt.legend(loc=1)
# plot running average
if avg == 1:
    xx = t[1::]
    yy = []
    aveRun = T[0]
    for i in range(0,indx-1):
        aveRun = aveRun + (T[i+1]-aveRun)/(i+1)
        yy.append(aveRun)
    plt.plot(xx,yy,'r')
    # compute standard deviation on last half of run
    stdDev = np.std(T[int(0.5*indx)::])
    print 'Average T =',aveRun, 'StdDev P =',stdDev
    plt.annotate('Avg = %1.0f' % aveRun, fontsize=13, xy=(0.14,0.86), color='r',xycoords='figure fraction')
    plt.annotate('Std = %1.2f' % stdDev, fontsize=13, xy=(0.14,0.55), color='r',xycoords='figure fraction')

# plot pressure vs time
plt.subplot(2,1,2)
plt.plot(t,P,label='Pressure')
plt.xlabel('Time (ns)')
plt.ylabel('Pressure (bar)')
plt.ticklabel_format(style='plain',useOffset=False,axis='both')
plt.legend(loc=1)
# plot running average
if avg == 1:
    xx = t[1::]
    yy = []
    aveRun = P[0]
    for i in range(0,indx-1):
        aveRun = aveRun + (P[i+1]-aveRun)/(i+1)
        yy.append(aveRun)
    plt.plot(xx,yy,'r')
    # compute standard deviation on last half of run
    stdDev = np.std(P[int(0.5*indx)::])
    print 'Average P =',aveRun, 'StdDev P =',stdDev
    plt.annotate('Avg = %1.0f' % aveRun, fontsize=13, xy=(0.14,0.42), color='r',xycoords='figure fraction')
    plt.annotate('Std = %1.2f' % stdDev, fontsize=13, xy=(0.14,0.12), color='r',xycoords='figure fraction')

# plot volumne vs time
plt.figure()
plt.subplot(2,1,1)
plt.plot(t,V,label='Volume')
plt.ylabel('Volume (nm$^3$)')
plt.ticklabel_format(style='plain',useOffset=False,axis='both')
plt.legend(loc=1)
# plot running average
if avg == 1:
    xx = t[1::]
    yy = []
    aveRun = V[0]
    for i in range(0,indx-1):
        aveRun = aveRun + (V[i+1]-aveRun)/(i+1)
        yy.append(aveRun)
    plt.plot(xx,yy,'r')
    # compute standard deviation on last half of run
    stdDev = np.std(V[int(0.5*indx)::])
    print 'Average V =',aveRun, 'StdDev V =',stdDev
    plt.annotate('Avg = %1.0f' % aveRun, fontsize=13, xy=(0.14,0.86), color='r',xycoords='figure fraction')
    plt.annotate('Std = %1.2f' % stdDev, fontsize=13, xy=(0.14,0.55), color='r',xycoords='figure fraction')


# plot density vs time
plt.subplot(2,1,2)
plt.plot(t,D,label='Density')
plt.xlabel('Time (ns)')
plt.ylabel('Density (kg/m$^3$)')
plt.ticklabel_format(style='plain',useOffset=False,axis='both')
plt.legend(loc=1)
# plot running average
if avg == 1:
    xx = t[1::]
    yy = []
    aveRun = D[0]
    for i in range(0,indx-1):
        aveRun = aveRun + (D[i+1]-aveRun)/(i+1)
        yy.append(aveRun)
    plt.plot(xx,yy,'r')
    # compute standard deviation on last half of run
    stdDev = np.std(D[int(0.5*indx)::])
    print 'Average D =',aveRun, 'StdDev D =',stdDev
    plt.annotate('Avg = %1.0f' % aveRun, fontsize=13, xy=(0.14,0.42), color='r',xycoords='figure fraction')
    plt.annotate('Std = %1.2f' % stdDev, fontsize=13, xy=(0.14,0.12), color='r',xycoords='figure fraction')
