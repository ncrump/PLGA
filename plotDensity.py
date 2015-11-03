import numpy as np
import matplotlib.pyplot as plt

# read data file
dX,densX = np.loadtxt('densityX_nvt0_h2o_bigbox3x.xvg',skiprows=19,unpack=True)
dY,densY = np.loadtxt('densityy_nvt0_h2o_bigbox3x.xvg',skiprows=19,unpack=True)
dZ,densZ = np.loadtxt('densityZ_nvt0_h2o_bigbox3x.xvg',skiprows=19,unpack=True)

# make plots
plt.figure()
plt.subplot(3,1,1)
plt.plot(dX,densX,'r')
plt.ylabel('X-Density (kg/m$^3$)')
plt.subplot(3,1,2)
plt.plot(dY,densY,'g')
plt.ylabel('Y-Density (kg/m$^3$)')
plt.subplot(3,1,3)
plt.plot(dZ,densZ,'b')
plt.ylabel('Z-Density (kg/m$^3$)')
plt.xlabel('Box Distance (nm)')

# print averages
aveX = np.average(densX)
aveY = np.average(densY)
aveZ = np.average(densZ)
print '\nAverages'
print 'Density X = ',aveX
print 'Density Y = ',aveY
print 'Density Z = ',aveZ