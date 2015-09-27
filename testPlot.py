"""
Created on Mon Nov 24 21:19:12 2014
"""

import numpy as np
import matplotlib.pyplot as plt

Temp = np.array([300, 400, 500, 600, 700, 800, 900])
Dens = np.array([1352.0, 1284.3, 1209.8, 1134.4, 1054.0, 966.8, 865.3])


# plot diff coeff vs temp & dens
D = np.array([9.00000000e-09,4.35000000e-07,3.80800000e-06,1.14930000e-05,2.43020000e-05,4.53830000e-05,8.82580000e-05])
# Quadratic fit to D vs Temp and D vs Dens
f1 = np.polyfit(Temp,D,2)
fit1 = f1[0]*Temp**2 + f1[1]*Temp + f1[2]
f2 = np.polyfit(Dens,D,2)
fit2 = f2[0]*Dens**2 + f2[1]*Dens + f2[2]
print f1
print f2
# plots
plt.figure()
plt.plot(Temp,D,'bo-',label='Data')
plt.plot(Temp,fit1,'ko--',label='Fit')
plt.xlabel('Temperature (K)')
plt.ylabel('Diffusion Coefficient (cm$^2$/s)')
plt.legend(loc=2)
plt.figure()
plt.plot(Dens,D,'ro-',label='Data')
plt.plot(Dens,fit2,'ko--',label='Fit')
plt.gca().invert_xaxis()
plt.xlabel('Density (kg/m$^3$)')
plt.ylabel('Diffusion Coefficient (cm$^2$/s)')
plt.legend(loc=2)
#
#
## plot gyrate vs temp
#t3,rg3,rgx3,rgy3,rgz3 = np.loadtxt('gyrate_T300_LGA.xvg',skiprows=22,unpack=True)
#t4,rg4,rgx4,rgy4,rgz4 = np.loadtxt('gyrate_T400_LGA.xvg',skiprows=22,unpack=True)
#t5,rg5,rgx5,rgy5,rgz5 = np.loadtxt('gyrate_T500_LGA.xvg',skiprows=22,unpack=True)
#t6,rg6,rgx6,rgy6,rgz6 = np.loadtxt('gyrate_T600_LGA.xvg',skiprows=22,unpack=True)
#t7,rg7,rgx7,rgy7,rgz7 = np.loadtxt('gyrate_T700_LGA.xvg',skiprows=22,unpack=True)
#t8,rg8,rgx8,rgy8,rgz8 = np.loadtxt('gyrate_T800_LGA.xvg',skiprows=22,unpack=True)
#t9,rg9,rgx9,rgy9,rgz9 = np.loadtxt('gyrate_T900_LGA.xvg',skiprows=22,unpack=True)
#fig = plt.figure()
#ax = plt.subplot(111)
#plt.plot(t3,rg3,label='300K')
#plt.plot(t3,rg4,label='400K')
#plt.plot(t3,rg5,label='500K')
#plt.plot(t3,rg6,label='600K')
#plt.plot(t3,rg7,label='700K')
#plt.plot(t3,rg8,label='800K')
#plt.plot(t3,rg9,label='900K')
#plt.ylabel('Radius of Gyration (nm)')
#plt.xlabel('Time (ns)')
#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
#plt.legend(loc=6,bbox_to_anchor=(1,0.5))
## plot average
#avgRg = [np.average(rg3),np.average(rg4),np.average(rg5),np.average(rg6),np.average(rg7),np.average(rg8),np.average(rg9)]
#plt.figure()
#plt.plot(Temp,avgRg,'bo-')
#plt.xlabel('Temperature (K)')
#plt.ylabel('Average Rg (nm)')
## print average
#print '300K', np.average(rg3)
#print '400K', np.average(rg4)
#print '500K', np.average(rg5)
#print '600K', np.average(rg6)
#print '700K', np.average(rg7)
#print '800K', np.average(rg8)
#print '900K', np.average(rg9)
#
#
## plot energy vs temp
#t3,p3 = np.loadtxt('potential_300K.xvg',skiprows=19,unpack=True)
#t4,p4 = np.loadtxt('potential_400K.xvg',skiprows=19,unpack=True)
#t5,p5 = np.loadtxt('potential_500K.xvg',skiprows=19,unpack=True)
#t6,p6 = np.loadtxt('potential_600K.xvg',skiprows=19,unpack=True)
#t7,p7 = np.loadtxt('potential_700K.xvg',skiprows=19,unpack=True)
#t8,p8 = np.loadtxt('potential_800K.xvg',skiprows=19,unpack=True)
#t9,p9 = np.loadtxt('potential_900K.xvg',skiprows=19,unpack=True)
#fig = plt.figure()
#ax = plt.subplot(111)
#plt.plot(t3,p3,label='300K')
#plt.plot(t3,p4,label='400K')
#plt.plot(t3,p5,label='500K')
#plt.plot(t3,p6,label='600K')
#plt.plot(t3,p7,label='700K')
#plt.plot(t3,p8,label='800K')
#plt.plot(t3,p9,label='900K')
#plt.ylabel('Potential (kJ/mol)')
#plt.xlabel('Time (ns)')
#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
#plt.legend(loc=6,bbox_to_anchor=(1,0.5))
## plot average and variance
#avgP = [np.average(p3),np.average(p4),np.average(p5),np.average(p6),np.average(p7),np.average(p8),np.average(p9)]
#varP = [np.var(p3),np.var(p4),np.var(p5),np.var(p6),np.var(p7),np.var(p8),np.var(p9)]
#plt.figure()
#plt.plot(Temp,avgP,'bo-')
#plt.xlabel('Temperature (K)')
#plt.ylabel('Average of Potential (kJ/mol)')
#plt.figure()
#plt.plot(Temp,varP,'ro-')
#plt.xlabel('Temperature (K)')
#plt.ylabel('Variance of Potential (kJ/mol)$^2$')
## print average and variance
#print '300K', np.average(p3), np.var(p3)
#print '400K', np.average(p4), np.var(p4)
#print '500K', np.average(p5), np.var(p5)
#print '600K', np.average(p6), np.var(p6)
#print '700K', np.average(p7), np.var(p7)
#print '800K', np.average(p8), np.var(p8)
#print '900K', np.average(p9), np.var(p9)