"""
Analyzes LGA+EA trajectory to get mean velocity of molecules in bins along x-axis
"""

import numpy as np
import matplotlib.pyplot as plt

# set input parameters
# -----------------------------------------------
fgro = 'nvt0_LGA+EA_0.5Agap_eabigbox3x_analysis_reimagedTRAJ.gro'
lgaMOL = 256      # number of LGA molecules per frame
lgaATM = 21       # number of LGA atoms per molecule
eeeMOL = 1620     # number of EA molecules per frame
eeeATM = 6        # number of EA atoms per molecule
frames = 101      # number of frames in trajectory
dtstep = 0.1      # time step between frames (ns)
xbox   = 20.342   # box length in x-direction (nm)
xbin   = 20       # number of bins in x-direction
intfce = 5.748    # interface distance from origin (nm)
# -----------------------------------------------

# read trajectory file (fixed column width)
print '\nreading file ...'
delimiter = [8,7,5,8,8,8,8,8,8]
nres0,natm0,num0,pxatm0,pyatm0,pzatm0,vxatm0,vyatm0,vzatm0 = np.genfromtxt(fgro,dtype=str,delimiter=delimiter,unpack=True)

# define atoms numbers
print 'analyzing trajectory ...'
LGAatms = lgaMOL*lgaATM
EEEatms = eeeMOL*eeeATM
TOTatms = LGAatms+EEEatms

# define dictionary of atom masses for LGA
mLGA    = 278.2166
massLGA = {'C0-L': 14.0270, 'H2-L': 01.0080, 'C7-L': 12.0110, 'O8-L': 15.9994,\
           'H1-L': 01.0080, 'OB-L': 15.9994, 'OA-L': 15.9994, 'C6-L': 13.0190,\
           'O7-L': 15.9994, 'C9-L': 12.011,  'C8-L': 14.0270, 'O9-L': 15.9994,\
           'O6-L': 15.9994, 'C3-L': 13.0190, 'C2-L': 15.035,  'O5-L': 15.9994,\
           'O4-L': 15.9994, 'O3-L': 15.9994, 'C1-L': 12.0110, 'C5-L': 15.0350,\
           'C4-L': 12.0110}

# define dictionary of atom masses for EA
mEEE    = 88.1068
massEEE = {'C4-E':15.0350,'C3-E':14.0270,'OB-E':15.9994,'C2-E':12.0110,\
           'OA-E':15.9994,'C1-E':15.0350}

# format arrays to strip trajectory header/footer lines
Ntot = len(natm0)
indx = np.arange(2,Ntot,TOTatms+3)
nres,natm,pxatm,pyatm,pzatm,vxatm,vyatm,vzatm = [],[],[],[],[],[],[],[]
# get molecule lines only
for i in indx:
    nres.append(nres0[i:i+TOTatms])
    natm.append(natm0[i:i+TOTatms])
    pxatm.append(pxatm0[i:i+TOTatms])
    pyatm.append(pyatm0[i:i+TOTatms])
    pzatm.append(pzatm0[i:i+TOTatms])
    vxatm.append(vxatm0[i:i+TOTatms])
    vyatm.append(vyatm0[i:i+TOTatms])
    vzatm.append(vzatm0[i:i+TOTatms])
# tag atom types for LGA or EA
nres  = np.char.strip(nres)
natm  = np.char.strip(natm)
for i in range(frames):
    for j in range(0,LGAatms):
        natm[i][j] = natm[i][j]+'-L'
    for j in range(LGAatms,TOTatms):
        natm[i][j] = natm[i][j]+'-E'
# reformat arrays
nres  = np.array(nres).flatten().astype(str)
natm  = np.array(natm).flatten().astype(str)
pxatm = np.array(pxatm).flatten().astype(float)
pyatm = np.array(pyatm).flatten().astype(float)
pzatm = np.array(pzatm).flatten().astype(float)
vxatm = np.array(vxatm).flatten().astype(float)
vyatm = np.array(vyatm).flatten().astype(float)
vzatm = np.array(vzatm).flatten().astype(float)
t     = np.linspace(0,(frames-1)*dtstep,frames)

# calculate center mass position and velocity of each molecule per frame
pxcm,pycm,pzcm,vxcm,vycm,vzcm = [],[],[],[],[],[]
indx = 0
# loop over frames
for i in range(frames):
    # loop over LGA molecules in frame
    for j in range(lgaMOL):
        pxcom,pycom,pzcom = 0,0,0
        vxcom,vycom,vzcom = 0,0,0
        # loop over LGA atoms in molecule
        for k in range(lgaATM):
            pxcom += massLGA[natm[indx]]*pxatm[indx]
            pycom += massLGA[natm[indx]]*pyatm[indx]
            pzcom += massLGA[natm[indx]]*pzatm[indx]
            vxcom += massLGA[natm[indx]]*vxatm[indx]
            vycom += massLGA[natm[indx]]*vyatm[indx]
            vzcom += massLGA[natm[indx]]*vzatm[indx]
            indx += 1
        pxcm.append(pxcom/mLGA)
        pycm.append(pycom/mLGA)
        pzcm.append(pzcom/mLGA)
        vxcm.append(vxcom/mLGA)
        vycm.append(vycom/mLGA)
        vzcm.append(vzcom/mLGA)
    # loop over EA molecules in frame
    for j in range(eeeMOL):
        pxcom,pycom,pzcom = 0,0,0
        vxcom,vycom,vzcom = 0,0,0
        # loop over EA atoms in molecule
        for k in range(eeeATM):
            pxcom += massEEE[natm[indx]]*pxatm[indx]
            pycom += massEEE[natm[indx]]*pyatm[indx]
            pzcom += massEEE[natm[indx]]*pzatm[indx]
            vxcom += massEEE[natm[indx]]*vxatm[indx]
            vycom += massEEE[natm[indx]]*vyatm[indx]
            vzcom += massEEE[natm[indx]]*vzatm[indx]
            indx += 1
        pxcm.append(pxcom/mEEE)
        pycm.append(pycom/mEEE)
        pzcm.append(pzcom/mEEE)
        vxcm.append(vxcom/mEEE)
        vycm.append(vycom/mEEE)
        vzcm.append(vzcom/mEEE)

# calculate mean velocity in each bin per frame
aveVx,aveVy,aveVz = [],[],[]
vx,vy,vz = np.zeros(xbin),np.zeros(xbin),np.zeros(xbin)
nx,ny,nz = np.zeros(xbin),np.zeros(xbin),np.zeros(xbin)
bins = np.linspace(0,xbox,xbin+1)
dx   = xbox/xbin
Nmol = lgaMOL+eeeMOL
indx = 0
# loop over frames
for i in range(frames):
    vx,vy,vz = vx*0,vy*0,vz*0
    nx,ny,nz = nx*0,ny*0,nz*0
    # loop over molecules in frame
    for j in range(Nmol):
        binNum = int(pxcm[indx]/dx)
        vx[binNum] += vxcm[indx]
        vy[binNum] += vycm[indx]
        vz[binNum] += vzcm[indx]
        nx[binNum] += 1
        ny[binNum] += 1
        nz[binNum] += 1
        indx += 1
    aveVx.append(vx/nx)
    aveVy.append(vy/ny)
    aveVz.append(vz/nz)

# make plots
print 'making plots ...'
plt.ioff()
cnt = 0
fontsize = 18
factor   = 1.20
vxmin = np.min(np.array(aveVx).flatten())*factor
vymin = np.min(np.array(aveVy).flatten())*factor
vzmin = np.min(np.array(aveVz).flatten())*factor
vxmax = np.max(np.array(aveVx).flatten())*factor
vymax = np.max(np.array(aveVy).flatten())*factor
vzmax = np.max(np.array(aveVz).flatten())*factor
for i in range(frames):
    plt.figure()
    plt.subplot(311)
    plt.plot(bins[:-1]+0.5*dx,aveVx[i],'ro-')
    plt.xlim(0,xbox)
    plt.ylim(vxmin,vxmax)
    plt.xticks([])
    plt.ylabel(r'$ \overline{Vx}$',fontsize=fontsize)
    plt.hlines(0,0,xbox,'k',linestyle='dashed')
    plt.vlines(bins,vxmin,vxmax)
    plt.vlines(intfce,vxmin,vxmax,'m',lw=2,linestyle='dashed')
    plt.subplot(312)
    plt.plot(bins[:-1]+0.5*dx,aveVy[i],'go-')
    plt.xlim(0,xbox)
    plt.ylim(vymin,vymax)
    plt.xticks([])
    plt.ylabel(r'$ \overline{Vy}$',fontsize=fontsize)
    plt.hlines(0,0,xbox,'k',linestyle='dashed')
    plt.vlines(bins,vymin,vymax)
    plt.vlines(intfce,vymin,vymax,'m',lw=2,linestyle='dashed')
    plt.subplot(313)
    plt.plot(bins[:-1]+0.5*dx,aveVz[i],'bo-')
    plt.xlim(0,xbox)
    plt.ylim(vzmin,vzmax)
    plt.ylabel(r'$ \overline{Vz}$',fontsize=fontsize)
    plt.xlabel('$x-bins$',fontsize=fontsize)
    plt.hlines(0,0,xbox,'k',linestyle='dashed')
    plt.vlines(bins,vzmin,vzmax)
    plt.vlines(intfce,vzmin,vzmax,'m',lw=2,linestyle='dashed')
    plt.annotate('f='+'%3i'%i,fontsize=fontsize,xy=(0.64,0.92),xycoords='figure fraction')
    plt.annotate('t='+'%4.1f ns'%t[i],fontsize=fontsize,xy=(0.76,0.92),xycoords='figure fraction')
    plt.annotate('$(nm/ps)$',fontsize=16,xy=(0.05,0.93),xycoords='figure fraction')
    plt.annotate('$(nm)$',fontsize=16,xy=(0.85,0.025),xycoords='figure fraction')
    plt.savefig('Images/'+str(cnt)+'velocitybins.png')
    plt.close()
    cnt += 1

print cnt,'images saved'
print 'done !\n'
plt.ion()