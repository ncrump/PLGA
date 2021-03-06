"""
Analyzes LGA+EA trajectory to get mean velocity of molecules in bins along x-axis
"""

import numpy as np
import matplotlib.pyplot as plt

# set input parameters
# -----------------------------------------------
path = 'Data/'
fgro = 'nvt1_LGA+EA_0.5Agap_eabigbox3x_analysis_originalTRAJ.gro'
lgaMOL = 256        # number of LGA molecules per frame
lgaATM = 21         # number of LGA atoms per molecule
eeeMOL = 1620       # number of EA molecules per frame
eeeATM = 6          # number of EA atoms per molecule
frames = 101        # number of frames in trajectory
dtstep = 0.1        # time step between frames (ns)
xbox   = 20.34200   # box length in x-direction (nm)
ybox   =  5.59744   # box length in y-direction (nm)
zbox   =  3.19734   # box length in z-direction (nm)
dxbin  = 2          # width of bin centered around interface (nm)
numbin = 3          # number of bins centered around interface
intfce = 5.748      # interface distance from origin (nm)
# -----------------------------------------------

# read trajectory file (fixed column width)
print '\nreading file ...'
delimiter = [8,7,5,8,8,8,8,8,8]
nres0,natm0,num0,pxatm0,pyatm0,pzatm0,vxatm0,vyatm0,vzatm0 = np.genfromtxt(path+fgro,dtype=str,delimiter=delimiter,unpack=True)

# define atom numbers
print 'analyzing trajectory ...'
LGAatms = lgaMOL*lgaATM
EEEatms = eeeMOL*eeeATM
TOTatms = LGAatms+EEEatms

# define dictionary of atom masses for LGA (in g/mol)
mLGA    = 278.2166
massLGA = {'C0-L': 14.0270, 'H2-L': 01.0080, 'C7-L': 12.0110, 'O8-L': 15.9994,\
           'H1-L': 01.0080, 'OB-L': 15.9994, 'OA-L': 15.9994, 'C6-L': 13.0190,\
           'O7-L': 15.9994, 'C9-L': 12.011,  'C8-L': 14.0270, 'O9-L': 15.9994,\
           'O6-L': 15.9994, 'C3-L': 13.0190, 'C2-L': 15.035,  'O5-L': 15.9994,\
           'O4-L': 15.9994, 'O3-L': 15.9994, 'C1-L': 12.0110, 'C5-L': 15.0350,\
           'C4-L': 12.0110}

# define dictionary of atom masses for EA (in g/mol)
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

# get mean velocity in each bin per frame per LGA/EA
VxLGA,VyLGA,VzLGA,nTotLGA = [],[],[],[]
VxEEE,VyEEE,VzEEE,nTotEEE = [],[],[],[]
vxL,vyL,vzL = np.zeros(numbin),np.zeros(numbin),np.zeros(numbin)
vxE,vyE,vzE = np.zeros(numbin),np.zeros(numbin),np.zeros(numbin)
nL,nE       = np.zeros(numbin),np.zeros(numbin)
binstrt = intfce - 0.5*numbin*dxbin
binfnsh = intfce + 0.5*numbin*dxbin
bins    = np.linspace(binstrt,binfnsh,numbin+1)
indx = 0
# loop over frames
for i in range(frames):
    vxL,vyL,vzL,nL = vxL*0,vyL*0,vzL*0,nL*0
    vxE,vyE,vzE,nE = vxE*0,vyE*0,vzE*0,nE*0
    # loop over LGA molecules in frame
    for j in range(lgaMOL):
        xpos = pxcm[indx]
        if xpos < binstrt or xpos > binfnsh:
            indx += 1
            continue
        else:
            binNum = int((xpos-binstrt)/dxbin)
            vxL[binNum] += vxcm[indx]
            vyL[binNum] += vycm[indx]
            vzL[binNum] += vzcm[indx]
            nL[binNum]  += 1
            indx += 1
    VxLGA.append(vxL/nL)
    VyLGA.append(vyL/nL)
    VzLGA.append(vzL/nL)
    nTotLGA.append(nL)
    # loop over EA molecules in frame
    for j in range(eeeMOL):
        xpos = pxcm[indx]
        if xpos < binstrt or xpos > binfnsh:
            indx += 1
            continue
        else:
            binNum = int((xpos-binstrt)/dxbin)
            vxE[binNum] += vxcm[indx]
            vyE[binNum] += vycm[indx]
            vzE[binNum] += vzcm[indx]
            nE[binNum]  += 1
            indx += 1
    VxEEE.append(vxE/nE)
    VyEEE.append(vyE/nE)
    VzEEE.append(vzE/nE)
    nTotEEE.append(nE)

# get mean velocities per frame in each bin for plotting
aveVxLGA,aveVyLGA,aveVzLGA = [],[],[]
aveVxEEE,aveVyEEE,aveVzEEE = [],[],[]
vxL,vyL,vzL = np.zeros(frames),np.zeros(frames),np.zeros(frames)
vxE,vyE,vzE = np.zeros(frames),np.zeros(frames),np.zeros(frames)
# loop over bins
for i in range(numbin):
    vxL,vyL,vzL = vxL*0,vyL*0,vzL*0
    vxE,vyE,vzE = vxE*0,vyE*0,vzE*0
    # loop over frames per bin
    for j in range(frames):
        # get LGA mols
        vxL[j] = VxLGA[j][i]
        vyL[j] = VyLGA[j][i]
        vzL[j] = VzLGA[j][i]
        # get EA mols
        vxE[j] = VxEEE[j][i]
        vyE[j] = VyEEE[j][i]
        vzE[j] = VzEEE[j][i]
    aveVxLGA.append(vxL)
    aveVyLGA.append(vyL)
    aveVzLGA.append(vzL)
    aveVxEEE.append(vxE)
    aveVyEEE.append(vyE)
    aveVzEEE.append(vzE)

# when no molecules exist in bin set it to zero
for i in range(numbin):
    aveVxLGA[i][np.where(np.isnan(aveVxLGA[i]) == True)] = 0
    aveVyLGA[i][np.where(np.isnan(aveVyLGA[i]) == True)] = 0
    aveVzLGA[i][np.where(np.isnan(aveVzLGA[i]) == True)] = 0
    aveVxEEE[i][np.where(np.isnan(aveVxEEE[i]) == True)] = 0
    aveVyEEE[i][np.where(np.isnan(aveVyEEE[i]) == True)] = 0
    aveVzEEE[i][np.where(np.isnan(aveVzEEE[i]) == True)] = 0

# get density per frame per bin (in kg/m^3)
avogad = 6.02214e23
mLGAkg = (mLGA/1000)/avogad
mEEEkg = (mEEE/1000)/avogad
mLGAu  = mLGA*np.array(nTotLGA)
mEEEu  = mEEE*np.array(nTotEEE)
volm3  = (dxbin*ybox*zbox)*(1e-27)
rhoLGA = (mLGAkg*np.array(nTotLGA))/volm3
rhoEEE = (mEEEkg*np.array(nTotEEE))/volm3


# make plots
print 'making plots ...'
# set plot parameters
plt.ioff()
cnt = 0
fontsize = 18
factor   = 1.20
color    = ['r','g','b']
# plot bin image
plt.figure(figsize=(8,5))
plt.title('Bins',fontsize=20)
plt.xlim(0,xbox)
plt.yticks([])
plt.xlabel('x-axis (nm)',fontsize=16)
plt.vlines(bins,0,1)
plt.vlines(intfce,0,1,'m',lw=2,linestyle='dashed')
for i in range(numbin):
    plt.text(bins[i]+0.20*dxbin,0.5,str(i+1),fontsize=15,color='b')
plt.savefig(path+'Plots/0binsvel.png')
plt.close()
# plot mols per bin per frame
plt.figure()
for i in range(numbin):
    plt.plot(t,np.array(nTotLGA)[:,i],color[i]+'-',label='bin '+str(i+1))
    plt.plot(t,np.array(nTotEEE)[:,i],color[i]+'.-')
plt.xlabel('Time (ns)',fontsize=15)
plt.ylabel('Molecules',fontsize=15)
lim = plt.ylim()
plt.ylim(lim[0]-lim[1]*0.01,lim[1]*factor)
plt.grid()
plt.legend(loc=1)
plt.savefig(path+'Plots/0molcountvel.png')
plt.close()
# plot mass per bin per frame
plt.figure()
for i in range(numbin):
    plt.plot(t,mLGAu[:,i],color[i]+'-',label='bin '+str(i+1))
    plt.plot(t,mEEEu[:,i],color[i]+'.-')
plt.xlabel('Time (ns)',fontsize=15)
plt.ylabel('Mass (g/mol)',fontsize=15)
lim = plt.ylim()
plt.ylim(lim[0]-lim[1]*0.01,lim[1]*factor)
plt.grid()
plt.legend(loc=1)
plt.savefig(path+'Plots/0massvel.png')
plt.close()
# plot density per bin per frame
plt.figure()
for i in range(numbin):
    plt.plot(t,rhoLGA[:,i],color[i]+'-',label='bin '+str(i+1))
    plt.plot(t,rhoEEE[:,i],color[i]+'.-')
plt.hlines(1340,0,t[-1],'k',linestyle='dashed',linewidth=2)
plt.hlines(897,0,t[-1],'k',linestyle='dashed',linewidth=2)
plt.xlabel('Time (ns)',fontsize=15)
plt.ylabel('Density (kg/m$^3$)',fontsize=15)
lim = plt.ylim()
plt.ylim(lim[0]-lim[1]*0.01,lim[1]*factor)
plt.grid()
plt.legend(loc=1)
plt.savefig(path+'Plots/0densityvel.png')
plt.close()
# plot velocities
for i in range(numbin):
    plt.figure()
    plt.subplot(311)
    plt.plot(t,aveVxLGA[i],'b-',label='LGA')
    plt.plot(t,aveVxEEE[i],'g-',label='EA')
    plt.legend(loc=9,bbox_to_anchor=(0.5,1.3),ncol=2)
    plt.grid()
    plt.ylabel(r'$ \overline{Vx}$',fontsize=fontsize)
    plt.hlines(0,0,t[-1],'k',linestyle='dashed')
    plt.subplot(312)
    plt.plot(t,aveVyLGA[i],'b-')
    plt.plot(t,aveVyEEE[i],'g-')
    plt.grid()
    plt.ylabel(r'$ \overline{Vy}$',fontsize=fontsize)
    plt.hlines(0,0,t[-1],'k',linestyle='dashed')
    plt.subplot(313)
    plt.plot(t,aveVzLGA[i],'b-')
    plt.plot(t,aveVzEEE[i],'g-')
    plt.grid()
    plt.ylabel(r'$ \overline{Vz}$',fontsize=fontsize)
    plt.xlabel('Time (ns)',fontsize=15)
    plt.hlines(0,0,t[-1],'k',linestyle='dashed')
    plt.annotate('Bin ='+'%2i'%(i+1),fontsize=fontsize,xy=(0.78,0.92),xycoords='figure fraction')
    plt.annotate('(nm/ps)',fontsize=14,xy=(0.05,0.93),xycoords='figure fraction')
    plt.savefig(path+'Plots/'+str(cnt+1)+'velocitytime.png')
    plt.close()
    cnt += 1

# print done
print cnt+4,'images saved'
print 'done !\n'
plt.ion()