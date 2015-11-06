"""
Analyzes LGA+H2O trajectory to get mean velocity of molecules in bins along x-axis
"""

import numpy as np
import matplotlib.pyplot as plt

# set input parameters
# -----------------------------------------------
path = 'C:/Users/Code 8000/Desktop/Grad Research/LGA+H2O/Merged/H2O BigBox2x/Analysis/Middle Run 10ns/'
fgro = 'nvt0_LGA+H2O_0.5Agap_bigbox3x_analysis_originalTRAJ.gro'
lgaMOL = 256        # number of LGA molecules per frame
lgaATM = 21         # number of LGA atoms per molecule
h2oMOL = 8442       # number of H2O molecules per frame
h2oATM = 3          # number of H2O atoms per molecule
frames = 101        # number of frames in trajectory
dtstep = 0.1        # time step between frames (ns)
xbox   = 20.20600   # box length in x-direction (nm)
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
H2Oatms = h2oMOL*h2oATM
TOTatms = LGAatms+H2Oatms

# define dictionary of atom masses for LGA (in g/mol)
mLGA    = 278.2166
massLGA = {'C0-L': 14.0270, 'H2-L': 01.0080, 'C7-L': 12.0110, 'O8-L': 15.9994,\
           'H1-L': 01.0080, 'OB-L': 15.9994, 'OA-L': 15.9994, 'C6-L': 13.0190,\
           'O7-L': 15.9994, 'C9-L': 12.011,  'C8-L': 14.0270, 'O9-L': 15.9994,\
           'O6-L': 15.9994, 'C3-L': 13.0190, 'C2-L': 15.035,  'O5-L': 15.9994,\
           'O4-L': 15.9994, 'O3-L': 15.9994, 'C1-L': 12.0110, 'C5-L': 15.0350,\
           'C4-L': 12.0110}

# define dictionary of atom masses for H2O
mH2O    = 18.0154
massH2O = {'OW-H':9.9514,'HW1-H':4.0320,'HW2-H':4.0320}

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
# tag atom types for LGA or H2O
nres  = np.char.strip(nres)
natm  = np.char.strip(natm)
for i in range(frames):
    for j in range(0,LGAatms):
        natm[i][j] = natm[i][j]+'-L'
    for j in range(LGAatms,TOTatms):
        natm[i][j] = natm[i][j]+'-H'
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
    # loop over H2O molecules in frame
    for j in range(h2oMOL):
        pxcom,pycom,pzcom = 0,0,0
        vxcom,vycom,vzcom = 0,0,0
        # loop over H2O atoms in molecule
        for k in range(h2oATM):
            pxcom += massH2O[natm[indx]]*pxatm[indx]
            pycom += massH2O[natm[indx]]*pyatm[indx]
            pzcom += massH2O[natm[indx]]*pzatm[indx]
            vxcom += massH2O[natm[indx]]*vxatm[indx]
            vycom += massH2O[natm[indx]]*vyatm[indx]
            vzcom += massH2O[natm[indx]]*vzatm[indx]
            indx += 1
        pxcm.append(pxcom/mH2O)
        pycm.append(pycom/mH2O)
        pzcm.append(pzcom/mH2O)
        vxcm.append(vxcom/mH2O)
        vycm.append(vycom/mH2O)
        vzcm.append(vzcom/mH2O)

# get mean velocity in each bin per frame per LGA/H2O
VxLGA,VyLGA,VzLGA,nTotLGA = [],[],[],[]
VxH2O,VyH2O,VzH2O,nTotH2O = [],[],[],[]
vxL,vyL,vzL = np.zeros(numbin),np.zeros(numbin),np.zeros(numbin)
vxH,vyH,vzH = np.zeros(numbin),np.zeros(numbin),np.zeros(numbin)
nL,nH       = np.zeros(numbin),np.zeros(numbin)
binstrt = intfce - 0.5*numbin*dxbin
binfnsh = intfce + 0.5*numbin*dxbin
bins    = np.linspace(binstrt,binfnsh,numbin+1)
indx = 0
# loop over frames
for i in range(frames):
    vxL,vyL,vzL,nL = vxL*0,vyL*0,vzL*0,nL*0
    vxH,vyH,vzH,nH = vxH*0,vyH*0,vzH*0,nH*0
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
    # loop over H2O molecules in frame
    for j in range(h2oMOL):
        xpos = pxcm[indx]
        if xpos < binstrt or xpos > binfnsh:
            indx += 1
            continue
        else:
            binNum = int((xpos-binstrt)/dxbin)
            vxH[binNum] += vxcm[indx]
            vyH[binNum] += vycm[indx]
            vzH[binNum] += vzcm[indx]
            nH[binNum]  += 1
            indx += 1
    VxH2O.append(vxH/nH)
    VyH2O.append(vyH/nH)
    VzH2O.append(vzH/nH)
    nTotH2O.append(nH)

# get mean velocities per frame in each bin for plotting
aveVxLGA,aveVyLGA,aveVzLGA = [],[],[]
aveVxH2O,aveVyH2O,aveVzH2O = [],[],[]
vxL,vyL,vzL = np.zeros(frames),np.zeros(frames),np.zeros(frames)
vxH,vyH,vzH = np.zeros(frames),np.zeros(frames),np.zeros(frames)
# loop over bins
for i in range(numbin):
    vxL,vyL,vzL = vxL*0,vyL*0,vzL*0
    vxH,vyH,vzH = vxH*0,vyH*0,vzH*0
    # loop over frames per bin
    for j in range(frames):
        # get LGA mols
        vxL[j] = VxLGA[j][i]
        vyL[j] = VyLGA[j][i]
        vzL[j] = VzLGA[j][i]
        # get H2O mols
        vxH[j] = VxH2O[j][i]
        vyH[j] = VyH2O[j][i]
        vzH[j] = VzH2O[j][i]
    aveVxLGA.append(vxL)
    aveVyLGA.append(vyL)
    aveVzLGA.append(vzL)
    aveVxH2O.append(vxH)
    aveVyH2O.append(vyH)
    aveVzH2O.append(vzH)

# when no molecules exist in bin set it to zero
for i in range(numbin):
    aveVxLGA[i][np.where(np.isnan(aveVxLGA[i]) == True)] = 0
    aveVyLGA[i][np.where(np.isnan(aveVyLGA[i]) == True)] = 0
    aveVzLGA[i][np.where(np.isnan(aveVzLGA[i]) == True)] = 0
    aveVxH2O[i][np.where(np.isnan(aveVxH2O[i]) == True)] = 0
    aveVyH2O[i][np.where(np.isnan(aveVyH2O[i]) == True)] = 0
    aveVzH2O[i][np.where(np.isnan(aveVzH2O[i]) == True)] = 0

# get density per frame per bin (in kg/m^3)
avogad = 6.02214e23
mLGAkg = (mLGA/1000)/avogad
mH2Okg = (mH2O/1000)/avogad
mLGAu  = mLGA*np.array(nTotLGA)
mH2Ou  = mH2O*np.array(nTotH2O)
volm3  = (dxbin*ybox*zbox)*(1e-27)
rhoLGA = (mLGAkg*np.array(nTotLGA))/volm3
rhoH2O = (mH2Okg*np.array(nTotH2O))/volm3


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
    plt.plot(t,np.array(nTotH2O)[:,i],color[i]+'.-')
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
    plt.plot(t,mH2Ou[:,i],color[i]+'.-')
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
    plt.plot(t,rhoH2O[:,i],color[i]+'.-')
plt.hlines(1340,0,t[-1],'k',linestyle='dashed',linewidth=2)
plt.hlines(997,0,t[-1],'k',linestyle='dashed',linewidth=2)
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
    plt.plot(t,aveVxH2O[i],'g-',label='H2O')
    plt.legend(loc=9,bbox_to_anchor=(0.5,1.3),ncol=2)
    plt.grid()
    plt.ylabel(r'$ \overline{Vx}$',fontsize=fontsize)
    plt.hlines(0,0,t[-1],'k',linestyle='dashed')
    plt.subplot(312)
    plt.plot(t,aveVyLGA[i],'b-')
    plt.plot(t,aveVyH2O[i],'g-')
    plt.grid()
    plt.ylabel(r'$ \overline{Vy}$',fontsize=fontsize)
    plt.hlines(0,0,t[-1],'k',linestyle='dashed')
    plt.subplot(313)
    plt.plot(t,aveVzLGA[i],'b-')
    plt.plot(t,aveVzH2O[i],'g-')
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
print cnt+3,'images saved'
print 'done !\n'
plt.ion()