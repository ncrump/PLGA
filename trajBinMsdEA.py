"""
Analyzes LGA+EA trajectory to get MSD of molecules in bins along x-axis
"""

import numpy as np
import matplotlib.pyplot as plt

# set input parameters
# -----------------------------------------------
path = 'C:/Users/Code 8000/Desktop/Grad Research/LGA+EA/Merged/EA BigBox3x/Analysis/1Middle Run 10ns/'
fgro = 'nvt1_LGA+EA_0.5Agap_eabigbox3x_analysis_originalTRAJ.gro'
lgaMOL  = 256        # number of LGA molecules per frame
lgaATM  = 21         # number of LGA atoms per molecule
eeeMOL  = 1620       # number of EA molecules per frame
eeeATM  = 6          # number of EA atoms per molecule
frames  = 101        # number of frames in trajectory
dtstep  = 0.1        # time step between frames (ns)
xbox    = 20.34200   # box length in x-direction (nm)
ybox    =  5.59744   # box length in y-direction (nm)
zbox    =  3.19734   # box length in z-direction (nm)
dxbin   = 2          # width of bin centered around interface (nm)
numbin  = 3          # number of bins centered around interface
intfce  = 5.748      # interface distance from origin (nm)
strmsd  = 10         # frame number to start msd calculation
strtime = 'yes'      # calculate msd from fixed start time or moving window
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

# calculate center mass position of each molecule per frame
res,pxcm,pycm,pzcm = [],[],[],[]
indx = 0
# loop over frames
for i in range(frames):
    # loop over LGA molecules in frame
    for j in range(lgaMOL):
        pxcom,pycom,pzcom = 0,0,0
        # loop over LGA atoms in molecule
        for k in range(lgaATM):
            pxcom += massLGA[natm[indx]]*pxatm[indx]
            pycom += massLGA[natm[indx]]*pyatm[indx]
            pzcom += massLGA[natm[indx]]*pzatm[indx]
            indx += 1
        res.append(nres[indx-1]+'-'+str(i))
        pxcm.append(pxcom/mLGA)
        pycm.append(pycom/mLGA)
        pzcm.append(pzcom/mLGA)
    # loop over EA molecules in frame
    for j in range(eeeMOL):
        pxcom,pycom,pzcom = 0,0,0
        # loop over EA atoms in molecule
        for k in range(eeeATM):
            pxcom += massEEE[natm[indx]]*pxatm[indx]
            pycom += massEEE[natm[indx]]*pyatm[indx]
            pzcom += massEEE[natm[indx]]*pzatm[indx]
            indx += 1
        res.append(nres[indx-1]+'-'+str(i))
        pxcm.append(pxcom/mEEE)
        pycm.append(pycom/mEEE)
        pzcm.append(pzcom/mEEE)

# convert lists to numpy arrays
res  = np.array(res)
pxcm = np.array(pxcm)
pycm = np.array(pycm)
pzcm = np.array(pzcm)

# get index list for sorting and binning along x
indx  = []
Nmol  = lgaMOL+eeeMOL
xsort = np.sort(pxcm[0:Nmol])
for i in range(Nmol):
    k = np.where(pxcm[0:Nmol] == xsort[i])[0][0]
    indx.append(k)

# define new arrays to hold re-ordered molecules
resf = np.chararray((frames,Nmol),itemsize=12)
xcmf = np.zeros((frames,Nmol))
ycmf = np.zeros((frames,Nmol))
zcmf = np.zeros((frames,Nmol))

# re-order molecules in each frame to match order of first frame sorted along x
for i in range(frames):
    for j in range(Nmol):
        resf[i,j] = res[indx[j]  + i*Nmol]
        xcmf[i,j] = pxcm[indx[j] + i*Nmol]
        ycmf[i,j] = pycm[indx[j] + i*Nmol]
        zcmf[i,j] = pzcm[indx[j] + i*Nmol]

# get MSD in each bin per frame per LGA/EA
msdXBinLGA,msdYBinLGA,msdZBinLGA,nTotLGA  = [],[],[],[]
msdXBinEEE,msdYBinEEE,msdZBinEEE,nTotEEE  = [],[],[],[]
msdxL,msdyL,msdzL = np.zeros(numbin),np.zeros(numbin),np.zeros(numbin)
msdxE,msdyE,msdzE = np.zeros(numbin),np.zeros(numbin),np.zeros(numbin)
nL,nE             = np.zeros(numbin),np.zeros(numbin)
binstrt = intfce - 0.5*numbin*dxbin
binfnsh = intfce + 0.5*numbin*dxbin
bins    = np.linspace(binstrt,binfnsh,numbin+1)
binNum  = ((xcmf[0]-binstrt)/dxbin).astype(int)
# loop over frames
for i in range(strmsd,frames-1):
    msdxL,msdyL,msdzL,nL = msdxL*0,msdyL*0,msdzL*0,nL*0
    msdxE,msdyE,msdzE,nE = msdxE*0,msdyE*0,msdzE*0,nE*0
    # loop over LGA molecules in frame
    for j in range(0,lgaMOL):
        xpos = xcmf[0,j]
        if xpos < binstrt or xpos > binfnsh:
            continue
        else:
            k      = binNum[j]
            nL[k] += 1
            # calculate msd from fixed start time
            if strtime == 'yes':
                msdxL[k] += (xcmf[i+1,j] - xcmf[strmsd,j])**2
                msdyL[k] += (ycmf[i+1,j] - ycmf[strmsd,j])**2
                msdzL[k] += (zcmf[i+1,j] - zcmf[strmsd,j])**2
            # otherwise calculate msd from moving window
            else:
                msdxL[k] += (xcmf[i+1,j] - xcmf[i,j])**2
                msdyL[k] += (ycmf[i+1,j] - ycmf[i,j])**2
                msdzL[k] += (zcmf[i+1,j] - zcmf[i,j])**2
    msdXBinLGA.append(msdxL/nL)
    msdYBinLGA.append(msdyL/nL)
    msdZBinLGA.append(msdzL/nL)
    nTotLGA.append(nL)
    # loop over EA molecules in frame
    for j in range(lgaMOL,Nmol):
        xpos = xcmf[0,j]
        if xpos < binstrt or xpos > binfnsh:
            continue
        else:
            k      = binNum[j]
            nE[k] += 1
            # calculate msd from fixed start time
            if strtime == 'yes':
                msdxE[k] += (xcmf[i+1,j] - xcmf[strmsd,j])**2
                msdyE[k] += (ycmf[i+1,j] - ycmf[strmsd,j])**2
                msdzE[k] += (zcmf[i+1,j] - zcmf[strmsd,j])**2
            # otherwise calculate msd from moving window
            else:
                msdxE[k] += (xcmf[i+1,j] - xcmf[i,j])**2
                msdyE[k] += (ycmf[i+1,j] - ycmf[i,j])**2
                msdzE[k] += (zcmf[i+1,j] - zcmf[i,j])**2
    msdXBinEEE.append(msdxE/nE)
    msdYBinEEE.append(msdyE/nE)
    msdZBinEEE.append(msdzE/nE)
    nTotEEE.append(nE)

# get MSD per frame in each bin for plotting
msdXLGA,msdYLGA,msdZLGA = [],[],[]
msdXEEE,msdYEEE,msdZEEE = [],[],[]
# loop over bins
for i in range(numbin):
    # get LGA mols
    msdXLGA.append(np.array(msdXBinLGA)[:,i])
    msdYLGA.append(np.array(msdYBinLGA)[:,i])
    msdZLGA.append(np.array(msdZBinLGA)[:,i])
    # get EA mols
    msdXEEE.append(np.array(msdXBinEEE)[:,i])
    msdYEEE.append(np.array(msdYBinEEE)[:,i])
    msdZEEE.append(np.array(msdZBinEEE)[:,i])

# get radial msd (y+z) in interface bin only
msdRLGA = msdYLGA[1]+msdZLGA[1]
msdREEE = msdYEEE[1]+msdZEEE[1]

# get fits to lines for diffusion coefficient
mlga,blga = np.polyfit(t[strmsd+1::],msdXLGA[1],1)
meee,beee = np.polyfit(t[strmsd+1::],msdXEEE[1],1)
# get diffusion coefficient (in m^2/s)
Dlga = (mlga/6.0)*10**-9
Deee = (meee/6.0)*10**-9

# when no molecules exist in bin set it to zero
for i in range(numbin):
    msdXLGA[i][np.where(np.isnan(msdXLGA[i]) == True)] = 0
    msdYLGA[i][np.where(np.isnan(msdYLGA[i]) == True)] = 0
    msdZLGA[i][np.where(np.isnan(msdZLGA[i]) == True)] = 0
    msdXEEE[i][np.where(np.isnan(msdXEEE[i]) == True)] = 0
    msdYEEE[i][np.where(np.isnan(msdYEEE[i]) == True)] = 0
    msdZEEE[i][np.where(np.isnan(msdZEEE[i]) == True)] = 0


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
plt.savefig(path+'Plots/0binsmsd.png')
plt.close()
# plot mols per bin at frame zero to check it is constant in each iteration
# note: this is because mols are always associated to their starting bin for msd
plt.figure()
for i in range(numbin):
    plt.plot(t[strmsd+1::],np.array(nTotLGA)[:,i],color[i]+'-',label='bin '+str(i+1))
    plt.plot(t[strmsd+1::],np.array(nTotEEE)[:,i],color[i]+'.-')
plt.xlabel('Time (ns)',fontsize=15)
plt.ylabel('Molecules',fontsize=15)
lim = plt.ylim()
plt.ylim(lim[0]-lim[1]*0.01,lim[1]*factor)
plt.grid()
plt.legend(loc=1)
plt.savefig(path+'Plots/0molcountmsd.png')
plt.close()
# plot msd
for i in range(numbin):
    plt.figure()
    plt.subplot(321)
    plt.plot(t[strmsd+1::],msdXLGA[i],'b-',label='LGA')
    plt.ylabel('MSD-x',fontsize=14)
    plt.legend(loc=2)
    plt.grid()
    plt.subplot(322)
    plt.plot(t[strmsd+1::],msdXEEE[i],'g-',label='EA')
    plt.legend(loc=2)
    plt.grid()
    plt.subplot(323)
    plt.plot(t[strmsd+1::],msdYLGA[i],'b-')
    plt.ylabel('MSD-y',fontsize=14)
    plt.grid()
    plt.subplot(324)
    plt.plot(t[strmsd+1::],msdYEEE[i],'g-')
    plt.grid()
    plt.subplot(325)
    plt.plot(t[strmsd+1::],msdZLGA[i],'b-')
    plt.ylabel('MSD-z',fontsize=14)
    plt.xlabel('Time (ns)',fontsize=14)
    plt.grid()
    plt.subplot(326)
    plt.plot(t[strmsd+1::],msdZEEE[i],'g-')
    plt.xlabel('Time (ns)',fontsize=14)
    plt.grid()
    plt.annotate('Bin ='+'%2i'%(i+1),fontsize=fontsize,xy=(0.78,0.92),xycoords='figure fraction')
    plt.annotate('(nm$^2$)',fontsize=14,xy=(0.05,0.93),xycoords='figure fraction')
    if strtime == 'yes':
        plt.savefig(path+'Plots/'+str(cnt+1)+'msdtime_f'+str(strmsd)+'.png')
    else:
        plt.savefig(path+'Plots/'+str(cnt+1)+'msdtime_fprev.png')
    plt.close()
    cnt += 1
# plot lateral and radial msd for interface bin only
plt.figure()
plt.subplot(221)
plt.plot(t[strmsd+1::],msdXLGA[1],'b-',label='LGA')
plt.plot(t[strmsd+1::],mlga*t[strmsd+1::]+blga,'r-')
plt.ylabel('MSD-x',fontsize=14)
plt.legend(loc=2)
plt.grid()
plt.subplot(222)
plt.plot(t[strmsd+1::],msdXEEE[1],'g-',label='EA')
plt.plot(t[strmsd+1::],meee*t[strmsd+1::]+beee,'r-')
plt.legend(loc=2)
plt.grid()
plt.subplot(223)
plt.plot(t[strmsd+1::],msdRLGA,'b-')
plt.ylabel('MSD-yz',fontsize=14)
plt.xlabel('Time (ns)',fontsize=14)
plt.grid()
plt.subplot(224)
plt.plot(t[strmsd+1::],msdREEE,'g-')
plt.xlabel('Time (ns)',fontsize=14)
plt.grid()
plt.annotate('Bin ='+'%2i'%(2),fontsize=fontsize,xy=(0.78,0.92),xycoords='figure fraction')
plt.annotate('(nm$^2$)',fontsize=14,xy=(0.05,0.93),xycoords='figure fraction')
plt.annotate('D=%6.3e'%Dlga+' (m$^2$/s)',fontsize=12,xy=(0.25,0.56),xycoords='figure fraction')
plt.annotate('D=%6.3e'%Deee+' (m$^2$/s)',fontsize=12,xy=(0.67,0.56),xycoords='figure fraction')
if strtime == 'yes':
    plt.savefig(path+'Plots/'+str(cnt+1)+'msdtime_f'+str(strmsd)+'.png')
else:
    plt.savefig(path+'Plots/'+str(cnt+1)+'msdtime_fprev.png')
plt.close()
cnt += 1

# print done
print cnt+2,'images saved'
print 'done !\n'
plt.ion()