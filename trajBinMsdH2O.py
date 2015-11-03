"""
Analyzes LGA+H2O trajectory to get MSD of molecules in bins along x-axis
"""

import numpy as np
import matplotlib.pyplot as plt

# set input parameters
# -----------------------------------------------
path = 'C:/Users/Code 8000/Desktop/Grad Research/LGA+H2O/Merged/H2O BigBox2x/Analysis/Middle Run 10ns/'
fgro = 'nvt0_LGA+H2O_0.5Agap_bigbox2x_analysis_reimagedTRAJ.gro'
lgaMOL  = 256        # number of LGA molecules per frame
lgaATM  = 21         # number of LGA atoms per molecule
h2oMOL  = 8442       # number of H2O molecules per frame
h2oATM  = 3          # number of H2O atoms per molecule
frames  = 101        # number of frames in trajectory
dtstep  = 0.1        # time step between frames (ns)
xbox    = 20.20600   # box length in x-direction (nm)
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
    # loop over H2O molecules in frame
    for j in range(h2oMOL):
        pxcom,pycom,pzcom = 0,0,0
        # loop over H2O atoms in molecule
        for k in range(h2oATM):
            pxcom += massH2O[natm[indx]]*pxatm[indx]
            pycom += massH2O[natm[indx]]*pyatm[indx]
            pzcom += massH2O[natm[indx]]*pzatm[indx]
            indx += 1
        res.append(nres[indx-1]+'-'+str(i))
        pxcm.append(pxcom/mH2O)
        pycm.append(pycom/mH2O)
        pzcm.append(pzcom/mH2O)

# convert lists to numpy arrays
res  = np.array(res)
pxcm = np.array(pxcm)
pycm = np.array(pycm)
pzcm = np.array(pzcm)

# get index list for sorting and binning along x
indx  = []
Nmol  = lgaMOL+h2oMOL
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

# get MSD in each bin per frame per LGA/H2O
msdXBinLGA,msdYBinLGA,msdZBinLGA,nTotLGA  = [],[],[],[]
msdXBinH2O,msdYBinH2O,msdZBinH2O,nTotH2O  = [],[],[],[]
msdxL,msdyL,msdzL = np.zeros(numbin),np.zeros(numbin),np.zeros(numbin)
msdxH,msdyH,msdzH = np.zeros(numbin),np.zeros(numbin),np.zeros(numbin)
nL,nH             = np.zeros(numbin),np.zeros(numbin)
binstrt = intfce - 0.5*numbin*dxbin
binfnsh = intfce + 0.5*numbin*dxbin
bins    = np.linspace(binstrt,binfnsh,numbin+1)
binNum  = ((xcmf[0]-binstrt)/dxbin).astype(int)
# loop over frames
for i in range(strmsd,frames-1):
    msdxL,msdyL,msdzL,nL = msdxL*0,msdyL*0,msdzL*0,nL*0
    msdxH,msdyH,msdzH,nH = msdxH*0,msdyH*0,msdzH*0,nH*0
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
    # loop over H2O molecules in frame
    for j in range(lgaMOL,Nmol):
        xpos = xcmf[0,j]
        if xpos < binstrt or xpos > binfnsh:
            continue
        else:
            k      = binNum[j]
            nH[k] += 1
            # calculate msd from fixed start time
            if strtime == 'yes':
                msdxH[k] += (xcmf[i+1,j] - xcmf[strmsd,j])**2
                msdyH[k] += (ycmf[i+1,j] - ycmf[strmsd,j])**2
                msdzH[k] += (zcmf[i+1,j] - zcmf[strmsd,j])**2
            # otherwise calculate msd from moving window
            else:
                msdxH[k] += (xcmf[i+1,j] - xcmf[i,j])**2
                msdyH[k] += (ycmf[i+1,j] - ycmf[i,j])**2
                msdzH[k] += (zcmf[i+1,j] - zcmf[i,j])**2
    msdXBinH2O.append(msdxH/nH)
    msdYBinH2O.append(msdyH/nH)
    msdZBinH2O.append(msdzH/nH)
    nTotH2O.append(nH)

# get MSD per frame in each bin for plotting
msdXLGA,msdYLGA,msdZLGA = [],[],[]
msdXH2O,msdYH2O,msdZH2O = [],[],[]
# loop over bins
for i in range(numbin):
    # get LGA mols
    msdXLGA.append(np.array(msdXBinLGA)[:,i])
    msdYLGA.append(np.array(msdYBinLGA)[:,i])
    msdZLGA.append(np.array(msdZBinLGA)[:,i])
    # get H2O mols
    msdXH2O.append(np.array(msdXBinH2O)[:,i])
    msdYH2O.append(np.array(msdYBinH2O)[:,i])
    msdZH2O.append(np.array(msdZBinH2O)[:,i])

# when no molecules exist in bin set it to zero
for i in range(numbin):
    msdXLGA[i][np.where(np.isnan(msdXLGA[i]) == True)] = 0
    msdYLGA[i][np.where(np.isnan(msdYLGA[i]) == True)] = 0
    msdZLGA[i][np.where(np.isnan(msdZLGA[i]) == True)] = 0
    msdXH2O[i][np.where(np.isnan(msdXH2O[i]) == True)] = 0
    msdYH2O[i][np.where(np.isnan(msdYH2O[i]) == True)] = 0
    msdZH2O[i][np.where(np.isnan(msdZH2O[i]) == True)] = 0


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
plt.savefig(path+'Images/0binsmsd.png')
plt.close()
# plot mols per bin at frame zero to check it is constant in each iteration
# note: this is because mols are always associated to their starting bin for msd
plt.figure()
for i in range(numbin):
    plt.plot(t[strmsd+1::],np.array(nTotLGA)[:,i],color[i]+'-',label='bin '+str(i+1))
    plt.plot(t[strmsd+1::],np.array(nTotH2O)[:,i],color[i]+'.-')
plt.xlabel('Time (ns)',fontsize=15)
plt.ylabel('Molecules',fontsize=15)
lim = plt.ylim()
plt.ylim(lim[0]-lim[1]*0.01,lim[1]*factor)
plt.grid()
plt.legend(loc=1)
plt.savefig(path+'Images/0molcountmsd.png')
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
    plt.plot(t[strmsd+1::],msdXH2O[i],'g-',label='H2O')
    plt.legend(loc=2)
    plt.grid()
    plt.subplot(323)
    plt.plot(t[strmsd+1::],msdYLGA[i],'b-')
    plt.ylabel('MSD-y',fontsize=14)
    plt.grid()
    plt.subplot(324)
    plt.plot(t[strmsd+1::],msdYH2O[i],'g-')
    plt.grid()
    plt.subplot(325)
    plt.plot(t[strmsd+1::],msdZLGA[i],'b-')
    plt.ylabel('MSD-z',fontsize=14)
    plt.xlabel('Time (ns)',fontsize=14)
    plt.grid()
    plt.subplot(326)
    plt.plot(t[strmsd+1::],msdZH2O[i],'g-')
    plt.xlabel('Time (ns)',fontsize=14)
    plt.grid()
    plt.annotate('Bin ='+'%2i'%(i+1),fontsize=fontsize,xy=(0.78,0.92),xycoords='figure fraction')
    plt.annotate('(nm$^2$)',fontsize=14,xy=(0.05,0.93),xycoords='figure fraction')
    if strtime == 'yes':
        plt.savefig(path+'Images/'+str(cnt+1)+'msdtime_f'+str(strmsd)+'.png')
    else:
        plt.savefig(path+'Images/'+str(cnt+1)+'msdtime_fprev.png')
    plt.close()
    cnt += 1

## alternate plot style
#for i in range(numbin):
#    plt.figure()
#    plt.subplot(311)
#    plt.plot(t[strmsd+1::],msdXLGA[i],'b-',label='LGA')
#    plt.plot(t[strmsd+1::],msdXH2O[i],'g-',label='H2O')
#    plt.legend(loc=9,bbox_to_anchor=(0.5,1.3),ncol=2)
#    plt.grid()
#    plt.ylabel('MSD-x',fontsize=14)
#    plt.subplot(312)
#    plt.plot(t[strmsd+1::],msdYLGA[i],'b-')
#    plt.plot(t[strmsd+1::],msdYH2O[i],'g-')
#    plt.grid()
#    plt.ylabel('MSD-y',fontsize=14)
#    plt.subplot(313)
#    plt.plot(t[strmsd+1::],msdZLGA[i],'b-')
#    plt.plot(t[strmsd+1::],msdZH2O[i],'g-')
#    plt.grid()
#    plt.ylabel('MSD-z',fontsize=14)
#    plt.xlabel('Time (ns)',fontsize=15)
#    plt.annotate('Bin ='+'%2i'%(i+1),fontsize=fontsize,xy=(0.78,0.92),xycoords='figure fraction')
#    plt.annotate('(nm$^2$)',fontsize=14,xy=(0.05,0.93),xycoords='figure fraction')
#    plt.savefig(path+'Images/'+str(cnt+1)+'msdtime.png')
#    plt.close()
#    cnt += 1

# print done
print cnt+2,'images saved'
print 'done !\n'
plt.ion()