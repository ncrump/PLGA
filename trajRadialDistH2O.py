"""
Analyzes LGA+H2O trajectory to get g(r) of molecules in bins along x-axis
"""

import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt

# set input parameters
# -----------------------------------------------
path = 'Data/'
fgro = 'nvt1_LGA+H2O_0.5Agap_bigbox3x_finalShort_originalTRAJ.gro'
lgaMOL  = 256        # number of LGA molecules per frame
lgaATM  = 21         # number of LGA atoms per molecule
h2oMOL  = 8442       # number of H2O molecules per frame
h2oATM  = 3          # number of H2O atoms per molecule
frames  = 201        # number of frames in trajectory
dtstep  = 0.001      # time step between frames (ns)
xbox    = 20.34200   # box length in x-direction (nm)
ybox    =  5.59744   # box length in y-direction (nm)
zbox    =  3.19734   # box length in z-direction (nm)
rcut    = 1.5        # cutoff radius (nm)
dxbin   = 2          # width of bin centered around interface (nm)
numbin  = 3          # number of bins centered around interface
intfce  = 5.748      # interface distance from origin (nm)
grbins  = 60         # number of bins for g(r)
# -----------------------------------------------

t0 = datetime.now()

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

# re-order molecules in arrays per frame
resfL = np.chararray((frames,lgaMOL),itemsize=12)
xcmfL = np.zeros((frames,lgaMOL))
ycmfL = np.zeros((frames,lgaMOL))
zcmfL = np.zeros((frames,lgaMOL))
resfE = np.chararray((frames,h2oMOL),itemsize=12)
xcmfE = np.zeros((frames,h2oMOL))
ycmfE = np.zeros((frames,h2oMOL))
zcmfE = np.zeros((frames,h2oMOL))
indx = 0
# loop over frames
for i in range(frames):
    # loop over LGA molecules
    for j in range(lgaMOL):
        resfL[i,j] = res[indx]
        xcmfL[i,j] = pxcm[indx]
        ycmfL[i,j] = pycm[indx]
        zcmfL[i,j] = pzcm[indx]
        indx += 1
    # loop over H2O molecules
    for j in range(h2oMOL):
        resfE[i,j] = res[indx]
        xcmfE[i,j] = pxcm[indx]
        ycmfE[i,j] = pycm[indx]
        zcmfE[i,j] = pzcm[indx]
        indx += 1

# get bins along x
binstrt = intfce - 0.5*numbin*dxbin
binfnsh = intfce + 0.5*numbin*dxbin
xbins   = np.linspace(binstrt,binfnsh,numbin+1)

# get g(r) in each bin per frame per LGA-H2O
deltaR = rcut/grbins
grRbin = np.linspace(deltaR,rcut,grbins)
xhlf,yhlf,zhlf = 0.5*xbox,0.5*ybox,0.5*zbox
grLL = np.zeros((2,grbins))
grEE = np.zeros((2,grbins))
grLE = np.zeros(grbins)
grEL = np.zeros(grbins)
grnL = np.zeros((2,frames))
grnE = np.zeros((2,frames))

# loop over frames
for f in range(frames):
    print 'frame ',f+1,' of ',frames,' ...'

    # loop over LGA
    for i in range(lgaMOL-1):
        binNum = int((xcmfL[f,i]-binstrt)/dxbin)+1
        xposi = xcmfL[f,i]
        yposi = ycmfL[f,i]
        zposi = zcmfL[f,i]
        # get LGA-LGA g(r)
        if binNum == 1 or binNum == 2:
            grnL[binNum-1,f] += 1
            for j in range(i+1,lgaMOL):
                xposj = xcmfL[f,j]
                yposj = ycmfL[f,j]
                zposj = zcmfL[f,j]
                # get distances
                dx = xposi - xposj
                dy = yposi - yposj
                dz = zposi - zposj
                # minimum image
                if dx > xhlf: dx = dx-xbox
                elif dx < -xhlf: dx = dx + xbox
                if dy > yhlf: dy = dy-ybox
                elif dy < -yhlf: dy = dy + ybox
                if dz > zhlf: dz = dz-zbox
                elif dz < -zhlf: dz = dz + zbox
                r = (dx*dx + dy*dy + dz*dz)**0.5
                if r <= rcut:
                    ndx = int(r/deltaR)
                    grLL[binNum-1,ndx] += 2
        # get LGA-H2O g(r)
        if binNum == 2:
            for j in range(h2oMOL):
                xposj = xcmfE[f,j]
                yposj = ycmfE[f,j]
                zposj = zcmfE[f,j]
                # get distances
                dx = xposi - xposj
                dy = yposi - yposj
                dz = zposi - zposj
                # minimum image
                if dx > xhlf: dx = dx-xbox
                elif dx < -xhlf: dx = dx + xbox
                if dy > yhlf: dy = dy-ybox
                elif dy < -yhlf: dy = dy + ybox
                if dz > zhlf: dz = dz-zbox
                elif dz < -zhlf: dz = dz + zbox
                r = (dx*dx + dy*dy + dz*dz)**0.5
                if r <= rcut:
                    ndx = int(r/deltaR)
                    grLE[ndx] += 1

    # loop over H2O
    for i in range(h2oMOL-1):
        binNum = int((xcmfE[f,i]-binstrt)/dxbin)+1
        # get H2O-H2O g(r)
        if binNum == 2 or binNum == 3:
            grnE[binNum-2,f] += 1
            xposi = xcmfE[f,i]
            yposi = ycmfE[f,i]
            zposi = zcmfE[f,i]
            for j in range(i+1,h2oMOL):
                xposj = xcmfE[f,j]
                yposj = ycmfE[f,j]
                zposj = zcmfE[f,j]
                # get distances
                dx = xposi - xposj
                dy = yposi - yposj
                dz = zposi - zposj
                # minimum image
                if dx > xhlf: dx = dx-xbox
                elif dx < -xhlf: dx = dx + xbox
                if dy > yhlf: dy = dy-ybox
                elif dy < -yhlf: dy = dy + ybox
                if dz > zhlf: dz = dz-zbox
                elif dz < -zhlf: dz = dz + zbox
                r = (dx*dx + dy*dy + dz*dz)**0.5
                if r <= rcut:
                    ndx = int(r/deltaR)
                    grEE[binNum-2,ndx] += 2
         # get H2O-LGA g(r)
        if binNum == 3:
            for j in range(lgaMOL):
                xposj = xcmfL[f,j]
                yposj = ycmfL[f,j]
                zposj = zcmfL[f,j]
                # get distances
                dx = xposi - xposj
                dy = yposi - yposj
                dz = zposi - zposj
                # minimum image
                if dx > xhlf: dx = dx-xbox
                elif dx < -xhlf: dx = dx + xbox
                if dy > yhlf: dy = dy-ybox
                elif dy < -yhlf: dy = dy + ybox
                if dz > zhlf: dz = dz-zbox
                elif dz < -zhlf: dz = dz + zbox
                r = (dx*dx + dy*dy + dz*dz)**0.5
                if r <= rcut:
                    ndx = int(r/deltaR)
                    grEL[ndx] += 1

# get average mols in bins
nL1,nL2 = np.average(grnL[0]),np.average(grnL[1])
nE2,nE3 = np.average(grnE[0]),np.average(grnE[1])
# get density in bins
rho1    = nL1/(dxbin*ybox*zbox)
rho2    = (nL2+nE2)/(dxbin*ybox*zbox)
rho3    = nE3/(dxbin*ybox*zbox)
# normalize g(r)
grLL1   = grLL[0]/(2*np.pi*frames*nL1*rho1*grRbin*grRbin*deltaR)
grLL2   = grLL[1]/(2*np.pi*frames*nL2*0.2*rho2*grRbin*grRbin*deltaR)
grLE2   = grLE/(2*np.pi*frames*nL2*rho2*grRbin*grRbin*deltaR)
grEL2   = grEL/(2*np.pi*frames*nE2*rho2*grRbin*grRbin*deltaR)
grEE2   = grEE[0]/(2*np.pi*frames*nE2*2*rho2*grRbin*grRbin*deltaR)
grEE3   = grEE[1]/(2*np.pi*frames*nE3*2*rho3*grRbin*grRbin*deltaR)


# make plots
print 'making plots ...'
# set plot parameters
plt.ioff()
cnt = 0
fontsize = 14
factor   = 1.20
# plot bin image
plt.figure(figsize=(8,5))
plt.title('Bins',fontsize=20)
plt.xlim(0,xbox)
plt.yticks([])
plt.xlabel('x-axis (nm)',fontsize=16)
plt.vlines(xbins,0,1)
plt.vlines(intfce,0,1,'m',lw=2,linestyle='dashed')
for i in range(numbin):
    plt.text(xbins[i]+0.20*dxbin,0.5,str(i+1),fontsize=15,color='b')
plt.savefig(path+'Plots/0binsradial.png')
plt.close()
# plot LGA-LGA g(r) in bin1
plt.figure()
plt.plot(grRbin,grLL1,'b-',label='LGA-LGA')
plt.xlabel('r (nm)',fontsize=fontsize)
plt.ylabel('g (r)',fontsize=fontsize)
plt.legend(loc=2)
plt.annotate('Bin = 1',fontsize=18,xy=(0.78,0.92),xycoords='figure fraction')
plt.savefig(path+'Plots/1radial.png')
plt.close()
# plot LGA-LGA g(r) in bin2
plt.figure()
plt.subplot(221)
plt.plot(grRbin,grLL2,'b-',label='LGA-LGA')
plt.ylabel('g (r)',fontsize=fontsize)
plt.legend(loc=1)
# plot LGA-H2O g(r) in bin2
plt.subplot(223)
plt.plot(grRbin,grLE2,'b-',label='LGA-H2O')
plt.xlabel('r (nm)',fontsize=fontsize)
plt.ylabel('g (r)',fontsize=fontsize)
plt.legend(loc=1)
# plot H2O-LGA g(r) in bin2
plt.subplot(222)
plt.plot(grRbin,grEL2,'g-',label='H2O-LGA')
plt.legend(loc=1)
# plot H2O-H2O g(r) in bin2
plt.subplot(224)
plt.plot(grRbin,grEE2,'g-',label='H2O-H2O')
plt.xlabel('r (nm)',fontsize=fontsize)
plt.ylabel('g (r)',fontsize=fontsize)
plt.legend(loc=1)
plt.annotate('Bin = 2',fontsize=18,xy=(0.78,0.92),xycoords='figure fraction')
plt.savefig(path+'Plots/2radial.png')
plt.close()
# plot H2O-H2O g(r) in bin3
plt.figure()
plt.plot(grRbin,grEE3,'g-',label='H2O-H2O')
plt.xlabel('r (nm)',fontsize=fontsize)
plt.ylabel('g (r)',fontsize=fontsize)
plt.legend(loc=1)
plt.annotate('Bin = 3',fontsize=18,xy=(0.78,0.92),xycoords='figure fraction')
plt.savefig(path+'Plots/3radial.png')
plt.close()


# print done
print 4,'images saved'
print 'done !\n'
plt.ion()

# print runtime
t1 = datetime.now()
print ''
print 'elapsed runtime:',t1-t0