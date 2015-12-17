"""
Analyzes LGA+H2O trajectory to count number of LGA that cross interface
"""

import numpy as np

# set input parameters
# -----------------------------------------------
path = 'Data/'
fgro = 'nvt1_LGA+H2O_0.5Agap_bigbox3x_analysis_originalTRAJ.gro'
lgaMOL = 256      # number of LGA molecules per frame
lgaATM = 21       # number of LGA atoms per molecule
h2oMOL = 8442     # number of H2O molecules per frame
h2oATM = 3        # number of H2O atoms per molecule
frames = 101      # number of frames in trajectory
dtstep = 0.1      # time step between frames (ns)
intfce = 5.748    # interface distance from origin (nm)
# -----------------------------------------------

# define dictionary of atom masses for LGA
Mtot = 278.2166
mass = {'C0': 14.0270, 'H2': 01.0080, 'C7': 12.0110, 'O8': 15.9994, 'H1': 1.0080,\
        'OB': 15.9994, 'OA': 15.9994, 'C6': 13.0190, 'O7': 15.9994, 'C9': 12.011,\
        'C8': 14.0270, 'O9': 15.9994, 'O6': 15.9994, 'C3': 13.0190, 'C2': 15.035,\
        'O5': 15.9994, 'O4': 15.9994, 'O3': 15.9994, 'C1': 12.0110, 'C5': 15.035,\
        'C4': 12.0110}

# read trajectory file (fixed column width)
print '\nreading file ...'
delimiter = [8,7,5,8,8,8,8,8,8]
nres0,natm0,num0,xatm0,yatm0,zatm0,vxatm0,vyatm0,vzatm0 = np.genfromtxt(path+fgro,dtype=str,delimiter=delimiter,unpack=True)

# format arrays to keep only LGA
print 'analyzing trajectory ...'
Ntot = len(nres0)
indx = np.arange(2,Ntot,lgaMOL*lgaATM+h2oMOL*h2oATM+3)
nres,natm,xatm,yatm,zatm = [],[],[],[],[]
# get LGA molecules only
for i in indx:
    nres.append(nres0[i:i+lgaMOL*lgaATM])
    natm.append(natm0[i:i+lgaMOL*lgaATM])
    xatm.append(xatm0[i:i+lgaMOL*lgaATM])
    yatm.append(yatm0[i:i+lgaMOL*lgaATM])
    zatm.append(zatm0[i:i+lgaMOL*lgaATM])
# reformat arrays
nres = np.array(nres).flatten().astype(str)
natm = np.array(natm).flatten().astype(str)
nres = np.char.strip(nres)
natm = np.char.strip(natm)
xatm = np.array(xatm).flatten().astype(float)
yatm = np.array(yatm).flatten().astype(float)
zatm = np.array(zatm).flatten().astype(float)
t    = np.linspace(0,(frames-1)*dtstep,frames)


# calculate center mass of each molecule per frame
frm,res,xcm,ycm,zcm = [],[],[],[],[]
indx = 0
for i in range(frames):
    for j in range(lgaMOL):
        xcom,ycom,zcom = 0,0,0
        for k in range(lgaATM):
            xcom += mass[natm[indx]]*xatm[indx]
            ycom += mass[natm[indx]]*yatm[indx]
            zcom += mass[natm[indx]]*zatm[indx]
            indx += 1
        frm.append(t[i])
        res.append(nres[indx-1])
        xcm.append(xcom/Mtot)
        ycm.append(ycom/Mtot)
        zcm.append(zcom/Mtot)

# count number of LGA on EA side of interface per frame
cnt,resList,posList = [],[],[]
indx = 0
for i in range(frames):
    count = 0
    resTemp = []
    posTemp = []
    for j in range(lgaMOL):
        if xcm[indx] > intfce:
            count += 1
            resTemp.append(res[indx])
            posTemp.append(round(xcm[indx],3))
        indx += 1
    cnt.append(count)
    resList.append(resTemp)
    posList.append(posTemp)

# print results
print 'printing results ...\n'
print '%5s %5s %3s %10s %10s' % ('frame','time','count','molecule','x-com')
for i in range(frames):
    print '%5i %5.1f %3i %10s %10s' % (i,t[i],cnt[i],resList[i],posList[i])
print 'average LGA across interface:'
print np.average(cnt)