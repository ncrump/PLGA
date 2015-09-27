"""
Generates trimmed system of atom coordinates from input .gro file
Removes molecules outside specified limits and outputs new .gro file
This is needed for merging two systems of different box sizes
"""

import numpy as np

# set input variables
# ---------------------------------------------------------------------------
infile = '2finalnpt_EA_resized.gro'
outfile = '3finalnpt_EA_trimmed.gro'
resname = 'EEE'
atms = 6      # number of atoms per molecule
posx = 0.0    # trim molecules by this much at +x (nm)
negx = 0.0    # trim molecules by this much at -x (nm)
posy = 0.48256    # trim molecules by this much at +y (nm)
negy = 0.0    # trim molecules by this much at -y (nm)
posz = 1.61966    # trim molecules by this much at +z (nm)
negz = 0.0    # trim molecules by this much at -z (nm)
# ---------------------------------------------------------------------------

# read in .gro file
res,atm,num,x,y,z = np.genfromtxt(infile,skiprows=2,skip_footer=1,dtype=str,unpack=True)
xbox,ybox,zbox = np.genfromtxt(infile,skiprows=int(num[-1])+2,usecols=(0,1,2),dtype=str,unpack=True)

# convert coords to floats
x = np.array([float(i) for i in x])
y = np.array([float(i) for i in y])
z = np.array([float(i) for i in z])

# get system box boundary
xlim = [min(x),max(x)]
ylim = [min(y),max(y)]
zlim = [min(z),max(z)]

# trim system
# find molecules outside limits
N = int(num[-1])
resTrim = []
for i in range(N):
    xi,yi,zi = x[i],y[i],z[i]
    if xi > xlim[1]-posx or xi < xlim[0]+negx: resTrim.append(res[i])
    if yi > ylim[1]-posy or yi < ylim[0]+negy: resTrim.append(res[i])
    if zi > zlim[1]-posz or zi < zlim[0]+negz: resTrim.append(res[i])
resTrim = list(set(resTrim))

# get index of molecules
indx = []
for i in resTrim:
    indx.append(np.where(res==i))
indx = np.array(indx).flatten()

# remove molecules
atmNew = np.delete(atm,indx)
xNew = np.delete(x,indx)
yNew = np.delete(y,indx)
zNew = np.delete(z,indx)

# generate sequential residue and atom lists
Nnew = len(xNew)
resNum = int(Nnew/atms)
numNew = range(1,Nnew+1)
resNew = range(1,resNum+1)*atms
resNew.sort()

# generate new .gro file of trimmed system
# open and write .gro header
f = open(outfile,'w')
f.write(resname+' TRIMMED\n')
f.write(str(Nnew)+'\n')

# write .gro body and atom coordinates
for i in range(Nnew):
    f.write('%8s%2s%5s%5s%8s%8s%8s\n' %
           (str(resNew[i])+resname,'',atmNew[i],numNew[i],xNew[i],yNew[i],zNew[i]))
f.write('%6s %6s %6s' % (xbox,ybox,zbox))
f.close()

# print output to screen
print '\n'
print 'new file generated:', outfile