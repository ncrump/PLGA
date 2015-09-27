"""
Generates merged system of atom coordinates from input .gro files
Combines equilibrated LGA and EA systems into single .gro file
Applies a shift to the atom numbers and chosen coordinate axis
Applies a shift to the combined system to center in new box
"""

import numpy as np

# set input variables
# ---------------------------------------------------------------------------
LGAfile = '1finalnvt_LGA_reimaged.gro'
EAfile =  '1finalnvt_EA_reimaged.gro'
outfile = '1startnvt_LGA+EA_1Agap.gro'
pad = 0.20  # pading along box edges  (nm)
dx  = 0.10  # interface gap size in x (nm)
dy  = 0.00  # interface gap size in y (nm)
dz  = 0.00  # interface gap size in z (nm)
# ---------------------------------------------------------------------------

# read in .gro files
res1,atm1,num1,x1,y1,z1 = np.genfromtxt(LGAfile,skiprows=2,skip_footer=1,dtype=str,unpack=True)
res2,atm2,num2,x2,y2,z2 = np.genfromtxt(EAfile,skiprows=2,skip_footer=1,dtype=str,unpack=True)

# shift atom numbers to be sequential
N1 = len(num1)
N2 = len(num2)
num2 = [str(int(i)+int(num1[-1])) for i in num2]
Ntot = num2[-1]

# convert strings to floats
x1 = np.array([float(i) for i in x1])
y1 = np.array([float(i) for i in y1])
z1 = np.array([float(i) for i in z1])
x2 = np.array([float(i) for i in x2])
y2 = np.array([float(i) for i in y2])
z2 = np.array([float(i) for i in z2])

# get distance for interface gap
if dx != 0: distx = abs(min(x2)) + max(x1)
else:       distx = 0
if dy != 0: disty = abs(min(y2)) + max(y1)
else:       disty = 0
if dz != 0: distz = abs(min(z2)) + max(z1)
else:       distz = 0

# shift atom coords to make interface gap
x2 = x2+distx+dx
y2 = y2+disty+dy
z2 = z2+distz+dz

# calculate new box vectors and box center
xmax = max(np.concatenate((x1,x2)))
xmin = min(np.concatenate((x1,x2)))
ymax = max(np.concatenate((y1,y2)))
ymin = min(np.concatenate((y1,y2)))
zmax = max(np.concatenate((z1,z2)))
zmin = min(np.concatenate((z1,z2)))
xbox = round(xmax-xmin+2*pad,3)
ybox = round(ymax-ymin+2*pad,3)
zbox = round(zmax-zmin+2*pad,3)
xctr = xbox/2
yctr = ybox/2
zctr = zbox/2

# shift system to center in new box
xshft = round(xbox-pad-xmax,3)
yshft = round(ybox-pad-ymax,3)
zshft = round(zbox-pad-zmax,3)
x1 = x1+xshft
y1 = y1+yshft
z1 = z1+zshft
x2 = x2+xshft
y2 = y2+yshft
z2 = z2+zshft

# calculate gap between outermost molecules
xsep = min(x2) - max(x1)
ysep = min(y2) - max(y1)
zsep = min(z2) - max(z1)
if xsep < 0: xsep = 0
if ysep < 0: ysep = 0
if zsep < 0: zsep = 0

# calculate point in interface plane
if dx != 0: xpt = max(x1)+0.5*dx
else:       xpt = 0
if dy != 0: ypt = max(y1)+0.5*dy
else:       ypt = 0
if dz != 0: zpt = max(z1)+0.5*dz
else:       zpt = 0

# generate new .gro file of merged system
# open and write .gro header
f = open(outfile,'w')
f.write('LGA+EA MERGED\n')
f.write(Ntot+'\n')

# write .gro body and atom coordinates
for i in range(N1):
    f.write('%8s%2s%5s%5s%8s%8s%8s\n' %
           (res1[i],'',atm1[i],num1[i],str(x1[i]),str(y1[i]),str(z1[i])))
for j in range(N2):
    f.write('%8s%2s%5s%5s%8s%8s%8s\n' %
           (res2[j],'',atm2[j],num2[j],str(x2[j]),str(y2[j]),str(z2[j])))
f.write('%7s %7s %7s' % (str(xbox),str(ybox),str(zbox)))
f.close()

# print output to screen
print '\n'
print '          gap size: %6.3f %6.3f %6.3f' % (xsep,ysep,zsep), 'nm'
print '   interface point: %6.3f %6.3f %6.3f' % (xpt,ypt,zpt), 'nm'
print ' padding box edges: %6.3f %6.3f %6.3f' % (pad,pad,pad), 'nm'
print '   new box vectors: %6.3f %6.3f %6.3f' % (xbox,ybox,zbox), 'nm'
print ' new system center: %6.3f %6.3f %6.3f' % (xctr,yctr,zctr), 'nm'
print 'new file generated:', outfile