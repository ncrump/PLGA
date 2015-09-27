"""
Generates LGA sandwiched between EA systems from input .gro files
Combines sandwiched system into single .gro file
Applies a shift to the atom numbers and chosen coordinate axis
Applies a shift to the combined system to center in new box
"""

import numpy as np

# set input variables
# ---------------------------------------------------------------------------
LGAfile = 'finalMD_noPBC_LGA_trimmed.gro'
EAfile =  'finalMD_noPBC_EA_trimmed.gro'
outfile = 'newbox_lga_ea_trimmed_1Asandwich_15Apad.gro'
pad = 1.50  # pading along box edges  (nm)
dx  = 0.10  # sandwich gap size on x sides (nm)
dy  = 0.00  # sandwich gap size on y sides (nm)
dz  = 0.00  # sandwich gap size on z sides (nm)
# ---------------------------------------------------------------------------

# read in .gro files
res1,atm1,num1,x1,y1,z1 = np.genfromtxt(LGAfile,skiprows=2,skip_footer=1,dtype=str,unpack=True)
res2,atm2,num2,x2,y2,z2 = np.genfromtxt(EAfile,skiprows=2,skip_footer=1,dtype=str,unpack=True)
res3,atm3,num3,x3,y3,z3 = np.genfromtxt(EAfile,skiprows=2,skip_footer=1,dtype=str,unpack=True)

# shift atom numbers to be sequential
N1 = len(num1)
N2 = len(num2)
N3 = len(num3)
num2 = [str(int(i)+int(num1[-1])) for i in num2]
num3 = [str(int(i)+int(num2[-1])) for i in num3]
Ntot = num3[-1]

# shift residue numbers to be sequential
nres = int(len(res2)/6)
res3 = []
for i in range(nres):
    for j in range(6):
        res3.append(str(nres+i+1)+'EEE')

# convert strings to floats
x1 = np.array([float(i) for i in x1])
y1 = np.array([float(i) for i in y1])
z1 = np.array([float(i) for i in z1])
x2 = np.array([float(i) for i in x2])
y2 = np.array([float(i) for i in y2])
z2 = np.array([float(i) for i in z2])
x3 = np.array([float(i) for i in x3])
y3 = np.array([float(i) for i in y3])
z3 = np.array([float(i) for i in z3])

# get distances for sandwich gaps
if dx != 0:
    distx1 = abs(min(x2)) + max(x1)
    distx2 = abs(min(x1)) + max(x3)
else:
    distx1 = 0
    distx2 = 0
if dy != 0:
    disty1 = abs(min(y2)) + max(y1)
    disty2 = abs(min(y1)) + max(y3)
else:
    disty1 = 0
    disty2 = 0
if dz != 0:
    distz1 = abs(min(z2)) + max(z1)
    distz2 = abs(min(y1)) + max(y3)
else:
    distz1 = 0
    distz2 = 0

# shift atom coords to make sandwich gaps
x2 = x2+distx1+dx
y2 = y2+disty1+dy
z2 = z2+distz1+dz
x3 = x3-distx2-dx
y3 = y3-disty2-dy
z3 = z3-distz2-dz

# calculate new box vectors and box center
xmax = max(np.concatenate((x1,x2,x3)))
xmin = min(np.concatenate((x1,x2,x3)))
ymax = max(np.concatenate((y1,y2,y3)))
ymin = min(np.concatenate((y1,y2,y3)))
zmax = max(np.concatenate((z1,z2,z3)))
zmin = min(np.concatenate((z1,z2,z3)))
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
x3 = x3+xshft
y3 = y3+yshft
z3 = z3+zshft

# calculate gap between outermost molecules
xsep1 = min(x2) - max(x1)
ysep1 = min(y2) - max(y1)
zsep1 = min(z2) - max(z1)
xsep2 = min(x1) - max(x3)
ysep2 = min(y1) - max(y3)
zsep2 = min(z1) - max(z3)
if xsep1 < 0: xsep1 = 0
if ysep1 < 0: ysep1 = 0
if zsep1 < 0: zsep1 = 0
if xsep2 < 0: xsep2 = 0
if ysep2 < 0: ysep2 = 0
if zsep2 < 0: zsep2 = 0

# calculate points in interface planes
if dx != 0:
    xpt1 = max(x1)+0.5*dx
    xpt2 = min(x1)-0.5*dx
else:
    xpt1 = 0
    xpt2 = 0
if dy != 0:
    ypt1 = max(y1)+0.5*dy
    ypt2 = min(y1)-0.5*dx
else:
    ypt1 = 0
    ypt2 = 0
if dz != 0:
    zpt1 = max(z1)+0.5*dz
    zpt2 = min(z1)-0.5*dx
else:
    zpt1 = 0
    zpt2 = 0

# generate new .gro file of merged system
# open and write .gro header
f = open(outfile,'w')
f.write('LGA+EA SANDWICH\n')
f.write(Ntot+'\n')

# write .gro body and atom coordinates
for i in range(N1):
    f.write('%8s%2s%5s%5s%8s%8s%8s\n' %
           (res1[i],'',atm1[i],num1[i],str(x1[i]),str(y1[i]),str(z1[i])))
for j in range(N2):
    f.write('%8s%2s%5s%5s%8s%8s%8s\n' %
           (res2[j],'',atm2[j],num2[j],str(x2[j]),str(y2[j]),str(z2[j])))
for k in range(N3):
    f.write('%8s%2s%5s%5s%8s%8s%8s\n' %
           (res3[k],'',atm3[k],num3[k],str(x3[k]),str(y3[k]),str(z3[k])))
f.write('%6s %6s %6s' % (str(xbox),str(ybox),str(zbox)))
f.close()

# print output to screen
print '\n'
print 'sandwich gap sizes: (%6.3f %6.3f %6.3f) (%6.3f %6.3f %6.3f)' % (xsep1,ysep1,zsep1,xsep2,ysep2,zsep2), 'nm'
print '   sandwich points: (%6.3f %6.3f %6.3f) (%6.3f %6.3f %6.3f)' % (xpt1,ypt1,zpt1,xpt2,ypt2,zpt2), 'nm'
print ' padding box edges: %6.3f %6.3f %6.3f' % (pad,pad,pad), 'nm'
print '   new box vectors: %6.3f %6.3f %6.3f' % (xbox,ybox,zbox), 'nm'
print ' new system center: %6.3f %6.3f %6.3f' % (xctr,yctr,zctr), 'nm'
print 'new file generated:', outfile