"""
Generates merged system of atom coordinates from input .gro files
Combines equilibrated LGA and EA systems into single .gro file
Applies a shift to the atom numbers and chosen coordinate axis
"""

import numpy as np

# set input variables
# ---------------------------------------------------------------------------
LGAfile = '1finalnvt_LGA_reimaged.gro'
EAfile =  '1finalnvt_EA_reimaged.gro'
outfile = '1startnvt_LGA+EA_n0.5Agap.gro'
dx   = -0.05     # interface gap size in x (nm)
xbox = 4.89938  # box size of LGA in x (nm)
ybox = 5.59744  # box size of LGA in y (nm)
zbox = 3.19734  # box size of LGA in z (nm)
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
distx = abs(min(x2)) + max(x1)

# shift atom coords to make interface gap
x2 = x2+distx+dx

# calculate shift for edge spacing
diff = distx-xbox
xshft = 0.8*diff

# calculate new box vectors and box center
#xbox = xbox + diff + 0.5*dx + xshft + 14.26400
xbox = max(x2)-min(x1)
xctr = xbox/2
yctr = ybox/2
zctr = zbox/2

# shift system to center in new box
x1 = x1+xshft
x2 = x2+xshft

# calculate gap between outermost molecules
xsep = min(x2) - max(x1)
ysep = min(y2) - max(y1)
zsep = min(z2) - max(z1)
if xsep < 0: xsep = 0
if ysep < 0: ysep = 0
if zsep < 0: zsep = 0

# calculate point in interface plane
xpt = max(x1)+0.5*dx
ypt = yctr
zpt = zctr

# generate new .gro file of merged system
# open and write .gro header
f = open(outfile,'w')
f.write('LGA+EA MERGED\n')
f.write(Ntot+'\n')

# write .gro body and atom coordinates
for i in range(N1):
    f.write('%8s%2s%5s%5s%8.3f%8.3f%8.3f\n' %
           (res1[i],'',atm1[i],num1[i],x1[i],y1[i],z1[i]))
for j in range(N2):
    f.write('%8s%2s%5s%5s%8.3f%8.3f%8.3f\n' %
           (res2[j],'',atm2[j],num2[j],x2[j],y2[j],z2[j]))
f.write('%8.5f %8.5f %8.5f' % (xbox,ybox,zbox))
f.close()

# print output to screen
print '\n'
print '          gap size: %6.3f %6.3f %6.3f' % (xsep,ysep,zsep), 'nm'
print '   interface point: %6.3f %6.3f %6.3f' % (xpt,ypt,zpt), 'nm'
print '   new box vectors: %6.3f %6.3f %6.3f' % (xbox,ybox,zbox), 'nm'
print ' new system center: %6.3f %6.3f %6.3f' % (xctr,yctr,zctr), 'nm'
print 'new file generated:', outfile