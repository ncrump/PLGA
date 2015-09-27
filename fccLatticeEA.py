"""
Generates FCC lattice of Ethyl Acetate molecules from input pdb
Places center of molecule at face centered orthorhombic lattice sites
"""

import numpy as np

# set input variables
# ---------------------------------------------------------------------------
filename = '1eee_opt_unitedatom'
nD = 0.006  # number density (molecules/cu angstrom)
pD = 0.90   # percentage of number density to use
dx = 7      # x distance between molecules (angstrom)
dy = 8      # y distance between molecules (angstrom)
dz = 8      # z distance between molecules (angstrom)
tB = 0      # number of times to translate unit cell
# ---------------------------------------------------------------------------

# get atom coords and info from input pdb file
ax,ay,az = np.loadtxt(filename+'.pdb',comments='END',skiprows=4,usecols=(5,6,7),unpack=True)
name = np.genfromtxt(filename+'.pdb',comments='END',skiprows=4,usecols=2,dtype=str,unpack=True)
res = np.genfromtxt(filename+'.pdb',comments='END',skiprows=4,usecols=3,dtype=str,unpack=True)
atm = np.genfromtxt(filename+'.pdb',comments='END',skiprows=4,usecols=10,dtype=str,unpack=True)

# get molecule parameters
natoms = len(ax)           # atoms per molecule
xL = max(ax)-min(ax)       # x molecule length
yL = max(ay)-min(ay)       # y molecule length
zL = max(az)-min(az)       # z molecule length

# get lattice parameters
N = 4*(tB+1)**3            # total number of molecules
nD = nD*pD                 # desired lattice density
a = 0.5*xL + dx            # x unit cell dimension
b = yL + dy                # y unit cell dimension
c = zL + dz                # z unit cell dimension
Lx = a*(tB+1)              # x lattice dimension
Ly = b*(tB+1)              # y lattice dimension
Lz = c*(tB+1)              # z lattice dimension
V = Lx*Ly*Lz               # total lattice volume
nDt = N/V                  # actual lattice density

# get initial unit cell positions
x = np.array([ax,ax+0.5*a,ax,ax+0.5*a]).flatten()
y = np.array([ay,ay+0.5*b,ay+0.5*b,ay]).flatten()
z = np.array([az,az,az+0.5*c,az+0.5*c]).flatten()

# get translation step of cube
xstep = np.linspace(a, a*tB, tB)
ystep = np.linspace(b, b*tB, tB)
zstep = np.linspace(c, c*tB, tB)

# translate in x direction
for xi in xstep:
    x = np.append(x,x[0:natoms*4]+xi)
    y = np.append(y,y[0:natoms*4])
    z = np.append(z,z[0:natoms*4])

# translate in y direction
for yi in ystep:
    x = np.append(x,x[0:natoms*4*(tB+1)])
    y = np.append(y,y[0:natoms*4*(tB+1)]+yi)
    z = np.append(z,z[0:natoms*4*(tB+1)])

# translate in z direction
for zi in zstep:
    x = np.append(x,x[0:natoms*4*(tB+1)**2])
    y = np.append(y,y[0:natoms*4*(tB+1)**2])
    z = np.append(z,z[0:natoms*4*(tB+1)**2]+zi)

# calculate atom-to-atom parameters
lx = max(x)-min(x)         # atom-atom x dimension
ly = max(y)-min(y)         # atom-atom y dimension
lz = max(z)-min(z)         # atom-atom z dimension
v = lx*ly*lz               # atom-atom volume
nDc = N/v                  # atom-atom density

# convert floats to strings for write
x = [np.str(i) for i in x]
y = [np.str(i) for i in y]
z = [np.str(i) for i in z]

# generate new pdb file with lattice positions
# open and write pdb header
fname = '1eee_'+str(N)+'Lattice.pdb'
f = open(fname,'w')
f.write('HEADER    PYTHON GENERATED PDB LATTICE DDMMMYYYY\n')
f.write('TITLE     1EEE\n')
f.write('COMPND    ETHYL ACETATE\n')
f.write('AUTHOR    N CRUMP\n')

# write pdb body and atom coordinates
indx = 0
anum = 0
for i in range(N):
    for j in range(natoms+1):
        anum += 1
        if j != natoms:
            indx += 1
            f.write('%-4s %6s %3s %4s %5s %11s %7s %7s %5s %5s %11s\n' %
                   ('ATOM',str(anum),name[j],res[j],str(i+1),x[indx-1],y[indx-1],z[indx-1],'1.00','0.00',atm[j]))
        if j == natoms:
            f.write('%-3s %7s %8s %5s\n' %
                   ('TER',str(anum),res[j-1],str(i+1)))

# write end line and close pdb
f.write('END')
f.close()

# print output to screen
print '\n'
print 'file generated:',fname
print N,'chains', '/', natoms,'atoms per chain', '/', natoms*N,'atoms'
print '\n'
print 'lattice parameters (angstroms)'
print '----------------------------------------------------------'
print 'input density = ', nD
print 'input distances = %i, %i, %i' % (dx,dy,dz)
print 'lattice dims: Lx Ly Lz V nDensity = %1.3f, %1.3f, %1.3f, %1.3f, %1.8f' % (Lx,Ly,Lz,V,nDt)
print 'atom-to-atom: Lx Ly Lz V nDensity = %1.3f, %1.3f, %1.3f, %1.3f, %1.8f' % (lx,ly,lz,v,nDc)