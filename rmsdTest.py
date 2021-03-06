import numpy as np

x1 = np.array([0.695,0.578,0.456,0.340,0.334,0.220])
y1 = np.array([0.236,0.332,0.258,0.328,0.451,0.235])
z1 = np.array([0.202,0.201,0.202,0.201,0.200,0.202])

x2 = np.array([0.696,0.581,0.456,0.339,0.328,0.221])
y2 = np.array([0.232,0.332,0.260,0.330,0.452,0.233])
z2 = np.array([0.202,0.202,0.202,0.201,0.200,0.202])

rmsd = ((sum((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2))/6.0)**0.5
print rmsd