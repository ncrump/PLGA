import numpy as np

# generates input text file to pdb2gmx scripted menu options

size = 2 + 2*540           # 2+2*(molecules in latttice)
filename = '1620EEE.txt'    # output file

x = np.zeros(size,dtype=int)

x[0] = 14                  # option 14 = GROMOS96 54a7 ff
x[1] = 3                   # option 3  = none for water model
x[2:size] = 2              # option 2  = no termini

np.savetxt(filename,x,fmt='%i')