"""
Reads LGA+EA trajectory file
"""

import numpy as np


# set input parameters
# -----------------------------------------------
fgro = 'nvt0_LGA+EA_0.5Agap_eabigbox3x_analysis_reimagedDUMP.gro'
# -----------------------------------------------

# read trajectory file (fixed column width)
delimiter = [8,7,5,8,8,8,8,8,8]
res,atm,num,x,y,z,vx,vy,vz = np.genfromtxt(fgro,dtype=str,delimiter=delimiter,skip_header=2,skip_footer=1,unpack=True)