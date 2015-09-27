import numpy as np
import matplotlib.pyplot as plt

# read data file
t,PE,KE,E,T,P = np.loadtxt('nvt0_15Apad_energy.xvg',skiprows=23,unpack=True)

# calculate running averages
PEave,KEave,Eave          = PE[0],KE[0],E[0]
Tave,Pave                 = T[0],P[0]
PEaveArr,KEaveArr,EaveArr = [PEave],[KEave],[Eave]
TaveArr,PaveArr           = [Tave],[Pave]
indx = len(t)
for i in range(1,indx):
    PEave = PEave + (PE[i] - PEave)/i
    KEave = KEave + (KE[i] - KEave)/i
    Eave  = Eave  +  (E[i] -  Eave)/i
    Tave  = Tave  +  (T[i] -  Tave)/i
    Pave  = Pave  +  (P[i] -  Pave)/i
    PEaveArr.append(PEave)
    KEaveArr.append(KEave)
    EaveArr.append(Eave)
    TaveArr.append(Tave)
    PaveArr.append(Pave)

# plot energies
plt.figure()
plt.subplot(3,1,1)
plt.plot(t,PE,'b-')
plt.plot(t,PEaveArr,'r-')
plt.ylabel('$<PE>$ (kJ/mol)')
plt.subplot(3,1,2)
plt.plot(t,KE,'b-')
plt.plot(t,KEaveArr,'r-')
plt.ylabel('$<KE>$ (kJ/mol)')
plt.subplot(3,1,3)
plt.plot(t,E,'b-')
plt.plot(t,EaveArr,'r-')
plt.ticklabel_format(style='plain',useOffset=False,axis='both')
plt.ylabel('$<E>$ (kJ/mol)')
plt.xlabel('Time (ps)')

# plot temp/pres
plt.figure()
plt.subplot(2,1,1)
plt.plot(t,T,'b-')
plt.plot(t,TaveArr,'r-')
plt.ylabel('$<T>$ (Kelvin)')
plt.subplot(2,1,2)
plt.plot(t,P,'b-')
plt.plot(t,PaveArr,'r-')
plt.ylabel('$<P>$ (bar)')
plt.xlabel('Time (ps)')