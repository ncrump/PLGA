import numpy as np
import matplotlib.pyplot as plt

# reads and plots xvg data file generated by GROMACS g_energy

def plotXVG(filename,avg=0):

    # read file to get header info
    f = open(filename)
    header = [f.next().strip() for i in range(19)]
    f.close()

    # parse header info for plot labels
    hxlabel = header[9][19:-1]
    hylabel = header[10][19:-1]
    hlabel = header[18][13:-1]

    # read file to get data values
    s = np.genfromtxt(filename,dtype=str,skiprows=19,unpack=True)

    # get data arrays
    N = len(s[0])
    x = [float(s[0][i]) for i in range(N)]
    y = [float(s[1][i]) for i in range(N)]

    # plot data
    plt.figure(figsize=(9,6),dpi=100)
    plt.plot(x,y,'b',label=hlabel)
    plt.title('Gromacs Plot: ' + hlabel)
    plt.xlabel(hxlabel)
    plt.ylabel(hlabel + ' ' + hylabel)
    plt.legend()

    # plot running average
    if avg == 1:
        xx = x[1::]
        yy = []
        aveRun = y[0]
        for i in range(0,N-1):
            aveRun = aveRun + (y[i+1]-aveRun)/(i+1)
            yy.append(aveRun)
        plt.plot(xx,yy,'r')
        print 'Average',hlabel,'=',aveRun

    # save file
    plt.savefig(filename[:-4]+'.png')