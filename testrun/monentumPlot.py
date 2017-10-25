#!/usr/bin/env python
#============================
#This code is to plot the energy spread
#Zhicong@ 03 May 2017
#plots are saved at '/post'

import os, sys
import matplotlib.pyplot as plt

def plot(arg):
    lineType = ['-','--','-.',':']
    picNum = len(arg) - 1
    plotPath = './post'
    if os.path.exists(plotPath) == False:
        os.makedirs(plotPath)
    for i in range(1,picNum+1):
        print("ploting From: " + arg[i])
        try:
            fin = open(arg[i],'r')
        except:
            print( "  ERRPR! Can't open file '" + arg[i] + "'")
            exit(-1)

        linesList  = fin.readlines()
        fin .close()
        linesList  = [line.split() for line in linesList ]
        x   = [float(xrt[0]) for xrt in linesList]
        y   = [float(xrt[4]) for xrt in linesList]
        plt.plot(x, y, linestyle=lineType[(i-1)%4], linewidth=2,label=arg[i])
        #print((i-1)%5,lineType)
    plt.xlabel('Z distance (m)')
    plt.ylabel('Rms momentum at Z (MeV) ')
    #plt.tick_params(axis="both", labelsize=48)
    plt.legend()


    plt.savefig(plotPath+'/'+'energySpread.pdf')
    plt.show()

    
if __name__ == '__main__':
    if len(sys.argv) == 1:
        arg=['ct','fort.26']
    else:
        arg= sys.argv
    plot(arg)
