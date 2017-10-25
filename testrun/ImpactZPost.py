#!/usr/bin/env python
#============================
#This code is to plot the result from ImpactZ
#Zhicong@21/10/2016
#Input : fort.xx
#Output: figures about beam size and emittance
# plots are saved at '/post'

import os
import matplotlib.pyplot as plt
import warnings


#===========Flag===========#
showFlag = True
saveFlag = True
saveFormat = 'pdf'
#===========Flag===========#

picNum = 4
fileList    = [[]*2]*picNum
saveName    = []
labelList   = [[]*2]*picNum
xdataList   = [[]*2]*picNum
ydataList   = [[]*2]*picNum
xyLabelList = [[]*2]*picNum

saveName.append('sizeX')
fileList[0]     = ['fort.24','fort.27']
labelList[0]    = ['rms.X','max.X']
xdataList[0]    = [1,1]
ydataList[0]    = [3,2]
xyLabelList[0]  = ['z drection (m)','rms and max beam size in X (m)']

saveName.append('sizeY')
fileList[1]     = ['fort.25','fort.27']
labelList[1]    = ['rms.Y','max.Y']
xdataList[1]    = [1,1]
ydataList[1]    = [3,4]
xyLabelList[1]  = ['z drection (m)','rms and max beam size in Y (m)']

saveName.append('sizeZ')
fileList[2]     = ['fort.26','fort.27']
labelList[2]    = ['rms.Z','max.Z']
xdataList[2]    = [1,1]
ydataList[2]    = [3,6]
xyLabelList[2]  = ['z drection (m)','rms and max phase size in Z (deg)']

saveName.append('emitXY')
fileList[3]     = ['fort.24','fort.25']
labelList[3]    = ['emit.nor.X','emit.nor.Y']
xdataList[3]    = [1,1]
ydataList[3]    = [7,7]
xyLabelList[3]  = ['z drection (m)','emittance at X and Y (m*rad)']

lineType = ['r-','b--']
plotPath = './post'
dataPath = plotPath + '/data'
if os.path.exists(plotPath) == False:
    os.makedirs(plotPath)
#if os.path.exists(dataPath) == False:
#    os.makedirs(dataPath)

f1 = plt.figure(0,figsize=(16,8))
#plt.tight_layout(pad=3, h_pad=2, w_pad=4)
plt.ioff()

#===============for saving====================#
if saveFlag:
    for i in range(0,picNum):
        print "ploting " + str(i+1) + '/' + str(picNum) +  " : " + saveName[i]
        for j in range(0,2):
            try:
                fin = open(fileList[i][j],'r')
            except:
                print "  ERRPR! Can't open file '" + fileList[i][j] + "'"
                exit(-1)
            linesList  = fin.readlines()
            fin .close()
            linesList  = [line.split() for line in linesList ]
            xId = xdataList[i][j]-1
            yId = ydataList[i][j]-1
            x   = [xrt[xId] for xrt in linesList]
            y   = [xrt[yId] for xrt in linesList]
            plt.plot(x, y, lineType[j], linewidth=2, label=labelList[i][j])
        plt.xlabel(xyLabelList[i][0])
        plt.ylabel(xyLabelList[i][1])
        plt.legend()
        plt.savefig(plotPath+'/'+saveName[i]+'.'+saveFormat)
        plt.clf()
    print ('Plots are save at "' + plotPath + '"')
else:
    print ('Plots are not saved, turn the "saveFlag" at line 10 of this python file to turn it on')


#===============for showing====================#
if showFlag:
    print ('Showing...')
    ax = []
    ax.append(plt.subplot(221))
    ax.append(plt.subplot(222))
    ax.append(plt.subplot(223))
    ax.append(plt.subplot(224))
    for i in range(0,picNum):
        for j in range(0,2):
            try:
                fin = open(fileList[i][j],'r')
            except:
                print "  ERRPR! Can't open file '" + fileList[i][j] + "'"
                exit(-1)
            linesList  = fin.readlines()
            fin .close()
            linesList  = [line.split() for line in linesList ]
            xId = xdataList[i][j]-1
            yId = ydataList[i][j]-1
            x   = [xrt[xId] for xrt in linesList]
            y   = [xrt[yId] for xrt in linesList]
            ax[i].plot(x, y, lineType[j], linewidth=2, label=labelList[i][j])
        ax[i].set_xlabel(xyLabelList[i][0])
        ax[i].set_ylabel(xyLabelList[i][1])
        box = ax[i].get_position()
        ax[i].set_position([box.x0, box.y0, box.width, box.height * 0.9])
        ax[i].legend(loc='upper center', bbox_to_anchor=(0.5, 1.17),fancybox=True, shadow=True, ncol=5)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        f1.show()
else:
    print ('Plots are not shown, turn the "showFlag" at line 9 of this python file to turn it on')

raw_input("Press [Enter] to finish")
