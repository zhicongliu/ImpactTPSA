#!/usr/bin/env python

#This script is run Impact-Z using a TraceWin lattice.
#It converts a TraceWin lattice file 
#to a Impact-Z input file and call Impact-Z.
#Some parameter need to be set as following.

import sys
import math
import subprocess
print 'Impact, From Python input file'

#-------------------------------------------------#
#--------------------parameter--------------------#
#-------------------------------------------------#

#1. Processor number for parallel
ProcessorY      = 4
ProcessorZ      = 4

#2. dimension
dimension       = 6
particleNumber  = 16000
integratorFlag  = 2
errorStudyFlag  = 0
outputFlag      = 1

#3. grid property for space charge calculation
gridNumberX     = 65
gridNumberY     = 65
gridNumberZ     = 129
boundryFlag     = 4
widthX          = 0.006
widthY          = 0.006
periodLength    = 7.64

#4. initial distribution
distributionType= 4
restartFlag     = 0
subcycleFlag    = 0
numberOfChargeState     = 4

#5. number of particle list for each charge state
particleNumberList      = [0]*numberOfChargeState
for i in range(numberOfChargeState):
  particleNumberList[i]=16000

#6. number of current list for each charge state
currentList     = [0.0]*numberOfChargeState
for i in range(numberOfChargeState):
  currentList[i]=0.0

#7. number of current list for each charge state
QtoM_List       = [0.0]*numberOfChargeState
for i in range(numberOfChargeState):
  QtoM_List[i]=4.06e-10

#8. beam size in X, Twiss, mismatch, offset
sigmaX          = '2.6052929781E-2'
lambdaX         = '3.6536731692E-4'
muX             = '-0.95073401'
mismatchX       = '1.0'
mismatchPx      = '1.0'
offsetX         = '0.0'
offsetPx        = '0.0'

#9. beam size in Y, Twiss, mismatch, offset
sigmaY          = '1.6921326621E-2'
lambdaY         = '2.3573085744E-4'
muY             = '-0.67283857'
mismatchY       = '1.0'
mismatchPy      = '1.0'
offsetY         = '0.0'
offsetPy        = '0.0'

#10. beam size in Z, Twiss, mismatch, offset
sigmaZ          = '6.1423272782E-2'
lambdaZ         = '1.7583221689E-4'
muZ             = '0.148895386'
mismatchZ       = '1.0'
mismatchPz      = '1.0'
offsetZ         = '0.0'
offsetPz        = '0.0'

#11. current and kinetic energy
currentAverage  = 0.070         # A
kineticE        = 120.0e6       # eV
particleMass    = 939.394e6     # eV
charge          = 1.0
frequency       = 352.2e6       # Hz
phaseInitial    = 0.0           #

#-----------------------------------------------#
#--------------------Lattice--------------------#
#-----------------------------------------------#




#==================Output Start======================#

adresssOut = 'test.in'
if os.path.exists(adresssOut):
    exit(-1)
else:
    fout = open(adresssOut,'w')

#1.
#number of procs in Y direction, number of procs in Z direction
#NOTE WELL: the product of these numbers must equal the number of
#           processors that you run on!
fout.write(str(ProcessorY).ljust(8)+' '+str(ProcessorZ)+'\n')

#2.
#"6" - 6D, "20" - # of particles for each charge state, "2" - nonlinear
#Lorentz integrator, "0" - no error study, "1" - standard output
fout.write(str(dimension)       .ljust(8)+' ')
fout.write(str(particleNumber)  .ljust(8)+' ')
fout.write(str(integratorFlag)  .ljust(8)+' ')
fout.write(str(errorStudyFlag)  .ljust(8)+' ')
fout.write(str(outputFlag)      .ljust(8)+' '+'\n')

#3.
#65x65x129 mesh grid, "4" - transverse round pipe, longitudinal periodic BC,
#"0.014" - x pipe width (m), y pipe width (m) and period length (m).
fout.write(str(gridNumberX)     .ljust(8)+' ')
fout.write(str(gridNumberY)     .ljust(8)+' ')
fout.write(str(gridNumberZ)     .ljust(8)+' ')
fout.write(str(boundryFlag)     .ljust(8)+' ')
fout.write(str(widthX)          .ljust(8)+' ')
fout.write(str(widthY)          .ljust(8)+' ')
fout.write(str(periodLength)    .ljust(8)+' '+'\n')

#4.
#oflagdist, orstartflg, oflagsbstp, onchrg
fout.write(str(distributionType).ljust(8)+' ')
fout.write(str(restartFlag)     .ljust(8)+' ')
fout.write(str(subcycleFlag)    .ljust(8)+' ')
fout.write(str(numberOfChargeState)     .ljust(8)+' ')

#5.
# particle list for each charge state
for i in particleNumberList:
  fout.write(str(i)+'  ')
fout.write('\n')

#6.
#  current list for each charge state
for i in currentList:
  fout.write(str(i)+'  ')
fout.write('\n')

#7.
#      q/m list for each charge state
for i in QtoM_List:
  fout.write(str(i)+'  ')
fout.write('\n')

#8
# twiss X
fout.write(str(sigmaX)          .ljust(20))
fout.write(str(lambdaX)         .ljust(20))
fout.write(str(muX)             .ljust(20))
fout.write(str(mismatchX)       .ljust(8)+' ')
fout.write(str(mismatchPx)      .ljust(8)+' ')
fout.write(str(offsetX)         .ljust(8)+' ')
fout.write(str(offsetPx)        .ljust(8)+' '+'\n')

#9
# twiss X
fout.write(str(sigmaY)          .ljust(20))
fout.write(str(lambdaY)         .ljust(20))
fout.write(str(muY)             .ljust(20))
fout.write(str(mismatchY)       .ljust(8)+' ')
fout.write(str(mismatchPy)      .ljust(8)+' ')
fout.write(str(offsetY)         .ljust(8)+' ')
fout.write(str(offsetPy)        .ljust(8)+' '+'\n')

#10
# twiss X
fout.write(str(sigmaZ)          .ljust(20))
fout.write(str(lambdaZ)         .ljust(20))
fout.write(str(muZ)             .ljust(20))
fout.write(str(mismatchZ)       .ljust(8)+' ')
fout.write(str(mismatchPz)      .ljust(8)+' ')
fout.write(str(offsetZ)         .ljust(8)+' ')
fout.write(str(offsetPz)        .ljust(8)+' '+'\n')

#11
# particle information
fout.write(str(currentAverage)  .ljust(20))
fout.write(str(kineticE)        .ljust(20))
fout.write(str(particleMass)    .ljust(20))
fout.write(str(charge)          .ljust(8)+' ')
fout.write(str(frequency)       .ljust(8)+' ')
fout.write(str(phaseInitial)    .ljust(8)+' '+'\n')

#=============Parameter Output End=======================#


fout.write('\n'+'!'+'-'*47+'!'+'\n')
fout.write('!'+'-'*20+'Lattice'+'-'*20+'!'+'\n')
fout.write('!'+'-'*47+'!'+'\n'*2)

#=============Lattice Output Start=======================#
ele = []
fin = open('lattice.txt','r')
failout = open('TWconvert.fail','w')
lines = fin.readlines()
for i in range(0,len(lines)-1):
  line=lines[i].expandtabs()
  if ';' in line:
    line=line[:line.index(';')]+'\n'
  ele = line.split()
  if len(ele) == 0 or ele[0][0] == ';':
    continue
  if ele[0].upper() == 'DRIFT':
    fout.write(str(float(ele[1])/1000)  .ljust(8)+' ')
    fout.write(str(4)                   .ljust(8)+' ')      #step
    fout.write(str(20)                  .ljust(8)+' ')      #map step
    fout.write(str(0)                   .ljust(8)+' ')      #drift
    fout.write(str(float(ele[2])/1000)  .ljust(8)+' ')
    fout.write('\n')
  elif ele[0].upper() == 'QUAD':
    fout.write(str(float(ele[1])/1000)  .ljust(8)+' ')
    fout.write(str(4)                   .ljust(8)+' ')
    fout.write(str(20)                  .ljust(8)+' ')
    fout.write(str(1)                   .ljust(8)+' ')      #quad
    fout.write(ele[2]                   .ljust(8)+' ')
    fout.write('0.'                     .ljust(8)+' ')
    fout.write(str(float(ele[3])/1000)  .ljust(8)+' ')
    fout.write('0.0'                    .ljust(8)+' ')      #x misalignment error
    fout.write('0.0'                    .ljust(8)+' ')      #y misalignment error
    fout.write('0.0'                    .ljust(8)+' ')      #x rotation error
    fout.write('0.0'                    .ljust(8)+' ')      #y rotation error
    fout.write('0.0'                    .ljust(8)+' ')      #z rotation error
    fout.write('\n')
  elif ele[0].upper() == 'BEND':
    fout.write(str(float(ele[2])/1000*float(ele[1]) / 360*2 * math.pi)  .ljust(8)+'  ')      #length
    fout.write(str(10)                  .ljust(8)+' ')      #step
    fout.write(str(20)                  .ljust(8)+' ')      #map step
    fout.write(str(4)                   .ljust(8)+' ')      #dipole
    fout.write(str(float(ele[1]) / 360*2 * math.pi).ljust(8)+' ')#angle at rad

    if i >= 1:
      for j in range(i-1,0,-1):
        lineEdge= lines[j].expandtabs()
        eleEdge = lineEdge.split()
        if len(eleEdge) == 0 or eleEdge[0][0] == ';':
          continue
        if eleEdge[0].upper() == 'EDGE':
          fout.write(str(float(eleEdge[4])).ljust(8)+' ') #k1
        else:
          fout.write(str(0.0).ljust(8)+' ')
        break
    else:
       fout.write(str(0.0).ljust(8)+' ')
       
    fout.write(str(150).ljust(8)+' ')                   #input switch ???
    fout.write(str(float(ele[2]) / 1000).ljust(8)+' ')  #radius
    
    if i >= 1:
      for j in range(i-1,0,-1):
        lineEdge= lines[j].expandtabs()
        eleEdge = lineEdge.split()
        if len(eleEdge) == 0 or eleEdge[0][0] == ';':
          continue
        if eleEdge[0].upper() == 'EDGE':
          fout.write(str(float(eleEdge[1])/360*2*math.pi).ljust(8)+' ') #enter angle
        else:
          fout.write(str(0.0).ljust(8)+' ')
        break
    else:
       fout.write(str(0.0).ljust(8)+' ')

    if i < len(lines)-1:
      for j in range(i+1,len(lines)-1):
        lineEdge= lines[j].expandtabs()
        eleEdge = lineEdge.split()
        if len(eleEdge) == 0 or eleEdge[0][0] == ';':
          continue
        if ele[0].upper() == 'EDGE':
          fout.write(str(float(eleEdge[1]-ele[1])/360*2*math.pi).ljust(8)+' ') #exit angle
        else:
          fout.write(str(0.0).ljust(8)+' ')
        break
    else:
       fout.write(str(0.0).ljust(8)+' ')
       
    if i >= 1:
      for j in range(i-1,0,-1):
        lineEdge= lines[j].expandtabs()
        eleEdge = lineEdge.split()
        if len(eleEdge) == 0 or eleEdge[0][0] == ';':
          continue
        if eleEdge[0].upper() == 'EDGE':
          fout.write(str(float(eleEdge[2])/1000).ljust(8)+' ') #curvature of enter face
        else:
          fout.write(str(0.0).ljust(8)+' ')
        break
    else:
       fout.write(str(0.0).ljust(8)+' ')

    if i < len(lines)-1:
      for j in range(i+1,len(lines)-1):
        lineEdge= lines[j].expandtabs()
        eleEdge = lineEdge.split()
        if len(eleEdge) == 0 or eleEdge[0][0] == ';':
          continue
        if eleEdge[0].upper() == 'EDGE':
          fout.write(str(float(eleEdge[2])/1000).ljust(8)+' ') #curvature of exit face
        else:
          fout.write(str(0.0).ljust(8)+' ')
        break
    else:
       fout.write(str(0.0).ljust(8)+' ')
    fout.write(str(0.0).ljust(8)+' ')           #intergrated fringe field???
    fout.write('\n')
  elif ele[0].upper() == 'EDGE':
    a = 'need to do nothing because EDGE had been taken into consideration at BEND'
#  elif ele[0].upper() == 'FIELD_MAP':
#    if ele[1] == '7700':
    
  else:
    failout.write(line+'\n')
fin.close( )
print 'convert finished'
print '='*50
cmd = 'mpiexec -n ' + str(ProcessorY*ProcessorZ) + './ImpactZexe'
print cmd
#subprocess.call(cmd,shell=True)

