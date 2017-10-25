#!/usr/bin/env python
import os
os.system('make')
print 'cd testrun'
os.chdir('./testrun')
os.system('mpirun -n 2 ./ImpactZexe')
os.chdir('..')
