#!/usr/bin/env python
import os
print 'cd testrun'
os.chdir('./testrun')
os.system('mpirun -n 1 ./ImpactZexe')
os.chdir('..')
