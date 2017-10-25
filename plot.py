#!/usr/bin/env python
import os
print 'Plot, cd testrun'
os.chdir('./testrun')
os.system('./ImpactZPost.py')
os.chdir('..')
