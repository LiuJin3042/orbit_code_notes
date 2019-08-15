# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 16:40:57 2019

@author: LJ
"""

import sys
import time
from commands import getstatusoutput as gso
from configuration import *

def monitor():
    # a program to monitor running status
    if sys.version[0] == '2':
        if nplot == 2:
            snap_time = 45
        if nplot in [1,3]:
            snap_time = 3
        else:
            snap_time = 15
        status, output = gso('qstat')
        while output != '':
            print(output)
            time.sleep(snap_time)
            status, output = gso('qstat')
            print('\n\n\n\n\n')
        print('job is done')

def pack(des_folder):
    # a program to package file 
    gso('cp {*.plt orbit.out configuration.py} ./%s/orbit_results'%(des_folder))
    gso('cp -r ./plot_functions ./%s'%(des_folder))