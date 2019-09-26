# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 15:59:25 2019

@author: LJ
"""
from commands import getstatusoutput as gso
import use
from yae_sakura import ying
import time

def sub_task(comment, pdist, numeric, submit, monitor):
    print(ying)
    # make file 
    print('making files, please wait')
    print('making equilibrium files')
    status, output = gso('make FC=pgf90 eqs')
    print(output)
    status, output = gso('./eqs')
    print(output)
    status, output = gso('make FC=pgf90')
    print(output)
    # let the user choose weather to submit the job
    # make sure it is 'y' when using batch test
    if submit == 'y':
        status, output = gso('qsub ./job.pbs')
        print(output)
        if monitor == 'y':
            use.monitor()
            use.pack(comment,pdist,numeric)
    else:
        print('not submitted')
