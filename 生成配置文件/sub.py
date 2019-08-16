# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 15:59:25 2019

@author: LJ
"""
from commands import getstatusoutput as gso
import use
from yae_sakura import *
import time

def sub_task(comment,pdist,numeric):
    print(ying)
    # make file 
    print('making files, please wait')
    status, output = gso('make FC=pgf90')
    print(output)
    # rename the desired file folder as 20010101-comment
    date = time.strftime('%Y%m%d',time.localtime(time.time()))
    # des_folder: destination of output files
    des_folder = date + '-' + comment
    # remove than creat the folder
    gso('rm -rf %s'%des_folder)
    gso('mkdir %s'%des_folder)
    # cp certain files based the value of pdist and numeric
    if pdist*numeric == 2:
        # numeric balance and distribution
        gso('cp {orbit,orbit.F,spdata,fbm_dist.dat,job.pbs} ./%s'%des_folder)
    elif pdist == 2:
        # numeric distribution
        gso('cp {orbit,orbit.F,fbm_dist.dat,job.pbs,spdata} ./%s'%des_folder)
    else:
        gso('cp {orbit,orbit.F,spdata,job.pbs} ./%s'%des_folder)
    # let the user choose weather to submit the job
    # submit = raw_input('configuration is completed, submit the mission?(y/n) ')
    submit = 'y'
    monitor = 'y'
    if submit == 'y':
        status, output = gso('qsub ./%s/job.pbs'%des_folder)
        print(output)
        # monitor = raw_input('monitor the results?(y/n) ')
        if monitor == 'y':
            use.monitor()
            use.pack(des_folder)
    else:
        print('not submitted')
