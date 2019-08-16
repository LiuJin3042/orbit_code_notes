# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 08:04:15 2019

@author: LJ
"""

# ekev from 20 to 80
# amp from 1e-5 to 1e-4

from __future__ import division
from configuration import *
import make

def linspace(start,stop,count):
    # generate a arithmetic progression, from start to stop, total number is count
    step = (stop-start)/(count-1) 
    ap = []
    for i in range(count):
        ap.append(start + i*step)
    return ap

l_ekev = linspace(20,77,6)
l_amp = linspace(1e-5,7e-4,5)

for iamp in l_amp:
    for iekev in l_ekev:
        amp = [iamp]
        ekev = iekev
        make.main(numeric,a,rmaj,rx,krip,q0,qed,qrx,modes,harm,nmod,mmod,omegv,alfv,amp,dele,a1,npert,polo,p1,p2,pchi,zprt,prot,ekev,bkg,ntor,nprt,nplot,pdist)
             

