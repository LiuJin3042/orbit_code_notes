# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 07:41:43 2019

@author: LJ
"""

from configuration import *
import numpy as np

def num2str(num):
    str = '%.6E'%(num)
    if num < 0:
        str += ' -'
    if num > 0:
        str += '  '
    return str

# generate R, this case we use normal dist
# creat sigma, ensure 99.7% samples is in this range
sigma = a/6
R = np.random.normal(loc=rmaj, scale=sigma, size=[nprt,1])

# generate Z
