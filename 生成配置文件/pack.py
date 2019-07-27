# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 22:04:41 2019

@author: LJ
"""

from commands import getstatusoutput as gso
import time

from configuration import comment

date = time.strftime('%Y%m%d',time.localtime(time.time()))
status, output = gso('tar -czvf %s-%s.tar.gz *.plt orbit.out configuration.py'%(comment,date))
print(output)
