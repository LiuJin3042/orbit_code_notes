# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 07:47:45 2019

@author: LJ
"""

import commands
(status, output) = commands.getstatusoutput('cat /proc/cpuinfo')
print status, output
os.popen('mkdir)