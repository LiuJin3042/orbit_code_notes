# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 08:40:28 2019

@author: LJ
"""

from configuration import *
import numpy as np 
import matplotlib.pyplot as plot 

"""
density funciton of R:
       7.59x + 0.2415   (0 < x < 0.15)
p(x) = 1.38             (0.15 < x < 0.55)
       -2.913x + 2.9823 (0.55 < x < 1)
       0                (others)
       
dsitribution function of R:
       3.795x^2 + 0.2415x           (0 < x < 0.15)
F(x) = -0.0854 + 1.38x              (0.15 < x < 0.55)
       -1.4565x^2 + 2.9823x - 0.5258   (0.55 < x < 1)

reversed distribution function:
          -0.0318 + 0.1318*sqrt(0.0583 +15.18x) (0 < x < 0.1216)
F(x)^-1 = 0.725x + 0.0619 (0.1216 < x < 0.674)
          1.0238 - 0.3433*sqrt(8.8941 - 5.826*(0.5258+x)) (0.674 < x < 1)
"""



def rdf(x):
    # reversed distribution function
    y = np.zeros(x.size)
    y[(0 <= x) * (x <= 0.1216)] = -0.0318 + 0.1318*np.sqrt(0.0583 +15.18*x[(0 <= x) * (x <= 0.1216)])
    y[(0.1216 < x) * (x <= 0.674)] = 0.725*x[(0.1216 < x) * (x <= 0.674)] + 0.0619
    y[(0.674 < x) * (x <= 1)] = 1.0238 - 0.3433*np.sqrt(8.8941 - 5.826*(0.5258+x[(0.674 < x) * (x <= 1)]))
    return y

def df(x):
    # density function
    y = np.zeros(x.size)
    y[(0 <= x) * (x <= 0.15)]  = 7.59*x[(0 <= x) * (x <= 0.15)] + 0.2415
    y[(0.15 < x) * (x <= 0.55)] = 1.38
    y[(0.55 < x) * (x <= 1)] = -2.913*x[(0.55 < x) * (x <= 1)] + 2.9823
    return y
    
gen_particle_number = nprt + 1
head = """
 5 %d
 N=   %d Emin= 0.000000E+00   rng_seed=        100001
 "AC" FILE 63887A01.DATA9  TIME =  5.8906E+00 +/-  1.2500E-02 sec.           
  tokamak: ???? NLSYM=F #= 1 A= 2.0 Z= 1.0 NEUTRAL BEAM                      
 host=login113                               date=29-Mar-2018 09:30:07
  R(cm)  Z(cm)  vpll/v  E(eV) -- 4(1x,1pe13.6)  (get_fbm)
"""
head = head%(gen_particle_number,gen_particle_number)
gen_particle_number = 30000

"""
generate R, the major radius of particle
"""
# generate y, whichi is normalized radius
# R = rmaj + a*y*cos(theta)
x = np.random.rand(gen_particle_number)
y = rdf(x)
# generate theta
theta = np.random.rand(gen_particle_number)*np.pi*2
# get R
gen_particle_R = rmaj + a*y*np.cos(theta)
# a brief view of the distribution of y
t = np.arange(0,1,0.001)
y0 = df(t)
y1 = rdf(t)
plot.plot(t, y0, 'r-', linewidth=1) 
plot.plot(t, y1, 'r-', linewidth=1) 
plot.hist(y, bins=80, normed=1, facecolor='green', alpha=0.75)
plot.show()
plot.hist(a*y*np.cos(theta), bins=80, normed=1, facecolor='green', alpha=0.75)
plot.show()
plot.hist(gen_particle_R, bins=10, normed=1, facecolor='green', alpha=0.75)
plot.show()

"""
generate pitch angle, which follows the normal distribution
"""
mu = 0.3
sigma = 0.1
ptch = np.random.normal(mu,sigma,gen_particle_number)











