# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 17:55:23 2019

Read the configuration file, read the file to be changed, and write the new file

@author: LJ
"""

from configuration import *
import sys

"""
read and rewrite eqs.f
output file should be in the same dir
"""
# read files - rewrite certain lines - write to new files
r_eqs = open('./source_file/eqs.f', 'r')
w_eqs = open('./eqs.f', 'w')
eqs = r_eqs.readlines()
eqs[19] = '      numeric = ' + str(numeric) + '\n'
if numeric == 0:
    eqs[50] = '        rmaj = ' + str(rmaj) + '\n'
    eqs[106] = '      eps = ' + str(a) + '.D0/rmaj\n'
    eqs[61] = '      krip = ' + str(krip) + '\n'
    eqs[109] = '      q0 = ' + str(q0) + '\n'
    eqs[110] = '      qed = ' + str(qed) + '\n'
    eqs[111] = '      rqx = ' + str(rx) + '\n'
    eqs[113] = '      qx = ' + str(float(qrx)) + '\n'
    qr3 = (qrx-q0+rx**2*q0-rx**2*qed)/(rx**2*(rx-1))
    qr2 = qed-q0-qr3
#    eqs[126] = '      qr2 = ' + str(qr2)[0:5] + '\n'
#    eqs[127] = '      qr3 = ' + str(qr3)[0:5] + '\n'
    
w_eqs.writelines(eqs)
r_eqs.close()
w_eqs.close()


"""
read and rewrite perturb.f
output file should be in the same dir
"""
r_ptrb = open('./source_file/perturb.f', 'r')
w_ptrb = open('./perturb.f', 'w')
ptrb = r_ptrb.readlines()
ptrb[17] = '      modes = ' + str(modes) + '\n'

# set mode params, e.g
#      harm(1) = 1
#      mmod(1) = 2
#      nmod(1) = 1
#      amp(1) = 5.0D-4
#      omegv(1) = 0.0*2.D3*pi/omeg0
#      alfv(1) = 1
mode_params = ''
for i in range(modes):
    j = i+1
    single_set =  '      harm(%d) = ' + str(harm[i]) + '\n      mmod(%d) = ' + \
    str(mmod[i]) + '\n      nmod(%d) = ' + str(nmod[i]) + '\n      amp(%d) = ' + \
    str(float(amp[i])) + '\n      omegv(%d) = ' + str(float(omegv[i])) + '*2.D3*pi/omeg0\n      alfv(%d) = '+ \
    str(alfv[i]) + '\n'
    mode_params += single_set%(j,j,j,j,j,j)
ptrb[18] = mode_params
for i in range(modes):
    if omegv[i] != 0:
        ptrb[19] = '      dele = ' + str(dele) + '\n'
        break
    
# set a1-alpha
#a1(j,md) = exp(-((xd-cnt(md))/wdt(md))**2)   ! gaussian
#a1(j,md) = a1(j,md)*(rm/rn - qdum)    !  gaussian MHD
#a1(j,md) = (eps*xd)**m*(1-n*qdum/m)/(gdum*qdum)    ! MHD
#a1(j,md) = (eps*xd)**m*(pw - px)    ! resistive
mode_type = ['exp(-((xd-cnt(md))/wdt(md))**2)',
             'a1(j,md)*(rm/rn-qdum)',
             '(eps*xd)**m*(1-n*qdum/m)/(gdum*qdum)',
             '(eps*xd)**m*(pw-px)']
if a1 == 2:
    ptrb[51] = '         a1(j,md) = ' + mode_type[0] + '\n         a1(j,md) = ' + \
    mode_type[1] + '\n'
else:
    ptrb[51] = '         a1(j,md) = ' + mode_type[a1-1] + '\n'
w_ptrb.writelines(ptrb)
r_ptrb.close()
w_ptrb.close()


"""
read and rewrite orbit.F
output file should be in the same dir
"""
r_orbit = open('./source_file/orbit.F', 'r')
w_orbit = open('./orbit.F', 'w')
orbit = r_orbit.readlines()
orbit[117] = '        npert = ' + str(npert) + '\n'
if pdist != 2:
    orbit[126] = '        polo = ' + str(polo) + '*pw\n'
    orbit[127] = '        p1 = ' + str(p1) + '*pw\n'
    orbit[128] = '        p2 = ' + str(p2) + '*pw\n'
    orbit[129] = '        pchi = ' + str(pchi) + '\n'
    orbit[133] = '      zprt = ' + str(zprt) + '.D0\n'
    orbit[134] = '      prot = ' + str(prot) + '.D0\n'
    orbit[135] = '      ekev = ' + str(ekev) + '\n'
    orbit[108] = '      bkg = ' + str(bkg) + '\n'
    orbit[94] = '        ntor = ' + str(ntor) + '\n'
    orbit[79] = '        nprt = ' + str(nprt) + '\n'
orbit[75] = '      nplot = ' + str(nplot) + '\n'
ndist = ['shelldep', 'sampledep', 'poindep']
orbit[242] = '        call ' + ndist[pdist-1] + '\n'

w_orbit.writelines(orbit)
r_orbit.close()
w_orbit.close()

def make_files():
    # make files
    cmd = '''
    from commands import getstatusoutput as gso
    status, output = gso('make FC=pgf90 eqs')
    print status, output
    status, output = gso('./eqs')
    print status, output
    status, output = gso('make FC=pgf90')
    print status, output
    # shell will take input 'n' as a variable rather than a string
    # so this step is necessary
    n = 'n'
    y = 'y'
    submit = input('configuration is completed, submit the mission?(y/n)')
    if submit == 'y':
        status, output = gso('qsub job.pbs')
        print status, output
    else:
        print 'not submitted'
    '''
    exec(cmd)

if sys.version[0] == '2':
    make_files()


























