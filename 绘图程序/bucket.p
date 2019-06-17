#device postencap :SY@: :OF@:   new.ps
#Read poincare data
erase
window 1 1 1 1
lweight 1
data bucket.out
read t 1
read f 2
read p 3
read e 4
read pol 5
ctype black
expand 2
limits t pol
lweight 2
ptype 0 0
connect t pol
ctype black
expand 2
relocate (2000 22000)
putlabel 5 f kHz
expand 1
ltype 0
relocate (2000 15000)
angle 0
expand 2
xlabel  t msec
expand 1.5
box
#quit
