        program eqs
      IMPLICIT NONE
      include 'orbcom'
      include 'o.cln'
      INTEGER netcdf   !Input file data format
      INTEGER kk
      LOGICAL file_exist
      open(6,file='eqout',status='unknown')
      open(20,file='spdata',status='unknown')
      open(21,file='data1',status='unknown')
cc   r. b. white   princeton, jan 1982
cc   files needed: o.cln, spline, bzio.f eqsub.f equilibrium mp0... and mp1...
ccc  compile subroutines eqsub.f with shesub
cc   shell she will run and print data
ccc
c --- S.ETHIER 03/14/2000  Here we choose the data format of the input file
c                          If netcdf=1 then the file is in the NetCDF format
c                          and its default name is "map01.cdf".
      netcdf = 1
      numeric = 1   ! 0 gives analytic equilibrium from tok0, 2=rfp
c
      if (numeric .ne. 0) then
cccc     numerical equilibrium choice
         if (netcdf .eq. 1) then
            mp0 = 'map01.cdf'
            mp1 = 'dummy'
            INQUIRE(FILE=mp0, EXIST=file_exist)
            if (.not.file_exist) then
               kk=LEN_TRIM(mp0)
               write(*,*) '*** Input file ',mp0(1:kk),' DOES NOT EXIST!'
               stop
            endif
         else
            mp0 = 'm0tok'
            mp1 = 'm1tok'
            INQUIRE(FILE=mp0, EXIST=file_exist)
            if (.not.file_exist) then
               kk=LEN_TRIM(mp0)
               write(*,*) '*** Input file ',mp0(1:kk),' DOES NOT EXIST!'
               stop
            endif
            INQUIRE(FILE=mp1, EXIST=file_exist)
            if (.not.file_exist) then
               kk=LEN_TRIM(mp1)
               write(*,*) '*** Input file ',mp1(1:kk),' DOES NOT EXIST!'
               stop
            endif
         endif
      endif
c
        rmaj = 165   ! magnetic axis in cm
ccc
         lsp = 201   ! b,x,z,giac,q,ripple- poloidal spline grid points
         lst = 141   ! b,x,z,giac,ripple, must be odd-theta spline grid points
         lemax = 4    ! g, I
         lrmax = 8              ! pol(r), r(pol), pressure
         write(6,26) lsp,lst
 26      format(' spline dimensions lsp,lst ',2i4)
ccc
cccc  ripple choice krip=1-TFTR, krip=2-Tore Supra, krip=3-ITER
cccc    krip=4-NSTX, krip=5-Ignitor
      krip = 0
cccc
      if(numeric.eq.0) call tok0
      if(numeric.eq.2) call rfp
      if(numeric.eq.1) call main0(rmaj,netcdf)     ! will change rmaj
      pw = psival(nosurf)   ! poloidal flux
ccc
      write(6,20)
 20   format(' code eqs.f, equilibrium for orbit.f')
      write(6,760) mp0
 760  format(a10)
      write(6,25) pw,ped,nosurf,mth
 25   format(' pw ped nosurf mth',1p2e12.4,2i6)
      call polrep
      write(6,22) xc,eps,bax,rmaj
 22   format(' xc,eps,bax,rmaj ',1p5e12.4)
      call rstor
      call wdat
ccc
ccc      call contin  ! checks continuity spline conditions
ccc   call wrtb
 999  continue
      stop
      end
cccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccc
      subroutine tok0
      IMPLICIT NONE
      include 'o.cln'
      include 'orbcom'
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER lmax,ier,i,j
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 deltap,delta,bet,qed,rqx,qx,shift,err0,pold,pz,rdum,
     > rfun,qdum,qfnr,del0,gfan,tdum,dd1,dd2,dd3,dum,psid
C============
      common /dblk/ deltap(500),delta(500)
      character*(80) wlabel
ccccccc-Shafranov-analytic-field
ccccc- can choose pressure profile, Shafranov shift, and q profile
      write(6,1)
 1    format('  analytic Shafranov equilibrium')
      mp0 = 'analytic'
      eps = 40.D0/rmaj  !  last flux surface  R=1 is axis
      bet = .0D0
cc- set q0, q(rqx) = qx, qed, -> qr2,qr3
      q0 =  .8   !  q = q0 + qr2*r**2 + qr3*r**3
      qed = 4.    ! q at wall
      rqx =10.D0/rmaj    ! q = qx radius
      if(rqx.gt.eps) stop
      qx = 2.
      shift = 0.1D0/rmaj   ! shift of plasma center, R = 1
ccc
      wk2(1) = qed - q0
      wk2(2) = qx - q0
      bmat(1) = eps**2
      bmat(3) = eps**3
      bmat(2) = rqx**2
      bmat(4) = rqx**3
      lmax = 2
cc      call gelg(wk2,bmat,lmax,1,err0,ier) !fit q with matrix inversion
      qr2 = wk2(1)
      qr3 = wk2(2)
      qr2 = 16*3.2
      qr3 = 0.
      write(6,4) q0,qr2,qr3
 4    format('  q0,qr2,qr3 ',1p5e10.2)
ccccc
      pi = 4.D0*atan(1.D0)
      pi2 = 2*pi
      pi2i = 1/pi2
      mth = lst
      do 15 i = 1,lst
         theval(i)  = (i-1)*pi2/lst
  15   continue
      write(6,77) eps
 77   format(' eps ',1p2e12.4)
         write(6,98) lsp
 98      format(' lsp ',i4)
         nosurf = lsp
         mth1 = lst
      pw = pold(eps)
      ped = pw
      rq1 = 0
      do 10 j = 1,lsp
      pz = (j-1)*pw/(lsp-1.D0) + 1.D-20    !  j=1 is axis
      psival(j) = pz    ! this is poloidal flux
      rdum = rfun(pz)   ! minor radius
      ps1(j) = .5D0*rdum**2   ! this is toroidal flux
      rp1(j) = rdum
      qdum = qfnr(rdum)
      if(qdum.lt.1.D0) rq1 = rdum
      qd1(j) = qdum
      pd1(j) = .5D0*bet*(.5D0 + .5D0*cos(pi*(rdum/eps)**1.5D0))
      delta(j) = del0(rdum,shift)
cccc functions ri and rip
      cur(j) = rdum*rdum/qdum
      rd1(j) = cur(j)
cccc
ccccc   function g
         gd1(j) = gfan(rdum)
         gd1(j) = 1.D0     !!!!!!!
cccc    b of pol,thet
      do 20 i = 1,lst
         tdum = theval(i) + rdum*sin(theval(i))
         tdum = theval(i) 
      x1(j,i) = 1 + rdum*cos(tdum) - delta(j)
      z1(j,i) = rdum*sin(tdum)
      b1(j,i) = sqrt((gd1(j)/x1(j,i))**2 + (rdum/qdum)**2)
      g1(j,i) = (gd1(j)*qdum + cur(j))/b1(j,i)**2
      pv(j) = x1(j,1)   !   for plot
      xv(j) = g1(j,1)   !   for plot
      x(i,nosurf) = x1(lsp,i)
      z(i,nosurf) = z1(lsp,i)
20    continue
 10   continue
      write(6,88) rq1,delta(lsp)
 88   format(' r(q=1),shafranov shift at edge, units rmaj ',1p2e12.4)
      dd1 = (1 - rp1(lsp) - delta(lsp))*rmaj
      dd2 = (1 + rp1(lsp) - delta(lsp))*rmaj
      dd3 = .5D0*(dd1+dd2)
      write(6,89) dd1,dd3,rmaj,dd2
 89   format(' edge,center,axis,edge ',1p4e12.4)
      dum = psid(pw)
      write(6,188) dum
 188  format('  toroidal psi edge ',1pe12.4)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccc
      subroutine rfp
      IMPLICIT NONE
      include 'o.cln'
      include 'orbcom'
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER lmax,ier,i,j
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 deltap,delta,bet,qed,rqx,qx,shift,err0,pold,pz,rdum,
     > rfun,qdum,qfnr,del0,gfan,tdum,dd1,dd2,dd3,dum,psid,th0,
     >  dm4,dm6,dm8,dm1
C============
      character*(80) wlabel
ccccccc-rfp-analytic-field
ccccc- can choose th0, eps, rmaj is set in main
      th0 = 1.5D0
      write(6,1) th0
 1    format(' RFP -polynomial function, th0=',1pe12.4)
      mp0 = 'analytic'
      eps = 70.D0/rmaj  !  last flux surface  R=1 is axis
ccccc
      pi = 4.D0*atan(1.D0)
      pi2 = 2*pi
      pi2i = 1/pi2
      mth = lst
      do 15 i = 1,lst
         theval(i)  = (i-1)*pi2/lst
  15   continue
      write(6,77) eps
 77   format(' eps ',1p2e12.4)
         write(6,98) lsp
 98      format(' lsp ',i4)
         nosurf = lsp
         mth1 = lst
      rq1 = 0
      do 10 j = 1,lsp
      rdum = (j-1)/(lsp-1.D0) + 1.D-20   ! normalized minor radius
      rp1(j) = rdum
ccccc   function g
      gd1(j) = 1.D0 - th0**2*rdum**2 + .5*th0**2*rdum**4
cccc functions ri and rip
      cur(j) = eps*(th0*rdum**2 - .5*th0**3*rdum**4
     a       + th0*(th0**2 - 1.D0)*rdum**6/3.)
      rd1(j) = cur(j)
      pv(j) = cur(j)/(eps*rdum)   ! B pol for plot
      qdum = eps**2*gd1(j)*rdum**2/rd1(j)
      qd1(j) = qdum
      psival(j) = eps*(th0*rdum**2/2. - th0**3*rdum**4/8.
     a    + th0*(th0**2 - 1.D0)*rdum**6/18.)    ! this is poloidal flux
      ps1(j) = eps**2*(.5D0*rdum**2 - th0**2*rdum**4/4.
     a   + th0**2*rdum**6/12.)      ! this is toroidal flux
      dm4 = rdum**4 - 1.D0
      dm6 = rdum**6 - 1.D0
      dm8 = rdum**8 - 1.D0
      dm1 = rdum**10 - 1.D0
ccc-pressure
      pd1(j) = .25*th0**4*dm4 + .5*th0**4*dm6
     a    - .5*th0**2*dm4 - .125*th0**4*dm8
     a    - 4*th0**2*(th0**2 - 1.D0)*dm6/9. 
     a    - th0**6*dm6/6. + 5.*th0**6*dm8/24.
     a    - 5.*th0**4*dm8/24. - th0**2*(th0**2 - 1.D0)**2*dm1/15.
cccc
cccc    b of pol,thet
      do 20 i = 1,lst
         tdum = theval(i)
      x1(j,i) = 1 + eps*rdum*cos(tdum)
      z1(j,i) = eps*rdum*sin(tdum)
      b1(j,i) = sqrt((gd1(j))**2 + (rd1(j)/rdum)**2)
      g1(j,i) = (gd1(j)*qdum + cur(j))/b1(j,i)**2
20    continue
 10   continue
      pw = psival(lsp)
      ped = pw
      dum = ps1(lsp)
      write(6,188) dum
 188  format('  toroidal psi edge ',1pe12.4)
      do 50 j = 1,lsp
      xv(j) = 2*pd1(j)/pv(lsp)**2   ! beta poloidal for plot
 50   continue
      return
      end
ccccccccccccccccccccccccccccccccccccccccccc
      function qfnr(dum)
      IMPLICIT NONE
      include 'orbcom'
ccc...function q(r) -can be changed
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 dum,qfnr
C============
      qfnr = q0 + qr2*dum**2 + qr3*dum**3
      return
      end
ccccccccccccc cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function qfnrp(dum)
      IMPLICIT NONE
      include 'orbcom'
cccc-this is dq/dr
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 dum,qfnrp
C============
      qfnrp = 2*qr2*dum + 3*qr3*dum**2
      return
      end
ccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccc
      function del0(rdum,shift)
      IMPLICIT NONE
      include 'orbcom'
ccc-shafranof shift function
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 rdum,shift,del0
C============
        del0 = (rdum/eps)**2*shift
      return
      end
ccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function ppfun(rd)
      IMPLICIT NONE
      include 'orbcom'
ccc-pressure profile, can be changed
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 rd,ppfun
C============
      ppfun = -beta*rd/eps**2
      return
      end
! 14Mar2000 fgtok -s r8_precision.sub "r8con.csh conversion"
