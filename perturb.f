ccccccccccccccccccccccccccc
      subroutine ampa
      IMPLICIT NONE
      include 'orbcom'
      INTEGER lptm,md,n,m,j
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 amp0,tau1,tau2,xd,gdum,qdum,px,rm,rn,qfun,gfun,rpol,anorm
      REAL*8 wdt,cnt,dpx,xx0,pp
      common /splmd5/ wdt(10),cnt(10)
C============
      common /magni/ amp0(5),tau1,tau2
cccc-  Analytic alpha perturbation
cccc       kHz  = omegv(md)*omeg0/(2.D3*pi)
cccc  amp(md) is the ideal MHD displacement in units of major radius
cccc   the ideal xi is amp*r**(m-1)
      nvalx = 2
      modes = 2
      harm(1) = 1
      mmod(1) = 2
      nmod(1) = 1
      amp(1) = 5.0D-4
      omegv(1) = 0.0*2.D3*pi/omeg0
      alfv(1) = 1
      harm(2) = 1
      mmod(2) = 3
      nmod(2) = 2
      amp(2) = 2.D-5
      omegv(2) =  10.0*2.D3*pi/omeg0
      dele = 10      !    necessary if omega not zero
      md1 = 1
      md2 = 1
      alfv(2) = 1
cccccccccccc-  restrict modes 
ccc      go to 10
      modes = 1
      nval = 1
      nvalx = 1
      harm(1) = 1
      md1 = 1
      md2 = 1
 10   continue
cccccccccccc-Choose mode structure
      lpt = 900
      lptm = lpt - 1
      dpx = pw/lptm
      wdt(1) = .3D0
      cnt(1) = .5D0
      wdt(2) = .2D0
      cnt(2) = .6D0
ccccc define a1
      do 25 md = md1,md2
      anorm = 0
      n = nmod(md)
      m = mmod(md)
      rn = n
      rm = m
      do 20 j = 1,lpt
         px = (j-1)*pw/lptm
         qdum = qfun(px)
         gdum = gfun(px)
         xd = rpol(px)/eps  ! xd is minor radius 0 < xd < 1
cc         a1(j,md) = exp(-((xd-cnt(md))/wdt(md))**2)   ! gaussian
ccc         a1(j,md) = a1(j,md)*(rm/rn - qdum)    !  gaussian MHD
ccc         a1(j,md) = (eps*xd)**m*(1-n*qdum/m)/(gdum*qdum)    ! MHD
             a1(j,md) = (eps*xd)**m*(pw - px)    ! resistive
cc              xx0= 0.6  ! xx0 is the position of O point of magnetic island
cc              pp = 4/3   ! pp is a fractional number for an approximate value of m*(1/xx0-1)
cc         a1(j,md)=(xd/xx0)**m*((1-xd)/(1-xx0))**pp
ccc      ifac = .6D0*(1.D0 - sign(1.D0,qdum-m))
ccc      a1(j,md) = a1(j,md)*ifac
      if(abs(a1(j,md)).gt.anorm) anorm = abs(a1(j,md))
 20      continue
ccccc-  Normalize spline to 1
         do 21 j = 1,lpt
           a1(j,md) = a1(j,md)/anorm
ccc  DIAGNOSTIC
ccc      write(myfile,77) md,j,a1(j,md)
 77   format('  md,j,a1 ',2i5,1pe12.4)
 21        continue
         a1(1,md) = 0
         a1(2,md) = 0  ! near axis zero
 25      continue
      call splna
      return
      end
ccccccccccccccccccccccccccc
ccccccccccccccccccccccccccc
      subroutine ampx
      IMPLICIT NONE
      include 'orbcom'
      INTEGER lptm,md,n,m,j
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 amp0,tau1,tau2,xd,gdum,qdum,px,rm,rn,qfun,gfun,rpol,anorm
      REAL*8 wdt,cnt,dpx,rifun,alnm,dum,dum1
      common /splmd5/ wdt(10),cnt(10)
C============
      common /magni/ amp0(5),tau1,tau2
cccc-  Analytic xi perturbation
cccc       kHz  = omegv(md)*omeg0/(2.D3*pi)
cccc  amp(md) is the ideal MHD displacement in units of major radius
cccc   the ideal xi is amp*r**(m-1)
      nvalx = 2
      modes = 2
      harm(1) = 1
      mmod(1) = 4
      nmod(1) = 3
      amp(1) = 5.D-3
      omegv(1) = 100*2.D3*pi/omeg0
      alfv(1) = 1
      harm(2) = 1
      mmod(2) = 3
      nmod(2) = 2
      amp(2) = 2.D-5
      omegv(2) =  10*2.D3*pi/omeg0
      dele = 10      !    necessary if omega not zero
      md1 = 1
      md2 = 1
      alfv(2) = 1
cccccccccccc-  restrict modes
ccc      go to 10
      modes = 1
      nval = 1
      nvalx = 1
      harm(1) = 1
      md1 = 1
      md2 = 1
 10   continue
cccccccccccc-Choose mode structure
      lpt = 900
      lptm = lpt - 1
      dpx = pw/lptm
      wdt(1) = .3D0
      cnt(1) = .5D0
      wdt(2) = .2D0
      cnt(2) = .6D0
ccccc define xi1
      do 25 md = md1,md2
      anorm = 0
      n = nmod(md)
      m = mmod(md)
      rn = n
      rm = m
      do 20 j = 1,lpt
         px = (j-1)*pw/lptm
         qdum = qfun(px)
         gdum = gfun(px)
      dum = m*gdum + n*rifun(px)
      dum1 = m - n*qdum
      xd = rpol(px)/eps  ! xd is minor radius 0 < xd < 1
      alnm = exp(-((xd-cnt(md))/wdt(md))**2)   ! gaussian alpha
ccc      alnm = (eps*xd)**m*(pw - px)    ! resistive alpha
         xi1(j,md) = dum*alnm/dum1
      if(abs(xi1(j,md)).gt.anorm) anorm = abs(xi1(j,md))
ccc  DIAGNOSTIC
ccc      write(92,77) md,j,xi1(j,md),alnm
 77   format('  md,j,xi1,alpha ',2i5,1p3e12.4)
 20      continue
ccccc-  Normalize spline to 1
         do 21 j = 1,lpt
           xi1(j,md) = xi1(j,md)/anorm
 21        continue
         xi1(1,md) = 0
         xi1(2,md) = 0  ! near axis zero
 25      continue
      call splnx
      return
      end
ccccccccccccccccccccccccccc
ccccccccccccccccccccccccccc
      subroutine ptrba(jpts)
      IMPLICIT NONE
      include 'orbcom'
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER jpts,k,md,n,m,jd,ndum
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 az,phnm,pdum,dpx,dp2,agg,qdum,fdum,qdp,tdum,dum5
      REAL*8 dbn,dbnt,dbnz,dbnp,dum,potnm,daa,dpp,dum1,dum2,dum3,dum4
      REAL*8 giac,gdum,sng,alnm,alnmp,cnm,snm
C============
ccc   perturbation calculation, given alpha, calculate potential
cc  Note that ideal modes must vanish at nq-m=0, but numerically they might 
cc  not.  Hence in the potential nq-m is given finite width, value sng.
cc  The value can be changed, depending on the accuracy of the modes.  Plot 
cc  the potential to make sure it is not singular. 
      sng = 1.D-2  !  sing width
      do 10 k = 1,nprt
      alp(k) = 0.D0
      dadp(k) = 0.D0
      dadt(k) = 0.D0
      dadz(k) = 0.D0
      padt(k) = 0.D0
      padt(k) = 0.D0
      pot(k) = 0.
      dptdt(k) = 0.
      dptdz(k) = 0.
      dptdp(k) = 0.
 10   continue
      do  312 md = md1,md2
      n = nmod(md)
      m = mmod(md)
      do 332 k = 1,jpts
         pdum = pol(k)
      jd = pol(k)*(lpt-1)/pw + 1
      jd = min(jd,lpt-1)
      jd = max(jd,1)
      dpx = pol(k) - (jd-1)*pw/(lpt-1)
      dp2 = dpx*dpx
      alnm = amp(md)*(a1(jd,md) + a2(jd,md)*dpx + a3(jd,md)*dp2)
      alnmp = amp(md)*( a2(jd,md) + 2*a3(jd,md)*dpx)
      agg = n*zet(k) - m*thet(k)
      cnm = cos(agg - phaz(k,md))
      snm = sin(agg - phaz(k,md))
      alp(k) = alp(k) + alnm*snm
      dadp(k) = dadp(k) + alnmp*snm
      dadt(k) = dadt(k) - m*alnm*cnm
      dadz(k) = dadz(k) + alnm*n*cnm
      padt(k) = padt(k) - omegv(md)*alnm*cnm
ccc      go to 332    !   No Potential
      dum = (n*q(k)-m)**2 + sng
      qdum = (n*q(k) - m)/dum   !  qdum is 1/(nq-m) with no singularity
      qdp = -n*qp(k)*(dum - 2*sng)/dum**2
            fdum= omegv(md)*(g(k)*q(k)+ri(k))*qdum
cccc- compute the potential for each particle
            pot(k)= pot(k) + alnm*snm*fdum
            dptdt(k) = dptdt(k) - m*alnm*cnm*fdum
            dptdz(k) = dptdz(k) + alnm*n*cnm*fdum
            dptdp(k) = dptdp(k) + alnmp*snm*fdum
     &  + alnm*snm*omegv(md)*(g(k)*qp(k)+gp(k)*q(k)+rip(k))*qdum
     &  + alnm*snm*omegv(md)*(g(k)*q(k)+ri(k))*qdp
cccccccccccccccccccccccccccccccccccccccccc
332   continue
 312  continue
      return
      end
ccccccccccccccccccccccccccc
ccccccccccccccccccccccccccc
      subroutine ptrbx(jpts)
      IMPLICIT NONE
      include 'orbcom'
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER jpts,k,md,n,m,jd,ndum
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 az,phnm,pdum,dpx,dp2,agg,qdum,qdp,tdum
      REAL*8 dbn,dbnt,dbnz,dbnp,dum,dpp,dum1,snm,cnm,alnm,alnmp
      REAL*8 giac,gdum,dump,xinm,xinmp,ptnm,ptnmp
C============
ccc   perturbation calculation, given xi, calculate alpha, potential
      do 10 k = 1,nprt
      alp(k) = 0.D0
      dadp(k) = 0.D0
      dadt(k) = 0.D0
      dadz(k) = 0.D0
      padt(k) = 0.D0
      padt(k) = 0.D0
      pot(k) = 0.
      dptdt(k) = 0.
      dptdz(k) = 0.
      dptdp(k) = 0.
 10   continue
      do  312 md = md1,md2
      n = nmod(md)
      m = mmod(md)
      do 332 k = 1,jpts
         pdum = pol(k)
      jd = pdum*(lpt-1)/pw + 1
      jd = min(jd,lpt-1)
      jd = max(jd,1)
      dpx = pdum - (jd-1)*pw/(lpt-1)
      dp2 = dpx*dpx
      xinm = amp(md)*(xi1(jd,md) + xi2(jd,md)*dpx + xi3(jd,md)*dp2)
      xinmp = amp(md)*(xi2(jd,md) + 2*xi3(jd,md)*dpx)
      dum = 1./(m*g(k) + n*ri(k))
      dump = -(m*gp(k) + n*rip(k))*dum**2
      dum1 = m - n*q(k)
      alnm = xinm*dum1*dum
      alnmp= xinmp*dum1*dum - n*qp(k)*xinm*dum + xinm*dum1*dump 
      agg = n*zet(k) - m*thet(k)
      cnm = cos(agg - phaz(k,md))
      snm = sin(agg - phaz(k,md))
cccc- compute alpha for each particle
      alp(k) = alp(k) + alnm*snm
      dadp(k) = dadp(k) + alnmp*snm
      dadt(k) = dadt(k) - m*alnm*cnm
      dadz(k) = dadz(k) + alnm*n*cnm
      padt(k) = padt(k) - omegv(md)*alnm*cnm
cccc- compute the potential for each particle
      ptnm = -(g(k)*q(k)+ri(k))*omegv(md)*xinm*dum
      pot(k) = pot(k) + ptnm*snm
      ptnmp = -(g(k)*q(k)+ri(k))*omegv(md)*xinmp*dum
     a   -(gp(k)*q(k) + g(k)*qp(k)+ rip(k))*omegv(md)*xinm*dum
     b   -(g(k)*q(k) + ri(k))*omegv(md)*xinm*dump
            dptdt(k) = dptdt(k) - m*ptnm*cnm
            dptdz(k) = dptdz(k) + ptnm*n*cnm
            dptdp(k) = dptdp(k) + ptnmp*snm
332   continue
 312  continue
      return
      end
ccccccccccccccccccccccccccc
      subroutine pt1(md,jpts)
      IMPLICIT NONE
      include 'orbcom'
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER jpts,md,k,jd,n,m
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 xinm,pdum,dpx,dp2,qdum,qdp,dum1,dum,dbn,alnm,ptnm
      REAL*8 gfun,rifun,qfun,gdum,rdum
C============
ccc   perturbation calculation for plot
      open(71,file='fix.plt',status='unknown')
      do 332 k = 1,jpts
         pdum = pol(k)
      n = nmod(md)
      m = mmod(md)
      jd = pol(k)*(lpt-1)/pw + 1
      jd = min(jd,lpt-1)
      jd = max(jd,1)
      dpx = pol(k) - (jd-1)*pw/(lpt-1)
      dp2 = dpx*dpx
      gdum = gfun(pdum)
      qdum = qfun(pdum)
      rdum = rifun(pdum)
      if(nptrbax.eq.1) then
      alnm = amp(md)*(a1(jd,md) + a2(jd,md)*dpx + a3(jd,md)*dp2)
      dum = m*gdum + n*rdum
      dum1 = m - n*qdum
      xinm = alnm*dum/dum1
      endif
      if(nptrbax.eq.2) then
      xinm = amp(md)*(xi1(jd,md) + xi2(jd,md)*dpx + xi3(jd,md)*dp2)
      dum = m*gdum + n*rdum
      dum1 = m - n*qdum
      alnm = xinm*dum1/dum
      endif
      ptnm = (gdum*qdum + rdum)*omegv(md)*alnm/(n*qdum - m)
      wk1(k,md) = alnm
      cur2(k,md) = alnm*mmod(md)/xx(k)
      cur3(k,md) = ptnm
      cur4(k,md) = xinm
      yy(k) = yy(k) + cur2(k,md)**2   !   mean square dB/B, all modes
      dum = pdum/pw
ccc       write(71,310) jd,pdum,xinm,dbn,alnm
 310     format(i4,1p20e12.4)
332   continue
      return
      end
cccccccccccccccccccc
ccccccccccccccccccccccccccc
      subroutine ptrb2(jpts)
      IMPLICIT NONE
      include 'orbcom'
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER jpts,k
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 dum
C============
cccccccccc-  Radial potential
      dum = pamp*engn/ekev
      do 10 k = 1,jpts
      pot(k) =  dum*(1 - pol(k)/pw)
      dptdp(k) = -dum/pw
      dptdt(k) = 0.D0
      dptdz(k) = 0.D0
 10   continue
      return
      end
cccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccc
      subroutine readptrbx
      IMPLICIT NONE
      include 'orbcom'
      REAL*8 dpx,dum,dum2,dbsmax,dnmax,omrat,ampdum,alimit
      REAL*8 thetd,zetd,sdum,agg,snmd,xx2,dp2,rdum,qmn
      REAL*8 qfun,qdum,pdum,rpol,giac,dum3
      INTEGER md,ndum,mdum,j,jm,jp,jpp,m,lptm,ldum,mload,k,jd,l,n
      INTEGER jdum,idum,kdum,mmin,mmax,nmd
C============
ccc-modified for file from Nikolai, 3/2012, read xi
      nval = 0
      modes = 0
      plabel = "Xin08w.1267E+01"
      open(61,file=plabel,status='unknown') 
      write(6,801) plabel
 801  format('  subroutine readptrbx, perturbation read=   ',A30)
      read(61,*) 
      read(61,*) 
      read(61,*) 
      read(61,*) lpt,nmd,mmin,mmax,dum,ndum
cc      write(6,*)  lpt,nmd,mmin,mmax,dum,ndum
      lptm = lpt - 1
      dpx = pw/lptm
      read (61,*) jdum,idum
ccc      write(6,*) jdum,idum
      do kdum = 1,2*idum
      read(61,*)
      enddo
      read (61,*) jdum,idum
      read (61,*) jdum,idum
      write(6,*) jdum,idum
ccccccccccccccccccccccccccccc
      read(61,*) ((xi1(j,md),j=1,lpt),md = 1,idum)
cccccccccccccccccccccccccccc
      modes = modes + mmax - mmin + 1
      nval = 1
      harm(nval) = mmax - mmin + 1
      do md = 1,harm(nval)
         alfv(md) = 1
         amp(md) = 1.D-3
         omegv(md) = 1.e-4
         nmod(md) = nmd
         mmod(md) = mmin - 1 + md
         enddo
 81   continue
      nvalx = nval
cccccccccccccccccccccccc
ccccc-  The perturbation harmonics are used only from md1 to md2
      md1 = 1
      md2 = modes
cccccccccccccccccccc  Select one mode
      if(nplot.eq.9.or.nplot.eq.8) then
      nvalx = 1
      nval = 1
      md1 = 1
      md2 = 1
      modes = md2 - md1 + 1
      endif
cccccccccc-renormalize 
      dum =  1.  ! amplitude renormalization
      dum2 = 1.  ! frequency renormalization
      write(6,57) dum,dum2,nval
 57   format('  change amp,freq, nval',1p2e12.4,i6)
            do 50 md = md1,md2
               amp(md) = amp(md)*dum   !   modify mode amplitude
               omegv(md) = omegv(md)*dum2     !  modify frequency
                dum3 = omegv(md)*omeg0/(2.D3*pi)
ccc               write(6,52) md,mmod(md),nmod(md),amp(md),dum3
 52               format(i4,' mode- m,n,amp, freq ',2i4,1p2e12.4)
ccc             write(6,121)(xi1(j,md),j=1,lpt)
 50            continue
               call splnx
               return
               end
ccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccc
      subroutine readptrba
      IMPLICIT NONE
      include 'orbcom'
      REAL*8 dpx,dum,dum2,dbsmax,dnmax,omrat,ampdum,alimit
      REAL*8 thetd,zetd,sdum,agg,snmd,xx2,dp2,rdum,qmn
      REAL*8 qfun,qdum,pdum,rpol,giac
      INTEGER md,ndum,mdum,j,jm,jp,jpp,m,lptm,ldum,mload,k,jd,l,n
C============
ccc-modified for file from Nikolai, read alpha
      plabel = "ptr1_sm_141711.dat"
      open(61,file=plabel,status='unknown') 
      write(6,801) plabel
 801  format('  subroutine readptrba, perturbation read=   ',A30)
      nval = 0
      read(61,*) 
      read(61,*) 
      read(61,*) lpt,mload
      lptm = lpt - 1
      dpx = pw/lptm
      nval = nval+1
      do  md=1,mload
         alfv(md) = nval
      read(61,*) 
      read(61,*) mmod(md),nmod(md),omrat,amp(md)
      omegv(md) = 13.69*omrat*6280/omeg0 
      read(61,111) (a1(j,md),j=1,lpt)
      enddo
      modes = modes + mload
      harm(nval) = mload
      open(62,file='ptr2_sm_141711.dat',status='unknown')
      read(62,*) 
      read(62,*) 
      read(62,*) lpt,mload
      lptm = lpt - 1
      dpx = pw/lptm
      nval = nval+1
      do  md = modes + 1,modes + mload
         alfv(md) = nval
      read(62,*) 
      read(62,*) mmod(md),nmod(md),omrat,amp(md)
      omegv(md) = 13.69*omrat*6280/omeg0 
      read(62,111) (a1(j,md),j=1,lpt)
      enddo
      modes = modes + mload
      harm(nval) = mload
  111 format(8e12.5)
      open(63,file='ptr3_sm_141711.dat',status='unknown')
      read(63,*) 
      read(63,*) 
      read(63,*) lpt,mload
      lptm = lpt - 1
      dpx = pw/lptm
      nval = nval+1
      do  md = modes + 1,modes + mload
         alfv(md) = nval
      read(63,*) 
      read(63,*) mmod(md),nmod(md),omrat,amp(md)
      omegv(md) = 13.69*omrat*6280/omeg0 
      read(63,111) (a1(j,md),j=1,lpt)
      enddo
      modes = modes + mload
      harm(nval) = mload
      open(64,file='ptr4_sm_141711.dat',status='unknown')
      read(64,*) 
      read(64,*) 
      read(64,*) lpt,mload
      lptm = lpt - 1
      dpx = pw/lptm
      nval = nval+1
      do  md = modes + 1,modes + mload
         alfv(md) = nval
      read(64,*) 
      read(64,*) mmod(md),nmod(md),omrat,amp(md)
      omegv(md) = 13.69*omrat*6280/omeg0 
      read(64,111) (a1(j,md),j=1,lpt)
      enddo
      modes = modes + mload
      harm(nval) = mload
      nvalx = nval
cccccccccccccccccccccccc
ccccc-  The perturbation harmonics are used only from md1 to md2
      md1 = 1
      md2 = modes
cccccccccccccccccccc  Select one mode
      if(nplot.eq.9.or.nplot.eq.8) then
      nvalx = 1
      nval = 1
      md1 = 1
      md2 = 14
      modes = md2 - md1 + 1
      endif
cccccccccc-renormalize 
      dum =  1.   ! mode renormalization
      dum2 = 1.
      write(6,57) dum,dum2,nval
 57   format('  change amp,freq, nval',1p3e12.4,i6)
            do 50 md = md1,md2
               amp(md) = amp(k)*dum   !   modify mode amplitude
               omegv(md) = omegv(k)*dum2          !  modify frequency
ccc               write(6,52) md,nmod(md),mmod(md),amp(md)
 52               format(' mode- n,m,amp ',3i4,1pe12.4)
 50            continue
ccccccccccccccc- now spline
               call splna
      return
      end
ccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine splna
      IMPLICIT NONE
      include 'orbcom'
C============
ccc-   Spline the alpha representation
ccc
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER lptm,md,j,jm,jp,jpp,m
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 dum,dpx
C============
cccc
      lptm = lpt - 1
         m = mmod(md)
      dpx = pw/lptm
      do 25 md = md1,md2
cccc  set a2(1,md)
            a2(1,md) = (10*a1(2,md) - 7*a1(1,md) - 3*a1(3,md))/(4*dpx)
            if(m.ne.1) a2(1,md) = 0
cccc
            do 100 j = 2,lptm
               jm = j - 1
               jp = j + 1
               jpp = min(j + 2,lpt)
               a2(j,md) = -a2(jm,md) + 2*(a1(j,md)-a1(jm,md))/dpx
cccc  smooth a1
               a1(jp,md) = .4D0*dpx*a2(j,md) + .3D0*a1(jpp,md)
     &               + .7D0*a1(j,md)
 100        continue
               do 105 j = 1,lptm
                  jp = j + 1
                  a3(j,md) = (a2(jp,md)-a2(j,md))/(2*dpx)
 105           continue
 25            continue
ccc
ccc-a1,a2,a3-finished
               nptrbax = 1  !  use alpha , calc potential
          return
          end
ccccccccccccccccccccccccccc
      subroutine splnx
      IMPLICIT NONE
      include 'orbcom'
      INTEGER lptm,md,n,m,j,jm,jp,jpp
      REAL*8 dum,dpx
      lptm = lpt - 1
      dpx = pw/lptm
ccccccccccccccc- now spline
      do 25 md = md1,md2
         m = mmod(md)
cccc
cccc  set xi2(1,md)
           xi2(1,md)=(10*xi1(2,md) - 7*xi1(1,md) - 3*xi1(3,md))/(4*dpx)
            if(m.ne.1) xi2(1,md) = 0
cccc
            do 100 j = 2,lptm
               jm = j - 1
               jp = j + 1
               jpp = min(j + 2,lpt)
               xi2(j,md) = -xi2(jm,md) + 2*(xi1(j,md)-xi1(jm,md))/dpx
cccc  smooth xi1
               xi1(jp,md) = .4D0*dpx*xi2(j,md) + .3D0*xi1(jpp,md)
     &               + .7D0*xi1(j,md)
 100        continue
               do 105 j = 1,lptm
                  jp = j + 1
                  xi3(j,md) = (xi2(jp,md)-xi2(j,md))/(2*dpx)
 105           continue
 25            continue
ccc
ccc-xi1,xi2,xi3-finished
               nptrbax = 2  !  use xi, calc alpha, potential
      return
      end
ccccccccccccccccccccccccccc
