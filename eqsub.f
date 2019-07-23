cccccccccccccccccccccccccccccccccccccc
      subroutine wdat
      IMPLICIT NONE
      include 'orbcom'
      include 'o.cln'
cccc writes data file to be read by ORBIT, version orbit.f
cccc- April 21, 1999
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER j,l
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 dum1,dum2
C============
      write(20,760) mp0,'version_05'
 760  format(a10,10x,a10)
      write(20,770) lsp,lst,lemax,lrmax
 770  format(6i4)
      write(20,777) pw,ped
 777        format(1p4e18.10)
         do 2 j = 1,lsp
            write(20,777) (b1(j,l),l=1,lst)
            write(20,777) (b2(j,l),l=1,lst)
            write(20,777) (b3(j,l),l=1,lst)
            write(20,777) (b4(j,l),l=1,lst)
            write(20,777) (b5(j,l),l=1,lst)
            write(20,777) (b6(j,l),l=1,lst)
            write(20,777) (b7(j,l),l=1,lst)
            write(20,777) (b8(j,l),l=1,lst)
            write(20,777) (b9(j,l),l=1,lst)
            write(20,777) (x1(j,l),l=1,lst)
            write(20,777) (x2(j,l),l=1,lst)
            write(20,777) (x3(j,l),l=1,lst)
            write(20,777) (x4(j,l),l=1,lst)
            write(20,777) (x5(j,l),l=1,lst)
            write(20,777) (x6(j,l),l=1,lst)
            write(20,777) (x7(j,l),l=1,lst)
            write(20,777) (x8(j,l),l=1,lst)
            write(20,777) (x9(j,l),l=1,lst)
            write(20,777) (z1(j,l),l=1,lst)
            write(20,777) (z2(j,l),l=1,lst)
            write(20,777) (z3(j,l),l=1,lst)
            write(20,777) (z4(j,l),l=1,lst)
            write(20,777) (z5(j,l),l=1,lst)
            write(20,777) (z6(j,l),l=1,lst)
            write(20,777) (z7(j,l),l=1,lst)
            write(20,777) (z8(j,l),l=1,lst)
            write(20,777) (z9(j,l),l=1,lst)
            write(20,777) (g1(j,l),l=1,lst)
            write(20,777) (g2(j,l),l=1,lst)
            write(20,777) (g3(j,l),l=1,lst)
            write(20,777) (g4(j,l),l=1,lst)
            write(20,777) (g5(j,l),l=1,lst)
            write(20,777) (g6(j,l),l=1,lst)
            write(20,777) (g7(j,l),l=1,lst)
            write(20,777) (g8(j,l),l=1,lst)
            write(20,777) (g9(j,l),l=1,lst)
            write(20,777) qd1(j),qd2(j),qd3(j)
            write(20,777) gd1(j),gd2(j),gd3(j)
            write(20,777) rd1(j),rd2(j),rd3(j)
            write(20,777) pd1(j),pd2(j),pd3(j)
            write(20,777) rp1(j),rp2(j),rp3(j)
            write(20,777) ps1(j),ps2(j),ps3(j)
 2       continue
         write(20,770) krip,nrip
         dum1 = wrip*rmaj/xc
         dum2 = xrip*rmaj/xc
      write(20,777) rmaj,d0,brip
      write(20,777) dum1,dum2
         do 16 j = 1,lsp
            write(20,777) (r1(j,l),l=1,lst)
            write(20,777) (r2(j,l),l=1,lst)
            write(20,777) (r3(j,l),l=1,lst)
            write(20,777) (r4(j,l),l=1,lst)
            write(20,777) (r5(j,l),l=1,lst)
            write(20,777) (r6(j,l),l=1,lst)
            write(20,777) (r7(j,l),l=1,lst)
            write(20,777) (r8(j,l),l=1,lst)
            write(20,777) (r9(j,l),l=1,lst)
 16      continue
         return
         end
cccccccccccccccccccccccccccccccccccccc
      subroutine polrep
      IMPLICIT NONE
      include 'orbcom'
      include 'o.cln'
clin  ispln=1: asymptotic match at pol=0. If ispln=0: original method.
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER j,k,lspm
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 px,qfun,gfun,rifun,pfun,psfun,b0,x0,g0,z0,tx,bfield,
     > xproj,zproj,giac,sdum,dpx
C============
      integer ispln
      ispln=1
 
      if(numeric.eq.0) go to 75
      if(numeric.eq.2) go to 75
clin   find q function first, q0=qd1(1)
      do j = 1,lsp
         px = (j-1)*pw/(lsp-1)
         qd1(j) = qfun(px)
      enddo
clin   find g function
      do j = 1,lsp
         px = (j-1)*pw/(lsp-1)
         gd1(j) = gfun(px)
      enddo
clin   find ri function
      do j = 1,lsp
         px = (j-1)*pw/(lsp-1)
         rd1(j) = rifun(px)
      enddo
clin   find p function
      do j = 1,lsp
         px = (j-1)*pw/(lsp-1)
         pd1(j) = pfun(px)
      enddo
clin   find psi function
      do j = 1,lsp
         px = (j-1)*pw/(lsp-1)
         ps1(j) = psfun(px)
      enddo
clin  initially use averaged z0,x0,b0,g0, modified in subroutine spln
      b0=0.0D0
      x0=0.0D0
      g0=0.0D0
      z0=0.0D0
      do k=1,lst
         tx = DBLE(k-1)*pi2/DBLE(lst)
         b0=b0+bfield(0.0D0,tx)
         x0=x0+xproj(0.0D0,tx)
         z0=z0+zproj(0.0D0,tx)
         g0=g0+giac(0.0D0,tx)
      enddo
      b0=b0/DBLE(lst)
      x0=x0/DBLE(lst)
      z0=z0/DBLE(lst)
      g0=g0/DBLE(lst)
 
clin  use right-handed coordinate, define z=z0+r*sin(theta)
      sdum=sign(1.0D0,zproj(0.1D0*pw,0.5D0*pi)-z0)
      z0=sdum*z0
      do j = 2,lsp
         do k = 1,lst
            px = (j-1)*pw/(lsp-1)
            tx = (k-1)*pi2/lst
            z1(j,k) = sdum*zproj(px,tx)
            x1(j,k) = xproj(px,tx)
            b1(j,k) = bfield(px,tx)
            g1(j,k) = giac(px,tx)
         enddo
      enddo
      go to 76
 
 75   continue
cccc  -analytic equilibrium
      ipos = 4
      lspm = lsp - 1
clin  magnetic axis quantities
      b0=b1(1,1)
      x0=x1(1,1)
      g0=g1(1,1)
      z0=z1(1,1)
 76   continue
 
clin  calculation of initial f2(1,k), modified in subroutine spln
clin  z,x,b,g: asymptotic expansion near the magnetic axis
      q0=qd1(1)
      dpx=pw/DBLE(lsp-1)
 
clin  match z function: z=z0+r*sin(theta)
      do k=1,lst
         tx = DBLE(k-1)*pi2/DBLE(lst)
         z1(1,k)=z0
         z2(1,k)=sqrt(2.0D0*q0/b0)*sin(tx)
clin match x function: x=x0+r*cos(theta)
         x1(1,k)=x0
         x2(1,k)=sqrt(2.0D0*q0/b0)*cos(tx)
clin  match b function: b*x=constant
         b1(1,k)=b0
         b2(1,k)=-sqrt(2.0D0*q0*b0)*cos(tx)/x0
      enddo
      call spln(ispln,z1,z2,z3,z4,z5,z6,z7,z8,z9)
      call spln(ispln,x1,x2,x3,x4,x5,x6,x7,x8,x9)
      call spln(ispln,b1,b2,b3,b4,b5,b6,b7,b8,b9)
ccc
 
clin  Jacobian uses old method
      do k=1,lst
         g1(1,k)=g0
         g2(1,k) = (10.D0*g1(2,k) - 7.D0*g1(1,k) - 3.D0*g1(3,k))/(4.D0*
     >     dpx)
      enddo
      call spln(0,g1,g2,g3,g4,g5,g6,g7,g8,g9)
      xc = x1(1,1)
      eps = x1(lsp,1)/xc - 1.D0
      bax = b1(1,1)
clin   find rad function along boozer angle zero
      do j = 1,lsp
         px = (j-1)*pw/(lsp-1)
         rp1(j) = sqrt((x1(j,1) - xc)**2 + (z1(j,1) - z1(1,1))**2)/xc
      enddo
      call splin1(qd1,qd2,qd3)
      call splin1(gd1,gd2,gd3)
      call splin1(rd1,rd2,rd3)
      call splin1(pd1,pd2,pd3)
      call splin1(rp1,rp2,rp3)
      call splin1(ps1,ps2,ps3)
ccccc
         return
         end
cccccccccccccccccccccccccccccccccccccc
ccccccccccccc
      subroutine wrtb
      IMPLICIT NONE
      include 'o.cln'
      include 'orbcom'
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER j,k
C============
      do 20 j = 1,lsp
         do 20 k = 1,lst
            write(6,21) j,k,b1(j,k)
 21         format('j,k,b1 ',2i4,1pe12.4)
            write(6,22) j,k,b2(j,k)
 22         format('j,k,b2 ',2i4,1pe12.4)
            write(6,23) j,k,b3(j,k)
 23         format('j,k,b3 ',2i4,1pe12.4)
            write(6,24) j,k,b4(j,k)
 24         format('j,k,b4 ',2i4,1pe12.4)
            write(6,25) j,k,b5(j,k)
 25         format('j,k,b5 ',2i4,1pe12.4)
            write(6,26) j,k,b6(j,k)
 26         format('j,k,b6 ',2i4,1pe12.4)
            write(6,27) j,k,b7(j,k)
 27         format('j,k,b7 ',2i4,1pe12.4)
            write(6,28) j,k,b8(j,k)
 28         format('j,k,b8 ',2i4,1pe12.4)
            write(6,29) j,k,b9(j,k)
 29         format('j,k,b9 ',2i4,1pe12.4)
 20      continue
         return
         end
ccccccccccccccccccccccccccccccccccccccccccc
      function pold(rdum)
      IMPLICIT NONE
      include 'orbcom'
ccc..poloidal flux at rdum
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER k
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 rdum,pold,dum,rd,dr,qd,qfnr
C============
      dum = 0
      rd = 0
      dr = .001D0*rdum
      do 10 k = 1,1000
       qd = qfnr(rd)
       dum = dum + rd*dr/qd
       rd = rd + dr
 10   continue
      pold = dum
      return
      end
ccccccccccccc cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccc
      function psid(pz)
      IMPLICIT NONE
      include 'orbcom'
ccc..toroidal flux at pz
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER k
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 pz,psid,dum,rd,dpz,qd,qfnr
C============
      dum = 0
      rd = 0
      dpz = .001D0*pz
      do 10 k = 1,1000
       qd = qfnr(rd)
       dum = dum + qd*dpz
       rd = sqrt(2*dum)
 10   continue
      psid = dum
      return
      end
ccccccccccccc cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccc
      function rfun(pz)
      IMPLICIT NONE
      include 'orbcom'
ccc...radius, function of poloidal flux
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER k
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 pz,rfun,dum,dpz,rd,qd,qfnr
C============
      dum = 0
      dpz = .001D0*pz
      do 10 k = 1,1000
       rd = sqrt(dum)
       qd = qfnr(rd)
      dum = dum + 2*qd*dpz
 10   continue
      rfun = sqrt(dum)
      return
      end
ccccccccccccc cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function gfan(rd)
      IMPLICIT NONE
      include 'orbcom'
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER k
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 rd,gfan,dr,rdum,dum,d1,ppfun,qd,qfnr,qpd,qfnrp,d2,d3
C============
      dr = .01D0*rd
      rdum = 0
      dum = 0
      do 10 k = 1,100
        rdum = rdum + dr
        d1 = -2*ppfun(rd)
        qd = qfnr(rdum)
        qpd = qfnrp(rdum)
        d2 = 2*rdum**2*qpd/qd**3
        d3 = -4*rdum/qd**2
        dum = dum + (d1 + d2 + d3)*dr
 10     continue
        gfan = sqrt(1 + dum)
      return
      end
ccccccccccccc
      subroutine splin1(f1,f2,f3)
      IMPLICIT NONE
      include 'orbcom'
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER lspm,lsp2,j,jm,jp,jpp
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 f2,f3,f1,dpx
C============
      dimension f1(idp),f2(idp),f3(idp)
ccc   perturbation spline calculation, dp = pol - pj,
ccc     grid given by    pj = (j-1)*pw/(lsp-1)
cccc-in domain j,  f(pol) = f1 + f2*dp + f3*dp**2
cccc  before calling splin1 load function f into f1(j)=f(pj)
cccc   after call,  f2(j),f3(j) will be determined
      if(lsp.gt.idp) write(*,22)
 22   format('  idp too small')
      if(lsp.gt.idp) stop
      lspm = lsp - 1
      lsp2 = lsp - 2
      dpx = pw/lspm
ccccc
cccc
cccc  set f2(1) to leave f1(2) unmoved by smoothing
            f2(1) = (10*f1(2) - 7*f1(1) - 3*f1(3))/(4*dpx)
 3       continue
cccc
            do 100 j = 2,lspm
               jm = j - 1
               jp = j + 1
               jpp = min(j + 2,lsp)
               f2(j) = -f2(jm) + 2*(f1(j)-f1(jm))/dpx
cccc  smooth f1
       if(jp.ne.lsp) f1(jp) =.4D0*dpx*f2(j)+.3D0*f1(jpp)+.7D0*f1(j) ! 4/23/99
 100        continue
              f2(lsp) = f2(lspm)
               do 105 j = 1,lspm
                  jp = j + 1
                  f3(j) = (f2(jp)-f2(j))/(2*dpx)
 105           continue
      f3(lspm) = (f1(lsp)-f1(lspm)-f2(lspm)*dpx)/dpx**2 ! 4/23/99
ccc-f1,f2,f3-finished
ccc
               return
               end
cccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccc
      subroutine spln(ispln,f1,f2,f3,f4,f5,f6,f7,f8,f9)
      IMPLICIT NONE
      include 'orbcom'
clin  ispln=1: asymptotic match at pol=0. If ispln=0: original method.
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER ld,lspm,k,j,jm,jp,jpp,lmax,l,lmax2,kp,ier,jd
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 f2,f3,f4,f5,f6,f7,f8,f9,f1,dpx,dtx,f0,err0
C============
      integer ispln
      dimension f1(idp,idt),f2(idp,idt),f3(idp,idt),f4(idp,idt)
     a  ,f5(idp,idt),f6(idp,idt),f7(idp,idt),f8(idp,idt),f9(idp,idt)
ccc   perturbation spline calculation, dp = pol - pj, dt = thet - tk
ccc     grid given by  tk = (k-1)*pi2/lst,  pj = (j-1)*pw/(lsp-1)
cccc-in domain j,k  f(pol,thet) = f1 + f2*dp + f3*dp**2 + f4*dt + f5*dt*dp
cccc  + f6*dt*dp**2 + f7*dt**2 + f8*dt**2*dp + f9*dt**2*dp**2
ccc-  gelg solves bmat(l,m)x(m) = wk2(m), m=1,ndim
ccc   wk2 contains x after the call
cccc  bmat stored as vector bmat(lm),lm = l + (m-1)*ndim
cccc  before calling spln load function f into f1(j,k)=f(pj,tk)
cccc   after call,  f2(j,k),f3(j,k) etc will be determined
      ld = mod(lst,2)
      if(ld.ne.1) stop          !  must be odd
      if(lsp.gt.idp) stop
      if(lst.gt.idt) write(6,22)
 22   format('  idt too small')
      if(lst.gt.idt) stop
      lspm = lsp - 1
      dpx = pw/lspm
      dtx = pi2/lst
ccccc
ccc
      if(ispln .eq. 1)then
clin     smooth sqrt(f)*df/dp
clin     first adjust f1(1,k), use average as on axis quantites.
         f0=0.0D0
         do k=1,lst
            f1(1,k)=2.0D0*f1(2,k)-f1(3,k)
     &           -(0.75D0-sqrt(1.0D0/32.0D0))*sqrt(dpx)*f2(1,k)
            f0=f0+f1(1,k)
         enddo
         f0=f0/DBLE(lst)
 
         do k=1,lst
            f1(1,k)=f0
clin    adjust f2(1,k) to enforce continuity of second order derivative
            f2(1,k)=(2.0D0*f1(2,k)-f1(1,k)-f1(3,k))
     &           *8.0D0/(6.0D0-sqrt(2.0D0))/sqrt(dpx)
         enddo
      else
         do k=1,lst
            f2(1,k) = (10.D0*f1(2,k) - 7.D0*f1(1,k) - 3.D0*f1(3,k))/
     >     (4.D0*dpx)
         enddo
      endif
 
      do k = 1,lst
         if(ispln .eq. 1)then
clin  match f2(2,k) at the first flux surface and smooth the second one
            f2(2,k)=-f2(1,k)*0.5D0/sqrt(dpx)+(f1(2,k)-f1(1,k))/dpx
            f1(3,k)=.4D0*dpx*f2(2,k)+.3D0*f1(4,k)+.7D0*f1(2,k)
         endif
         do j = 2+ispln,lspm
            jm = j - 1
            jp = j + 1
            jpp = min(j + 2,lsp)
            f2(j,k) = -f2(jm,k) + 2*(f1(j,k)-f1(jm,k))/dpx
cccc  smooth f1
         if(jp.ne.lsp)    ! 4/23/99
     a      f1(jp,k) = .4D0*dpx*f2(j,k) + .3D0*f1(jpp,k) + .7D0*f1(j,k)
         enddo
      enddo
      do k = 1,lst
         f2(lsp,k) = f2(lspm,k)
clin  match f3(1,k) at the first flux surface
         if(ispln .eq. 1)f3(1,k)=f2(2,k)-f2(1,k)*0.5D0/sqrt(dpx)
         do j = 1+ispln,lspm
            jp = j + 1
            f3(j,k) = (f2(jp,k)-f2(j,k))/(2*dpx)
         enddo
      f3(lspm,k) = (f1(lsp,k)-f1(lspm,k)-f2(lspm,k)*dpx)/dpx**2 ! 4/23/9
      enddo
ccc-f1,f2,f3-finished
ccc
ccc
ccc   find matrix for f4
            lmax = lst
         do 300 j = 1,lspm
ccc
cccc-clear wk2
            do 51 l = 1,lmax
            wk2(l) = 0
 51      continue
         lmax2 = lmax*lmax
ccc
ccccc-right hand side
         do 110 k = 1,lst
            kp = k + 1
            if(k.eq.lst) kp = 1
            wk2(k) = 2*f1(j,kp) - 2*f1(j,k)
 110     continue
cccc            write(6,177) (wk2(jd),jd=1,lmax)
 177         format('f4 rt ',1p5e10.2)
         call ldbmat
         do 179 ld = 1,lmax
cccc            write(6,178)(bmat(lm),lm=ld,lmax2,lmax)
 178         format( 'f4-matrix ',1p8e10.2)
 179      continue
                 err0 = 1.D-6
               call gelg(wk2,bmat,lmax,1,err0,ier)
               do 210 k = 1,lst
                  f4(j,k) = wk2(k)
 210           continue
               do 211 k = 1,lst
                  kp = k + 1
                  if(k.eq.lst) kp = 1
                  f7(j,k) = (f4(j,kp) - f4(j,k))/(2*dtx)
 211           continue
 300        continue
ccc
ccc-f4-f7-finished
cccccc
ccc   find matrix for f5
         do 400 j = 1,lspm
ccc
cccc-clear wk2,bmat
            do 52 l = 1,lmax
            wk2(l) = 0
 52      continue
ccc
ccccc-right hand side
         do 120 k = 1,lst
            kp = k + 1
            if(k.eq.lst) kp = 1
            wk2(k) = 2*f2(j,kp) - 2*f2(j,k)
 120     continue
cccc            write(6,277) (wk2(jd),jd=1,lmax)
 277         format('f5 rt ',1p5e10.2)
         call ldbmat
                  err0 = 1.D-6
               call gelg(wk2,bmat,lmax,1,err0,ier)
               do 220 k = 1,lst
                  f5(j,k) = wk2(k)
 220           continue
               do 221 k = 1,lst
                  kp = k + 1
                  if(k.eq.lst) kp = 1
                  f8(j,k) = (f5(j,kp) - f5(j,k))/(2*dtx)
 221           continue
 400        continue
ccc
ccc-f5-f8-finished
ccc
ccc   find matrix for f6
         do 500 j = 1,lspm
ccc
cccc-clear wk2,bmat
            do 53 l = 1,lmax
            wk2(l) = 0
 53      continue
ccc
ccccc-right hand side
         do 130 k = 1,lst
            kp = k + 1
            if(k.eq.lst) kp = 1
            wk2(k) = 2*f3(j,kp) - 2*f3(j,k)
 130     continue
cccc            write(6,377) (wk2(jd),jd=1,lmax)
 377         format('f6 rt ',1p5e10.2)
         call ldbmat
                  err0 = 1.D-6
               call gelg(wk2,bmat,lmax,1,err0,ier)
               do 230 k = 1,lst
                  f6(j,k) = wk2(k)
 230           continue
               do 231 k = 1,lst
                  kp = k + 1
                  if(k.eq.lst) kp = 1
                  f9(j,k) = (f6(j,kp) - f6(j,k))/(2*dtx)
 231           continue
 500        continue
ccc
ccc-f6-f9-finished
ccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ldbmat
      IMPLICIT NONE
      include 'orbcom'
cccc  bmat stored as vector bmat(lm),lm = l + (m-1)*lmax
cccc
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER lmax,lmax2,lm,l,m,mm,lmm
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 dtx
C============
      dtx = pi2/lst
            lmax = lst
ccc
cccc-clear bmat
         lmax2 = lmax*lmax
         do 61 lm = 1,lmax2
            bmat(lm) = 0
 61      continue
ccc
cccc  now matrix-stored-column by -column
         do 115 l = 1,lmax
            m = l
            mm = m + 1
            if(m.eq.lmax ) mm = 1
           lm = l + (m-1)*lmax
           lmm = l + (mm-1)*lmax
                  bmat(lm) = dtx
                  bmat(lmm) = dtx
 115            continue
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccc
      function xspl(px,tx)
      IMPLICIT NONE
      include 'orbcom'
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER jd,idum,kd
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 px,tx,xspl,dpx,dp2,tdum,dtx,dt2,dum
C============
      jd = px*(lsp-1)/pw + 1
      jd = min(jd,lsp-1)
      dpx = px - (jd-1)*pw/(lsp-1)
      dp2 = dpx*dpx
      idum = tx/pi2
      tdum = tx - pi2*(idum-1)
      idum = tdum *.1591549431D0
      tdum = tdum - pi2*idum
      kd = tdum*lst/pi2 + 1
      dtx = tdum - (kd-1)*pi2/lst
      dt2 = dtx*dtx
      dum = x1(jd,kd) + x2(jd,kd)*dpx + x3(jd,kd)*dp2
     a  + x4(jd,kd)*dtx + x5(jd,kd)*dpx*dtx + x6(jd,kd)*dtx*dp2
     b  + x7(jd,kd)*dt2 + x8(jd,kd)*dt2*dpx + x9(jd,kd)*dt2*dp2
ccc      write(6,3) px,tx,jd,kd,x1(jd,kd),x2(jd,kd),x3(jd,kd)
 3    format('px,tx,j,k,x1,x2,x3 ',1p2e12.4,2i4,1p3e12.4)
      xspl = dum
      end
cccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine defolt(netcdf)
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER istat,fshell,netcdf
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 pi,pi2,pi2i
C============
      character*50 string
      character*7 fort30,fort31
      include 'o.cln'
      character*60 mp0, mp1
      common / data / mp0, mp1
      pi = 4.D0*atan(1.D0)
      pi2 = 2*pi
      pi2i = 1/pi2
      twopi = 2.0D0*pi
      mth1 = mth + 1
      mth2 = mth + 2
      mth3 = mth + 3
      mth4 = mth + 4
c
c S.ETHIER 05/11/00  Use fort.30 and fort.31 only if using BZIO
      if (netcdf .ne. 1) then
ccc      copy equilibrium files to fort.30 fort.31
         fort30 = "fort.30"
         fort31 = "fort.31"
         write(string,115) mp0,fort30
         istat=fshell (string)
         write(string,115) mp1,fort31
 115     format('cp ',a10,a10)
         istat=fshell (string)
      endif
      indat = 15
      outdat = 6
      iomode = 30
      outmap1 = 31
      iotty = 59
      return
      end
ccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccc
      subroutine dskin(rmaj)
      IMPLICIT NONE
      include 'o.cln'
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER ndsk,nmap1,nmpdsk,nadres,i,lgivup,length,j
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 rmaj,ss,psi00
      character*60 mp0, mp1
      common / data / mp0, mp1
C============
       ndsk = 1
      call zop ( outmap1, mp1, nmap1, ndsk, ss,  999 )
      ndsk = 3
      call zop ( iomode, mp0, nmpdsk, ndsk, ss, 999 )
      nadres = 1
      call zrd ( outmap1,ntitle8(1),42,nadres,1,999 )
      ntitle = ntitle8
      nx = nx8
      nz = nz8
      nosurf = nosurf8
      mth = mth8
      lpless = lpless8
ccc      write(outdat,50) lj, mj, nj
   50 format ( 1x, "lj/mj/nj = ",i5,"/",i5,"/",i5 )
      mth1 = mth + 1
      mth2 = mth + 2
      mth3 = mth + 3
      mth4 = mth + 4
      dth = twopi / mth
      do 60 i = 1, mth1
   60 theval(i) = (i-1) * dth
      rx = r*r
      lgivup = 1
      nadres = 50 + 2*mth2
      call zrd(iomode,psival(1),nosurf,nadres,lgivup,999)
      nadres = nadres + nosurf
      length = nths*nsf
      call zrd(iomode,grpssq(1,1),length,nadres,lgivup,999)
      nadres = nadres + 8*length
      call zrd(iomode,xjacob(1,1),length,nadres,lgivup,999)
      nadres = 50
      call zrd(outmap1,pdat(1),nosurf,nadres,lgivup,999)
      nadres = nadres + nosurf
      call zrd(outmap1,ppdat(1),nosurf,nadres,lgivup,999)
      nadres = nadres + nosurf
      call zrd(outmap1,qdat(1),nosurf,nadres,lgivup,999)
      nadres = nadres + nosurf
      call zrd(outmap1,qpdat(1),nosurf,nadres,lgivup,999)
      nadres = nadres + nosurf
      call zrd(outmap1,gdat(1),nosurf,nadres,lgivup,999)
      nadres = nadres + 4*nosurf
      call zrd(outmap1,x(1,1),length,nadres,lgivup,999)
      nadres = nadres + length
      call zrd(outmap1,z(1,1),length,nadres,lgivup,999)
      nadres = nadres + length
      call zrd(outmap1,xdth(1,1),length,nadres,lgivup,999)
      nadres = nadres + length
      call zrd(outmap1,zdth(1,1),length,nadres,lgivup,999)
      nadres = nadres + length
      call zrd(outmap1,xpsi(1,1),length,nadres,lgivup,999)
      nadres = nadres + length
      call zrd(outmap1,zpsi(1,1),length,nadres,lgivup,999)
      psi00 = psival(1) / twopi
      do 200 j = 1, nosurf
      psival(j) = psival(j) / twopi - psi00
      gdat(j) = r * gdat(j)
  200 continue
      delpsi = psival(nosurf) - psival(1)
      dpsi = delpsi / (nosurf-1)
      rmaj = 100*xma
      return
  999 write ( outdat,'(5hdskin)' )
      end
ccccccccccccccccccccccccccccccccccccccccccc
      subroutine funcij
      IMPLICIT NONE
      include 'o.cln'
      include 'orbcom'
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER j,i
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 bsq
C============
      do 100 j = 1, nosurf
      cur(j) = 0.0D0
      do 80 i = 1, mth1
      bsq = ( grpssq(i,j) + gdat(j)**2 ) / x(i,j)**2
      bdat(i,j) = sqrt(bsq)
      if ( i .eq. mth1 ) go to 80
      cur(j) = cur(j) + xjacob(i,j)*grpssq(i,j)/x(i,j)**2
   80 continue
      cur(j) = cur(j) * dth / twopi
      bdat(mth2,j) = bdat(2,j)
      bdat(mth3,j) = bdat(3,j)
      bdat(mth4,j) = bdat(4,j)
 100  continue
      pi = 4.D0*atan(1.D0)
      pi2 = 2*pi
      pi2i = 1/pi2
      ped = psival(nosurf)
      return
      end
cccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccc
      subroutine contin
      IMPLICIT NONE
      include 'o.cln'
      include 'orbcom'
ccc   check continuity of spline representation of b
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER lspm,j,jm,k,km
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 dpx,dtx,dum
C============
      lspm = lsp-1
      dpx = pw/lspm
      dtx = pi2/lst
      write(6,5)
 5    format(' b and dbdt to left ')
      write(6,6)
 6    format(' b1(j,k) - b1(jm,k) - b2(jm,k)*dpx - b3(jm,k)*dpx*dpx ')
      do 20 j = 2,lsp
         jm = j - 1
         do 20 k = 1,lst
            km = k - 1
            if(k.eq.1) km = lst
            dum = b1(j,k) - b1(jm,k) - b2(jm,k)*dpx - b3(jm,k)*dpx*dpx
            if(abs(dum).gt.1.D-10) write(6,600) j,k,dum
 600        format(' j,k ',2i4,1pe12.4)
 20      continue
      write(6,7)
 7    format(' b4(j,k) - b4(jm,k) - b5(jm,k)*dpx - b6(jm,k)*dpx*dpx ')
      do 21 j = 2,lspm
         jm = j - 1
         do 21 k = 1,lst
            km = k - 1
            if(k.eq.1) km = lst
            dum =  b4(j,k) - b4(jm,k) - b5(jm,k)*dpx - b6(jm,k)*dpx*dpx
            if(dum.gt.1.D-10) write(6,600) j,k,dum
 21      continue
      write(6,8)
 8    format(' b7(j,k) - b7(jm,k) - b8(jm,k)*dpx - b9(jm,k)*dpx*dpx ')
      do 22 j = 2,lspm
         jm = j - 1
         do 22 k = 1,lst
            km = k - 1
            if(k.eq.1) km = lst
            dum = b7(j,k) - b7(jm,k) - b8(jm,k)*dpx - b9(jm,k)*dpx*dpx
            if(abs(dum).gt.1.D-10) write(6,600) j,k,dum
 22      continue
      write(6,9)
 9    format(' b and dbdp down ')
      write(6,10)
 10   format(' b1(j,k) - b1(j,km) - b4(j,km)*dtx - b7(j,km)*dtx*dtx ')
      do 30 j = 1,lspm
         jm = j - 1
         do 30 k = 1,lst
            km = k - 1
            if(k.eq .1D0) km = lst
            dum =  b1(j,k) - b1(j,km) - b4(j,km)*dtx - b7(j,km)*dtx*dtx
            if(abs(dum).gt.1.D-10) write(6,600) j,k,dum
 30      continue
      write(6,11)
 11   format(' b2(j,k) - b2(j,km) - b5(j,km)*dtx - b8(j,km)*dtx*dtx ')
      do 31 j = 1,lspm
         jm = j - 1
         do 31 k = 1,lst
            km = k - 1
            if(k.eq.1) km = lst
            dum =  b2(j,k) - b2(j,km) - b5(j,km)*dtx - b8(j,km)*dtx*dtx
            if(abs(dum).gt.1.D-10) write(6,600) j,k,dum
 31      continue
      write(6,12)
 12   format(' b3(j,k) - b3(j,km) - b6(j,km)*dtx - b9(j,km)*dtx*dtx ')
      do 32 j = 1,lspm
         jm = j - 1
         do 32 k = 1,lst
            km = k - 1
            if(k.eq.1) km = lst
            dum = b3(j,k) - b3(j,km) - b6(j,km)*dtx - b9(j,km)*dtx*dtx
            if(abs(dum).gt.1.D-10) write(6,600) j,k,dum
 32      continue
      write(6,13)
 13   format(' dbdp left ')
      write(6,14)
 14   format(' b2(j,k) - b2(jm,k) - 2*b3(jm,k)*dpx ')
      do 40 j = 2,lspm
         jm = j - 1
         do 40 k = 1,lst
            km = k - 1
            if(k.eq.1) km = lst
            dum =  b2(j,k) - b2(jm,k) - 2*b3(jm,k)*dpx
            if(abs(dum).gt.1.D-10) write(6,600) j,k,dum
 40      continue
      write(6,15)
 15   format(' b5(j,k) - b5(jm,k) - 2*b6(jm,k)*dpx ')
      do 41 j = 2,lspm
         jm = j - 1
         do 41 k = 1,lst
            km = k - 1
            if(k.eq.1) km = lst
            dum =   b5(j,k) - b5(jm,k) - 2*b6(jm,k)*dpx
            if(abs(dum).gt.1.D-10) write(6,600) j,k,dum
 41      continue
      write(6,16)
 16   format(' b8(j,k) - b8(jm,k) - 2*b9(jm,k)*dpx ')
      do 42 j = 2,lspm
         jm = j - 1
         do 42 k = 1,lst
            km = k - 1
            if(k.eq.1) km = lst
            dum =  b8(j,k) - b8(jm,k) - 2*b9(jm,k)*dpx
            if(abs(dum).gt.1.D-10) write(6,600) j,k,dum
 42      continue
      write(6,17)
 17   format(' dbdt down ')
      write(6,18)
 18   format(' b4(j,k) - b4(j,km) - 2*b7(j,km)*dtx ')
      do 50 j = 1,lspm
         jm = j - 1
         do 50 k = 1,lst
            km = k - 1
            if(k.eq.1) km = lst
            dum =   b4(j,k) - b4(j,km) - 2*b7(j,km)*dtx
            if(abs(dum).gt.1.D-10) write(6,600) j,k,dum
 50      continue
      write(6,19)
 19   format(' b5(j,k) - b5(j,km) - 2*b8(j,km)*dtx ')
      do 51 j = 1,lspm
         jm = j - 1
         do 51 k = 1,lst
            km = k - 1
            if(k.eq.1) km = lst
            dum =  b5(j,k) - b5(j,km) - 2*b8(j,km)*dtx
            if(abs(dum).gt.1.D-10) write(6,600) j,k,dum
 51      continue
      write(6,27)
 27   format(' b6(j,k) - b6(j,km) - 2*b9(j,km)*dtx ')
      do 52 j = 1,lspm
         jm = j - 1
         do 52 k = 1,lst
            km = k - 1
            if(k.eq.1) km = lst
            dum =   b6(j,k) - b6(j,km) - 2*b9(j,km)*dtx
            if(abs(dum).gt.1.D-10) write(6,600) j,k,dum
 52      continue
      return
      end
ccccccccccccc
cccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine lag3 ( xa,fa,n,h,x,f,fp,iop)
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER n,iop,i
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 fa,h,x,f,fp,xa,p,pm,pp,p11,pmh,pph
C============
      dimension xa(1), fa(1)
      i = ( x-xa(1) + .49999D0*h )/h + 1
      i = max(i,2)
      i = min(n-1,i)
      p = ( x-xa(i) ) / h
      pm = p - 1.0D0
      pp = p + 1.0D0
      p11 = pp*pm
      if ( iop .eq. 1 ) go to 20
      f = p*pm * fa(i-1)*.5D0 - p11 * fa(i) + p*pp * fa(i+1)*.5D0
   20 continue
      if ( iop .eq. 0 ) go to 30
      pmh = p - 0.5D0
      pph = p + 0.5D0
      fp = ( pmh*fa(i-1) - 2.0D0*p*fa(i) + pph*fa(i+1) ) / h
   30 continue
      return
      end
ccccccccccccccccccccccccccccccccccccc
      subroutine main0(rmaj,netcdf)
c.....calculates various quantities as a function of psi and theta
c     using output from equilibrium-mapping codes.   this uses
c     three point lagrange in theta and psi.  other
c     methods may be more suitable.
      IMPLICIT NONE
      include 'o.cln'
C============
      INTEGER netcdf
C idecl:  explicitize implicit REAL declarations:
      REAL*8 rmaj
C============
      call defolt(netcdf)
      if (netcdf .eq. 1) then
         call dskinCDF(rmaj)
      else
         call dskin(rmaj)
      endif
      call funcij
      return
      end
ccccccccccccccccccccccccccccccccccccc
      subroutine rstor
      IMPLICIT NONE
      include 'orbcom'
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER lspm,j,k,jz,jx,jxp,jzp
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 b0d,tau,xdum,zdum,dum,pdum,px,tx,xd,dum1,dum2,dxd,
     > dzd,dumx,dumz
C============
      dimension b0d(33,65)
      pi2 = 2*acos(-1.0D0)
cccc store ripple data on arrays in pol and theta
cccc the ripple is given by rdat(pol,thet)*sin(nrip*zet), fit with splines
      lspm = lsp - 1
      if(krip.eq.1) go to 10
      if(krip.eq.2) go to 20
      if(krip.eq.4) go to 40
      if(krip.eq.5) go to 50
      do 13 j = 1,lsp
      do 13 k = 1,lst
      r1(j,k) = -1.D8
 13   continue
      return
 10   continue
ccccc TFTR ripple (Scott 12/10/85 xrip=244cm, wrip=19cm, d0=4.3e-5, brip=1)
ccccc Redi:  xrip=223 cm, wrip=18.3 cm , d0 = 1.4e-5, brip = 1.1
      xrip = 223*xc/rmaj
      wrip = 18.3D0*xc/rmaj
      brip = 1.1D0
      nrip = 20
      d0 = 1.4D-5
      do 12 k = 1,lst
      do 12 j = 1,lsp
      tau = sqrt((x1(j,k)-xrip)**2 + brip*(z1(j,k))**2)
      r1(j,k) = d0*b1(j,k)*exp(tau/wrip)
      r1(j,k) = log(r1(j,k))
 12   continue
      write(6,5)
 5    format(' TFTR ripple entered')
      go to 100
 20   continue
ccccccc    ITER ripple
      nrip = 20
ccccc Redi/Miller  xrip=xrip(z) cm, wrip=53.5 cm , d0 = 3.75e-6, brip = .268
      wrip = 53.5D0*xc/rmaj
      brip = .268D0
      nrip = 20
      d0 = 3.75D-6
      do 32 k = 1,lst
         do 32 j = 1,lsp
            xrip = (xc/rmaj)*(675.D0-.00034D0*(z1(j,k)*rmaj/xc)**2)
      tau = sqrt((x1(j,k)-xrip)**2 + brip*(z1(j,k))**2)
      r1(j,k) = d0*b1(j,k)*exp(tau/wrip)
      r1(j,k) = log(r1(j,k))
 32      continue
      write(6,7)
 7    format('  ITER ripple entered')
      go to 100
 40   continue
ccccccc    NSTX ripple
      nrip = 16
      d0 = 7.07e-9
      do 42 k = 1,lst
         do 42 j = 1,lsp
cccc  xd is R in meters
ccc   zd is z in meters
            xd = x1(j,k)*rmaj*.01/xc
            dum = xd/.1136
            r1(j,k) = d0*b1(j,k)*exp(dum)
            r1(j,k) = log(r1(j,k))
 42      continue
      write(6,47)
 47    format('  NSTX ripple entered')
       go to 100
ccccc
 50   continue
ccccccc   Ignitor ripple
      open(51,file='ripple.dat',status='unknown')
ccc      x = 74 + (jx-1)*3.625     74 < x < 190  cm
ccc      z = -120 + (jz-1)*3.75   -120 < z < 120  cm
      read(51,*)
      read(51,*)
      read(51,*)
      read(51,*)
      do 51 jz = 1,65
      do 51 jx = 1,33
        read(51,*) dum1,dum2,b0d(jx,jz)
ccc        write(6,58) dum1,dum2,b0d(jx,jz)
 58     format(' x,z,dB ',1p3e12.4)
 51     continue
      nrip = 24
      write(6,189) rmaj,xc
 189  format(' rmaj xc ',1p2e12.4)
      do 52 k = 1,lst
         do 52 j = 1,lsp
      px = (j-1)*pw/lspm
       tx = (k-1)*pi2/lst
cccc  xd is R in cm
ccc   zd is z in cm
            xd = x1(j,k)*rmaj/xc
            zd = z1(j,k)*rmaj/xc
            jx = (xd - 74)/3.625D0 + 1
            jz = (zd + 120)/3.75D0 + 1
            jxp = jx + 1
             jzp = jz + 1
ccc             write(6,89) jx,jz,xd,zd
 89          format(' jx jz x z ',2i4,1p2e12.4)
            if(jxp.gt.33) stop
            if(jzp.gt.65) stop
            dxd = xd - (74 + (jx-1)*3.625D0)
            dzd = zd - (-120 + (jz-1)*3.75D0)
            dumx = (b0d(jxp,jz) - b0d(jx,jz))*dxd/3.625D0
            dumz = (b0d(jx,jzp) - b0d(jx,jz))*dzd/3.75D0
            dum = b0d(jx,jz) + dumx + dumz
            r1(j,k) = max(1.D-6,dum)
ccc            write(6,59) xd,zd,r1(j,k)
 59         format(' x z db ',1p3e12.4)
            r1(j,k) = log(r1(j,k))
 52      continue
      write(6,57)
 57    format('  Ignitor ripple entered')
ccccc
 100  continue
            call spln(0,r1,r2,r3,r4,r5,r6,r7,r8,r9)
      return
      end
ccccccccccccccccccccccccccc
      function xproj( polx,thetx )
      IMPLICIT NONE
      include 'o.cln'
      include 'orbcom'
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER idum,i0,i,il,j
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 polx,thetx,xproj,bj,bs,dbs,tdum,pt,ptp,ptm,pt11
C============
      dimension bj(nsf), bs(3), dbs(3)
      idum = thetx*.1591549431D0
      tdum = thetx - twopi*(idum-1)
      idum = tdum*.1591549431D0
      tdum = tdum - twopi*idum
      i0 = tdum/dth + 1.49999D0
      do 100 i = 1, 3
      il = i0 - 2 + i
      if ( il .eq. 0 ) il = mth
      do 80 j = 1, nosurf
   80 bj(j) = x(il,j)
      call lag3(psival,bj,nosurf,dpsi,polx,bs(i),dbs(i),0)
  100 continue
      pt = ( tdum-theval(i0) ) / dth
      ptp = pt + 1
      ptm = pt - 1
      pt11 = ptp*ptm
      xproj = 0.5D0*pt*ptm*bs(1) - pt11*bs(2) + 0.5D0*pt*ptp*bs(3)
      return
      end
cccccccccccccccccccccccccccccccccccc
      function bfield(polx,thx)
      IMPLICIT NONE
      include 'o.cln'
      include 'orbcom'
ccccc   model dependent functions
ccc   dependence on pol is fit using 3 point lagrange polynomials
ccc   dependence on thet using 5 point lagrange polynomials
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER k
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 polx,thx,bfield,dpi,dthi
C============
      integer jl(idm),il(idm),id(idm),i0(idm)
      real*8 th0(idm),ps0(idm),gd(idm,3),qd(idm,3),cr(idm,3)
      real*8 al(idm,5),bl(idm,5),bj(idm,5),pdm(idm,3)
      real*8 bs(idm,5),dbs(idm,5),ag(idm)
      real*8 q1(idm),q2(idm),q3(idm),q4(idm)

      dpi = 1.D0/dpsi
      dthi = 1.D0/dth
      k = 1
      id(k) = (polx - psival(1))*dpi + 1.499999D0
      id(k) = max(id(k),2)
      id(k) = min(nosurf-1,id(k))
      ps0(k) = psival(id(k))
      q1(k) = ( polx-ps0(k) ) *dpi
      q2(k) = q1(k)*q1(k)
      al(k,1) = .5D0*q1(k)*(q1(k)-1.D0)
      al(k,2) = 1 - q2(k)
      al(k,3) = .5D0*q1(k)*(q1(k) + 1.D0)
cccc non vectoriqable
      jl(k) = thx*.1591549431D0
      ag(k) = thx - twopi*(jl(k)-1)
      jl(k) = ag(k)*.1591549431D0
      ag(k) = ag(k) - twopi*jl(k)
      i0(k) = ag(k)*dthi + 1.49999D0
ccccc     i = 1
      il(k) = i0(k) - 2
      jl(k) = il(k)/(il(k) - 1.D-12)
      il(k) = mth*(1 - jl(k)) + il(k)
      bj(k,1) = bdat(il(k),id(k)-1)
      bj(k,2) = bdat(il(k),id(k))
      bj(k,3) = bdat(il(k),id(k)+1)
cccc    copy of lag3
      bs(k,1) = al(k,1)*bj(k,1) + al(k,2)*bj(k,2) + al(k,3)*bj(k,3)
cccc  end of lag3
ccccc     i = 2
      il(k) = i0(k) - 1
      jl(k) = il(k)/(il(k) - 1.D-12)
      il(k) = mth*(1 - jl(k)) + il(k)
17    continue
cccc non vectoriqable
      bj(k,1) = bdat(il(k),id(k)-1)
      bj(k,2) = bdat(il(k),id(k))
      bj(k,3) = bdat(il(k),id(k)+1)
cccc    copy of lag3
      bs(k,2) = al(k,1)*bj(k,1) + al(k,2)*bj(k,2) + al(k,3)*bj(k,3)
cccc  end of lag3
cccc   i = 3
      il(k) = i0(k)
      jl(k) = il(k)/(il(k) - 1.D-12)
      il(k) = mth*(1 - jl(k)) + il(k)
ccc   non vectorizable
      bj(k,1) = bdat(il(k),id(k)-1)
      bj(k,2) = bdat(il(k),id(k))
      bj(k,3) = bdat(il(k),id(k)+1)
cccc    copy of lag3
      bs(k,3) = al(k,1)*bj(k,1) + al(k,2)*bj(k,2) + al(k,3)*bj(k,3)
cccc  end of lag3
cccc   i = 4
      il(k) = i0(k) + 1
ccc   non vectorizable
      bj(k,1) = bdat(il(k),id(k)-1)
      bj(k,2) = bdat(il(k),id(k))
      bj(k,3) = bdat(il(k),id(k)+1)
cccc    copy of lag3
      bs(k,4) = al(k,1)*bj(k,1) + al(k,2)*bj(k,2) + al(k,3)*bj(k,3)
cccc  end of lag3
cccc   i = 5
      il(k) = i0(k) + 2
ccc   non vectorizable
      bj(k,1) = bdat(il(k),id(k)-1)
      bj(k,2) = bdat(il(k),id(k))
      bj(k,3) = bdat(il(k),id(k)+1)
      th0(k) = theval(i0(k))
cccc    copy of lag3
      bs(k,5) = al(k,1)*bj(k,1) + al(k,2)*bj(k,2) + al(k,3)*bj(k,3)
cccc  end of lag3
      q1(k) = ( ag(k) - th0(k) ) *dthi
      q2(k) = q1(k)*q1(k)
      q3(k) = q1(k)*q2(k)
      q4(k) = q1(k)*q3(k)
      al(k,1) = (q4(k)-2*q3(k)-q2(k)+2*q1(k))*.041666666666D0
      al(k,2) = (-q4(k)+q3(k)+4*q2(k)-4*q1(k))*.1666666666666666D0
      al(k,3) = (q4(k)-5*q2(k)+4.D0)*.25D0
      al(k,4) = (-q4(k)-q3(k)+4*q2(k)+4*q1(k))*.16666666666666666D0
      al(k,5) = (q4(k)+2*q3(k)-q2(k)-2*q1(k))*.041666666666666666D0
      bfield = al(k,1)*bs(k,1) + al(k,2)*bs(k,2) + al(k,3)*bs(k,3)
     1      + al(k,4)*bs(k,4) + al(k,5)*bs(k,5)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccc
      function zproj ( polx,thetx )
      IMPLICIT NONE
      include 'o.cln'
      include 'orbcom'
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER idum,i0,i,il,j
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 polx,thetx,zproj,bj,bs,dbs,tdum,pt,ptp,ptm,pt11,dum
C============
      dimension bj(nsf), bs(3), dbs(3)
      idum = thetx*.1591549431D0
      tdum = thetx - twopi*(idum-1)
      idum = tdum*.1591549431D0
      tdum = tdum - twopi*idum
      i0 = tdum/dth + 1.49999D0
      do 100 i = 1, 3
      il = i0 - 2 + i
      if ( il .eq. 0 ) il = mth
      do 80 j = 1, nosurf
   80 bj(j) = z(il,j)
      call lag3(psival,bj,nosurf,dpsi,polx,bs(i),dbs(i),0)
  100 continue
      pt = ( tdum-theval(i0) ) / dth
      ptp = pt + 1
      ptm = pt - 1
      pt11 = ptp*ptm
      dum = 0.5D0*pt*ptm*bs(1) - pt11*bs(2) + 0.5D0*pt*ptp*bs(3)
      if(krip.eq.3) dum = -dum  ! invert the ITER equilibria
      zproj = dum
      return
      end
ccccccccccccccccccccccccccccccccc
      function giac(px,tx)
      IMPLICIT NONE
      include 'o.cln'
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER ntor,nprt,id,jl,i0,il
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 px,tx,giac,q0,qw,pw,rq1,rmaj,ekev,bkg,zprt,prot,engn,
     > tran,xc,bs,al,bl,bj,pi,pi2,dpi,dthi,ps0,z1,z2,ag,th0,z3,z4
C============
      common /anal/ q0,qw,pw,rq1
      common /prm1/ rmaj,ekev,bkg,zprt,prot,engn,tran,xc,ntor,nprt
ccccc   model dependent functions
ccc   dependence on px is fit using 3 point lagrange pxynomials
ccc   dependence on thet using 5 point lagrange pxynomials
      dimension bs(5),al(5),bl(5),bj(5)
      pi = 4*atan(1.D0)
      pi2 = 2*pi
      dpi = 1.D0/dpsi
      dthi = 1.D0/dth
cccc non vectorizable operations
      id = (px - psival(1))*dpi + 1.499999D0
      id = max(id,2)
      id = min(nosurf-1,id)
      ps0 = psival(id)
      z1 = ( px-ps0 ) *dpi
      z2 = z1*z1
      al(1) = .5D0*z1*(z1-1.D0)
      al(2) = 1 - z2
      al(3) = .5D0*z1*(z1 + 1.D0)
cccc    b of px,thet and db/dpx,  db/dth
ccc       del of px, thet
      jl = tx*.1591549431D0
      ag = tx - twopi*(jl-1)
      jl = ag*.1591549431D0
      ag = ag - twopi*jl
      i0 = ag*dthi + 1.49999D0
ccccc     i = 1
      il = i0 - 2
      jl = il/(il - 1.D-12)
      il = mth*(1 - jl) + il
cccc non vectorizable
      bj(1) = xjacob(il,id-1)
      bj(2) = xjacob(il,id)
      bj(3) = xjacob(il,id+1)
cccc    copy of lag3
      bs(1) = al(1)*bj(1) + al(2)*bj(2) + al(3)*bj(3)
cccc  end of lag3
ccccc     i = 2
      il = i0 - 1
      jl = il/(il - 1.D-12)
      il = mth*(1 - jl) + il
      bj(1) = xjacob(il,id-1)
      bj(2) = xjacob(il,id)
      bj(3) = xjacob(il,id+1)
cccc    copy of lag3
      bs(2) = al(1)*bj(1) + al(2)*bj(2) + al(3)*bj(3)
cccc  end of lag3
cccc   i = 3
      il = i0
      jl = il/(il - 1.D-12)
      il = mth*(1 - jl) + il
      bj(1) = xjacob(il,id-1)
      bj(2) = xjacob(il,id)
      bj(3) = xjacob(il,id+1)
cccc    copy of lag3
      bs(3) = al(1)*bj(1) + al(2)*bj(2) + al(3)*bj(3)
cccc  end of lag3
cccc   i = 4
      il = i0 + 1
      bj(1) = xjacob(il,id-1)
      bj(2) = xjacob(il,id)
      bj(3) = xjacob(il,id+1)
cccc    copy of lag3
      bs(4) = al(1)*bj(1) + al(2)*bj(2) + al(3)*bj(3)
cccc  end of lag3
cccc   i = 5
      il = i0 + 2
      bj(1) = xjacob(il,id-1)
      bj(2) = xjacob(il,id)
      bj(3) = xjacob(il,id+1)
      th0 = theval(i0)
cccc    copy of lag3
      bs(5) = al(1)*bj(1) + al(2)*bj(2) + al(3)*bj(3)
cccc  end of lag3
      z1 = ( ag - th0 ) *dthi
      z2 = z1*z1
      z3 = z1*z2
      z4 = z1*z3
      al(1) = (z4-2*z3-z2+2*z1)*.04166666666666666D0
      al(2) = (-z4+z3+4*z2-4*z1)*.1666666666666666D0
      al(3) = (z4-5*z2+4.D0)*.25D0
      al(4) = (-z4-z3+4*z2+4*z1)*.16666666666666666D0
      al(5) = (z4+2*z3-z2-2*z1)*.041666666666666666D0
      giac = al(1)*bs(1) + al(2)*bs(2) + al(3)*bs(3)
     1      + al(4)*bs(4) + al(5)*bs(5)
      return
      end
cccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccc
      function qfun(px)
      IMPLICIT NONE
      include 'o.cln'
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER ntor,nprt,id
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 px,qfun,q0,qw,pw,rq1,rmaj,ekev,bkg,zprt,prot,engn,tran,
     > xc,al,bj,dpi,ps0,z1,z2
C============
      common /anal/ q0,qw,pw,rq1
      common /prm1/ rmaj,ekev,bkg,zprt,prot,engn,tran,xc,ntor,nprt
ccccc   model dependent functions
ccc   dependence on px is fit using 3 point lagrange pxynomials
      dimension al(5),bj(5)
      dpi = 1.D0/dpsi
cccc non vectorizable operations
      id = (px - psival(1))*dpi + 1.499999D0
      id = max(id,2)
      id = min(nosurf-1,id)
      ps0 = psival(id)
      z1 = ( px-ps0 ) *dpi
      z2 = z1*z1
      al(1) = .5D0*z1*(z1-1.D0)
      al(2) = 1 - z2
      al(3) = .5D0*z1*(z1 + 1.D0)
cccc    q of px
cccc non vectorizable
      bj(1) = qdat(id-1)
      bj(2) = qdat(id)
      bj(3) = qdat(id+1)
cccc    copy of lag3
      qfun = al(1)*bj(1) + al(2)*bj(2) + al(3)*bj(3)
      return
      end
cccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccc
      function psfun(px)
      IMPLICIT NONE
      include 'o.cln'
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER ntor,nprt,k
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 px,psfun,q0,qw,pw,rq1,rmaj,ekev,bkg,zprt,prot,engn,
     > tran,xc,dum,pdum,dpol,qfun
C============
      common /anal/ q0,qw,pw,rq1
      common /prm1/ rmaj,ekev,bkg,zprt,prot,engn,tran,xc,ntor,nprt
ccccc   find toroidal flux psi of pol
      dum = 0
      pdum = 0
      dpol = .001D0*px
      do 10 k = 1,1000
        pdum = pdum + dpol
        dum = dum + dpol*qfun(pdum)
 10     continue
        psfun = dum
      return
      end
cccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccc
      function gfun(px)
      IMPLICIT NONE
      include 'o.cln'
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER ntor,nprt,id
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 px,gfun,q0,qw,pw,rq1,rmaj,ekev,bkg,zprt,prot,engn,tran,
     > xc,al,bj,dpi,ps0,z1,z2
C============
      common /anal/ q0,qw,pw,rq1
      common /prm1/ rmaj,ekev,bkg,zprt,prot,engn,tran,xc,ntor,nprt
ccccc   model dependent functions
ccc   dependence on px is fit using 3 point lagrange pxynomials
      dimension al(5),bj(5)
      dpi = 1.D0/dpsi
cccc non vectorizable operations
      id = (px - psival(1))*dpi + 1.499999D0
      id = max(id,2)
      id = min(nosurf-1,id)
      ps0 = psival(id)
      z1 = ( px-ps0 ) *dpi
      z2 = z1*z1
      al(1) = .5D0*z1*(z1-1.D0)
      al(2) = 1 - z2
      al(3) = .5D0*z1*(z1 + 1.D0)
cccc    q of px
cccc non vectorizable
      bj(1) = gdat(id-1)
      bj(2) = gdat(id)
      bj(3) = gdat(id+1)
cccc    copy of lag3
      gfun = al(1)*bj(1) + al(2)*bj(2) + al(3)*bj(3)
      return
      end
cccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccc
      function rifun(px)
      IMPLICIT NONE
      include 'o.cln'
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER ntor,nprt,id
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 px,rifun,q0,qw,pw,rq1,rmaj,ekev,bkg,zprt,prot,engn,
     > tran,xc,al,bj,dpi,ps0,z1,z2
C============
      common /anal/ q0,qw,pw,rq1
      common /prm1/ rmaj,ekev,bkg,zprt,prot,engn,tran,xc,ntor,nprt
ccccc   model dependent functions
ccc   dependence on px is fit using 3 point lagrange pxynomials
      dimension al(5),bj(5)
      dpi = 1.D0/dpsi
cccc non vectorizable operations
      id = (px - psival(1))*dpi + 1.499999D0
      id = max(id,2)
      id = min(nosurf-1,id)
      ps0 = psival(id)
      z1 = ( px-ps0 ) *dpi
      z2 = z1*z1
      al(1) = .5D0*z1*(z1-1.D0)
      al(2) = 1 - z2
      al(3) = .5D0*z1*(z1 + 1.D0)
cccc    q of px
cccc non vectorizable
      bj(1) = cur(id-1)
      bj(2) = cur(id)
      bj(3) = cur(id+1)
cccc    copy of lag3
      rifun = al(1)*bj(1) + al(2)*bj(2) + al(3)*bj(3)
      return
      end
cccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccc
      function pfun(px)
      IMPLICIT NONE
      include 'o.cln'
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER ntor,nprt,id
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 px,pfun,q0,qw,pw,rq1,rmaj,ekev,bkg,zprt,prot,engn,tran,
     > xc,al,bj,dpi,ps0,z1,z2
C============
      common /anal/ q0,qw,pw,rq1
      common /prm1/ rmaj,ekev,bkg,zprt,prot,engn,tran,xc,ntor,nprt
ccccc   model dependent functions
ccc   dependence on px is fit using 3 point lagrange pxynomials
      dimension al(5),bj(5)
      dpi = 1.D0/dpsi
cccc non vectorizable operations
      id = (px - psival(1))*dpi + 1.499999D0
      id = max(id,2)
      id = min(nosurf-1,id)
      ps0 = psival(id)
      z1 = ( px-ps0 ) *dpi
      z2 = z1*z1
      al(1) = .5D0*z1*(z1-1.D0)
      al(2) = 1 - z2
      al(3) = .5D0*z1*(z1 + 1.D0)
cccc    q of px
cccc non vectorizable
      bj(1) = pdat(id-1)
      bj(2) = pdat(id)
      bj(3) = pdat(id+1)
cccc    copy of lag3
      pfun = al(1)*bj(1) + al(2)*bj(2) + al(3)*bj(3)
      return
      end
cccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C        SUBROUTINE GELG                                                GELG  40
C                                                                       GELG  50
C        PURPOSE                                                        GELG  60
C           TO SOLVE A GENERAL SYSTEM OF SIMULTANEOUS LINEAR EQUATIONS. GELG  70
C                                                                       GELG  80
C        USAGE                                                          GELG  90
C           CALL GELG(R,A,M,N,EPS,IER)                                  GELG 100
C                                                                       GELG 110
C        DESCRIPTION OF PARAMETERS                                      GELG 120
C           R      - THE M BY N MATRIX OF RIGHT HAND SIDES.  (DESTROYED)GELG 130
C                    ON RETURN R CONTAINS THE SOLUTION OF THE EQUATIONS.GELG 140
C           A      - THE M BY M COEFFICIENT MATRIX.  (DESTROYED)        GELG 150
C           M      - THE NUMBER OF EQUATIONS IN THE SYSTEM.             GELG 160
C           N      - THE NUMBER OF RIGHT HAND SIDE VECTORS.             GELG 170
C           EPS    - AN INPUT CONSTANT WHICH IS USED AS RELATIVE        GELG 180
C                    TOLERANCE FOR TEST ON LOSS OF SIGNIFICANCE.        GELG 190
C           IER    - RESULTING ERROR PARAMETER CODED AS FOLLOWS         GELG 200
C                    IER=0  - NO ERROR,                                 GELG 210
C                    IER=-1 - NO RESULT BECAUSE OF M LESS THAN 1 OR     GELG 220
C                             PIVOT ELEMENT AT ANY ELIMINATION STEP     GELG 230
C                             EQUAL TO 0,                               GELG 240
C                    IER=K  - WARNING DUE TO POSSIBLE LOSS OF SIGNIFI-  GELG 250
C                             CANCE INDICATED AT ELIMINATION STEP K+1,  GELG 260
C                             WHERE PIVOT ELEMENT WAS LESS THAN OR      GELG 270
C                             EQUAL TO THE INTERNAL TOLERANCE EPS TIMES GELG 280
C                             ABSOLUTELY GREATEST ELEMENT OF MATRIX A.  GELG 290
C                                                                       GELG 300
C        REMARKS                                                        GELG 310
C           INPUT MATRICES R AND A ARE ASSUMED TO BE STORED COLUMNWISE  GELG 320
C           IN M*N RESP. M*M SUCCESSIVE STORAGE LOCATIONS. ON RETURN    GELG 330
C           SOLUTION MATRIX R IS STORED COLUMNWISE TOO.                 GELG 340
C           THE PROCEDURE GIVES RESULTS IF THE NUMBER OF EQUATIONS M IS GELG 350
C           GREATER THAN 0 AND PIVOT ELEMENTS AT ALL ELIMINATION STEPS  GELG 360
C           ARE DIFFERENT FROM 0. HOWEVER WARNING IER=K - IF GIVEN -    GELG 370
C           INDICATES POSSIBLE LOSS OF SIGNIFICANCE. IN CASE OF A WELL  GELG 380
C           SCALED MATRIX A AND APPROPRIATE TOLERANCE EPS, IER=K MAY BE GELG 390
C           INTERPRETED THAT MATRIX A HAS THE RANK K. NO WARNING IS     GELG 400
C           GIVEN IN CASE M=1.                                          GELG 410
C                                                                       GELG 420
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  GELG 430
C           NONE                                                        GELG 440
C                                                                       GELG 450
C        METHOD                                                         GELG 460
C           SOLUTION IS DONE BY MEANS OF GAUSS-ELIMINATION WITH         GELG 470
C           COMPLETE PIVOTING.                                          GELG 480
C                                                                       GELG 490
C     ..................................................................GELG 500
C                                                                       GELG 510
      SUBROUTINE GELG(R,A,M,N,EPS,IER)
C                                                                       GELG 530
C                                                                       GELG 540
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER m,n,ier,mm,nm,l,i,lst,k,j,ll,lend,ii,ist
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 a,eps,r,piv,tb,tol,pivi
C============
      DIMENSION A(1),R(1)
      IF(M)23,23,1
C                                                                       GELG 570
C     SEARCH FOR GREATEST ELEMENT IN MATRIX A                           GELG 580
    1 IER=0
      PIV=0.D0
      MM=M*M
      NM=N*M
      DO 3 L=1,MM
      TB=ABS(A(L))
      IF(TB-PIV)3,3,2
    2 PIV=TB
      I=L
    3 CONTINUE
      TOL=EPS*PIV
C     A(I) IS PIVOT ELEMENT. PIV CONTAINS THE ABSOLUTE VALUE OF A(I).   GELG 700
C                                                                       GELG 710
C                                                                       GELG 720
C     START ELIMINATION LOOP                                            GELG 730
      LST=1
      DO 17 K=1,M
C                                                                       GELG 760
C     TEST ON SINGULARITY                                               GELG 770
      IF(PIV)23,23,4
    4 IF(IER)7,5,7
    5 IF(PIV-TOL)6,6,7
    6 IER=K-1
    7 PIVI=1.D0/A(I)
      J=(I-1)/M
      I=I-J*M-K
      J=J+1-K
C     I+K IS ROW-INDEX, J+K COLUMN-INDEX OF PIVOT ELEMENT               GELG 860
C                                                                       GELG 870
C     PIVOT ROW REDUCTION AND ROW INTERCHANGE IN RIGHT HAND SIDE R      GELG 880
      DO 8 L=K,NM,M
      LL=L+I
      TB=PIVI*R(LL)
      R(LL)=R(L)
    8 R(L)=TB
C                                                                       GELG 940
C     IS ELIMINATION TERMINATED                                         GELG 950
      IF(K-M)9,18,18
C                                                                       GELG 970
C     COLUMN INTERCHANGE IN MATRIX A                                    GELG 980
    9 LEND=LST+M-K
      IF(J)12,12,10
   10 II=J*M
      DO 11 L=LST,LEND
      TB=A(L)
      LL=L+II
      A(L)=A(LL)
   11 A(LL)=TB
C                                                                       GELG1070
C     ROW INTERCHANGE AND PIVOT ROW REDUCTION IN MATRIX A               GELG1080
   12 DO 13 L=LST,MM,M
      LL=L+I
      TB=PIVI*A(LL)
      A(LL)=A(L)
   13 A(L)=TB
C                                                                       GELG1140
C     SAVE COLUMN INTERCHANGE INFORMATION                               GELG1150
      A(LST)=J
C                                                                       GELG1170
C     ELEMENT REDUCTION AND NEXT PIVOT SEARCH                           GELG1180
      PIV=0.D0
      LST=LST+1
      J=0
      DO 16 II=LST,LEND
      PIVI=-A(II)
      IST=II+M
      J=J+1
      DO 15 L=IST,MM,M
      LL=L-J
      A(L)=A(L)+PIVI*A(LL)
      TB=ABS(A(L))
      IF(TB-PIV)15,15,14
   14 PIV=TB
      I=L
   15 CONTINUE
      DO 16 L=K,NM,M
      LL=L+J
   16 R(LL)=R(LL)+PIVI*R(L)
   17 LST=LST+M
C     END OF ELIMINATION LOOP                                           GELG1380
C                                                                       GELG1390
C                                                                       GELG1400
C     BACK SUBSTITUTION AND BACK INTERCHANGE                            GELG1410
   18 IF(M-1)23,22,19
   19 IST=MM+M
      LST=M+1
      DO 21 I=2,M
      II=LST-I
      IST=IST-LST
      L=IST-M
      L=A(L)+.5D0
      DO 21 J=II,NM,M
      TB=R(J)
      LL=J
      DO 20 K=IST,MM,M
      LL=LL+1
   20 TB=TB-A(K)*R(LL)
      K=J+L
      R(J)=R(K)
   21 R(K)=TB
   22 RETURN
C                                                                       GELG1600
C                                                                       GELG1610
C     ERROR RETURN                                                      GELG1620
   23 IER=-1
      RETURN
      END
 
 
 
! 14Mar2000 fgtok -s r8_precision.sub "r8con.csh conversion"
! 14Mar2000 fgtok
