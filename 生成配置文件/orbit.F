      program main
      IMPLICIT NONE
      include 'orbcom'
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER nspin,kripsav,ndum,name0,nterm,k,nend,kp,nlim
C============
C idecl:  explicitize implicit REAL declarations:

      REAL*8 dele2,dum2,qdum,qfun,pzdum,rmdum,tdum,dum,dum3
      REAL*8 starttime,endtime
C============
      integer  iseed,iws_id,nran

      common /randomg/ iseed,nran

      call cpu_time(starttime)
 
      myfile = 6

      nran = 0
cl3   master process opens other files 
      open(6,file='orbit.out',status='unknown')
      open(20,file='spdata',status='unknown')
      open(55,file='evst',status='unknown')
      open(92,file='diag.out',status='unknown')    
      open(34,file='lost.plt',status='unknown')
      open(97,file='poincare.plt',status='unknown')
      open(128,file='modstp.plt',status='unknown')
cc   r. b. white   princeton, jan 1982
cc   files needed:common blocks orbcom, equilibrium spdata written by eqs.f
cc   nplot =1 gives single particle orbit data,
cc       data recorded at intervals dt1, stored in trajectory
cc   nplot=2 many particles, for diffusion and loss study,stored in diffusion
cc   nplot=3  Poincare, data=poincare
cc   nplot=4  mean particle distribution, see call rcrd4, data=distave
cc       dist(pol) for modification of distribution due to modes
ccc      also construct    fdist(E,P,mu)  time averaged
cc   nplot = 5 Local Bootstrap current
cc   nplot = 6  Tent diffusion determination
cc   nplot = 7  flight distribution determination
cc   nplot = 8  mode excitation
ccc   nplot=9  kinetic Poincare
cc   nplot = 11  stoch1d and FLI
cc   nplot = 12  stochmu- Adjust #steps in stochmu, limits in rcrdrot, # modes
cc   nplot = 13  stoche- Adjust #steps in stoche, limits in rcrdrot, # modes
cc   nplot = 14  anneal
cc   nplot = 15  improve statistics
cc   nplot = 16  frequency
cc   nplot = 17  E, Pz to omega map
cc
cc--CODE PARAMETERS
cc   col = the collision frequency,
ccc- ncol=1, routine scat gives only energy**(-3/2) dependence
ccc- ncol=2, routine scatr gives T, n profiles and error function dependence
cc   ekev is particle energy in kev,
cc   bkg is b at the magnetic axis in kgauss
ccc  field is b(pol,thet) + rpl(pol,thet)sin(N*zeta), rpl not usual ripple!
cc   zprt is the charge, and prot the mass in proton units
cc   rmaj is the major axis in cm- given in equilibrium file
cc   dele is the allowed fractional energy change per step
cc   the time step dt is adjusted accordingly,
cc   to run at dt0 set dele > 1.
cc   nprt is the number of particles.
cc   trun is total run time. need ten steps per transit time, tran
cc    nploteq=1 plots some equilibrium functions, 2 more, 0 none
ccc      subroutine gc_xz2pt converts from x,z to pol,theta
ccc--computing time(Cray A): alpha particles in 272cm 50kg TFTR equilibrium
ccc --A run with 1000 particles 100 transits
ccc--with-ripple takes 38 sec
ccc--time-is 53% field, 31% onestep
cccccccc-  Perturbations - nptrbav = 1, calc potential from alpha
cccccccc-  Perturbations - nptrbav = 2, calc alpha and potential from xi
cccc
cccc  CHOICE OF RUN, PLOTS
      nplot = 2
      nploteq = 1
cc
ccccc NUMBER OF PARTICLES-  
        nprt = 1000
	if(nplot.eq.3) nprt = 30  ! poincare 
	if(nplot.eq.9) nprt = 40    ! kinetic poincare
	if(nplot.eq.11) nprt = 4000    ! stoch and FLI
	if(nplot.eq.14) nprt = 999900    ! anneal
	if(nplot.eq.16) nprt = 500    ! freq
	if(nplot.eq.17) nprt = 1000    ! E, Pz - freq map
      if(nplot.eq.7) open(8,file='flight.plt',status='unknown')
	if(nplot.eq.7) write(8,88)
88	format('#  flights   dpol,  dzet')
c -- Seed to initialize the random number generator. If iseed < 0, then
c    the rng is initialized with the date and time.
      iseed = 0
ccc
ccccc SIMULATION RUN TIME, in on-axis toroidal transits
    ntor = 2500
	if(nplot.eq.9) ntor = 200   !   kinetic poincare
	if(nplot.eq.8) ntor = 100   ! mode excitation
	if(nplot.eq.11) ntor = 50   ! rotation
	if(nplot.eq.15) ntor = 40   ! statistics
	if(nplot.eq.16) ntor = 1000   !  freq
	if(nplot.eq.17) ntor = 200   !  freq
	tim4 = 2.e-3   !   rcrd4 data dump
ccccccccc
ccc   READ EQUILIBRIUM DATA IN SPDATA
      call readvec        ! reads numerical equilibrium  
    
ccc      call drev    ! reverses drift direction
cccc  FIELD STRENGTH
      bkg = 23
cccc
cccc  numerical equilibrium spdata read
cccc  can set ripple zero with krip=0
      krip = 0  !  comment out  for ripple
      kripsav = krip
ccc
cccc  PERTURBATION CHOICE, 0=none, 1=delta B, 2=potential, 3=both
ccc    4= read delta B perturbation
        npert = 1
	if(nplot.eq.10) npert = 1
	if(nplot.eq.12) npert = 1
	if(nplot.eq.8) npert = 1
        if(nplot.eq.3) npert = 1
	if(nplot.eq.17) npert = 0
      pamp =1000    ! potential amplitude, Kev
ccccc
ccc   PARTICLE DISTRIBUTION, see shelldep, poindep
        polo = 0.5*pw
    p1 = 0.01*pw
    p2 = 0.99*pw
        pchi = 0.7

cccccc
ccccc  PARTICLE PARAMETERS
      zprt = 1.D0
      prot = 1.D0
      ekev = 93
cccccccccc
cccc-  DISTRIBUTION PARAMETERS -for stochastic, either stochmu or stoche
	mubk = 10    !   mu*B0 in keV
	mubmin = 0
	mubmax = 110
	emink = 40
	emaxk = 60
	dek = 1
	dmubk = 1
	nstoche = 30
	nstochmu = 50
	nstochp = 30
	nco = 1          !  co 1, counter -1, both 0  particles
        if(nplot.eq.3.or.nplot.eq.10) then
        emink = ekev
        mubk = 0.
        emaxk = ekev
        endif
cccccccccccccccccccccccccccccccccccc
      if(nplot.eq.3) ekev = 40

cccccc
ccccc   STEP CONTROL, USING ENERGY CONSERVATION 
        dele = 5.D-8
	dele = 5
	if(nplot.eq.2) dele = 5.
        dele2 = .1D0*dele
	if(nplot.eq.4) dele = 5.
	if(nplot.eq.5) dele = 5.
	if(nplot.eq.6) dele = 5.
	if(nplot.eq.8) dele = 5.
	if(nplot.eq.9) dele = 5.
ccc

         call set1

cccccc   BEAM BIRTH PROFILE
cccccc      open(89,file="142111M18_birth_4_nogt_rho_1.sv",status='unknown')
cccccc      read(89,*) ndum,nlim
cccccc	ndum = ndum + 3000    !  to move through the beam data
cccccc      do 5 j = 1,ndum
cccccc        read(89,*)
cccccc 5      continue
cc
ccc   PERTURBATION DATA
	if(npert.eq.1) call ampa  !  analytic alpha
cc	if(npert.eq.1) call ampx  ! analytic  xi
cc	if(npert.eq.4) call readptrba  !  read alpha
	if(npert.eq.4) call readptrbx  !  read xi
	n4start = tim4*omeg0/dt0
ccccc
cccccc  COLLISIONS , drag -units of on-axis gyro freq
cc   routines in collisions.f
      ncol = 0  ! ncol=1, energy dep, ncol=2 full profiles
      if(nplot.eq.1) ncol = 0
      if(nplot.eq.3) ncol = 0
      if(nplot.eq.9) ncol = 0
      if(nplot.eq.10) ncol = 0
      if(nplot.eq.11) ncol = 0
      if(nplot.eq.12) ncol = 0
      if(nplot.eq.14) ncol = 0
      if(nplot.eq.16) ncol = 0
      if(nplot.eq.17) ncol = 0
ccc-for ncol=2 see scatr and functions denb, deni, tempi,tempe
      massb = 2     ! background plasma mass
      chgb = 1     ! background plasma charge
      imp = 0    ! imp=1 impurity species, 0=none
      massi = 5  ! impurity mass
      chgi = 5   ! impurity charge
      eion = 1.1     ! ion energy in kev for plots
cccc- assume totally ionized
	col = 1./(.010*omeg0)   !  10 msec
	drag = 0./(.051*omeg0)  !  51 msec
ccc
cccc  LOSS CONDITION
      nwal = 0 ! nwal=0 gives last flux surface, nwal=1 gives TFTR wall
       call walset
ccc
cccc--MODIFY RIPPLE MJAGNITUDE, COIL NUMBER
ccc      call setrip(20,2.)
ccccc
cccc-   PLOT ROUTINES
cc  routines in orbplot.f
ccc these routines insert particles, so cannot be used during a run
cccc harmon plots field perturbations
      if(npert.ne.0) call harmon
      if(ncol.eq.2) call scatrplot(eion,imp)
ccccc  plotf plots equilibrium functions
      if(nploteq.gt.0) call plotf
	call mupplane
	call peplane
ccccc
cccccc-  POSTPROCESSING ROUTINES
ccc     call detector   !   Find detector space in P, mu
cc       call fixworm    !  postprocessor for nplot=12, worm
cc        call wormplot   !  postprocessor for nplot=9, poinkdep
cc	call fliworm    !  postprocessor for nplot=11, FLI and worm
cc       go to 501
ccc

cccccc    PARTICLE INSERTION
cc   routines in deposit.f
ccc shelldep gives a general uniform distribution.
ccc sampledep reads from TRANSP data
ccc    poincare,  alphas,  fulldep,  read distf, btail
	tdum = -2.6983
    call shelldep
cc        if(nplot.eq.1) call shelldep(tdum)
cc        if(nplot.eq.2) call shelldep(tdum)
cc        if(nplot.eq.2) call sampledep   !  read from TRANSP
ccc	if(nplot.eq.2) call fulldepe   !  fixed E
ccc	if(nplot.eq.2) call fulldepmu   !  fixed mu
ccc	if(nplot.eq.2) call distread  !  Continuation using dumped dist
ccc	if(nplot.eq.2) call btail     !  bump on tail
	if(nplot.eq.3) call poindep   !  Poincare
ccc	if(nplot.eq.2) call alphdep   !  Alphas
	if(nplot.eq.5) call bootdep   ! bootstrap
	if(nplot.eq.6) call tent      !  local diffusion
        if(nplot.eq.7) call shelldep(tdum)
	if(nplot.eq.8) call fulldep2   !  all confined
	if(nplot.eq.9) call poinkdep  !  
	if(nplot.eq.10) call poinkdep
	if(nplot.eq.11) call poinkrot
        if(nplot.eq.16) call freqdep
        if(nplot.eq.17) call epfdep
        call trapfrac    ! calculation of trapped fraction
ccccccccccccccccccccccccccccc
cccc  end input
ccc-SAVE INITIAL PARTICLE DISTRIBUTION 
        call rcrd0
ccc-  LABEL PARTICLES ACCORDING TO TYPE
	call domain 
cccc    pdist(0) plots full particle distribution
       if(nplot.ne.1) call pdist(0)
CCC
ccc AUXILLIARY ROUTINES 
cc
ccc	call resonance  !  resonance locations
ccc      call modb2    !   mod b contours
ccc      call modrip    !  ripple contours
cc	call dfplot   !  changes in P,mu,E  see plot9
cc	call pediff   !   diffusion in E,P
cc	call fcheck
	if(nplot.eq.12) then
	call stochmu
	go to 501
	endif
	if(nplot.eq.13) then
	call stoche
	go to 501
	endif
	if(nplot.eq.14) then
        call sampledep
	call trapfrac
	call anneal2d
	go to 501
	endif
ccc         go to 501    ! to exit after auxilliary use
ccccc
ccc
cccccc-first run nplot=2, find lost particle, then use runone to plot
  
      name0 = 0
         if(name0.eq.0) go to 66
      nstart =  00              ! start recording at nstep=nstart 
         call runone(name0)     ! runs only particle name
 66   continue
  
ccccc
      qdum = qfun(pol(1)) 
      pzdum = (g(1)*rho(1) - pol(1))/pw
      rmdum =   rmu(1)*bax/en(1)
      if(nplot.eq.1) write(myfile,102) pol(1),
     & thet(1),ptch(1),pzdum,rmdum
 102  format(' deposit-pol,thet,ptch,pz,rmu',1p5e10.2)

ccc      go to 501
cccccccccccccccccccccc
        nstep = 0
        npts = 0
        kplt = 0
        nterm = 300*ntor
      epart = 0
      do 106 k = 1,nprt
      epart = epart + en(k)
 106  continue
ccccccccccccccccccccccccccccccc
ccc      if(nplot.eq.2) write(myfile,101)
 101  format(' -loss#  pol0,thet,zet,x,z,ptch,t')
ccccccccccccccccccccccccccccccccccccccc
	if(nplot.eq.4) call prep4
CCC
ccc   TIME STEP
200   continue
      ndum = mod(nstep,1000)
	dum = nstep
	dum2 = 1000*time(1)/omeg0
      if(ndum.eq.0) write(0,888) dum,dum2
 888  format('begin step,  t msec',1p6e12.4)
      nstep = nstep + 1
	if(nplot.eq.2) call rcrd2(nend)
      if(nplot.eq.3) call rcrd3
      if(nplot.eq.4) call rcrd4
      if(nplot.eq.5) call rcrd5
      if(nplot.eq.6) call rcrd6
      if(nplot.eq.7) call rcrd7
      if(nplot.eq.8) call rcrd8
      if(nplot.eq.9) call rcrd9
      if(nplot.eq.10) call rcrd10
      if(nplot.eq.11) call rcrdrot
      if(nplot.eq.12) call rcrdrot
      if(nplot.eq.13) call rcrdrot
      if(nplot.eq.15) call statbin
      if(nplot.eq.17) call rcrdfreq
ccc	call lowe    !   replace low energy and lost
 77   continue
cccccccccccccccccccccccccccccccccccc
      call onestep
cccc	if(nplot.eq.2) call modestep  !  advance mode frequencies
cc  pol,thet,zet,rho  variables stepped
cc    update all subsidiary variables for next step
      call update(nend)
      if(nplot.eq.2) call reduce(1)
ccccccccccccccccccccccccccccccccccccccc
      if(kplt.gt.10000) go to 999
      if(amp(1).gt.1.) go to 999   ! bail out
      if(nstep.gt.nterm) go to 999   ! bail out
      if(nend.eq.1) go to 999
      if(dele.lt.1.D0)  call stepfix
      if(ncol.eq.1) call scat
      if(ncol.eq.1) call slow
      if(ncol.eq.2) call scatr(imp)
      if(nplot.eq.1) call rcrd1
      go to 200
ccc    END OF RUN
 999  continue
ccc--OUTPUT DATA AT END OF RUN
      krip = kripsav
      if(nplot.ne.1) call pdist(1)   !  confined particle distribution
      if(nplot.ne.1) call pdist(2)    !  lost particle distribution
      if(nplot.eq.2)  call plot2
      if(nplot.eq.5)  call plot5
      if(nplot.eq.6)  call plot6
      if(nplot.eq.8)  call plot8
      if(nplot.eq.9)  call statout
      if(nplot.eq.16)  call freqwrite
      if(nplot.eq.17)  call epfwrite
ccc

501   continue
      
      call wrt6
      call phys6
      call cpu_time(endtime)
	dum = endtime - starttime
      write(myfile,*) 'COMPUTE TIME = ',dum,' sec'
      stop
      end

