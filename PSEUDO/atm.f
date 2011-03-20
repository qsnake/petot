      program atm
c               ______________
c              /              \  
c             /   Main - ATM   \
c             \     START      /           
c              \______________/
c                      |
c                ______V_____
c               |            |          __________          _______      ______
c               |   Zesec    |         |          |        |       |    /      \
c               |____________|    |--->|  Charge  |-error->|  Ext  |-->|  STOP  |
c                      |          |    |__________|        |_______|    \______/
c                ______V______    | 
c               /             /<--|   __________
c              /             /       |          |
c             /    Input    /<------>|  Zedate  |
c            /             /         |__________|
c           /_____________/------|         _______       _____
c                    |           |        |       |     /      \
c   |---->---------->|           |-error->|  Ext  |--->|  STOP  |
c   |                |                    |_______|     \______/
c   |               _V_                                        ____________
c   |             /     \        _________      _______       /            \
c   |            / More  \      |         |    |       |     /  Main - ATM  \
c   |           <         >-no->|  Zesec  |--->|  Ext  |--->|                |
c   A            \ Data? /      |_________|    |_______|     \     STOP     /
c   |             \ ___ /                                     \____________/
c   |                | yes
c   |               _V_
c   |             /     \
c   A            /Config-\
c   |           < uration >yes->|
c   |            \ Test? /      |
c   |             \ ___ /       |
c   |                | no       |
c   |           _____V_____     |
c   A          |           |    |
c   |          |  Vionic   |    |        __________
c   |          |___________|    |       |          |
c   |                |          V  |--->|  Splift  |
c   |                |<---------|  |    |__________|
c   |           _____V_____        |
c   |          |           |<------|   _______       ______
c   A          |           |          |       |     /      \
c   |          |  Velect   |--error-->|  Ext  |--->|  STOP  |
c   |          |           |          |_______|     \______/
c   |          |           |            __________            _____________
c   |          |           |           |          |          |             |
c   |          |           |<--------->|  atomxc  |<-------->|  XC-Package |
c   |          |           |           |__________|          |_____________|
c   |          |           |         
c   |          |___________|<------|     _________
c   |                |             |    |         |
c   |                |             |--->|  Spliq  |
c   A   |--->------->|                  |_________|    
c   |   A            |                                   __________
c   |   |           _V_                                 |          |
c   |   |         / 1'st\          __________     |---->|  Tridib  |
c   |   |        /or 2'nd\        |          |<---|     |__________|
c   |   |       <  Itera- >--yes->|  Dsolv1  |           __________
c   |   |        \ tion? /        |__________|<---|     |          |
c   A   |         \ ___ /               |         |---->|  Tinvit  |
c   |   A            |no                |               |__________|
c   |   |            |                  V
c   |   |            |                  |--->------>----->----->---->----->------|
c   |   |            V               ___                                         |
c   |   |            |             /     \          __________         _______   |
c   |   |            |            / Rela- \        |          |       |       |  |
c   |   |       _____V____   |--><  tivis- ><-yes->|  Difrel  |-error>|  Ext  |  |
c   |   A      |          |<-|    \  tic? /        |__________|       |_______|  V
c   |   |      |  Dsolv2  |        \ ___ /                                |      |
c   A   |      |__________|<-|        A                                 __V___   |
c   |   |            |       |        | no                             /      \  |
c   |   |            |       |        |                               |  STOP  | |
c   |   |            |       |        |                                \______/  V
c   |   |            |       |   _____V____          _______      ______         |
c   |   A            V       |  |          |        |       |    /      \        |
c   |   |            |       |  |  Difnrl  |-error->|  Ext  |-->|  STOP  |       |
c   |   |            |       |  |__________|        |_______|    \______/        |
c   A   |            |       |      ___                                          V
c   |   |            |       |    /     \          _________                     |
c   |   |            V       |   / Con-  \        |         |                    |
c   |   |            |       |-><   verg  ><-yes->|  Orban  |                    |
c   |   A            |           \  ed?  /        |_________|                    |
c   |   |            |            \ ___ /                                        V
c   |   |            |                                                           |
c   |   |            |<------<-----<-------<-------<-------<------<------<-----<-|
c   |   |            |                   __________
c   |   |            |                  |          |
c   |   A            |             |--->|  Splift  |
c   |   |       _____V_____        |    |__________|
c   |   |      |           |<------|   _______       ______
c   A   |      |           |          |       |     /      \
c   |   |      |  Velect   |--error-->|  Ext  |--->|  STOP  |
c   |   |      |           |          |_______|     \______/
c   |   A      |___________|<------|     _________
c   |   |            |             |    |         |
c   |   |            |             |--->|  Spliq  |
c   |   |            |                  |_________|    
c   A   |           _V_                            ___
c   |   A         /     \       __________       /Val- \         __________
c   |   |        / Con-  \     |          |     / ence  \       |          |
c   |   |       <   verg  >--->|  Etotal  |--->< Modify? >-yes->|  Vionic  |
c   |   |        \  ed?  /     |__________|     \       /       |__________|
c   |   |         \ ___ /                        \ ___ /              |
c   A   A            | no                           |no               |
c   |   |        ____V____                          V                 V
c   |   |       |         |            |<------------<--------------<-|
c   |   |       |  Dmixp  |            |           __________________________
c   |   |       |_________|           _V_         |              0)Pseudo    |
c   |   |            |              /     \       |              1)Pseudk    |
c   |   A           _V_            /Pseudo \      |  Pseudo-     2)Pseudt    |
c   |   |         /Pass \         < Generate>-yes>|  Potential   3)Pseudv    |
c   A   |        / Max.  \         \   ?   /      |  Generation  4)Datout    |
c   |   |<--no--< Itera-  >         \ ___ /       |  Block       5)Pseudb    |
c   |            \ tion? /             |no        |              6)Pseud2    |
c   |             \ ___ /              V          |__________________________|
c   |                | yes             |---<----------<-------|
c   |             ___V___             _V_
c   A            |       |          /     \        __________
c   |            |  Ext  |         /Config-\      |          |
c   |            |_______|        < uration >-yes>|  Prdiff  |
c   |                |             \ Test? /      |__________|
c   |             ___V__            \ ___ /            |
c   |            /      \              |no             |
c   |           |  STOP  |             |               V
c   A            \______/              |<---------<----|
c   |                                  V
c   |---------<------------<-----------|
c
c
c  *************************************************************
c  *     Program for atomic calculations                       *
c  *     Copyright Norman J. Troullier Jr &                    *
c  *     Jose Luis Martins, S. Froyen et al.                   *
c  *     Written by Norman J. Troullier Jr., Sept 89           *
c  *     while at U of Minn, from a Sverre Froyen              *
c  *     UC Berkeley code.  Program input/output is            *
c  *     compatible with earlier Berkeley versions.            *
c  *                                                           *
c  *     Send comments/suggestions/bug reports to:             *
c  *     troullie@cs.umn.edu                                   *
c  *     Version 5.06, Dated Oct. 19, 1990                     *
c  *     Jose.L.Martins@inesc.pt                               *
c  *                                                           *
c  *     Version 5.06, Dated Oct. 19, 1990                     *
c  *     Version 5.51, Dated June 1,1995                       *
c  *     Version 5.60, Dated January 29, 1997                  *
c  *                                                           *
c  *************************************************************
c
c    Some parameters are set inside the program,
c  the most important ones are:
c  1)the tolerance for selfconsistency in the screening 
c    potential (set in the main-atm program-named tol),
c  2)the accuracy in the eigenvalue for a given potential
c    (set in difnrl-named tol or difrel-named tol),
c  3)the dimensions of the work space used: nrmax,
c    norbmx, lmax(needs to be set only in main-atm),
c  4)the machine precision - MACHEP, for use in the 
c    eispack group of fuctions: tinvit, and tridib.
c    (The current value is ok for this application.)
c  5)the machine precision exp(-2*expzer), set in difnrl 
c    and difrel for the predictor-corrector methods 
c    (both named expzer), 
c  
c    For traditional applications the default values 
c  should be enough for 1-4.  For 5, expzer should be
c  found for the system in use.
c    NOTE: that for most VAX systems the value of expzer
c  can be very very small !!
c
c    The subroutine orban is called once for each orbital
c  after selfconsistency has been achieved.
c  Any additional analysis of the orbitals should therefore
c  be placed in orban.  Note that arrays ar and br have
c  different meaning for non-relativistic (ar-wave, 
c  br-d(wave)/dj) and relativistic (ar-major, br-minor)
c  calculations.
c
c    There are six ways to generate the pseudopotential :
c  ikerk = 6 Improved Troullier and Martins
c  ikerk = 5 Bachelet, Hamann, and Schluter
c  ikerk = 4 generates data file another pseudopotential
c   generation program.
c  ikerk = 3 Vanderbilt
c  ikerk = 2 Troullier and Martins 
c  ikerk = 1 Kerker
c  ikerk = 0 Hamann Schluter and Chiang
c
c      This main - atm routine has had extremly major
c    modifications with respect to the Berkeley version.
c    However, all input and output files are still compatible 
c    with earlier Berkeley versions of this program.
c   
c    1)Machine dependent timing calls were placed
c      in the program as dummy subroutines.
c      The user will either have to change
c      these calls for his machine or ask the authors
c      if they have subroutines for a specific computer.
c    2)The plot.dat file is now opened as a formatted file,
c      this is user/plot method dependent.  The atom.job
c      file is no longer used.  Note that the Apollo
c      system does not use standard Fortran methods to
c      open unformatted files.
c    3)The charge density startup is scaled with
c      an empirical formula that reduces the
c      number of iterations needed for the screening
c      potential convergence.
c    4)The screening potential mixing parameter is
c      an empirical function of the nuclear charge.
c      Larger atoms require a slower convergence
c      then smaller atoms.
c    5)The screening potential is intially mixed with a
c      percentage of old and new for the first itsm 
c      iterations. This brings the potential to a stable
c      region after which an Anderson's extrapolation scheme
c      is used.
c    6)The files pseudo.dat and plot.dat files are closed
c      and new ones are opened before the generation of a 
c      pseudopotential.  This allows the generation of 
c      multiple pseudopotentials(up to 99).
c    7)The pseudopotentail generation scheme of Troullier
c      and Martins is added - pseudt.  The pseudopotential
c      scheme of Vanderbilt has been added - pseudv. 
c      The improved pseudopotential scheme of Troullier
c      and Martins has been added - pseud2. 
c      The datout routine generates a data file for use 
c      in external pseudopotential generation programs.
c      The user may wish to alter for his own use or ignore.
c    8)Only the major modifications(not programming style) 
c      to the algorithms of each subroutine are commented
c      in that routine. 
c    9)Cray(and other machines) conversions are indicated 
c      at the begining of each routine.
c   10)The difrel and difnrl allow for the calculation of
c      a nonbound state(zero eigenvalue).  These states
c      can only be used in the pseudt, pseudk and
c      pseud2 generation 
c      routines.  The pseudo, pseudb and pseudv will fail with
c      a zero eigenvalue due to the generation method.
c      The user should be very careful in using a nonbound
c      state, and should always  compare the resulting pseudopotential
c      to a bound state pseudopotential calculation.
c   11)What comes free comes with no guarantee!!
c 
c  tolerance for self-consistency
c
c
      implicit real*8 (a-h,o-z)
c
      parameter (tol=1.d-8)
c
      parameter (zero=0.0d0,one=1.0d0)
c
      parameter(lmax=5,nrmax=9022,norbmx=40)
c
      dimension r(nrmax),rab(nrmax),no(norbmx),lo(norbmx),
     1 so(norbmx),zo(norbmx),cdd(nrmax),cdu(nrmax),cdc(nrmax),
     2 viod(lmax,nrmax),viou(lmax,nrmax),vid(nrmax),viu(nrmax),
     3 vod(nrmax),vou(nrmax),vn1d(nrmax),vn1u(nrmax),
     4 vn11d(nrmax),vn11u(nrmax),vn2d(nrmax),vn2u(nrmax),
     5 vn22d(nrmax),vn22u(nrmax),ev(norbmx),evi(norbmx),ek(norbmx),
     6 ep(norbmx),wk1(nrmax),wk2(nrmax),wk3(nrmax),wk4(nrmax),
     7 wk5(nrmax),wk6(nrmax),wk7(nrmax),wk8(nrmax),wk9(nrmax),
     8 wkb(7*nrmax)
c
      dimension econf(100),etot(10)
c
      character*1 ispp
      character*2 naold,icold,icorr,nameat
      character*10 plotfile
      character*12 pseudofile
c
c  njtj  ***  machine call  ***
c    Call to machine dependent cpu time routine.
c    User may want to comment out timing calls -
c    here and at exit of main - atm
c
      call zesec(t1)
c 
c  njtj  ***  machine call  ***
c
c  Startup values for doing multiple input data sets.
c
      naold = '  '
      icold = '  '
      zsold = zero
      nconf = 0
c
c      open files
c
      open(unit=1,file='pseudo.dat',form='unformatted',
     1  status='unknown')
      rewind(1)
c
c  njtj  ***  modification  start ***
c    The plot.dat file is now opened as a formatted file.
c    This is user/plot method dependent.  The atom.job
c    file is no longer used.
c
      open(unit=3,file='plot.dat',status='unknown',form='formatted')
      open(unit=5,file='atom.input',status='old',form='formatted')
      open(unit=6,file='atom.out',status='unknown',form='formatted')
      rewind(3)
      rewind(5)
      rewind(6)
c
c  njtj  ***  modification end  ***
c
c   Start of loop over configuration.
c   Read the input data.
c

 20   nr = nrmax
      norb = norbmx


      call input(itype,ikerk,icorr,ispp,zsh,rsh,
     1 nr,a,b,r,rab,nameat,norb,ncore,no,lo,so,zo,
     2 znuc,zel,evi,nval_orig)

c
c  njtj  ***  machine call  ***
c  Stop - no more data in input file,
c  Find time taken for total calculation.
c  zesec - machine dependent routine
c
      if (itype .lt. 0) then
        call zesec(t2)
        write(6,2000)t2-t1
 2000 format(//,' The total time for the calculation is ',
     1 f12.5,' seconds')
        call ext(0)
      endif
c 
c  njtj  *** machine call ***
c
c  Jump charge density 'set up' and ionic data input if
c  configuration test.
c
      itsm=znuc/9+3
      if (zsold .eq. zsh .and. naold .eq. nameat 
     1 .and. itype .lt. 1 ) then
      else
        if (itype .lt. 4) then
c
c  Set up the initial charge density.
c  cdd and cdu  =  (4 pi r**2 rho(r))/2
c
c  njtj  ***  modification  ***
c    The charge density setup (aa) is scaled with
c    an empirical formula that reduces the
c    number of iterations needed for the screening
c    potential convergence.
c
          aa = sqrt(sqrt(znuc))/2+one
          do 30 i=1,nr
            cdd(i) = zel*aa**3*exp(-aa*r(i))*r(i)**2/4
            cdu(i) = cdd(i)
 30       continue
        endif                
c  
c  njtj  ***  modification end  ***
c
c  set up ionic potentials
c
        call vionic(ispp,itype,icorr,ifcore,zsh,rsh,
     1   lmax,nr,a,b,r,rab,nameat,ncore,znuc,
     2   cdd,cdu,cdc,viod,viou)
      endif
c
c   Set up the electronic potential.
c
      call velect(0,0,icorr,ispp,ifcore,
     1 nr,r,rab,zel,cdd,cdu,cdc,vod,vou,etot,wk1,wk2,
     2 wk3,wk4,wk5,wkb)
c
      do 50 i=1,nr
        vid(i) = vod(i)
        viu(i) = vou(i)
 50   continue
c
c   Start the iteration loop for electronic convergence.
c
      iconv = 0
      icon2 = 0
      maxit = 100
c
c  njtj  ***  modification start  ***
c    The screening potential mixing parameter is
c    an empirical function of the nuclear charge.
c    Larger atoms require a slower convergence
c    then smaller atoms.
c
      xmixo = one/log(znuc+7*one)
c
c  njtj  ***  modifications end  ***
c
      do 100 iter=1,maxit
        if (iter .eq. maxit) iconv=1
c
c  compute orbitals
c
        if (icon2 .lt. 2) then
          call dsolv1(lmax,nr,a,b,r,rab,norb,ncore,
     1     no,lo,so,zo,cdd,cdu,viod,viou,vid,viu,ev,
     2     wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wkb)
        else
          call dsolv2(iter,iconv,ispp,ifcore,lmax,nr,
     1     a,b,r,rab,norb,ncore,no,lo,so,zo,znuc,cdd,
     2     cdu,cdc,viod,viou,vid,viu,ev,ek,ep,wk1,wk2,
     3     wk3,wk4,wk5,wk6,wk7,evi)
        endif
c
c  set up output electronic potential from charge density
c
        call velect(iter,iconv,icorr,ispp,ifcore,
     1   nr,r,rab,zel,cdd,cdu,cdc,vod,vou,etot,wk1,wk2,
     2   wk3,wk4,wk5,wkb)
c
c  check for convergence
c
        if (iconv .gt. 0) goto 120
        dvmax = zero
        do 60 i=1,nr
          dv = (vod(i)-vid(i))/(1.D0+vod(i)+vou(i))
          if (abs(dv) .gt. dvmax) dvmax=abs(dv)
          dv = (vou(i)-viu(i))/(1.D0+vou(i)+vod(i))
          if (abs(dv) .gt. dvmax) dvmax=abs(dv)
 60     continue
        icon2 = icon2+1
        if (dvmax .le. tol) iconv=1
c
c  Mix the input and output electronic potentials.
c
c  njtj  ***  modification  start  ***
c    The screening potential is initially mixed with a
c    percentage of old and new for itsm iterations.
c    This brings the potential to a stable region
c    after which an Anderson's extrapolation scheme
c    is used.
c
        if (iter .lt. itsm) then
          iiter=2
        else
          iiter=iter-itsm+3
        endif
        call dmixp(vod,vid,xmixo,iiter,3,nr,wk1,wk2,
     1   vn1d,vn11d,vn2d,vn22d)
        call dmixp(vou,viu,xmixo,iiter,3,nr,wk1,wk2,
     1   vn1u,vn11u,vn2u,vn22u)
 100  continue
c
c   End of iteration of electronic convergence loop.
c
      write(6,110) dvmax,xmixo
 110  format(/,' potential not converged - dvmax =',e10.4,
     1 ' xmixo =',f5.3)
      call ext(1)
c
c  njtj  ***  modification end  ***
c
c  Find the total energy.
c
 120  write(6,121)icon2
 121  format(/,'Total number of iterations needed for',
     1 ' electron screening potential is ',i2,/)
      call etotal(itype,zsh,nameat,norb,
     1 no,lo,so,zo,etot,ev,ek,ep)
c
c   Replace the valence charge density.
c
      if (itype .eq. 5) call vionic(ispp,6,icorr,
     1 ifcore,zsh,rsh,lmax,nr,a,b,r,rab,nameat,
     2 ncore,znuc,cdd,cdu,cdc,viod,viou)
c
c  Pseudopotential generation.
c                     
c  njtj  ***  modification  start  ***
c    Current pseudo.dat and plot.dat files are closed
c    and new ones are opened.  This allows the 
c    generation of multiple pseudopotentials(up to 99).
c
      if (itype .ge.1 .and. itype .le. 3) then
        if (ikerk .ne. 4 ) then
          close(unit=1)
          close(unit=3)
          if (nconf .le. 8) then
            write(pseudofile,8000)nconf+1
            write(plotfile,8002)nconf+1
          else
            write(pseudofile,8001)nconf+1
            write(plotfile,8003)nconf+1
          endif                        
          write(6,8004)nconf+1
 8000 format('pseudo.dat0',i1)  
 8001 format('pseudo.dat',i2)
 8002 format('plot.dat0',i1)
 8003 format('plot.dat',i2)   
 8004 format(//,' Pseudopotential generation file number ',i2)  
c
          open(unit=1,file=pseudofile,form='unformatted',
     1     status='unknown')
          open(unit=3,file=plotfile,status='unknown',
     1     form='formatted')
           rewind(1)
           rewind(3)
        endif
c
c  njtj  ***  modification  end  ***
c
c  njtj  ***  modification start  ***
c    The pseudopotentail generation scheme of Troullier 
c    and Martins is added - pseudt.  The pseudopotential
c    scheme of Vanderbilt is added - pseudv.  The 
c    pseudopotential scheme of BHS is added - pseudb.
c    The improved pseudopotential scheme of Troullier 
c    and Martins is added - pseud2.  The 
c    dataout routine generates a data file for use in 
c    external pseudopotential generation programs. 
c    The user may wish to alter for their own use or ignore.
c
        if(ikerk.eq.0) then
          call pseudo(itype,icorr,ispp,lmax,nr,a,b,r,rab,
     1     nameat,norb,ncore,no,lo,so,zo,znuc,zel,cdd,cdu,cdc,
     2     viod,viou,vid,viu,vod,vou,etot,ev,ek,ep,wk1,wk2,
     3     wk3,wk4,wk5,wk6,wk7,wk8,wk9,vn1d,vn1u,vn2d,vn2u,
     4     vn11d,wkb,evi,nval_orig)
        elseif (ikerk .eq. 1) then
          call pseudk(itype,icorr,ispp,lmax,nr,a,b,r,rab,
     1     nameat,norb,ncore,no,lo,so,zo,znuc,zel,cdd,cdu,cdc,
     2     viod,viou,vid,viu,vod,vou,etot,ev,ek,ep,wk1,wk2,
     3     wk3,wk4,wk5,wk6,wk7,wk8,wk9,vn1d,vn1u,wkb,evi,
     4     nval_orig)
        elseif (ikerk .eq. 2) then
          call pseudt(itype,icorr,ispp,lmax,nr,a,b,r,rab,
     1     nameat,norb,ncore,no,lo,so,zo,znuc,zel,cdd,cdu,cdc,
     2     viod,viou,vid,viu,vod,vou,etot,ev,ek,ep,wk1,wk2,
     3     wk3,wk4,wk5,wk6,wk7,wk8,wk9,vn1d,vn1u,wkb,evi,
     4     nval_orig)
        elseif (ikerk .eq. 3) then
          call pseudv(itype,icorr,ispp,lmax,nr,a,b,r,rab,
     1     nameat,norb,ncore,no,lo,so,zo,znuc,zel,cdd,cdu,cdc,
     2     viod,viou,vid,viu,vod,vou,etot,ev,ek,ep,wk1,wk2,
     3     wk3,wk4,wk5,wk6,wk7,wk8,wk9,vn1d,vn1u,vn2d,vn2u,
     4     vn11d,wkb,evi,nval_orig)
        elseif (ikerk .eq. 4) then
          call datout(itype,icorr,ispp,lmax,nr,a,b,r,rab,
     1     nameat,norb,ncore,no,lo,so,zo,znuc,zel,cdc,
     2     viod,viou,vid,viu,ev) 
        elseif (ikerk .eq. 5) then
          call pseudb(itype,icorr,ispp,lmax,nr,a,b,r,rab,
     1     nameat,norb,ncore,no,lo,so,zo,znuc,zel,cdd,cdu,cdc,
     2     viod,viou,vid,viu,vod,vou,etot,ev,ek,ep,wk1,wk2,
     3     wk3,wk4,wk5,wk6,wk7,wk8,wk9,vn1d,vn1u,vn2d,vn2u,
     4     vn11d,wkb,evi,nval_orig)
        elseif (ikerk .eq. 6) then
          call pseud2(itype,icorr,ispp,lmax,nr,a,b,r,rab,
     1     nameat,norb,ncore,no,lo,so,zo,znuc,zel,cdd,cdu,cdc,
     2     viod,viou,vid,viu,vod,vou,etot,ev,ek,ep,wk1,wk2,
     3     wk3,wk4,wk5,wk6,wk7,wk8,wk9,vn1d,vn1u,wkb,evi,
     4     nval_orig)
        endif
      endif 
c 
c  njtj   ***  modification end  ***
c
c  printout difference from first configuration
c
      nconf = nconf + 1
      econf(nconf) = etot(10)
      if(naold .eq. nameat .and. icold .eq. icorr .and. nconf .ne. 1
     1 .and. (itype .lt. 1 .or. itype .gt. 3)) then
        call prdiff(nconf,econf)
        write(6,130) etot(10)-econf(1)
      endif
      write(6,135)
 130  format(//,' excitation energy         =',f18.8,/)
 135  format(//,60('%'),//)
      naold = nameat
      icold = icorr
      zsold = zsh
c
c   End loop of configuration.
c
cccc      goto 20
      end
