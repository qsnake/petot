c
c Copyright (C) 2002 Vanderbilt group
c This file is distributed under the terms of the GNU General Public
c License as described in the file 'License' in the current directory.
c
      program runatom
c
c----------------------------------------------------------------------------
c
c     program for obtaining self-consistent solutions of the
c     kohn-sham equations for atoms and generating ultra soft
c     separable pseudopotentials
c
c     written by dominic king-smith from programs due to david vanderbilt
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
      parameter( idim1 = 1000      , idim2 = 26         )
c     be careful if idim1 is changed, then also the dimension of 
c     qfplot in pssubs.f and ndm in funz must be changed

      parameter( idim3 = 10        , idim4 = 5          )
      parameter( idim5 = 4         , idim6 = 2*idim5-1  )
      parameter( idim7 = 4         , idim8 = 20         )
c
c     idim1  .ge.  no. of points in the radial mesh
c     idim2  .ge.  no. of shells in the atom
c     idim3  .ge.  no. of beta functions in vanderbilt potential
c     idim4  .ge.  no. of valence states for pseudization
c     idim5  .ge.  lmax + 1, lmax is maximum l value of potential
c     idim7  .ge.  no. of reference states for each angular momentum value
c     idim8  .ge.  no. of terms in taylor expansion of pseudized q_i,j
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1),sqr(idim1)
c.....potentials for the all electron calculation
      dimension ruae(idim1),runuc(idim1)
c.....hartree, exchange and correlation potentials
      dimension ruhar(idim1),reexch(idim1),ruexch(idim1)
c.....charge densities
      dimension rscore(idim1),rsvale(idim1),rsatom(idim1),rspsco(idim1)
      dimension rspscoef(idim8)
c.....wavefunctions
      dimension snl(idim1,idim2),snlo(idim1),yy(idim1,2)
c.....shell labels, energies, occupancies and real-space cutoff index
      dimension nnlz(idim2),ee(idim2),wwnl(idim2),nkkk(idim2),
     +  krel(idim2)
c.,...work arrays for the tcm routine call
      dimension bndmat(idim1,5),npivot(idim1)
c.....pseudo quantum numbers, energies occupancies and cutoff radii
      dimension nnlzps(idim4),eeps(idim4),wwnlps(idim4),rc(idim5),
     +rinner(idim6)
c.....beta and q functions and q pseudization coefficients
      dimension beta(idim1,idim3),qfunc(idim1,idim3,idim3),
     +qfcoef(idim8,idim6,idim3,idim3),vloc(idim1),vloc0(idim1)
c.....extra fixed points to avoid negative densities when pseudizing q_ij
      dimension ibfix(idim3),rfix(idim3),qfix(idim3)
c.....work space for solution to inhomogeneous equation and overlaps
      dimension eta(idim1,idim7),alpha(idim3,idim2)
c.....angular momenta, reference energy, qij and dij of vanderbilt scheme
      dimension lll(idim3),eee(idim3),iptype(idim3),qqq(idim3,idim3),
     +ddd0(idim3,idim3),ddd(idim3,idim3),sss(idim3,idim3),
     +ttt(idim3,idim3),bbb(idim3,idim3),bbbi(idim3,idim3),iwork(idim3)
c.....indexing for the beta functions
      dimension nbl0(idim5),nbl(idim5),nbl0d(idim5),nbld(idim5)
c.....psi,phi and chi functions for generation of vanderbilt potential
      dimension psi(idim1,idim3),phi(idim1,idim3),chi(idim1,idim3)
c.....version number and date
      dimension iver(3),idmy(3)
c.....type of wave function
      dimension ikeyee(idim3)
      dimension refwf(idim1,idim3)
c.....scratch arrays
      dimension xwa(idim1,5)

c.....this is used only in trou.f
c
c.....title, and exchange correlation type
      character*20 title,xctype
c.....file name array
      character*40 flname(6)
      character*40 filewfref
c
c      external idate
c
c=======================================================================
c     c h a n g e s  m a d e  t o  c r e a t e  v e r s i o n  3.0.0
c=======================================================================
c
c     [1] added nodless pseudization of wavefunctions using exponential
c         form.  controlled by iptype switch: 
c                if  iptype = 0    original polynomial pseudization
c                if  iptype = 1    exponential pseudization
c
c     [2] added new pseudization of the local potential to get log
c         derivative correct in one particular angular momentum channel
c         at one particular energy with local potential alone.  
c         controlled by lloc and eloc switches:
c                if  lloc  = -1    original pseudization scheme
c                if  lloc  = +ve   get l'deriv correct in channel lloc
c                at the energy specified by eloc
c
c     [3] added new pseudization of q_ij function using rappe, rabe,
c         kaxiras and joannopoulos type approach.  controlled by the
c         switches ifqopt, nqf and qtryc switches:
c                if  ifqopt = 0    original polynomial pseudization
c                if  ifqopt = 1    new scheme with taylor expansion
c                of nqf terms optimising for cutoff wavevector qtryc
c
c     [4] added version number and date to program.  version number
c         should greatly ease the maintanance of backward compatibilty
c         between potential and car-parrinello program in future 
c         upgrades
c
c     [5] split the variable nvales into two new variables nvales and
c         nang.  nang is now lmax + 1, where lmax is the maximum angular
c         momentum of the beta functions.  nvales is the number of 
c         states to be treated as valence states.  this allows greater
c         flexibilty in the use of the programs and for example 
c         affords the oportunity to treat 3s,3p,4s and 3d states in
c         titanium as valence states using only up to d non-locality
c         in the projector functions.
c
c     [6] general improvement of print outs and fourier transforms of
c         q_ij and the pseudo wavefunctions
c
c     [7] added sweeping of energy on first call to schgef on
c         pseudopotential transferability tests to aid convergence.
c
c=======================================================================
c     c h a n g e s  m a d e  t o  c r e a t e  v e r s i o n  4.0.0
c=======================================================================
c
c     [8] added soft pseudization of the wavefunctions.  controlled by 
c         iptype switch: 
c                if  iptype = 2    soft pseudization
c         cutoff and number of terms in taylor series set by nqf and
c         qtryc parameters.
c         ************************************************************
c         note: the program takes the nqf and qtryc parameters that
c         are specified for the q-functions, and uses them also for
c         the wavefunctions.  this seems rather artificial and I
c         would like to change it in a future version.  -dv 7/97
c         ************************************************************
c
c     [9] added check to ensure that there are no negative densities
c         introduced in the process of pseudizing q_ij function. 
c
c    [10] added source code of nag routines for generalized eigenvalue
c         problems (instead of eispak).
c
c    [11] added printout to check tails of potentials are precisely
c         zv / r
c
c    [12] added routines for relativistic corrections.  controlled
c         by the switch irel.  irel = 1 is for dirac calculations
c         (all electron only).  irel = 2 is for koelling and harmon
c         scalar relativistic equation (all electron and pseudopotential
c         tests and generation).
c
c=======================================================================
c     c h a n g e s  m a d e  t o  c r e a t e  v e r s i o n  4.1.0
c=======================================================================
c
c    [13] added extra constraint to force charge density to be positive
c         when pseudizing q_ij (ifqopt=2).
c
c=======================================================================
c     c h a n g e s  m a d e  t o  c r e a t e  v e r s i o n  5.0.0
c=======================================================================
c
c    [14] when ifqopt=2, program now pseudizes q_ij correctly
c         for each l = l_min,l_max,2 .  also there are cosmetic
c         improvements, eg in printouts.
c
c=======================================================================
c     c h a n g e s  m a d e  t o  c r e a t e  v e r s i o n  5.0.1
c=======================================================================
c
c    [15] treatment of tail of atom is improved, so that the "alpha"
c         term used to correct the ewald energy is less sensitive
c         to numerical errors:
c         (a) on last self-consistent iteration, charge is updated
c             without damping
c         (b) a polynomial interpolation is used to take exc and uxc
c             smoothly to zero for densities less than 1.0d-15.
c         (c) program checks convergence of alpha and stops if need be.
c
c=======================================================================
c     c h a n g e s  m a d e  t o  c r e a t e  v e r s i o n  5.1.0
c=======================================================================
c
c    [16] allows variable r_inner as a function of l for pseudization
c         of the q_ij
c         additional plot of all-electron and L-dependent pseudized
c         q_ij function
c
c
c=======================================================================
c     c h a n g e s  m a d e  t o  c r e a t e  v e r s i o n  6.0.0
c=======================================================================
c
c      Kurt Stokbro: stokbro@bohr.sissa.it
c    [17] calculation of eigenstates in besselfunction basisset
c
c    [18] prints also the wavefunctions to the pseudo file
c
c    [19] keyee =3 produces norm conserved wavefunctions
c         used for constructing mixed normconserved/ultrasoft pseudo-
c         potentials, see stokbro PRB 53, 6869 (1996)
c
c    [20] one may readin a second reference state, and use that to
c         construct pseudo wavefunctions. 
c
c=======================================================================
c     c h a n g e s  m a d e  t o  c r e a t e  v e r s i o n  7.0.0
c=======================================================================
c
c      Kurt Stokbro: stokbro@mic.dtu.dk
c    [21]   Included gga exchange correlation
c         if (exfact.eq. 1.) xctype=' C-A + B88gx + LYPgc'
c         if (exfact.eq. 2.) xctype=' C-A + B88gx        '
c         if (exfact.eq. 3.) xctype=' C-A + B88gx + P86gc'
c         if (exfact.eq. 4.) xctype='       PW(91)       '
c
c    [22] one may write out a variable number of atomic wave functions
c
c    [23] psqf2 has been changed to include 5 constraints:
c         area of r^l Q, Q(r), d1 Q(r), d2 Q(r), d3 Q(r)
c         *********************************************************
c         NOTE: WITH THIS CHANGE, THE PROGRAM IS NO LONGER BACKWARD
c         COMPATIBLE!!  IE, IT GENERATES A DIFFERENT PSEUDOPOTENTIAL
c         THAN BEFORE WHEN ACTING ON AN IDENTICAL INPUT FILE (WHEN
c         IFQOPT=2).
c         *********************************************************
c         basically, when ifqopt=2, the program now uses 5 constraints
c         instead of 3 when choosing the polynomial for the q-function
c         pseudization; the matching of second and third derivatives
c         was added.  presumably this is because in the gga schemes the
c         matching of higher derivatives of rho is needed if the xc
c         potential is to remain smooth.
c
c    [24] The nonlinear core correction has been included, note
c         the change in the inputfile and the pseudopotential file
c
c=======================================================================
c     c h a n g e s  m a d e  t o  c r e a t e  v e r s i o n  7.0.1
c=======================================================================
c
c         No functional changes, but calls to essl subroutines replaced
c         by calls to nag subroutines
c
c         Note that improvements in going from 5.1.0 to 5.1.3 are not
c         all incorporated.
c
c=======================================================================
c     c h a n g e s  m a d e  t o  c r e a t e  v e r s i o n  7.0.2
c=======================================================================
c
c    [25] Modified so as to use eispack routines, so that all proprietary
c         (essl, nag) routines are gone
c
c=======================================================================
c     c h a n g e s  m a d e  t o  c r e a t e  v e r s i o n  7.0.2_NG
c=======================================================================
c
c    [26] Add option for PBE(96) GGA:
c         if (exfact.eq. 5.) xctype='       PBE-GGA      '
c
c         I requested to Keith to incorporate the improvements in
c         going up to a5.1.3 into this code.  I think he said he
c         did this but I haven't checked.  Also, he reports some
c         bugs.
c
c=======================================================================
c     c h a n g e s  m a d e  t o  c r e a t e  v e r s i o n  7.0.3
c=======================================================================
c
c         Miscellaneous bug fixes and cosmetic improvements.
c
c=======================================================================
c     c h a n g e s  m a d e  t o  c r e a t e  v e r s i o n  7.1.0
c=======================================================================
c
c    [27] Restore backwards compatibility with a5.1.4 regarding
c         pseudization of the q-functions.  The choices are now:
c
c           ifqopt  subrout  action
c           ------  -------  ------------------------------------
c           0       psqf     original polynomial expansion
c           1       psqf1    optimize to nqf,qtryc; 3 constraints
c           2       psqf2    optimize to nqf,qtryc; 3 constraints
c                              plus optional nfix,ifix
c           3       psqf2    optimize to nqf,qtryc; 5 constraints
c                              plus optional nfix,ifix
c
c         Note that input files intended for gga runs have to be
c         modified to set ifqopt=3 instead of 2.
c
c         Note also that outputs of a5.1.3 and a5.1.4 have minor
c         differences for irel=2 (see note [19] in a5.1.4).  This
c         version should match a5.1.4, not a5.1.3, in that case.
c
c=======================================================================
c     c h a n g e s  m a d e  t o  c r e a t e  v e r s i o n  7.2.0
c=======================================================================
c
c    [28] Provide separate input parameters npf,ptryc for pseudizing
c         the wavefunctions, so that they do not automatically have
c         to be the same as nqf and qtryc.
c
c     *** This destroys compatibility with input files of previous
c     *** versions, but the program will issue a warning and
c     *** terminate if it thinks it is reading an old-style
c     *** input file.
c
c=======================================================================
c     c h a n g e s  m a d e  t o  c r e a t e  v e r s i o n  7.2.1
c=======================================================================
c
c    [29] Fix bug in the printout of the pseudopotential report
c
c    [30] Fix bug in aesubs that affects xc energy in case of pbe
c
c    [31] Fix bug in pswf1: affects pseudo wf when iptype=1; and
c         also affects local potential except when lloc=-1 (***)
c
c    [32] Some cosmetic improvements
c
c    *** Egad, this bug has been present since a3.0 !
c    *** This bug fix also destroys backward compatibility: this
c    *** version run on old input file can generate a different psp
c
c=======================================================================
c     c h a n g e s  m a d e  t o  c r e a t e  v e r s i o n  7.3
c=======================================================================
c
c    [33] Additions of C. Pickard 2/98 to handle f-electrons
c
c    [34] The initialization of the inhomogenous solution 'eta' of
c         the radial schroedinger equation has been changed in
c         subroutine schgef: it is now set to have zero value and
c         slope at the origin.  This seems to prevent numerical
c         problems that arose with f-states in some cases, and
c         otherwise has little practical effect.
c
c    [35] Improved stability of convergence on eigenvalue ecur in
c         subroutines schsol, dirsol, koesol
c
c    [36] Some minor cosmetic improvements
c
c=======================================================================
c     c h a n g e s  m a d e  t o  c r e a t e  v e r s i o n  7.3.1
c=======================================================================
c
c    [37] Add support in makefile for ibm compiler
c
c    [38] Add ability to generate plots of beta and chi functions
c
c    [39] Add example input files for generating potentials needed
c         for CUSP tutorials
c
c=======================================================================
c     c h a n g e s  m a d e  t o  c r e a t e  v e r s i o n  7.3.2
c=======================================================================
c
c    [40] Fixed minor problem that made reported "target and achieved
c         values of dqfunc" appear to differ more than they should.
c
c    [41] Overhauled 'ifprt' print control throughout the program
c         (see INPUT_GEN documentation for details).
c
c=======================================================================
c     c h a n g e s  f o r  c h a n g e o v e r  t o  cusp-733
c=======================================================================
c
c    [42] Added GNU GPL license information to each subroutine.
c
c    [43] Fixed write statement formats in psqf1 and pswf2.
c
c    [44] Replaced call to dsqrt by call to sqrt in corpbe.
c
c    [45] Fixed minor bug in reform.f.  This bug arises:
c
c         ONLY if you use 'reform.f', which compiles to become
c         'reform.x'.  This is a standalone utility program.  Note
c         that the main executable 'runatom.x' is *NOT* affected.
c
c         ONLY if the number of pseudo valence wavefunctions is > 4.
c
c=======================================================================


c
c.....variable for file = 0
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c.....set the version number
      data (iver(i),i=1,3) /7,3,3/
c
c     unit stderr=0 automatically connected to stderr
c
c     set up the files common block
      stderr = 0
      input  = 10 
      iout   = 6
      ioae   = 3
      iplot  = 15
      iologd = 16
      iops   = 14
      iwfref = 18
c     
c     get filenames from command argument string
      do 10 ifile=1,6
        call getarg(ifile,flname(ifile))
   10 continue
      if (iargc() .gt. 6) then
c write extra ref filename
         call getarg(7,filewfref)
         open( unit = iwfref , file = filewfref , status = 'old' )
         read(iwfref,*) igrid,nwf
c always save potential,wavefunction
c first read the all electron wavefunctions         
         read(iwfref,*) ((refwf(i,3+j),i=1,igrid),j=1,nwf)
c next read the all electron potential, screning potential of ref system, and last
c the current screening potential.
         read(iwfref,*) ((refwf(i,j),i=1,igrid),j=1,3)
         close(iwfref)
      endif
c
c     get the date of current run
c      call idate(idmy)
       idmy=0
c
c     open the file for standard out put of the program
c
      open( unit = iout , file = flname(2) , status = 'unknown' )
c
c     print a header line in output file
      
      write(iout,30) (iver(i),i=1,3),idmy(2),idmy(1),idmy(3)
   30 format(72('='),/,'pseudopotential program version',i3,'.',
     +  i1,'.',i1,'   date:',i3,' - ',i2,' - ',i4,/,72('='))
c
c     and also notify terminal of startup
      write (stderr,20) (iver(i),i=1,3)
   20 format('beginning execution pseudopotential program version',i3,
     +  '.',i1,'.',i1)
c
c     now do everything
c
      call atom(r,rab,sqr,ruae,runuc,ruhar,reexch,ruexch,
     +  rscore,rsvale,rsatom,rspsco,rspscoef,snl,snlo,nnlz,ee,wwnl,
     +  nkkk,krel,yy,xwa,title,xctype,flname,bndmat,npivot,
     +  nnlzps,eeps,wwnlps,rc,lll,eee,iptype,qqq,ddd0,ddd,
     +  sss,ttt,bbb,bbbi,iwork,alpha,eta,beta,qfunc,qfcoef,
     +  vloc,vloc0,nbl0,nbl,nbl0d,nbld,psi,phi,chi,iver,idmy,
     +  ibfix,rfix,qfix,rinner,
     +  ikeyee,ikeyeeloc,refwf,
     +  idim1,idim2,idim3,idim4,idim5,idim6,idim7,idim8)
c
      end
c
c----------------------------------------------------------------------------
c
      subroutine atom(r,rab,sqr,ruae,runuc,ruhar,reexch,ruexch,
     +  rscore,rsvale,rsatom,rspsco,rspscoef,snl,snlo,nnlz,ee,wwnl,
     +  nkkk,krel,yy,xwa,title,xctype,flname,bndmat,npivot,
     +  nnlzps,eeps,wwnlps,rc,lll,eee,iptype,qqq,ddd0,ddd,
     +  sss,ttt,bbb,bbbi,iwork,alpha,eta,beta,qfunc,qfcoef,
     +  vloc,vloc0,nbl0,nbl,nbl0d,nbld,psi,phi,chi,iver,idmy,
     +  ibfix,rfix,qfix,rinner,
     +  ikeyee,ikeyeeloc,refwf,
     +  idim1,idim2,idim3,idim4,idim5,idim6,idim7,idim8)
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1),sqr(idim1)
c.....potentials for the all electron calculation
      dimension ruae(idim1),runuc(idim1)
c.....hartree, exchange and correlation potentials
      dimension ruhar(idim1),reexch(idim1),ruexch(idim1)
c.....charge densities
      dimension rscore(idim1),rsvale(idim1),rsatom(idim1),rspsco(idim1)
      dimension rspscoef(idim8)
c.....wavefunctions
      dimension snl(idim1,idim2),snlo(idim1),yy(idim1,2)
c.....shell labels, energies, occupancies and real-space cutoff index
      dimension nnlz(idim2),ee(idim2),wwnl(idim2),nkkk(idim2),
     +  krel(idim2)
c.,...work arrays for the tcm routine call and relativity routines
      dimension bndmat(idim1,5),npivot(idim1)
c.....pseudo quantum numbers, energies occupancies and cutoff radii
      dimension nnlzps(idim4),eeps(idim4),wwnlps(idim4),rc(idim5)
c.....l-dependent pseudization qfunction cutoff radii
      dimension rinner(idim6)
c.....beta and q functions and q pseudization coefficients
      dimension beta(idim1,idim3),qfunc(idim1,idim3,idim3),
     +qfcoef(idim8,idim6,idim3,idim3),vloc(idim1),vloc0(idim1)
c.....extra fixed points to avoid negative densities when pseudizing q_ij
      dimension ibfix(idim3),rfix(idim3),qfix(idim3)
c.....work space for solution to inhomogeneous equation and overlaps
      dimension eta(idim1,idim7),alpha(idim3,idim2)
c.....angular momenta, reference energy, qij and dij of vanderbilt scheme
      dimension lll(idim3),eee(idim3),iptype(idim3),qqq(idim3,idim3),
     +ddd0(idim3,idim3),ddd(idim3,idim3),sss(idim3,idim3),
     +ttt(idim3,idim3),bbb(idim3,idim3),bbbi(idim3,idim3),iwork(idim3)
c.....indexing for the beta functions
      dimension nbl0(idim5),nbl(idim5),nbl0d(idim5),nbld(idim5)
c.....psi,phi and chi functions for generation of vanderbilt potential
      dimension psi(idim1,idim3),phi(idim1,idim3),chi(idim1,idim3)
c.....version number and date
      dimension iver(3),idmy(3)
c.....logic to determine whether to scan in scheqg
      logical lsweep
c.....type of wave function
      dimension ikeyee(idim3)
      dimension refwf(idim1,idim3)
c.....scratch arrays
      dimension xwa(idim1,5)
c
c.....title, and exchange correlation type
      character*20 title,xctype
c.....file name array
      character*40 flname(6)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c     ------------------------------------------------------
c     note: in this program:
c     potentials, e.g. rucore, are really r*v(r)
c     wave funcs, e.g. snl, are really proportional to r*psi(r)
c       and are normalized so int dr (snl**2) = 1
c     thus psi(r-vec)=(1/r)*snl(r)*y_lm(theta,phi)
c     conventions carry over to beta, etc
c     charge dens, e.g. rscore, really 4*pi*r**2*rho
c     ------------------------------------------------------
c
c     rydberg units are used throughout this program
c
c     ------------------------------------------------------
c
c     ******************************************
c     r o u t i n e  i n i t i a l i s a t i o n
c     ******************************************
c
c     read in the basic parameters determining the run
      call readin( + 1 ,flname,ifae,ifpsp,ifprt,ifplw,ilogd,
     +  rlogd,emin,emax,nnt,thresh,tol,damp,maxit,title,z,xion,
     +  exfact,xctype,rmax,aasf,bbsf,ncspvs,irel,nnlz,wwnl,ee,
     +  ncores,nvales,nang,keyps,lmaxps,ifpcor,rpcor,rinner,nbeta,
     +  rcloc,rc,lll,eee,iptype,npf,ptryc,nbl,nbl0,zv,lloc,eloc,
     +  ifqopt,nqf,qtryc,iploctype,besrmax,besemin,besemax,besde,
     +  ikeyee,ikeyeeloc,
     +  ibfix,rfix,qfix,nfix,idim2,idim3,idim4,idim5,idim6,idim7,
     +  idim8)
c
c     decide where to get the data for the all electron calculation
      if ( ifae .eq. 0 ) then
c
c       read in the data from a previous run
        call rwae( + 1 ,z,xion,exfact,mesh,irel,aasf,bbsf,rmax,
     +    ncspvs,nnlz,ee,wwnl,nkkk,snl,ruae,flname,idim1,idim2)
c
      elseif ( ifae .eq. 1 ) then
c
c       read in controls for a new all electron calculation
        call readin( +2 ,flname,ifae,ifpsp,ifprt,ifplw,ilogd,
     +    rlogd,emin,emax,nnt,thresh,tol,damp,maxit,title,z,xion,
     +    exfact,xctype,rmax,aasf,bbsf,ncspvs,irel,nnlz,wwnl,ee,
     +    ncores,nvales,nang,keyps,lmaxps,ifpcor,rpcor,rinner,nbeta,
     +    rcloc,rc,lll,eee,iptype,npf,ptryc,nbl,nbl0,zv,lloc,eloc,
     +    ifqopt,nqf,qtryc,iploctype,besrmax,besemin,besemax,besde,
     +    ikeyee,ikeyeeloc,
     +    ibfix,rfix,qfix,nfix,idim2,idim3,idim4,idim5,idim6,idim7,
     +    idim8)
c
      endif
c
      if ( ifpsp .gt. 0 ) then
c
c       read in controls for pseudopotential test or generation
        call readin( +3 ,flname,ifae,ifpsp,ifprt,ifplw,ilogd,
     +    rlogd,emin,emax,nnt,thresh,tol,damp,maxit,title,z,xion,
     +    exfact,xctype,rmax,aasf,bbsf,ncspvs,irel,nnlz,wwnl,ee,
     +    ncores,nvales,nang,keyps,lmaxps,ifpcor,rpcor,rinner,nbeta,
     +    rcloc,rc,lll,eee,iptype,npf,ptryc,nbl,nbl0,zv,lloc,eloc,
     +    ifqopt,nqf,qtryc,iploctype,besrmax,besemin,besemax,besde,
     +    ikeyee,ikeyeeloc,
     +    ibfix,rfix,qfix,nfix,idim2,idim3,idim4,idim5,idim6,idim7,
     +    idim8)
c
c
      endif
c
c
c     generate the logarithmic mesh
      call rinit(r,rab,sqr,rmax,aasf,bbsf,a,b,z,rlogd,klogd,
     +  mesh,ifprt,idim1)
c
c     establish the starting nuclear charge configuration
      call setnuc(z,r,runuc,sigma,mesh,idim1,ifprt)
c
c     possible generation of an analytic starting potential
      if ( ifae .eq. 1 ) then
c
        ipass = 1
        call startv(ipass,z,xion,mesh,r,runuc,ruae,
     +    wwnl,wwnlps,nvales,nvalps,ncores,zv,vloc,vloc0,
     +    idim1,idim2,idim4)
c
       endif
c
c     *********************************************
c     a l l  e l e c t r o n  c a l c u l a t i o n
c     *********************************************
c
      write (stderr,*) 'beginning the all electron calculation'
c
      if ( irel .eq. 0 ) then
        write(iout,9998) irel,'schroedinger'
      elseif ( irel .eq. 1 ) then
        write(iout,9998) irel,'dirac'
      elseif ( irel .eq. 2 ) then
        write(iout,9998) irel,'koelling-harmon'
      endif
 9998 format(' irel =',i2,' so all electron calculations',
     +  ' use ',a,' equations')
c
c     (mainly for gga): when starting from scratch, start with
c       lda and switch to true xc choice when close to convergence
      exfacttru = exfact
      if (ifae.eq.1) exfact = 0.
c
c     iterative loop to converge the potential of all electron calculation
      do 100 iter = 1,maxit
c
c
        if ( irel .eq. 0 ) then
c
c         solve the schroedinger equation
          call schsol(ruae,r,rab,sqr,a,b,xwa(1,1),z,ee,nnlz,wwnl,nkkk,
     +      snl,snlo,ncspvs,mesh,thresh,ifprt,idim1,idim2)
c
        elseif ( irel .eq. 1 ) then
c
c         solve the dirac equation
          call dirsol(ruae,r,rab,bndmat,ee,nnlz,wwnl,nkkk,krel,
     +      snl,yy,ncspvs,mesh,thresh,ifprt,idim1,idim2)
c
        elseif ( irel .eq. 2 ) then
c
          call koesol(ruae,r,rab,bndmat,ee,nnlz,wwnl,nkkk,
     +      snl,yy,xwa,ncspvs,mesh,thresh,ifprt,idim1,idim2)
c
        endif
c
c       compute the charge density from the wave functions
        call rsae(snl,wwnl,nkkk,nnlz,rscore,rsvale,ncores,ncspvs,
     +    mesh,alpha,nbl,nbl0,rab,1,0,kkbeta,qfunc,xwa(1,1),
     +    ifprt,idim1,idim2,idim3,idim5)
c
c       compute the exchange and correlation potential
        call vhxc( + 1 ,ifpcor,ruexch,ruhar,reexch,rscore,
     +    rsvale,rsatom,rspsco,xwa(1,1),exfact,r,rab,
     +    xwa(1,2),xwa(1,3),xwa(1,4),ehar,emvxc,mesh,idim1)
c
c       test for self-consistency and mix old and new potentials
        call mixer(ruae,runuc,ruhar,ruexch,xwa(1,1),damp,
     +    tol,delta,iselfc,mesh,idim1)
c
        write(stderr,9990) iter,delta
        if (delta .le. 1.d-1  .and. exfact .ne. exfacttru) then
           exfact = exfacttru
        endif
c       check whether self-consistency has been achieved
        if ( iselfc .eq. + 1 ) goto 200
c
c       check whether iterations complete but still not consistent
        if ( iselfc .eq. - 1 .and. iter .eq. maxit ) then
          write(iout,*) '***error in subroutine atom'
          write(iout,*) 'potential could not be converged'
          call exit(1)
        endif
c
 9990   format('completed self-consistent cycle',i4,' delta =',d15.7)
c
  100 continue
c
  200 continue
c
c     compute the total energy of the atom and print results
      call eigprt(delta,nnlz,ee,wwnl,nkkk,ncspvs,ehar,emvxc,a,b,
     +  r,ruae,rsatom,rscore,xwa(1,1),mesh,ifprt,idim1,idim2)
c
c     write the results of the all electron calculation to file
      if ( ifae .eq. 1 ) then
        call rwae( - 1 ,z,xion,exfact,mesh,irel,aasf,bbsf,rmax,
     +    ncspvs,nnlz,ee,wwnl,nkkk,snl,ruae,flname,idim1,idim2)
      endif
c
c     set up for print outs of wavefunctions
      ipass = 1
      call prwf(ifprt,ncores,nvales,ncspvs,nnlz,ee,nkkk,
     +  r,rab,a,b,rsvale,snl,ruae,ipass,flname,ifplw,mesh,idim1,idim2)
c
c     compute logarithmic derivatives
      zz = 0.0d0
      vzero = ruae(2) / r(2)
      ipass = 1
c
c     no logderivatives if irel is 1 (dirac equation)
      if ( irel .ne. 1 )
     +  call lderiv(ifprt,ipass,keyps,ilogd,klogd,irel,emax,emin,nnt,
     +  zz,vzero,b,r,rab,sqr,ruae,xwa(1,1),snlo,yy,flname,mesh,
     +  vloc,beta,nbl,nbl0,xwa(1,2),alpha,eta,kkbeta,ddd,qqq,
     +  bndmat,idim1,idim3,idim5,idim7)
c
      write(stderr,*) 'all electron calculation completed'
c
      if ( ifpsp .eq. 0 ) return
c
c     **************************************
c     test of the ultra soft pseudopotential
c     **************************************
c
      if ( ifpsp .eq. 1 ) then
c
        write(stderr,*) 'reading the pseudopotential from file'
c
        call rwps( + 1 ,title,zps,zvps,exftps,nvalps,irel,
     +    z,zv,exfact,nvales,mesh,etot,nnlzps,wwnlps,eeps,
     +    nnlz,wwnl,ee,ncores,keyps,ifpcor,rpcor,rinner,rc,nbeta,
     +    kkbeta,lll,eee,iptype,npw,ptryc,
     +    beta,ddd0,ddd,qqq,qfunc,qfcoef,
     +    rcloc,vloc0,vloc,snl,lloc,eloc,rsatom,rspsco,flname,
     +    ifqopt,nqf,qtryc,nang,nbl0,nbl,r,rab,iver,idmy,
     +    idim1,idim2,idim3,idim4,idim5,idim6,idim8)
c
c       for pseudopotential test ipass parameter is + 2
        ipass = 2
c
c       generate a starting guess for the potential
        call startv(ipass,z,xion,mesh,r,runuc,ruae,
     +    wwnl,wwnlps,nvales,nvalps,ncores,zv,vloc,vloc0,
     +    idim1,idim2,idim4)
c
c       generate ddd for this starting potential
        keydir = + 1
        call dddscr(keydir,r,rab,vloc,qfunc,xwa(1,1),kkbeta,nbeta,
     +    lll,ddd,ddd0,mesh,ifprt,idim1,idim3)
c
        write(stderr,*) 'beginning pseudopotential test'
c
c       iterative loop to converge the test pseudopotential configuration
        do 300 iter = 1,maxit
c
c         scan energy in scheqg if first cycle
          if ( iter .eq. 1 ) then
            lsweep = .true.
          else
            lsweep = .false.
          endif
c
c         solve the schroedinger equation
          call scheqg(vloc,r,rab,sqr,a,b,xwa(1,1),z,ee,rcloc,rc,
     +      nnlz,wwnl,nkkk,snl,snlo,ncores,ncspvs,mesh,
     +      thresh,beta,nbl,nbl0,xwa(1,2),alpha,eta,klogd,kkbeta,
     +      ddd,qqq,ifprt,lsweep,idim1,idim2,idim3,idim5,idim7,1)
c
c         compute the charge density from the wave functions
          call rsae(snl,wwnl,nkkk,nnlz,rscore,rsvale,ncores,ncspvs,
     +      mesh,alpha,nbl,nbl0,rab,ipass,keyps,kkbeta,qfunc,xwa(1,1),
     +      ifprt,idim1,idim2,idim3,idim5)
c
c         compute the exchange and correlation potential
          call vhxc(ipass,ifpcor,ruexch,ruhar,reexch,rscore,
     +      rsvale,rsatom,rspsco,xwa(1,1),exfact,r,rab,
     +      xwa(1,2),xwa(1,3),xwa(1,4),ehar,emvxc,mesh,idim1)
c
c         test for self-consistency and mix old and new potentials
          call mixer(vloc,vloc0,ruhar,ruexch,xwa(1,1),damp,
     +      tol,delta,iselfc,mesh,idim1)
c
c         progress report
          write(stderr,9989) iter,delta
 9989     format('completed self-consistent cycle',i4,' delta =',d15.7)
c
c         check whether self-consistency has been achieved
          if ( iselfc .eq. + 1 ) goto 400
c
c         update the ddd matrix
          keydir = + 1
          call dddscr(keydir,r,rab,vloc,qfunc,xwa(1,1),kkbeta,nbeta,
     +      lll,ddd,ddd0,mesh,ifprt,idim1,idim3)
c
c         check whether iterations complete but still not consistent
          if ( iselfc .eq. - 1 .and. iter .eq. maxit ) then
            write(iout,*) '***error in subroutine atom'
            write(iout,*) 'potential could not be converged'
            call exit(1)
          endif
c
  300   continue
  400   continue
          call scheqg(vloc,r,rab,sqr,a,b,xwa(1,1),z,ee,rcloc,rc,
     +      nnlz,wwnl,nkkk,snl,snlo,ncores,ncspvs,mesh,
     +      thresh,beta,nbl,nbl0,xwa(1,2),alpha,eta,klogd,kkbeta,
     +      ddd,qqq,ifprt,lsweep,idim1,idim2,idim3,idim5,idim7,0)

c
        write(stderr,*) 'pseudopotential test converged'
        if ( ifplw .eq. 2 ) then
           write(iplot,*) (ruhar(i)+ruexch(i),i=1,mesh)
        endif

c
      endif
c
c     *************************************
c     ultra soft pseudopotential generation
c     *************************************
c
      if ( ifpsp .eq. 2 .and. keyps .eq. 3 ) then
c
c       generate the vanderbilt soft core pseudopotential
        call scpgef(r,rab,sqr,a,b,rlogd,klogd,irel,ruae,
     +    zz,vzero,mesh,ifprt,lmaxps,nbeta,rcloc,rc,rinner,
     +    lll,eee,iptype,npf,ptryc,iploctype,nbl,nbl0,
     +    nang,ddd,ddd0,kkbeta,beta,qqq,qfunc,qfcoef,ifqopt,nqf,qtryc,
     +    lloc,eloc,sss,ttt,bbb,bbbi,iwork,ibfix,rfix,qfix,nfix,
     +    psi,phi,chi,vloc,yy,bndmat,npivot,
     +    xwa(1,2),xwa(1,3),xwa(1,4),
     +    ikeyee,ikeyeeloc,refwf,
     +    idim1,idim3,idim5,idim6,idim7,idim8)
c
c       print the real-space qfunctions
        call prqf(ifprt,qfunc,kkbeta,nbeta,r,a,b,idim1,idim3)
c
        write(stderr,*) 'fourier transforming the qfunctions'
c
c       fourier transform and print reciprocal-space q functions
        call fanalq(r,rab,qfunc,nbeta,xwa(1,1),mesh,ifprt,
     +    idim1,idim3)

c       calculate the core charge
        if (ifpcor .gt. 0) then
           do i=1,mesh
              rspsco(i) = rscore(i)
           enddo
c test ?
c          write(6,*)rspsco,rspscoef,xwa(1,2),rpcor,mesh,r,rab
c          write(6,*)a,b,ifprt,nqf,qtryc,idim1,idim8
c
           call pspcor(rspsco,rspscoef,xwa(1,2),rpcor,mesh,r,rab,
     +          a,b,ifprt,nqf,qtryc,idim1,idim8)
c
c       print the real-space core charge
           call prcor(ifprt,rspsco,rscore,rsvale,kkbeta,r,a,b,idim1)
c
           write(stderr,*) 'fourier transforming the core charge'
c
c       fourier transform and print reciprocal-space core charge
           call fanalcor(r,rab,rspsco,rscore,xwa(1,1),mesh,ifprt,
     +          idim1)
        endif


c       solve the schroedinger equation
c
c       no scanning of log derivatives so lsweep = .false.
        lsweep = .false.
c
        write(stderr,*) 'solving the schroedinger equation'
c
        call scheqg(vloc,r,rab,sqr,a,b,xwa(1,1),z,ee,rcloc,rc,
     +    nnlz,wwnl,nkkk,snl,snlo,ncores,ncspvs,mesh,
     +    thresh,beta,nbl,nbl0,xwa(1,2),alpha,eta,klogd,kkbeta,
     +    ddd,qqq,ifprt,lsweep,idim1,idim2,idim3,idim5,idim7,0)
c
c       correct the charge density
        call rsae(snl,wwnl,nkkk,nnlz,rscore,rsvale,ncores,ncspvs,
     +    mesh,alpha,nbl,nbl0,rab,2,keyps,kkbeta,qfunc,xwa(1,1),
     +    ifprt,idim1,idim2,idim3,idim5)
c
        write(stderr,*) 'descreening the potential'
c
c       descreen both the vloc and ddd in the vanderbilt scheme
        call descrn(vloc,vloc0,ddd,ddd0,ifpcor,mesh,
     +       idim5,ikeyee,ikeyeeloc,refwf,ruae,
     +    ruexch,ruhar,reexch,rscore,rsvale,rsatom,rspsco,
     +    exfact,r,rab,xwa,etot,wwnl,ee,ncores,ncspvs,zv,
     +    keyps,nbeta,kkbeta,lll,qfunc,ifprt,idim1,idim2,idim3)

        if ( ifplw .eq. 2 ) then
           write(iplot,*) (ruhar(i)+ruexch(i),i=1,mesh)
        endif
c
        write(stderr,*) 'writing the pseudopotential to file'
c
        call rwps( - 1 ,title,z,zv,exfact,nvales,irel,
     +    z,zv,exfact,nvales,mesh,etot,
     +    nnlz(ncores+1),wwnl(ncores+1),ee(ncores+1),
     +    nnlz,wwnl,ee,ncores,keyps,ifpcor,rpcor,rinner,rc,nbeta,
     +    kkbeta,lll,eee,iptype,npf,ptryc,
     +    beta,ddd0,ddd,qqq,qfunc,qfcoef,
     +    rcloc,vloc0,vloc,snl,lloc,eloc,rsatom,rspsco,flname,
     +    ifqopt,nqf,qtryc,nang,nbl0,nbl,r,rab,iver,idmy,
     +    idim1,idim2,idim3,idim4,idim5,idim6,idim8)
c
        if ( ifprt .ge. 2 ) then
c
          write (stderr,*) 'computing logarithmic derivatives of vloc'
          write (iout,'(/a)')
     +      ' calling lderiv: log derivs of LOCAL potential only'
c
          do i = 1,nang
            nbl0d(i) = 0
            nbld(i) = 0
          end do
c
c         compute logarithmic derivatives
          zz = 0.0d0
          vzero = ruae(2) / r(2)
          ipass = 2
          call lderiv(ifprt,ipass,keyps,ilogd,klogd,irel,emax,emin,nnt,
     +      zz,vzero,b,r,rab,sqr,ruae,xwa(1,1),snlo,yy,flname,mesh,
     +      vloc,beta,nbld,nbl0d,xwa(1,2),alpha,eta,kkbeta,ddd,qqq,
     +      bndmat,idim1,idim3,idim5,idim7)
c
        endif
c
      endif
c
c     set up for print outs of wavefunctions
      ipass = 2
      call prwf(ifprt,ncores,nvales,ncspvs,nnlz,ee,nkkk,
     +  r,rab,a,b,rsvale,snl,ruae,ipass,flname,ifplw,mesh,idim1,idim2)
c
c     fourier analyse the pseudo wavefunctions
      write(stderr,*) 'fourier analysing pseudowavefunctions'
      call fanal(r,rab,nnlz,wwnl,nkkk,snl,ncores,
     +  ncspvs,xwa(1,1),mesh,ifprt,idim1,idim2)
c
c     compute logarithmic derivatives
      write (stderr,*) 'computing the logarithmic derivatives'
      if ( ifprt .ge. 1) write (iout,'(/a)')
     +    ' calling lderiv: log derivatives of pseudopotential'
c
      zz = 0.0d0
      vzero = ruae(2) / r(2)
      ipass = 2
      call lderiv(ifprt,ipass,keyps,ilogd,klogd,irel,emax,emin,nnt,
     +  zz,vzero,b,r,rab,sqr,ruae,xwa(1,1),snlo,yy,flname,mesh,
     +  vloc,beta,nbl,nbl0,xwa(1,2),alpha,eta,kkbeta,ddd,qqq,
     +  bndmat,idim1,idim3,idim5,idim7)
c
      write (stderr,*) 'beginning search for possible ghost states'
c
c     find energy eigenvalues of nonlocal part by itself
      if ( ifprt .ge. 1) then
        write (iout,'(/1x,a/)') 'eigenvalues of nonlocal potential '
        do lp=1,nang
          nn=nbl(lp)
          ns=nbl0(lp)+1
          if (nn.gt.0) then
            call qdiag1(lp,nn,bbbi(ns,ns),qqq(ns,ns),sss(ns,ns),
     +        ttt(ns,ns),ddd(ns,ns),beta(1,ns),kkbeta,idim1,idim3)
          endif
        end do
      endif
c
c     solve the ultrasoft secular equation in a bessel-function basis
c
c     begin loop over bescut values
      bescut = besemin
 600  if (bescut .gt. besemax) goto 700
c
      if ( ifprt .ge. 0) then
c
        write (stderr,'(1x,a,f6.2)')
     +    'eigensolution in bessel basis with cutoff',bescut
        write (iout,'(/1x,49(1h-)/1x,a,f6.2/1x,49(1h-))')
     +    'eigensolution in bessel basis with cutoff',bescut
c
        nebes = 4
        call besdiag(ifprt,iout,besrmax,bescut,nang-1,mesh,
     +     r,rab,vloc,kkbeta,idim1,idim3,nbl0,nbl,ddd,qqq,beta,nebes)
c
      endif
c
c     end loop over bescut values
      bescut = bescut+besde
      goto 600
 700  continue

      return
      end
