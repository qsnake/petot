c
c Copyright (C) 2002 Vanderbilt group
c This file is distributed under the terms of the GNU General Public
c License as described in the file 'License' in the current directory.
c
c----------------------------------------------------------------------------
c
      subroutine rwae(iread,z,xion,exfact,mesh,irel,aasf,bbsf,rmax,
     +  ncspvs,nnlz,ee,wwnl,nkkk,snl,ruae,flname,idim1,idim2)
c
c     suboutine for reading and writing the all electron data
c     to file.
c              iread = + 1 ......... read ae data in
c              iread = - 1 ......... write ae data out
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c
c.....potentials for the all electron calculation
      dimension ruae(idim1)
c.....wavefunctions
      dimension snl(idim1,idim2)
c.....shell labels, energies, occupancies and real-space cutoff index
      dimension nnlz(idim2),ee(idim2),wwnl(idim2),nkkk(idim2)
c
c.....file name array
      character*40 flname(6)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
      if ( iread .eq. + 1 ) then
c
c       read in the data from a previous calculation
c
        open( unit = ioae , file = flname(3) , status = 'old' ,
     +        form = 'unformatted' )
        read (ioae) z,xion,exfact,mesh,irel,aasf,bbsf,rmax
        read (ioae) ncspvs
        do 542 k=1,ncspvs
  542   read (ioae) nnlz(k),wwnl(k),ee(k)
        read (ioae) (ruae(i),i=1,mesh)
        do 544 k=1,ncspvs
        read (ioae) kkk,(snl(i,k),i=1,kkk)
  544   nkkk(k)=kkk
        close (ioae)
c
      elseif ( iread .eq. - 1 ) then
c
c       write out the data from the all electron calculation
c
c       write data for later pseudo generation on unit ioae
        open( unit = ioae , file = flname(3) , status = 'unknown' , 
     +        form = 'unformatted' )
        write (ioae) z,xion,exfact,mesh,irel,aasf,bbsf,rmax
        write (ioae) ncspvs
        do 552 k=1,ncspvs
  552   write (ioae) nnlz(k),wwnl(k),ee(k)
        write (ioae) (ruae(i),i=1,mesh)
        do 554 k=1,ncspvs
        kkk=nkkk(k)
  554   write (ioae) kkk,(snl(i,k),i=1,kkk)
c
      else
c
        write(iout,*) '***error in subroutine rwae'
        write(iout,*) 'iread = ',iread,' is not allowed'
        call exit(1)
c
      endif
c
      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine readin(igoto,flname,ifae,ifpsp,ifprt,ifplw,ilogd,
     +  rlogd,emin,emax,nnt,thresh,tol,damp,maxit,title,z,xion,
     +  exfact,xctype,rmax,aasf,bbsf,ncspvs,irel,nnlz,wwnl,ee,
     +  ncores,nvales,nang,keyps,lmaxps,ifpcor,rpcor,rinner,nbeta,
     +  rcloc,rc,lll,eee,iptype,npf,ptryc,nbl,nbl0,zv,lloc,eloc,
     +  ifqopt,nqf,qtryc,iploctype,besrmax,besemin,besemax,besde,
     +  ikeyee,ikeyeeloc,
     +  ibfix,rfix,qfix,nfix,idim2,idim3,idim4,idim5,idim6,idim7,
     +  idim8)
c
c     routine for general reading in of controls for 
c     all electron and pseudopotential calculations
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c.....shell labels, energies, occupancies and real-space cutoff index
      dimension nnlz(idim2),ee(idim2),wwnl(idim2)
c.....pseudo cutoff radii
      dimension rc(idim5),rinner(idim6)
c.....angular momenta, reference energy, qij and dij of vanderbilt scheme
      dimension lll(idim3),eee(idim3),iptype(idim3)
c.....extra fixed points to avoid negative densities when pseudizing q_ij
      dimension ibfix(idim3),rfix(idim3),qfix(idim3)
c.....indexing for the beta functions
      dimension nbl0(idim5),nbl(idim5)
      dimension ikeyee(idim3)
c
c.....test whether any iptype=2
      logical ttype2
c
c.....title, and exchange correlation type
      character*20 title,xctype
c.....file name array
      character*40 flname(6)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
      if ( igoto .lt. 1 .or. igoto .gt. 3 ) then
        write(iout,*) '***error in subroutine readin'
        write(iout,*) 'igoto =',igoto,' is out of range'
        call exit(1) 
      endif
c
c     igoto = 1  ; read fundamental switches to determine run type
c     igoto = 2  ; read control cards for the all electron calculation
c     igoto = 3  ; read controls for vanderbilt potential
c
      goto (100,200,300) igoto
c
c--------------------------------------------------------------------------
c
  100 continue
c
c     r e a d  i n  b a s i c  p a r a m e t e r s  f o r  r u n
c
      open ( unit = input , file = flname(1) , status = 'old' )
c
c     ------------------------------------------------------
c
      read (input,110) ifae,ifpsp,ifprt,ifplw,ilogd
  110 format (6i5)
      write (iout,120) ifae,ifpsp,ifprt,ifplw,ilogd
  120 format (/' ifae =',i4,4x,'ifpsp =',i4,4x,'ifprt =',i4,4x,
     1   'ifplw =',i4,4x,'ilogd =',i4)
c
c     ------------------------------------------------------
c     ifae:    0  read all-electron ae data file
c              1  do ae from scratch and write ae data file
c     ------------------------------------------------------
c     ifpsp:   0  stop after ae part
c              1  read pseudo data file and compare to ae
c              2  generate pseudopotential and compare to ae
c     ------------------------------------------------------
c     ifprt:  -3  no prints or graphs
c             -2  prints but no graphs
c             -1  more or less standard
c              1  everything
c     ------------------------------------------------------
c     ifplw  if generate data file for wf plots
c     ------------------------------------------------------
c     ilogd  num l values for which log derivs calc'd and output
c     ------------------------------------------------------
c
      read (input,130) rlogd,emin,emax,nnt
  130 format (3f10.5,i5)
      write (iout,140) rlogd,emin,emax,nnt
  140 format (7h rlogd=,f10.5,5x,5hemin=,f10.5,5x,
     1   5hemax=,f10.5,5x,4hnnt=,i5)
c
c     rlogd is radius at which log derivs are calc'd
c     energy mesh spans (emin,emax) with nnt intervals
c
      read (input,150) thresh,tol,damp,maxit
  150 format (2e10.1,f10.5,i5)
      if ( thresh .eq. 0.0d0 ) thresh = 1.0d-06
      if ( tol    .eq. 0.0d0 )    tol = 1.0d-05
      if ( damp   .eq. 0.0d0 )   damp = 5.0d-01
      if ( maxit  .eq. 0     )  maxit = 250
      write (iout,160) thresh,tol,damp,maxit
  160 format (' thresh, tol =',1p2e12.3,4x,'damp =',0pf6.3,
     1   4x,'maxit =',i4)
      read (input,170) title
  170 format (a20)
c
c     end read in of fundamental switches for the calculation
      return
c
c--------------------------------------------------------------------------
c
  200 continue
c
c     c o n t r o l  c a r d s  f o r  a l l  e l e c t r o n  c a s e
c
c     (ifae = 1)
c
c     run a self-consistency check
      if ( ifae .ne. 1 ) then
        write(iout,*) '***error in subroutine readin'
        write(iout,*) 'reading all electron controls but ifae=',ifae
        call exit(1)
      endif
c
c     get the nuclear charge, net charge on atom and 
c     type of exchange and correlation
      read (input,210) z,xion,exfact
  210 format (f5.0,2f10.5)
c
c     get the parameters for generating the logarithmic mesh
      read (input,215) rmax,aasf,bbsf
c     write(iout,*) 'rmax=',rmax,'aasf=',aasf,'bbsf=',bbsf
  215 format (3f10.5)
      if ( rmax .eq. 0.0d0 ) rmax = 80.0d0
      if ( aasf .eq. 0.0d0 ) aasf =  6.0d0
      if ( bbsf .eq. 0.0d0 ) bbsf = 40.0d0
c
c     number core states plus valence states and whether relativistic
      read (input,220) ncspvs,irel
  220 format (2i5)
c     write(iout,*) 'ncspvs=',ncspvs,'irel=',irel
c
c     check that ncspvs is within range
      if ( ncspvs .gt. idim2 ) then
        write(iout,*) '***error in subroutine readin'
        write(iout,*) 'ncspvs=',ncspvs,' gt idim2=',idim2
        call exit(1)
        endif
c     check that irel is zero or one
      if ( irel .lt. 0 .or. irel .gt. 2 ) then
        write(iout,*) '***error in subroutine readin'
        write(iout,*) 'irel =',irel,' is out of range'
        call exit(1)
        endif
c
c     find quantum numbers, occupancies and estimated energies
      read (input,230) (nnlz(i),wwnl(i),ee(i),i=1,ncspvs)
  230 format(i4,f7.3,f14.6)
c
c       wwnl < 0 used to specify that orbital should be skipped,
c         e.g., if it is not bound
c
      do 240 i=1,ncspvs
  240   if ( wwnl(i) .lt. 0.0d0 ) wwnl(i) = -1.d-20
c
c     perform a self-consistency check
      www = 0.0d0
      do 243 i = 1,ncspvs
  243   www = www + wwnl(i)
      if ( abs(z-www-xion) .gt. 1.e-3 ) then
        write (iout,*) '***error in subroutine readin'
        write (iout,245) www,xion,z,ncores,nvales,ncspvs
  245   format(6h www= ,f4.0,7h xion= ,f4.0,6h   z= ,f4.0,
     +    11h   ncores= ,i4,11h   nvales= ,i4,11h   ncspvs= ,i4/
     +    25h control cards incorrect? )
        call exit(1) 
      endif
c
c     note: for consistency with ifpsp>0, "valence" orbitals
c           should be listed last in order s or s,p or s,p,d !
c
c     exfact:  0  ceperley-alder (perdew-zunger param)
c             -1  wigner ex-corr
c             -2  hedin-lundquist ex-corr
c             -3  gunnarson-lundquist ex-corr
c              1  LDA (CA) + BLYP GC xc 
c              2  LDA (CA) + Becke GC x 
c              3  LDA (CA) + BP GC xc 
c              4         PW(91)
      
      if (exfact.eq. 0.) xctype='      ceperley-alder'
      if (exfact.eq.-1.) xctype='              wigner'
      if (exfact.eq.-2.) xctype='     hedin-lundqvist'
      if (exfact.eq.-3.) xctype=' gunnarson-lundqvist'
      if (exfact.eq. 1.) xctype=' C-A + B88gx + LYPgc'
      if (exfact.eq. 2.) xctype=' C-A + B88gx        '
      if (exfact.eq. 3.) xctype=' C-A + B88gx + P86gc'
      if (exfact.eq. 4.) xctype=' Perdew Wang 1991   '
      if (exfact.eq. 5.) xctype=' PBE - GGA          '
c
c     print out of all electron control parameters
      write (iout,250) title,xctype
  250 format (/10x,61(1h*)/10x,2a20,' exchange-correlation'
     1        /10x,61(1h*)/)
      write (iout,260) z,xion,exfact,irel
  260 format (' z =',f5.0,4x,'xion =',f5.2,4x,
     1   'exfact =',f10.5,4x,'irel=',i4)
      write (iout,270) ncspvs
  270 format (' ncspvs =',i5)
      write (iout,280) rmax,aasf,bbsf 
  280 format (' rmax =',f10.5,4x,'aasf =',f10.5,4x,'bbsf =',f10.5)
c
c     set the number of core and valence states
      ncores = ncspvs
      nvales = 0
c
c     read in of all electron control cards completed
      return
c--------------------------------------------------------------------------
c
  300 continue
c
c     r e a d  i n  o f  p s e u d o p o t e n t i a l  c o n t r o l s
c
c     ( ifpsp = 1 or ifpsp = 2 )
c
c     run a self-consistency check
      if ( ifpsp .lt. 1 .or. ifpsp .gt. 2 ) then
        write(iout,*) '***error in subroutine readin'
        write(iout,*) 'reading for pseudopotential but ifpsp=',ifpsp
        call exit(1)
      endif
c
c     check that this is not a dirac calculation
      if ( irel .eq. 1 ) then
        write(iout,*) '***error in subroutine readin'
        write(iout,*) 'no dirac pseudopotential test or generation'
        call exit(1)
      endif
c
c     get the number of core and valence states
      read (input,310) ncores,nvales,nang
  310 format (3i5)
      write (iout,320) ncores,nvales,nang
  320 format (' ncores =',i5,4x,'nvales =',i5,4x,'nang =',i5)
c
      if ( ncores + nvales .ne. ncspvs ) then
        write(iout,*) '***error in subroutine readin'
        write(iout,*) 'ncores+nvales=',ncores+nvales,
     +    ' but ncspvs =',ncspvs
        call exit(1)
      endif
c
      if ( nvales .gt. idim4 ) then
        write(iout,*) '***error in subroutine readin'
        write(iout,*) 'nvales =',nvales,' but idim4 =',idim4
        call exit(1)
      endif
c
      if ( nang .gt. idim5 ) then
        write(iout,*) '***error in subroutine readin'
        write(iout,*) 'nang =',nang,' but idim5 =',idim5
        call exit(1)
      endif
c
c     compute the valence charge
      wwc = 0.0d0
      do 430 i = 1,ncores
        wwc = wwc + wwnl(i)
  430 continue
      zv = z - wwc
      write(iout,440) wwc,zv
  440 format(' core charge =',f9.5,' and valence charge =',f9.5)
c
c   get information for the besselfunction calculation
        read(input,515) besrmax,besemin,besemax,besde
 515    format(4f10.5)
        write(iout,517) besrmax,besemin,besemax,besde
 517    format(' besrmax =',f6.2,4x,'besemin =',f6.2,4x,'besemax =',
     +  f6.2,4x,'besde=',f6.2)

c     possible return if pseudopotential test
      if ( ifpsp .eq. 1 ) return
c
      read (input,330) keyps,ifpcor,rinner1,rpcor
  330 format (2i5,2f10.5)
      write (iout,340) keyps,ifpcor,rinner1,rpcor
  340 format (' keyps =',i5,3x,'ifpcor =',i5,3x,'rinner =',f10.5, 
     +     3x,'rpcor =',f10.5)
      if (rinner1 .eq. 0.0) then
        read (input,335) (rinner(i),i=1,nang*2-1)
  335   format (5f10.5)
        write (iout,345) (rinner(i),i=1,nang*2-1)
  345   format ('   rinner for L=0,1,2... =',5f10.5)
      else
        do i=1,nang*2-1
        rinner(i)=rinner1
        enddo
      endif
c
c     ------------------------------------------------------------
c     keyps = 0 --> standard hsc pseudopotential with exponent 4.0
c     keyps = 1 --> standard hsc pseudopotential with exponent 3.5
c     keyps = 2 --> vanderbilt modifications using defaults
c     keyps = 3 --> new generalized eigenvalue pseudopotentials
c     keyps = 4 --> frozen core all-electron case
c     ifpcor = 1 if "partial core correction" of louie, froyen,
c       & cohen to be used; 0 otherwise
c     rinner = L-dependent radius at which to cut off partial core or q_ij
c     for true frozen core case, use keyps=4, ifpcor=1, rinner=0.
c     rpcor the cutoff radius of the partial core. The partial core is
c          pseudized the same way as the qfunctions
c     -----------------------------------------------------------
c
      if ( keyps .ne. 3 ) then
        write(iout,*) '***error in subroutine readin'
        write(iout,*) 'keyps =',keyps,' not programed'
        call exit(1)
      endif
c
c      if ( ifpcor .ne. 0 ) then
c        write(iout,*) '***error in subroutine readin'
c        write(iout,*) 'ifpcor =',ifpcor,' not programed'
c        call exit(1)
c      endif
c
      do i = 1, nang*2-1
      if ( rinner(i) .gt. rlogd ) then
        write(iout,*) '***error in subroutine readin'
        write(iout,*) 'rinner =',rinner(i),' but rlogd =',rlogd
        call exit(1)
      endif
      enddo
c
c     maximum angular momentum of potential is  s  with nang = 1
c     maximum angular momentum of potential is  p  with nang = 2
c     maximum angular momentum of potential is  d  with nang = 3
c
      if ( keyps .eq. 3 ) then
c
c       ---------------------------------------------------------------
c       nbeta = number of beta's (no of l-epsilon combinations)
c       rcloc = core radius used to pseudize the local potential
c       rc are cut-off radii for pseudo gen
c       lll = angular momentum l
c       keyee = 0 : read in reference energy epsilon
c       keyee > 0 : set reference energy epsilon equal to ae eigenvalue
c       ---------------------------------------------------------------
c
        read (input,350) nbeta,rcloc
  350   format (i5,f10.5)
c
        if ( nbeta .gt. idim3 ) then
          write(iout,*) '***error in subroutine readin'
          write(iout,*) 'nbeta =',nbeta,' but idim3',idim3
          call exit(1)
        endif
        write (iout,360) nbeta,rcloc
  360   format (' nbeta =',i3,5x,'rcloc =',f8.3)
c
        read (input,370) (rc(i),i=1,nang)
  370   format (4f10.5)
        write (iout,380) (rc(i),i=1,nang)
  380   format (' rc for s, p, d, f =',4(f10.5,3x))
c
        do 385 i = 1,idim5
          nbl(i) = 0
  385   continue
        lmaxps = 0
        ttype2=.false.
c
        do 410 ib = 1,nbeta
c
          read (input,390) ll,keyee,eeread,ipread
  390     format (2i5,f10.5,i5)
c
c         check that ll s have been specified in correct order
c         must be input first s, then p, then d beta's
          if ( ll .lt. lmaxps ) then
            write(iout,*) '***error in subroutine readin'
            write(iout,*) 'll =',ll,' is .lt. lmaxps =',lmaxps
            call exit(1)
          endif 
c
          lmaxps=ll
c
          lll(ib)=ll
          nbl(ll+1)=nbl(ll+1)+1
c
c         check that idim7 is sufficient
          if ( nbl(ll+1) .gt. idim7 ) then
            write(iout,*) '***error in subroutine readin'
            write(iout,*) 'nbl(ll+1)=',nbl(ll+1),' .gt. idim7=',idim7
            call exit(1)
          endif
c
          ikeyee(ib) = keyee
          if (keyee.le.0) then
c           read in reference energy explicitly
            eee(ib)=eeread
          elseif (keyee .le. nvales) then
c           set equal to one of the ae eigenvalues
            eee(ib)=ee(ncores+keyee)
          else
            write(iout,*) '***error in subroutine readin'
            write(iout,*) 'keyee .gt. nvales, keyee =',keyee
            call exit(1)
          endif
c
c         determine pseudisation type 0=polynomial 1=exponential
          if ( ipread .lt. 0 .or. ipread .gt. 4 ) then
            write(iout,*) '***error subroutine readin'
            write(iout,*) 'ipread =',ipread,' not allowed'
            call exit(1)
          endif
c
          iptype(ib) = ipread
c
          if (ipread.eq.2) ttype2=.true.
c
          write (iout,400) ib,lll(ib),eee(ib),iptype(ib)
  400     format (' ib, lll, eee, iptype =',2i4,f8.3,i4)
c
  410   continue
c
c       nang is by definition lmaxps + 1 so run a check
        if ( nang .ne. lmaxps + 1 ) then
          write(iout,*) '***error in subroutine readin'
          write(iout,*) ' nang .ne. lmaxps + 1, nang,lmaxps',nang,lmaxps
          call exit(1)
        endif
c
c       beta's for l+1=lp go from nbl0(lp)+1 to nbl0(lp)+nbl(lp)
        nbl0(1)=0
        do 420 i = 1,lmaxps
          nbl0(i+1) = nbl0(i) + nbl(i)
  420   continue
c
c       if npf,ptryc are needed (for iptype=2) then read them in
        if (ttype2) then
          read (input,490) npf,ptryc
  490     format(i5,f10.5)
          write (iout,495) npf,ptryc
  495     format(' npf =',i3,5x,' ptryc =',f10.5)
          if (npf.le.4) then
c           this looks more like an lloc value than an npf value
c           ie, user forgot to add the needed input line
            write (iout,*) '***error in subroutine readin'
            write (iout,*) 'npf.le.4, npf =',npf
            write (iout,'(a/a/a)')
     +       ' probably a sign that old-style (pre-7.2.0) input file',
     +       ' has been provided.  as of 7.2.0, if any iptype=2,',
     +       ' add a line to the input file specifying npf and ptryc'
            call exit(1)
          endif
          if ( idim8 .lt. npf ) then
            write(iout,*) '***error subroutine readin'
            write(iout,*) 'idim8 .lt. npf with iptype=2, idim8=',idim8
            call exit(1)
          endif
        endif
c
c       get information about where to get the local potential from
        read(input,500) lloc,keyee,eloc,iploctype
        ikeyeeloc = keyee
  500   format(i5,i5,f10.5,i5)
        if (keyee.gt.0) then
           if (keyee .le. nvales) then
c           set equal to one of the ae eigenvalues
              eloc=ee(ncores+keyee)
           else
              write(iout,*) '***error in subroutine readin'
              write(iout,*) 'keyee .gt. nvales, keyee =',keyee
              call exit(1)
           endif
           if ((.not.ttype2).and.(iploctype.eq.2)) then
              write(iout,*) '***error in subroutine readin'
              write(iout,*) 'npf and ptryc not read in'
              write(iout,*) 'they are needed if iploctype=2'
              call exit(1)
           endif
        endif
c
        write(iout,510) lloc,eloc,iploctype
  510   format(' lloc =',i3,5x,'eloc =',f10.5,5x,'iploctype=',i3)
c
c       get information about how to pseudize the q_ij function
        read(input,520) ifqopt,nqf,qtryc,nfix
  520   format(2i5,f10.5,i5)
c
        write(iout,530) ifqopt,nqf,qtryc,nfix
  530   format(' ifqopt =',i3,5x,'nqf =',i3,5x,'qtryc =',f10.5,
     +    5x,'nfix =',i3)
c
        if ( ifqopt .eq. 0 ) then
          if ( idim8 .lt. 3 ) then
            write(iout,*) '***error subroutine readin'
            write(iout,*) 'idim8 .lt. 3 with ifqopt=0, idim8=',idim8
            call exit(1)
          endif
        elseif ( ifqopt .ge. 1 .and. ifqopt .le. 3 ) then
          if ( idim8 .lt. nqf ) then
            write(iout,*) '***error subroutine readin'
            write(iout,*) 'idim8 .lt. nqf with ifqopt=1-3, idim8=',idim8
            call exit(1)
          endif
        else
          write(iout,*) '***error in subroutine readin'
          write(iout,*) 'ifqopt out of range, ifqopt=',ifqopt
          call exit(1)
        endif
c
        if ( ifqopt .eq. 0 ) then
          if ( nqf .ne. 3 ) then
            write(iout,*) '*** warning subroutine readin'
            write(iout,*) 'iqopt=0 requires nqf=3, but nqf=',nqf
c           write(iout,*) 'resetting nqf'
c           nqf = 3
          endif
        else
          if ( nqf .le. 3 ) then
            write(iout,*) '*** error in subroutine readin'
            write(iout,*) 'ifqopt=1,2 requires nqf .gt. 3 but nqf=',nqf
            call exit(1)
          endif
        endif
c
        if ( ifqopt .lt. 2 ) then
          if ( nfix .ne. 0 ) then
            write(iout,*) '***error in subroutine readin'
            write(iout,*) 'ifqopt =0,1 requires nfix=0 but nfix=',nfix
            call exit(1)
          endif
        else
          if ( nfix .lt. 0 .or. nfix .gt. nbeta ) then
            write(iout,*) '***error in subroutine readin'
            write(iout,*) 'need 0 .le. nfix .le. nbeta but nfix=',nfix
            call exit(1)
          endif
        endif
c
        if ( ifqopt .lt. 3 ) then
          if ( exfact .ge. 1. ) then
            write(iout,*) '***error in subroutine readin'
            write(iout,*) 'need ifqopt=3 for gga, exfact=',exfact
            call exit(1)
          endif
        endif
c
c       read in the extra fixed constraints for diagonal of qfunc
        do 540 ifix = 1,nfix
          read(input,550) ibfix(ifix),rfix(ifix),qfix(ifix)
          write(iout,560) ifix,ibfix(ifix),rfix(ifix),qfix(ifix)
c
  550     format(i5,2f10.5)
  560     format(' ifix',i3,' ibeta',i3,' with rfix,qfix =',2f12.7)
c
          do i=1,nang*2-1
          if ( rfix(ifix) .gt. rinner(i) ) then
            write(iout,*) '***error in subroutine readin'
            write(iout,*) 'rfix must .lt. rinner',rfix(ifix),rinner(i)
            call exit(1)
          endif
          enddo
          if ( ifix .eq. 1 ) then
            if ( ibfix(ifix) .lt. 0 ) then
              write(iout,*) '***error in subroutine readin'
              write(iout,*) 'ibfix .lt. 0 not allowed',ibfix(ifix)
              call exit(1)
            endif
          else
            if ( ibfix(ifix) .le. ibfix(ifix-1) ) then
              write(iout,*) '***error in subroutine readin'
              write(iout,*) 'must have increasing ibfix'
              call exit(1)
            endif
          endif
          if ( ibfix(ifix) .gt. nbeta ) then
            write(iout,*) '***error in subroutine readin'
            write(iout,*) 'require ibfix.le.nbeta, ibfix=',ibfix(ifix)
            call exit(1)
          endif
c
  540   continue
c
c
      endif
c
c     end pseudopotential readin
      return
c
      end
c
c----------------------------------------------------------------------------
c
      subroutine rwps(iread,title,zps,zvps,exftps,nvalps,irel,
     +  z,zv,exfact,nvales,mesh,etot,nnlzps,wwnlps,eeps,
     +  nnlz,wwnl,ee,ncores,keyps,ifpcor,rpcor,rinner,rc,nbeta,
     +  kkbeta,lll,eee,iptype,npf,ptryc,
     +  beta,ddd0,ddd,qqq,qfunc,qfcoef,
     +  rcloc,vloc0,vloc,snl,lloc,eloc,rsatom,rspsco,flname,
     +  ifqopt,nqf,qtryc,nang,nbl0,nbl,r,rab,iver,idmy,
     +  idim1,idim2,idim3,idim4,idim5,idim6,idim8)
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c.....radial mesh information
      dimension r(idim1),rab(idim1),rinner(idim6)
c.....charge densities
      dimension rsatom(idim1),rspsco(idim1)
c.....shell labels, energies, occupancies and real-space cutoff index
      dimension nnlz(idim2),ee(idim2),wwnl(idim2)
c.....pseudo quantum numbers, energies occupancies and cutoff radii
      dimension nnlzps(idim4),eeps(idim4),wwnlps(idim4),rc(idim5)
c.....beta and q functions and q pseudization coefficients
      dimension beta(idim1,idim3),qfunc(idim1,idim3,idim3),
     +qfcoef(idim8,idim6,idim3,idim3),vloc(idim1),vloc0(idim1)
c.....indexing for the beta functions
      dimension nbl0(idim5),nbl(idim5)
c.....angular momenta, reference energy, qij and dij of vanderbilt scheme
      dimension lll(idim3),eee(idim3),iptype(idim3),qqq(idim3,idim3),
     +ddd0(idim3,idim3),ddd(idim3,idim3)
c.....version number and date
      dimension iver(3),idmy(3)
c.....logic for logarithmic mesh types
      logical tlog
c.....wave functions
      dimension snl(idim1,idim2)
c
c.....title
      character*20 title,xctype
c.....file name array
      character*40 flname(6)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
      if ( iread .eq. + 1 ) then
c
c--------------------------------------------------------------------------
c
c       r e a d  p s e u d o p o t e n t i a l  f r o m  f i l e
c
c       save data so compatibility checks can be run later
        mold  = mesh
        r2    = r(2)
        rmesh = r(mesh)
c
        open( unit = iops , file = flname(6) , status = 'old' ,
     +    form = 'unformatted' )
c
        read (iops) (iver(i),i=1,3),(idmy(i),i=1,3)
c
        if ( iver(1) .gt. 99 .or. iver(1) .lt. 1 .or. 
     +       iver(2) .gt. 9  .or. iver(2) .lt. 0 .or.
     +       iver(3) .gt. 9  .or. iver(3) .lt. 0       ) then
c
c         handle errors caused by pre-version number files
c
          write(stderr,*) '***error reading pseudopotential file'
          write(stderr,*) 'file does not contain version number'
          write(stderr,*) 'taking corrective action - will assume:'
c
          iver(1) = 2
          iver(2) = 1
          iver(3) = 0
c
          idmy(1) = 0
          idmy(2) = 0
          idmy(3) = 0
c
          write(stderr,90) (iver(i),i=1,3),idmy(2),idmy(1),idmy(3)
   90     format('pseudopotential program version',i3,'.',
     +      i1,'.',i1,'   date:',i3,' - ',i2,' - ',i4)
c
          rewind (iops)
c
        endif
c
c       set the mesh type
        if ( iver(1) .eq. 1 ) then
          tlog = .false.
        else
          tlog = .true.
        endif
c
c
        read (iops) title,zps,zvps,exftps,nvalps,mesh,etot
        read (iops) (nnlzps(i),wwnlps(i),eeps(i),i=1,nvalps)
        read (iops) keyps,ifpcor,rinner(1)
c
        if ( keyps .ne. 3 ) then
          write(iout,*) '***error in subroutine rwps'
          write(iout,*) 'keyps =',keyps,' out of programed range'
          call exit(1)
        endif
c
        if ( iver(1) .ge. 3 ) then
          read(iops) nang,lloc,eloc,ifqopt,nqf,qtryc
        endif
c
        if (10*iver(1)+iver(2).ge.51) then
          read (iops) (rinner(i),i=1,nang*2-1)
        else
          if (nang.gt.1) then
            do i=2,2*nang-1
              rinner(i)=rinner(1)
            end do
          endif
        endif
c
        if ( iver(1) .ge. 4 ) then
          read(iops) irelps
        endif
c
c
c       set the number of angular momentum terms in q_ij to read in
        if ( iver(1) .eq. 1 ) then
c         no distinction between nang and nvalps
          nang = nvalps
c         no optimisation of q_ij so 3 term taylor series
          nqf = 3
          nlc = 5
        elseif ( iver(1) .eq. 2 ) then
c         no distinction between nang and nvalps
          nang = nvalps
c         no optimisation of q_ij so 3 term taylor series
          nqf = 3
          nlc = 2 * nvalps - 1
        else
          nlc = 2 * nang - 1
        endif
c
c
        read (iops) (rc(i),i=1,nang)
        read (iops) nbeta,kkbeta
        do 100 j=1,nbeta
          read (iops) lll(j),eee(j),(beta(i,j),i=1,kkbeta)
          do 100 k=j,nbeta
            read (iops) ddd0(j,k),ddd(j,k),qqq(j,k),
     +      (qfunc(i,j,k),i=1,kkbeta),
     +      ((qfcoef(i,lp,j,k),i=1,nqf),lp=1,nlc)
            ddd0(k,j)=ddd0(j,k)
            ddd (k,j)=ddd (j,k)
            qqq (k,j)=qqq (j,k)
            do 100 i=1,kkbeta
              qfunc(i,k,j)=qfunc(i,j,k)
  100   continue
        do 110 i=1,nqf
        do 110 lp = 1,nlc
          qfcoef(i,lp,k,j)=qfcoef(i,lp,j,k)
  110   continue
c
        if (10*iver(1)+iver(2).ge.72) then
          read (iops) (iptype(j),j=1,nbeta),npf,ptryc
        endif
c
        read (iops) rcloc,(vloc0(i),i=1,mesh)
c
c       set index arrays nbl and nbl0
        do 150 lp = 1,nang
          nbl(lp) = 0
  150   continue
        do 160 ib = 1,nbeta
          lp = lll(ib) + 1
          nbl(lp) = nbl(lp) + 1
  160   continue
        lmaxps = lll(nbeta)
        nbl0(1) = 0
        do 170 i = 1,lmaxps
          nbl0(i+1) = nbl0(i) + nbl(i)
  170   continue
c
c       possible readin of the frozen core correction
        if (ifpcor.gt.0) then
           read(iops) rpcor
           read (iops) (rspsco(i),i=1,mesh)
        endif
c
        read (iops) (vloc(i),i=1,mesh)
        read (iops) (rsatom(i),i=1,mesh)
c
c       with the logarithmic mesh potentials grid information included
        if ( tlog ) then
          read (iops) (r(i),i=1,mesh)
          read (iops) (rab(i),i=1,mesh)
        endif
        if (iver(1) .ge. 6) then
           nchi = nvales
           if (iver(1) .ge. 7)  read (iops) nchi
           read (iops) ((snl(i,j), i=1,mesh),j=ncores+1,ncores+nchi)
        endif
c
        close (iops)
c
        write (iout,'(/a)') (' pseudopotential has been read in')
c
c       checks for compatibility
c
c       check all no of mesh points in ae and ps cases agree
        if ( mold .ne. mesh ) then
          write(iout,*) '***error in subroutine rwps'
          write(iout,*) 'mesh =',mesh,' but mold =',mold
          call exit(1)
        endif
c
c       chesk all electron mesh and pseudopotential mesh agree
        if ( abs(r2/r(2)-1.0d0) .gt. 1.0d-09 .or.
     +             abs(rmesh/r(mesh)-1.0d0) .gt. 1.0d-09 ) then
          write(iout,*) '***error in subroutine rwps'
          write(iout,*) 'r2,rmesh =',r2,rmesh
          write(iout,*) 'r(2),r(mesh) =',r(2),r(mesh)
          call exit(1)
        endif
c
        if (z.ne.zps.or.zv.ne.zvps.or.exfact.ne.exftps) then
          write(iout,*) '***error in subroutine rwps'
          write(iout,*) '     z  ,zv  ,exfact =',z,zv,exfact
          write(iout,*) '.ne. zps,zvps,exftps =',zps,zvps,exftps
          call exit(1)
        endif
c
c       nvales can be smaller, e.g. psp has s,p,d but testing
c       neutral config so only want to test s,p
        nvalesmin= nvales
        if ( nvales .gt. nvalps ) then
           nvalesmin= nvalps
c          write(iout,*) '***error in subroutine rwps'
c          write(iout,*) 'nvales=',nvales,'.gt. nvalps=',nvalps
c          call exit(1)
        endif
c
c       i=1,nvales must label s, p, d respectively
        do 200 i=1,nvalesmin
          if ( nnlz(ncores+i) .ne. nnlzps(i) ) then
            write(iout,*) '***error in subroutine rwps'
            write(iout,*) 'nnlz',nnlz(ncores+i),'.ne. nnlzps',nnlzps(i)
            call exit(1)
          endif
  200   continue
c
c       check that relativity switches compatible
        if ( iver(1) .lt. 4 ) then
c         potentials generated prior to version 4 are never relativistic
          if ( irel .ne. 0 ) then
            write(iout,*) '***error in subroutine rwps'
            write(iout,*) 'irel=',irel,' .and. iver=',iver(1),' illegal'
            call exit(1)
          endif
        else
c         only proceed if relativistic switches of ae and ps agree
          if ( irelps .ne. irel ) then
            write(iout,*) '***error in subroutine rwps'
            write(iout,*) 'irel =',irel,' .ne. irelps =',irelps
            call exit(1)
          endif
        endif
c
c       read in complete
c
      elseif ( iread .eq. -1 ) then
c
c--------------------------------------------------------------------------
c
c       w r i t e  p s e u d o p o t e n t i a l  t o  f i l e
c
        open( unit = iops , file=flname(6) , status = 'unknown' ,
     +    form = 'unformatted')
c
        write (iops) (iver(i),i=1,3),(idmy(i),i=1,3)
c
c       set the mesh type
        if ( iver(1) .eq. 1 ) then
          tlog = .false.
        else
          tlog = .true.
        endif
c
        write (iops) title,zps,zvps,exftps,nvalps,mesh,etot
        write (iops) (nnlzps(i),wwnlps(i),eeps(i),i=1,nvalps)
        write (iops) keyps,ifpcor,rinner(1)
c
        if ( iver(1) .ge. 3 ) then
          write(iops) nang,lloc,eloc,ifqopt,nqf,qtryc
        endif
c
        if (10*iver(1)+iver(2).ge.51) then
          write(iops) (rinner(i),i=1,nang*2-1)
        endif
c
        if ( iver(1) .ge. 4 ) then
          write(iops) irel
        endif
c
c       set the number of angular momentum terms in q_ij to write out
        if ( iver(1) .eq. 1 ) then
c         no distinction between nang and nvalps
          nang = nvalps
c         no optimisation of q_ij so 3 term taylor series
          nqf = 3
          nlc = 5
        elseif ( iver(1) .eq. 2 ) then
c         no distinction between nang and nvalps
          nang = nvalps
c         no optimisation of q_ij so 3 term taylor series
          nqf = 3
          nlc = 2 * nvalps - 1
        else
          nlc = 2 * nang - 1
        endif
c
        write (iops) (rc(i),i=1,nang)
        write (iops) nbeta,kkbeta
        do 300 j=1,nbeta
          write (iops) lll(j),eee(j),(beta(i,j),i=1,kkbeta)
          do 300 k=j,nbeta
            write (iops) ddd0(j,k),ddd(j,k),qqq(j,k),
     +      (qfunc(i,j,k),i=1,kkbeta),
     +      ((qfcoef(i,lp,j,k),i=1,nqf),lp=1,nlc)
  300   continue
c
        if (10*iver(1)+iver(2).ge.72) then
          write (iops) (iptype(j),j=1,nbeta),npf,ptryc
        endif
c
        write (iops) rcloc,(vloc0(i),i=1,mesh)
        if (ifpcor.gt.0)then
            write (iops) rpcor
           write (iops) (rspsco(i),i=1,mesh)
        endif
        write (iops) (vloc(i),i=1,mesh)
        write (iops) (rsatom(i),i=1,mesh)
c
c       with the logarithmic mesh potentials grid information included
        if ( tlog ) then
          write (iops) (r(i),i=1,mesh)
          write (iops) (rab(i),i=1,mesh)
        endif
        nchi = nvales
        write (iops) nchi
        write (iops) ((snl(i,j), i=1,mesh), j=ncores+1, ncores+nchi)
c
        close (iops)
c
        write (iout,'(/a)') (' pseudopotential has been written out')
c
      else
c
        write(iout,*) '***error in subroutine rwps'
        write(iout,*) 'iread =',iread,' is out of range'
        call exit(1)
c
      endif
c
c-----------------------------------------------------------------------
c
c              p s e u d o p o t e n t i a l  r e p o r t
c
      if (exftps.eq. 0.) xctype = '      ceperley-alder'
      if (exftps.eq.-1.) xctype = '              wigner'
      if (exftps.eq.-2.) xctype = '     hedin-lundqvist'
      if (exftps.eq.-3.) xctype = ' gunnarson-lundqvist'
      if (exftps.eq. 1.) xctype = ' C-A + B88gx + LYPgc'
      if (exftps.eq. 2.) xctype = ' C-A + B88gx        '
      if (exftps.eq. 3.) xctype = ' C-A + B88gx + P86gc'
      if (exftps.eq. 4.) xctype = ' Perdew Wang 1991   '
      if (exftps.eq. 5.) xctype = ' PBE - GGA          '
c
      write (iout,1000) (iver(i),i=1,3),idmy(2),idmy(1),idmy(3)
 1000 format (/4x,60(1h=)/4x,'|  pseudopotential report: version',
     +   i3,'.',i1,'.',i1,' date',i3,'-',i2,'-',i4,2x,'|',
     +   /4x,60(1h-))
      write (iout,1010) title,xctype
 1010 format (4x,'|  ',2a20,' exchange-corr  |')
      write (iout,1020) zps,zvps,exftps
 1020 format (4x,'|  z =',f5.0,4x,'zv =',f5.0,4x,'exfact =',f10.5,
     +   13x,'|')
      write (iout,1030) etot
 1030 format (4x,'|     ',9x,'    ',9x,' etot  =',f10.5,13x,'|')
      write (iout,1040)
 1040 format(4x,'|  index    orbital      occupation    energy',14x,'|')
      write (iout,1050) (i,nnlzps(i),wwnlps(i),eeps(i),i=1,nvalps)
 1050 format(4x,'|',i5,i11,5x,f10.2,f12.2,15x,'|')
      if (10*iver(1)+iver(2).ge.51) then
        write (iout,1061) keyps,ifpcor,rpcor
 1061   format(4x,'|  keyps =',i2,5x,'ifpcor =',i2,5x,'rpcor =',
     +       f10.5,10x,'|')
        do 1064 i = 1,2*nang-1
          write (iout,1065) rinner(i),i
 1064   continue
 1065   format(4x,'|  rinner =',f10.2,5x,'for L=',i5,22x,'|')
      else
        write (iout,1060) keyps,ifpcor,rinner(1)
 1060   format(4x,'|  keyps =',i2,5x,'ifpcor =',i2,5x,'rinner =',
     +     f10.4,9x,'|')
      endif
c
      if (keyps.le.2) then
c
c       standard hsc or with vanderbilt smoothness
c       write (iout,1070)
c1070   format(4x,'|  hsc-type pseudo generation scheme:',22x,'|')
c       write (iout,1080) aa
c1080   format(4x,'|    aa = ',f8.3,41x,'|')
c       write (iout,1090) key
c1090   format(4x,'|    key =',3i8,25x,'|')
c       write (iout,1100) aaa
c1100   format(4x,'|    aaa =',3f8.3,25x,'|'
c    +   /4x,'|         l     rc',41x,'|')
c       write (iout,1110) (i,rc(i),i=1,nvalps)
c1110   format(4x,'|',4x,i6,f8.3,40x,'|')
c
      else if (keyps.eq.3) then
c
c       new scheme
        write (iout,1120)
 1120   format(4x,'|    new generation scheme:',32x,'|')
        write (iout,1130) nbeta,kkbeta,rcloc
 1130   format(4x,'|    nbeta = ',i2,5x,'kkbeta =',i5,5x,
     +                 'rcloc =',f10.4,4x,'|'/
     +         4x,'|    ibeta      l    epsilon   rcut iptype',17x,'|')
c
        do 1140 ib=1,nbeta
          lp=lll(ib)+1
          if (10*iver(1)+iver(2).lt.72) iptype(ib)=-1
c           -1 means iptype was not recorded and is unknown
          write (iout,1150) ib,lll(ib),eee(ib),rc(lp),iptype(ib)
 1140   continue
 1150   format(4x,'|',6x,i2,5x,i2,4x,2f7.2,6x,i1,18x,'|')
c
        ifip2=0
        do ib=1,nbeta
          if (iptype(ib).eq.2) ifip2=1
        end do
        if (ifip2.eq.1) then
          write (iout,1155) npf,ptryc
 1155     format(4x,'|  npf    =',i2,'  ptryc =',f8.3,29x,'|')
        endif
c
        if ( iver(1) .gt. 2 ) then
          write(iout,1160) lloc,eloc
 1160     format(4x,'|  lloc   =',i2,'  eloc   =',f8.3,28x,'|')
          write(iout,1170) ifqopt,nqf,qtryc
 1170     format(4x,'|  ifqopt =',i2,'  nqf    =',i2,'  qtryc =',f8.3,
     +    17x,'|')
        endif
c
        if ( iver(1) .gt. 3 .and. irel .eq. 2 ) then
          write(iout,1180) 'koelling-harmon equation'
        else
          write(iout,1180) 'schroedinger equation   '
        endif
 1180   format(4x,'|',2x,'all electron calculation used ',a24,2x,'|')
c
        if ( .not. tlog ) then
          write(iout,2000) 'herman skillman mesh'
        else
          write(iout,2000) '**logarithmic mesh**'
        endif
c
 2000   format (4x,'|',9x,10('*'),a20,10('*'),9x,'|')
c
      else if (keyps.eq.4) then
c
c       frozen core
c       write (iout,2500)
c2500   format(4x,'|    frozen core all-electron case',25x,'|')
c
      endif
c
      write (iout,3000)
 3000 format (4x,60(1h=))
c
      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine eigprt(delta,nnlz,ee,wwnl,nkkk,ncspvs,ehar,emvxc,a,b,
     +  r,ruae,rsatom,rscore,rudif,mesh,ifprt,idim1,idim2)
c
c     routine computes the total enegy of the atom and prints out
c     results. Charge densities and potentials also shown if ifprt.ge.2
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c     parameter to determine the size of the plot out arrays
      parameter ( nind = 55 )
c
c.....logarithmic radial mesh information
      dimension r(idim1)
c.....potentials for the all electron calculation
      dimension ruae(idim1)
c.....charge densities
      dimension rscore(idim1),rsatom(idim1)
c.....shell labels, energies, occupancies
      dimension nnlz(idim2),ee(idim2),wwnl(idim2),nkkk(idim2)
c.....scratch arrays
      dimension rudif(idim1)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c     compute the total energy of the atom
c
c     compute the eigen value sum
      etot=0.0d0
      do 100 i=1,ncspvs
        etot=etot+wwnl(i)*ee(i)
  100 continue
      ebsr=etot
c     correct for the over counting hartree and exchange
c     and correlation contributions
      etot=etot-ehar+emvxc
c
c     print out of the results
c
      write (iout,200) delta
  200 format(//,t26,29hthe system has converged with,/,t34,6hdelta=,
     +  d14.7,///,t30,19hthe eigenvalues are)
      write (iout,210)
      write (iout,220) (nnlz(i),wwnl(i),ee(i),i=1,ncspvs)
  210 format(//,10x,7horbital,10x,10hoccupation,10x,6henergy)
  220 format(10x,i4,16x,f7.2,11x,f11.5)
      etev=etot*13.605
      write (iout,230) etot,etev,ebsr,ehar,emvxc
  230 format (6h etot=,f13.5,7h ryd ; ,3x,f13.5,3h ev,12x,
     +   16hebsr,ehar,emvxc=,3f13.5)
c
c     possible further print outs of densities and potentials
      if ( ifprt .ge. 2 ) then
c
c       print out a header
        write(iout,240)
c
c       find maximum value in nkkk
        imax = 0
        do 300 i = 1,ncspvs
          imax = max ( imax , nkkk(i) )
  300   continue
c
c       set the step size
        rinc = r(imax) / dble( nind + 1 )
        rcur = rinc / 2.0d0
c
        do 250 j = 1,nind
          i = int( dlog( rcur / a + 1.0d0 ) / b ) + 1
          write (iout,260) i,r(i),ruae(i),rsatom(i),rscore(i),rudif(i)
          rcur = rcur + rinc
  250   continue
c
      endif
c
  240 format(///,t14,1hr,t31,3hr*v,t44,18h4pi*r*r*rho(total),
     +t66,17h4pi*r*r*rho(core),t88,21hdelta(difference r*v))
  260 format(i5,t8,f14.7,t26,d14.7,t44,d14.7,t66,d14.7,t88,d14.7)
c
      return
      end
c
c----------------------------------------------------------------------------
c     ***************************************************************
      subroutine prmat(a,n,c,idim3)
c     ***************************************************************
c ---------------------------------------------------------------
c
c     routine for printing out a square matrix
c
      implicit double precision (a-h,o-z)
c ---------------------------------------------------------------
c
      character c
      dimension a(idim3,idim3)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
      write (iout,10) c
   10 format (/4x,'matrix ',a1,'(i,j)'/)
      do 40 j=1,n
        write (iout,20) (a(i,j),i=1,n)
   20   format (7f14.5)
   40 continue
c
      return
      end
c
c----------------------------------------------------------------------------
c     ***************************************************************
      subroutine prwf(ifprt,ncores,nvales,ncspvs,nnlz,ee,nkkk,
     +  r,rab,a,b,rsvale,snl,vpot,ipass,flname,ifplw,mesh,idim1,idim2)
c     ***************************************************************
c
c     routine for printing comparisons between pseudo and
c     all electron energies, wavefunctions and densities
c
c     modified 9/9/91 to print out wavefunctions in order s,p,d
c
c------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c     parameter to determine the size of the plot out arrays
      parameter ( nind = 55 )
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1),vpot(idim1)
c.....charge densities
      dimension rsvale(idim1)
c.....wavefunctions
      dimension snl(idim1,idim2)
c.....shell labels, energies, occupancies and real-space cutoff index
      dimension nnlz(idim2),ee(idim2),nkkk(idim2)
c.....file names
      character*40 flname(6)
c.....print out arrays
      dimension index(nind),plot(nind,14),eeprt(12),rhostr(nind,2)
c.....symbols for the line printer plot out
      dimension isymbl(13)
c.....character arrays for print outs
      character*1 spd(4),orblet(12)
c
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c     initialise the symbol array
      data isymbl /1h!,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9,1hA,1hB,1hC/
c
c     intialise the strings in spd
      data spd / 's' , 'p' , 'd', 'f' /
c
      save index, orblet, eeprt, plot, rhostr
c-----------------------------------------------------------------------
c
c            i n i t i a l i s a t i o n 
c
      if ( ipass .lt. 1 .or. ipass .gt. 2 ) then
        write(iout,*) '***error in subroutine prwf'
        write(iout,*) 'ipass =',ipass,' is out of range'
        call exit(1)
      endif
c
      if ( nvales .gt. 6 ) then
        write(iout,*) '***error in subroutine prwf'
        write(iout,*) 'nvales =',nvales,' will over run local arrays'
        write(stderr,*) '***premature exit of prwf - nvales too large'
        return
      endif
c
c     routine only plots out valence wavefunctions
c     so return if nvales is equal to zero
      if ( nvales .eq. 0 ) return
c
c     possible print outs to wavefunction plot file
      if ( ifplw .eq. 1 ) then
        write(iout,'(/a)')
     1    ' subroutine prwf: writing wavefunctions to unit iplot'
c
        if ( ipass .eq. 1 ) then
c
c         first call to prwf so open plot out file
          open( unit = iplot , file = flname(4) , status = 'unknown' , 
     +      form = 'formatted' )
c
        endif
c
        do iorb = ncores+1,ncspvs
           do ir=1,mesh
              write(iplot,*) r(ir),snl(ir,iorb)
           enddo
           write(iplot,*) '&'
        enddo
c
      endif

      if ( ifplw .eq. 2 ) then
        if ( ipass .eq. 1 ) then
c
        write(iout,*) 'writing out wavefunctions for refstate input'
c
c         first call to prwf so open plot out file
          open( unit = iplot , file = flname(4) , status = 'unknown' , 
     +      form = 'formatted' )
c
          
          write(iplot,*) mesh,nvales
          do  iorb = ncores+1,ncspvs
             write(iplot,*) (snl(ir,iorb),ir=1,mesh)
          enddo
          write(iplot,*) (vpot(ir),ir=1,mesh)
c
        endif
c
c
      endif
c
      if ( ipass .eq. 1 ) then
c
        do 90 i = 1,12
          orblet(i) = ' '
   90   continue
c
c       zero the plot out array
        do 100 j = 1,14
        do 100 i = 1,nind
          plot(i,j) = 0.0d0
  100   continue
        do 110 i = 1,12
          eeprt(i) = 0.0d0
  110   continue
c
c       find the points for printing out
        nkmax = 0
        do 120 i = ncores+1,ncspvs
          nkmax = max(nkmax,nkkk(i))
  120   continue
c
c       in practice much of the tail is uninteresting so divide by 3
        rinc = r(nkmax) / dble( 3 * ( nind + 1 ) )
        rcur = rinc / 2.0d0
c
        do 130 i = 1,nind
          index(i) = int( dlog( rcur / a + 1.0d0 ) / b ) + 1
          rcur = rcur + rinc
  130   continue
c
        do 140 i = 1,nind
          plot(i,1) = r(index(i))
  140   continue
c
      endif
c
c-----------------------------------------------------------------------
c
c          l o a d  t h e  p r i n t  o u t  a r r a y s
c
c     first set up the offset parameters
      if ( ipass .eq. 1 ) then
        ioff = 2
      elseif ( ipass .eq. 2 ) then
        ioff = 8
      endif
c
c     now load the energy
      do 200 i = 1,nvales
        iorb = ncores + i
        eeprt( 2*(i-1) + ipass ) = ee(iorb)
  200 continue
c
c     load the wavfunctions
      do 210 j = 1,nvales
c
c       compute orbital number of current state
        iorb = ncores + j
c       compute angular momentum of current state
        lcur = nnlz(iorb) / 10 - 10 * ( nnlz(iorb) / 100 )
        lp = lcur + 1
c
c       ensure maximum value of snl is +ve
        snlmax = 0.0d0
        do 220 i = 1,nind
          if ( abs(snl(index(i),iorb)) .gt. snlmax ) then
            snlmax = abs(snl(index(i),iorb))
            imax = i
          endif
  220   continue
c
        if ( snlmax .gt. 1.0d-10 ) then
          sig = snl(index(imax),iorb) / abs(snl(index(imax),iorb))
        else
          sig = 0.0d0
        endif
c
        do 230 i = 1,nind
          plot(i,ioff+j) = sig * snl(index(i),iorb)
  230   continue
c
c       set the orbital letter
        if ( lp .gt. 4 ) then
          write(iout,*) '***error in subroutine prwf'
          write(iout,*) 'lp .gt. 4 too large, lp =',lp
          call exit(1)
        endif
c
        orblet(j+nvales*(ipass-1)) = spd(lp)
c
  210 continue
c
c     load the valence charge density
      do 240 i = 1,nind
        rhostr(i,ipass) = rsvale(index(i))
  240 continue
c
c---------------------------------------------------------------------
c
c               p r i n t  o u t  o f  r e s u l t s
c
      if ( ipass .eq. 2 ) then
c
        write (iout,400) (nnlz(mp),mp=ncores+1,ncspvs)
  400   format (/,'   all electron and pseudo energies',//,12x,
     +    4(11x,i3,11x))
        write (iout,410) ('   full     pseudo  ',i=1,nvales)
  410   format (/15x,4(2x,a20,3x)/)
        write (iout,420) 'energy',(eeprt(i),i=1,2*nvales)
  420   format (4x,a6,2x,4(2x,f10.7,1x,f10.7,2x))
c
        if ( ifprt .ge. 1 ) then
          write (iout,430) (orblet(i),i=1,2*nvales)
  430     format(/,'   valence all-electron and pseudo wavefunctions',
     +      //,15x,8(5x,a1,6x))
          write (iout,440) (' full ',i=1,nvales),('pseudo',i=1,nvales)
  440     format(8x,'r',6x,8(3x,a6,3x))
          do 445 in = 1,nind
            write (iout,450) plot(in,1),(plot(in,2+i),i=1,nvales),
     +        (plot(in,8+j),j=1,nvales)
  445     continue
  450     format(9(2x,f10.5))
        endif
c
      endif
c
c     print and line-printer plot of wavefunctions
      if ( ifprt .ge. 1 .and. ipass .eq. 2 ) then
        write (iout,'(/3x,a)')
     +    'chart: valence all-electron and pseudo wavefunctions'
        call lpplot(2,plot,nind,14,nind,0,nind,14,isymbl,
     +    -0.6,1.4)
      endif
c
c     print and line-printer plot of charge densities
      if ( ifprt .ge. 2 .and. ipass .eq. 2 ) then
  500   format(/,'   valence pseudo and all electron charge densities',
     +    //,8x,'r',9x,' full ',6x,'pseudo')
        write(iout,510) (plot(in,1),rhostr(in,1),rhostr(in,2),
     +    in=1,nind)
  510   format(3(2x,f10.5))
c
      endif
c
      return
      end
c
c-----------------------------------------------------------------------
c     ***************************************************************
      subroutine prbeta(ifprt,name,beta,kkbeta,nbeta,lll,r,a,b,
     +   idim1,idim3)
c     ***************************************************************
c
c     routine for producing line printer plotout of beta (or chi)
c
c------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c     parameter to determine the size of the plot out arrays
      parameter ( nind = 55 )
c
c.....logarithmic radial mesh information
      dimension r(idim1)
c.....beta or chi function
      dimension beta(idim1,idim3),lll(idim3)
c.....print out arrays
      dimension index(nind),plot(nind,12)
c.....symbols for the line printer plot out
      dimension isymbl(13)
c.....beta or chi ?
      character*4 name
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c     initialise the symbol array
      data isymbl /1h!,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9,1hA,1hB,1hC/
c
c-----------------------------------------------------------------------
c
c            i n i t i a l i s a t i o n 
c
c
      if ( nbeta .gt. 10) then
        write(iout,*) '***error in subroutine prqf'
        write(iout,*) 'nbeta =',nbeta,' will over run local arrays'
        call exit(1)
      endif
c
      if ( ifprt .ge. 3) write (iout,'(/3x,a)') 'subroutine prbeta'
c
c     zero the plot out array
      do 100 j = 1,12
      do 100 i = 1,nind
        plot(i,j) = 0.0d0
  100 continue
c
      rinc = r(kkbeta) / dble( nind - 1 )
      rcur = rinc * 0.5d0
c
      do 130 i = 1,nind
        index(i) = int( dlog( rcur / a + 1.0d0 ) / b ) + 1
        rcur = rcur + rinc
  130 continue
c
      do 140 i = 1,nind
        plot(i,1) = r(index(i))
  140 continue
c
c-----------------------------------------------------------------------
c
c          l o a d  t h e  p r i n t  o u t  a r r a y s
c
c
c     load the beta functions into plot
      do 210 j = 1,nbeta
c
        do 230 i = 1,nind
          plot(i,2+j) = beta(index(i),j)
  230   continue
c
  210 continue
c
c
c---------------------------------------------------------------------
c
c     simple tabular listing of beta or chi functions
      if ( ifprt .ge.  3 ) then
        write(iout,'(/1x,a4,10i11)') name,(lll(ib),ib=1,nbeta)
        do i = 1,nind
          write (iout,'(f8.3,10f11.5)') plot(i,1),
     +       (plot(i,2+ib),ib=1,nbeta)
        end do
      endif
c
c     line printer plot out of beta or chi functions
      if ( ifprt .ge.  3 ) then
        write(iout,'(/1x,a4,a)') name,' functions in real-space'
        call lpplot(0,plot,nind,nbeta+2,nind,0,nind,12,isymbl,0.0,0.0)
      endif
c
      return
      end

c-----------------------------------------------------------------------
c     ***************************************************************
      subroutine prqf(ifprt,qfunc,kkbeta,nbeta,r,a,b,idim1,idim3)
c     ***************************************************************
c
c     routine for producing line printer plotout of qfunc
c
c------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c     parameter to determine the size of the plot out arrays
      parameter ( nind = 55 )
c
c.....logarithmic radial mesh information
      dimension r(idim1)
c.....qfunction
      dimension qfunc(idim1,idim3,idim3)
c.....print out arrays
      dimension index(nind),plot(nind,10)
c.....symbols for the line printer plot out
      dimension isymbl(13)
c
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c     initialise the symbol array
      data isymbl /1h!,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9,1hA,1hB,1hC/
c
c-----------------------------------------------------------------------
c
c            i n i t i a l i s a t i o n 
c
c
      if ( nbeta .gt. 8 ) then
        write(iout,*) '***error in subroutine prqf'
        write(iout,*) 'nbeta =',nbeta,' will over run local arrays'
        call exit(1)
      endif
c
c
c     zero the plot out array
      do 100 j = 1,10
      do 100 i = 1,nind
        plot(i,j) = 0.0d0
  100 continue
c
      rinc = r(kkbeta) / dble( nind - 1 )
      rcur = rinc * 0.5d0
c
      do 130 i = 1,nind
        index(i) = int( dlog( rcur / a + 1.0d0 ) / b ) + 1
        rcur = rcur + rinc
  130 continue
c
      do 140 i = 1,nind
        plot(i,1) = r(index(i))
  140 continue
c
c-----------------------------------------------------------------------
c
c          l o a d  t h e  p r i n t  o u t  a r r a y s
c
c
c     load the qfunctions into plot
      do 210 j = 1,nbeta
c
        do 230 i = 1,nind
          plot(i,2+j) = qfunc(index(i),j,j)
  230   continue
c
  210 continue
c
c
c---------------------------------------------------------------------
c
c               p r i n t  o u t  o f  r e s u l t s
c
c
c     line printer plot out of real-space qfunctions
      if ( ifprt .ge. 1 ) then
        write (iout,'(3x,a)') 'chart: qfunction Q_{ib,ib}(r)'
        call lpplot(0,plot,nind,nbeta+2,nind,0,nind,10,isymbl,0.0,0.0)
      endif
c
      return
      end

c------------------------------------------------------------------------
      subroutine prcor(ifprt,rspsco,rscore,rsvale,kkbeta,r,a,b,idim1)
c     ***************************************************************
c
c     routine for producing line printer plotout of core charge
c
c------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c     parameter to determine the size of the plot out arrays
      parameter ( nind = 55 )
c
c.....logarithmic radial mesh information
      dimension r(idim1)
c.....qfunction
      dimension rspsco(idim1),rscore(idim1),rsvale(idim1)
c.....print out arrays
      dimension index(nind),plot(nind,5)
c.....symbols for the line printer plot out
      dimension isymbl(13)
c
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c     initialise the symbol array
      data isymbl /1h!,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9,1hA,1hB,1hC/
c
c-----------------------------------------------------------------------
c
c            i n i t i a l i s a t i o n 
c
c
c
c
c     zero the plot out array
      do 100 j = 1,5
      do 100 i = 1,nind
        plot(i,j) = 0.0d0
  100 continue
c
      rinc = r(kkbeta) / dble( nind - 1 )
      rcur = rinc * 0.5d0
c
      do 130 i = 1,nind
        index(i) = int( dlog( rcur / a + 1.0d0 ) / b ) + 1
        rcur = rcur + rinc
  130 continue
c
      do 140 i = 1,nind
        plot(i,1) = r(index(i))
  140 continue
c
c-----------------------------------------------------------------------
c
c          l o a d  t h e  p r i n t  o u t  a r r a y s
c
c
c     load the qfunctions into plot
c
      do 230 i = 1,nind
          plot(i,3) = rscore(index(i))
          plot(i,4) = rspsco(index(i))
          plot(i,5) = rsvale(index(i))
  230   continue
c
c
c
c---------------------------------------------------------------------
c
c               p r i n t  o u t  o f  r e s u l t s
c
c
c     line printer plot out of real-space qfunctions
      if ( ifprt .gt. -3 ) then
        write(iout,500)
  500   format(//,' core charge densities and valence density',
     +  //,8x,'r',9x,' full core',6x,'pseudo core',6x,'full valence')
        do 200 ig = 1,nind
          write(iout,9998) plot(ig,1),(plot(ig,2+i),i=1,3)
  200   continue
 9998   format(1x,f10.5,3f16.5)
c
        write(iout,9997) 
 9997   format(/,' real space core charge and valence densities',/)
c
        call lpplot(0,plot,nind,5,nind,0,nind,5,isymbl,0.0,0.0)
      endif
c
      return
      end

c
c-----------------------------------------------------------------------
c     ***************************************************************
      subroutine prqf1(ifprt,qfunc,kkbeta,nltotp,r,a,b,idim1,idim3)
c     ***************************************************************
c
c     routine for producing line printer plotout of qfunc
c
c------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c     parameter to determine the size of the plot out arrays
      parameter ( nind = 55 ,idm6=7)
c
c.....logarithmic radial mesh information
      dimension r(idim1)
c.....qfunction
      dimension qfunc(idim1,idm6)
c.....print out arrays
      dimension index(nind),plot(nind,10)
c.....symbols for the line printer plot out
      dimension isymbl(13)
c
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c     initialise the symbol array
      data isymbl /1h!,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9,1hA,1hB,1hC/
c
c-----------------------------------------------------------------------
c
c            c h e c k   i f p r t   f o r   r e t u r n
c
      if ( ifprt .lt. 4 ) return
c
c-----------------------------------------------------------------------
c
c            i n i t i a l i s a t i o n 
c
c
      if ( nltotp .gt. 6 ) then
        write(iout,*) '***error in subroutine prqf1'
        write(iout,*) 'nltotp =',nltotp,' will over run local arrays'
        call exit(1)
      endif
c
c
c     zero the plot out array
      do 100 j = 1,10
      do 100 i = 1,nind
        plot(i,j) = 0.0d0
  100 continue
c
      rinc = r(kkbeta) / dble( nind - 1 )
      rcur = rinc * 0.5d0
c
      do 130 i = 1,nind
        index(i) = int( dlog( rcur / a + 1.0d0 ) / b ) + 1
        rcur = rcur + rinc
  130 continue
c
      do 140 i = 1,nind
        plot(i,1) = r(index(i))
  140 continue
c
c-----------------------------------------------------------------------
c
c          l o a d  t h e  p r i n t  o u t  a r r a y s
c
c
c     load the qfunctions into plot
      do 210 j = 1,nltotp
c
        do 230 i = 1,nind
          plot(i,2+j) = qfunc(index(i),j)
  230   continue
c
  210 continue
c
c
c---------------------------------------------------------------------
c
c               p r i n t  o u t  o f  r e s u l t s
c
c
c     line printer plot out of real-space qfunctions
      if ( ifprt .gt. -3 ) then
        write(iout,*) 'l-dep. qfunction in real-space (diagonal only)'
        call lpplot(0,plot,nind,nltotp+2,nind,0,nind,10,isymbl,0.0,0.0)
      endif
c
      return
      end
c
c-----------------------------------------------------------------------
c     ***************************************************************
      subroutine prvloc(ifprt,vloc,rcloc,r,a,b,idim1)
c     ***************************************************************
c
c     routine for producing line printer plotout of vloc
c
c------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c     parameter to determine the size of the plot out arrays
      parameter ( nind = 55 )
c
c.....logarithmic radial mesh information
      dimension r(idim1)
c.....local potential
      dimension vloc(idim1)
c.....print out arrays
      dimension index(nind),plot(nind,3)
c.....symbols for the line printer plot out
      dimension isymbl(13)
c
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c     initialise the symbol array
      data isymbl /1h!,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9,1hA,1hB,1hC/
c
c-----------------------------------------------------------------------
c
c            i n i t i a l i s a t i o n 
c
c
c
c     zero the plot out array
      do 100 j = 1,3
      do 100 i = 1,nind
        plot(i,j) = 0.0d0
  100 continue
c
      rinc = 3.0d0 * rcloc / dble( nind - 1 )
      rcur = rinc * 0.5d0
c
      do 130 i = 1,nind
        index(i) = int( dlog( rcur / a + 1.0d0 ) / b ) + 1
        rcur = rcur + rinc
  130 continue
c
      do 140 i = 1,nind
        plot(i,1) = r(index(i))
  140 continue
c
c-----------------------------------------------------------------------
c
c          l o a d  t h e  p r i n t  o u t  a r r a y s
c
c
c
      do 200 i = 1,nind
        plot(i,3) = vloc(index(i))
  200 continue
c
c
c
c---------------------------------------------------------------------
c
c               p r i n t  o u t  o f  r e s u l t s
c
c
c     line printer plot out of real-space qfunctions
      if ( ifprt .ge. 1 ) then
        write (iout,'(/2x,a)') 'chart: pseudized local potential'
        call lpplot(0,plot,nind,3,nind,0,nind,3,isymbl,0.0,0.0)
      endif
c
      return
      end
c
c----------------------------------------------------------------------
c
c     ****************************************************************
      subroutine lpplot(nplot, array, nvals, nvars, nlines, norder,
     1  idim, jdim,isymbl,yminc,ymaxc)
c     ****************************************************************
c
c              l i n e  p r i n t e r  p l o t  o u t s
c
c     nplot - chart number for current plot
c     array - data to be plotted
c             array(*,1) - base axis scale
c             array(*,3 --> nvars ) - data for plotting
c     nvals - number of base variables to be plotted
c     nvars - no of data sets for plotting + 2
c     nlines ??
c     norder - possible sorting of data (norder .le. 0 -- > no sort)
c     idim1 - first dimension of array
c     idim2 - second dimension of array
c     isymbol - symbols for plot out
c     yminc - minimum value on cross varibale axis
c     ymaxc - maximum value on cross varibale axis
c
c----------------------------------------------------------------------
c
      integer blank, isymbl, out, i, j, k, nll, np1, line, nsize
      double precision array
      double precision yprint, f, ymax, ymin, xbase, xpmax,
     1   xprint, xscale, yscale, yscalf
      dimension out(101), yprint(11), array(idim, jdim)
c
      dimension isymbl(nvars-1)
c     note: isymbl must be passed with dimension nvars-1 or larger
c           to avoid reading past end of array
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
      data blank/1h /
c
      nstart = 1
c
c     sort base variable into ascending order
      if (norder .gt. 0) then
        do i = 1, nvals
        do j = i, nvals
          if (array(i, 1) .gt. array(j, 1)) then
            do k = 1, nvars
              f = array(i, k)
              array(i, k) = array(j, k)
              array(j, k) = f
            end do
          endif
        end do
        end do
      endif
c
c     find scale for base variable.
      nll = nlines
      xbase = array(nstart, 1)
      xscale = (array(nvals, 1) - xbase) / (nll - 1)
c
      if (nll .le. 0) nll = 50
c     find scale factor for cross-variables.
      ymin = array(nstart, 2)
      ymax = ymin
      do j = 2, nvars
         do i = nstart, nvals
            f = array(i, j)
            if (f .gt. ymax) ymax = f
            if (f .lt. ymin) ymin = f
         end do
      end do
      if(yminc .ne. 0.) ymin=yminc
      if(ymaxc .ne. 0.) ymax=ymaxc
      yscalf = 100.0e0 / (ymax - ymin)
c
c  find independent (base) variable print position.
      i = nstart
      np1 = nvars + 1
      do 80 line = nstart, nll
         xprint = xbase + (line - 1) * xscale
         nsize = 2
         if (xprint .ge. (- 99999.e0) .and. xprint .le. 999999.e0)
     1  nsize = 1
         xpmax = xprint + 0.5e0 * xscale
         if (array(i, 1) .gt. xpmax) go to 70
c
c  initialize print line to blanks.
         do jout = 1, 101
            out(jout) = blank
         end do
c
c  find dependent (cross) variables.
   57    do 60 j = 2, nvars
            jj = np1 - j
            jout = (array(i, jj + 1) - ymin) * yscalf + 1.5e0
      if(jout .le. 0 .or. jout .ge. 102) go to 60
            out(jout) = isymbl(jj)
   60    continue
c
c  should another set of values be plotted on this same line?
         if (i .ge. nvals) go to 61
         i = i + 1
         if (array(i, 1) .le. xpmax) go to 57
c
c  print line and clear, or skip.
   61    if (nsize .eq. 2) go to 63
         write (iout, 1002) xprint, (out(jout), jout = 1, 101)
         go to 80
   63    write (iout, 1004) xprint, (out(jout), jout = 1, 101)
         go to 80
   70    write (iout, 1003)
   80 continue
c
c  print marks for independent (cross) variables
      write (iout, 1007)
      yscale = 10.0e0 / yscalf
      do 90 i = 1, 10
   90    yprint(i) = ymin + (i - 1) * yscale
      yprint(11) = ymax
      write (iout, 1008) (yprint(i), i = 1, 11, 2)
      write (iout, 1009) (yprint(i), i = 2, 10, 2)
      write (iout,*)
      return
c
 1001 format(1h1, 60x, 6h chart, i4//)
 1002 format(1h , f11.4, 5x, 101a1)
 1003 format(1h )
 1004 format(1h , 1pe15.8, 1x, 101a1)
 1007 format(1h , 16x, 1hi,10(10h         i))
 1008 format(1h , 1x, 6(8x, g12.4))
 1009 format(1h , 11x, 5(8x, g12.4))
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine fanal(r,rab,nnlz,wwnl,nkkk,snl,ncores,
     +  ncspvs,rho,mesh,ifprt,idim1,idim2)
c
c     routine for fourier analysing pseudo wavefunctions
c
      implicit double precision(a-h,o-z)
c
      parameter( nind = 76 )
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1)
c.....wavefunctions
      dimension snl(idim1,idim2)
c.....shell labels, energies, occupancies and real space cutoff
      dimension nnlz(idim2),wwnl(idim2),nkkk(idim2)
c.....plot out arrays
      dimension plot(nind,8),isymbl(13)
c.....scratch
      dimension rho(idim1)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c     set the maximum wavevector for fourier analysis
      data gmax / 15.0d0 /
c
c     initialise the symbol array
      data isymbl /1h!,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9,1hA,1hB,1hC/
c
c-----------------------------------------------------------------------
c
c                i n i t i a l i z a t i o n
c
      if ( ifprt .ge. 1) write (iout,'(/1x,a)')
     +   'subroutine fanal: fourier transform pseudo wavefunctions'
c
c     check that the dimensions of plot are sufficient
      nvales = ncspvs - ncores
      if ( nvales .gt. 6 ) then
        write(iout,*) '***error in subroutine fanal'
        write(iout,*) 'nvales =',nvales,' .gt. 6'
        call exit(1)
      endif
c
c     compute the wavevector increment
      ginc = gmax / dble( nind - 1 )
c
      gcur = - ginc
      do 10 ig = 1,nind
        gcur = gcur + ginc
        plot(ig,1) = gcur
        plot(ig,2) = 0.0d0
   10 continue
c
c-----------------------------------------------------------------------
c
c               f o u r i e r  a n a l y s i s
c
      do 100 iorb = ncores+1,ncspvs
c
c       set the index for the plot out array
        index = 2 + iorb - ncores
c
c       compute angular momentum of this state
        ncur = nnlz(iorb) / 100
        lcur = nnlz(iorb) / 10 - 10 * ncur
c
c       check whether to skip fourier analysis on this state
        if ( wwnl(iorb) .lt. 0.0d0 ) then
          do 110 ig = 1,nind
            plot(ig,index) = 0.0d0
  110     continue
          goto 100
        endif
c
c       loop over wavevectors
        gcur = - ginc
c
        do 120 ig = 1,nind
c
c         compute the current wavevector
          gcur = gcur + ginc
c
c         remeber in performing transform that snl = r * psi
          do 130 ir = 1,nkkk(iorb)
            rho(ir) = r(ir) * snl(ir,iorb) * bessel(gcur*r(ir),lcur)
  130     continue
c
          do 140 ir = nkkk(iorb)+1,mesh
            rho(ir) = 0.0d0
  140     continue
c
          call radin(mesh,rho,rab,asum,idim1)
          plot(ig,index) = gcur * asum
c
  120   continue
c
c     next band
  100 continue
c
c     print out and line printer plot out
      if ( ifprt .ge. 1 ) then
c
        write (iout,'(/3x,a)') 'q * fourier transform of wfns'
        write (iout,'(/2x,a,13x,i3,10x,i3,10x,i3,10x,i3/)')
     +    ' wavevector',(nnlz(iorb),iorb=ncores+1,ncspvs)
c
        do ig = 1,nind
          write(iout,9998) plot(ig,1),(plot(ig,2+i),i=1,nvales)
        end do
 9998   format(4x,f9.5,7x,f9.5,4x,f9.5,4x,f9.5,4x,f9.5,4x,f9.5,4x,f9.5)
c
        write (iout,'(/2x,a)') 'chart: q * fourier transform of wfns'
c
        call lpplot(0,plot,nind,2+nvales,nind,0,nind,8,
     +    isymbl,0.0,0.0)
c
      endif
c
      if ( ifprt .ge. 1) write (iout,'(/1x,a)') 'leaving fanal'
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine fanalb(ifprt,name,beta,kkbeta,nbeta,lll,r,rab,rho,
     +   idim1,idim3)
c
c     routine for fourier analysing beta or chi functions
c
      implicit double precision(a-h,o-z)
c
      parameter( nind = 101 )
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1)
c.....scratch
      dimension rho(idim1)
c.....beta functions and their angular momentum labels
      dimension beta(idim1,idim3),lll(idim3)
c.....plot out arrays
      dimension plot(nind,12),isymbl(13)
c.....beta or chi ?
      character*4 name
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c     set the maximum wavevector for fourier analysis
      data gmax / 20.0d0 /
c
c     initialise the symbol array
      data isymbl /1h!,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9,1hA,1hB,1hC/
c
c-----------------------------------------------------------------------
c
c                i n i t i a l i z a t i o n
c
      if ( ifprt .ge. 3) write (iout,'(/3x,a)') 'subroutine fanalb'
c
c     check that the dimensions of plot are sufficient
      if ( nbeta .gt. 10 ) then
        write(iout,*) '***error in subroutine fanalb'
        write(iout,*) 'nbeta =',nbeta,' .gt. 10'
        call exit(1)
      endif
c
c     compute the wavevector increment
      ginc = gmax / dble( nind - 1 )
c
      gcur = - ginc
      do ig = 1,nind
        gcur = gcur + ginc
        plot(ig,1) = gcur
        plot(ig,2) = 0.0d0
      end do
c
c-----------------------------------------------------------------------
c
c               f o u r i e r  a n a l y s i s
c
      do j=1,nbeta
c
c       set the index for the plot out array
        index = 2 + j
c
c       compute angular momentum of this state
        lcur = lll(j)
c
c       loop over wavevectors
        gcur = - ginc
c
        do ig = 1,nind
c
c         compute the current wavevector
          gcur = gcur + ginc
c
c         remeber in performing transform that snl = r * psi
          do ir = 1,kkbeta
            rho(ir) = r(ir) * beta(ir,j) * bessel(gcur*r(ir),lcur)
          end do
c
          call radin(kkbeta,rho,rab,asum,idim1)
          plot(ig,index) = asum
c
        end do
c
      end do
c
c     simple tabular listing of ft of beta or chi functions
      if ( ifprt .ge.  3 ) then
        write(iout,9999) name,(lll(j),j=1,nbeta)
 9999   format(/,' fourier transform of ',a4,' functions',//,
     +    ' wavevector',10(10x,i3))
        do ig = 1,nind
          write(iout,9998) plot(ig,1),(plot(ig,2+j),j=1,nbeta)
        end do
 9998   format(1x,f9.5,7x,6(f9.5,4x))
      endif
c
c     line printer plot out of ft of beta or chi functions
      if ( ifprt .ge.  3 ) then
        write(iout,9997)  name
 9997   format(/,' fourier transform of ',a4,' functions')
        call lpplot(0,plot,nind,2+nbeta,nind,0,nind,12,
     +    isymbl,0.0,0.0)
      endif
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine fanalq(r,rab,qfunc,nbeta,rho,mesh,ifprt,
     +  idim1,idim3)
c
c     routine for fourier analysing q_ij function
c
      implicit double precision(a-h,o-z)
c
      parameter( nind = 51 )
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1)
c.....pseudized qij(r)
      dimension qfunc(idim1,idim3,idim3)
c.....plot out arrays
      dimension plot(nind,10),isymbl(13)
c.....scratch
      dimension rho(idim1)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c     set the maximum wavevector for fourier analysis
      data gmax / 20.0d0 /
c
c     initialise the symbol array
      data isymbl /1h!,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9,1hA,1hB,1hC/
c
c-----------------------------------------------------------------------
c
c                i n i t i a l i z a t i o n
c
c     check that the dimensions of plot are sufficient
      if ( nbeta .gt. 8 ) then
        write(iout,*) '***error in subroutine fanalq'
        write(iout,*) 'nbeta =',nbeta,' .gt. 8'
        call exit(1)
      endif
c
c     compute the wavevector increment
      ginc = gmax / dble( nind - 1 )
c
      gcur = - ginc
      do 10 ig = 1,nind
        gcur = gcur + ginc
        plot(ig,1) = gcur
        plot(ig,2) = 0.0d0
   10 continue
c
c-----------------------------------------------------------------------
c
c               f o u r i e r  a n a l y s i s
c
      do 100 ibeta = 1,nbeta
c
c       set the index for the plot out array
        index = 2 + ibeta
c
c       loop over wavevectors
        gcur = - ginc
c
        do 120 ig = 1,nind
c
c         compute the current wavevector
          gcur = gcur + ginc
c
c         remeber in performing transform that qfunc = r * r * psi * psi
          do 130 ir = 1,mesh
            rho(ir) = qfunc(ir,ibeta,ibeta) * bessel(gcur*r(ir),0)
  130     continue
c
          call radin(mesh,rho,rab,asum,idim1)
          plot(ig,index) = gcur * asum
c
  120   continue
c
c     next beta
  100 continue
c
c     print out and line printer plot out
      if ( ifprt .ge. 1 ) then
c
        write (iout,'(/3x,a)')
     +     'q * fourier transform of pseudized qfunctions'
        write (iout,'(/7x,a1,5x,8(4x,i1,1x,i1,3x))')
     +     'q',(ibeta,ibeta,ibeta=1,nbeta)
c
        do ig = 1,nind
          write(iout,'(1x,9f10.5)') plot(ig,1),(plot(ig,2+i),i=1,nbeta)
        end do
c
        write (iout,'(/3x,a)') 'chart: qfunction q*Q_{ib,ib}(q)'
c
        call lpplot(0,plot,nind,2+nbeta,nind,0,nind,10,
     +    isymbl,0.0,0.0)
c
      endif
c
      return
      end

      subroutine fanalcor(r,rab,rspsco,rscore,rho,mesh,ifprt,
     +  idim1)
c
c     routine for fourier analysing core charge function
c
      implicit double precision(a-h,o-z)
c
      parameter( nind = 51 )
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1)
c.....pseudized qij(r)
      dimension rspsco(idim1), rscore(idim1)
c.....plot out arrays
      dimension plot(nind,4),isymbl(13)
c.....scratch
      dimension rho(idim1)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c     set the maximum wavevector for fourier analysis
      data gmax / 20.0d0 /
c
c     initialise the symbol array
      data isymbl /1h!,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9,1hA,1hB,1hC/
c
c-----------------------------------------------------------------------
c
c                i n i t i a l i z a t i o n
c
c     check that the dimensions of plot are sufficient
c
c     compute the wavevector increment
      ginc = gmax / dble( nind - 1 )
c
      gcur = - ginc
      do 10 ig = 1,nind
        gcur = gcur + ginc
        plot(ig,1) = gcur
        plot(ig,2) = 0.0d0
   10 continue
c
c-----------------------------------------------------------------------
c
c               f o u r i e r  a n a l y s i s
c

c
c       set the index for the plot out array

c
c       loop over wavevectors
        gcur = - ginc
c
        do  ig = 1,nind
c
c         compute the current wavevector
          gcur = gcur + ginc
c
c         remeber in performing transform that qfunc = r * r * psi * psi
          do  ir = 1,mesh
            rho(ir) = rscore(ir) * bessel(gcur*r(ir),0)
         enddo
c
          call radin(mesh,rho,rab,asum,idim1)
          plot(ig,3) = gcur * asum
c
       enddo


c       loop over wavevectors
        gcur = - ginc
c
        do  ig = 1,nind
c
c         compute the current wavevector
          gcur = gcur + ginc
c
c         remeber in performing transform that qfunc = r * r * psi * psi
          do  ir = 1,mesh
            rho(ir) = rspsco(ir) * bessel(gcur*r(ir),0)
         enddo
c
          call radin(mesh,rho,rab,asum,idim1)
          plot(ig,4) = gcur * asum
c
       enddo
c
c
c     print out and line printer plot out
      if ( ifprt .ge. 0 ) then
c
        write(iout,9999) 
 9999   format(/,' q * fourier transform of the core charge',//,
     +  //,8x,'q',9x,' full core',6x,'pseudo core')
c
        do 200 ig = 1,nind
          write(iout,9998) plot(ig,1),(plot(ig,2+i),i=1,2)
  200   continue
 9998   format(1x,f10.5,2f16.5)
c
        write(iout,9997) 
 9997   format(/,' q *  fourier transform of pseudized core charge',/)
c
        call lpplot(0,plot,nind,4,nind,0,nind,4,
     +    isymbl,0.0,0.0)
c
      endif
c
      return
      end

c
c-----------------------------------------------------------------------
c
      subroutine fanalq1(r,rab,qfunc,nltotp,mesh,ifprt,
     +  idim1,idim3)
c
c     routine for fourier analysing q_ij function
c
      implicit double precision(a-h,o-z)
c
      parameter( nind = 51 ,idm6=7)
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1)
c.....pseudized qij(r)
      dimension qfunc(idim1,idm6)
c.....plot out arrays
      dimension plot(nind,10),isymbl(13)
c.....scratch
      dimension rho(1000)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c     set the maximum wavevector for fourier analysis
      data gmax / 20.0d0 /
c
c     initialise the symbol array
      data isymbl /1h!,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9,1hA,1hB,1hC/
c
c-----------------------------------------------------------------------
c
c            c h e c k   i f p r t   f o r   r e t u r n
c
      if ( ifprt .lt. 4 ) return
c
c-----------------------------------------------------------------------
c
c                i n i t i a l i z a t i o n
c
c     check that the dimensions of plot are sufficient
      if ( nltotp .gt. 6 ) then
        write(iout,*) '***error in subroutine fanalq'
        write(iout,*) 'nltotp =',nltotp,' .gt. 6'
        call exit(1)
      endif
c
c     compute the wavevector increment
      ginc = gmax / dble( nind - 1 )
c
      gcur = - ginc
      do 10 ig = 1,nind
        gcur = gcur + ginc
        plot(ig,1) = gcur
        plot(ig,2) = 0.0d0
   10 continue
c
c-----------------------------------------------------------------------
c
c               f o u r i e r  a n a l y s i s
c
      do 100 iltot = 1,nltotp
c
c       set the index for the plot out array
        index = 2 + iltot
c
c       loop over wavevectors
        gcur = - ginc
c
        do 120 ig = 1,nind
c
c         compute the current wavevector
          gcur = gcur + ginc
c
c         remeber in performing transform that qfunc = r * r * psi * psi
          do 130 ir = 1,mesh
            rho(ir) = qfunc(ir,iltot) * bessel(gcur*r(ir),0)
  130     continue
c
          call radin(mesh,rho,rab,asum,idim1)
          plot(ig,index) = gcur * asum
c
  120   continue
c
c     next ltot
  100 continue
c
c     print out and line printer plot out
      if ( ifprt .ge. 5 ) then
c
        write(iout,9999) (iltot,iltot,iltot=1,nltotp)
 9999   format(/,' q * fourier transform of the qfunctions',//,
     +    7x,'q',5x,8(4x,i1,1x,i1,3x))
c
        do 200 ig = 1,nind
          write(iout,9998) plot(ig,1),(plot(ig,2+i),i=1,nltotp)
  200   continue
 9998   format(1x,9f10.5)
c
      endif
c
      if ( ifprt .ge. 4 ) then
c
        write(iout,9997) 
 9997   format(/,' q *  fourier transform of pseudized qfunctions',/)
c
        call lpplot(0,plot,nind,2+nltotp,nind,0,nind,10,
     +    isymbl,0.0,0.0)
c
      endif
c
      return
      end
c
c--------------------------------------------------------------------------
c
      double precision function bessel(qr,lcur)
c
c     double precision function returns value of spherical bessel function
c
c     Modified: Chris J. Pickard Feb. 1998 l=5,6
c
      implicit double precision (a-h,o-z)
c
      dimension xsplit(0:6)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c     data to determine cross over between taylor and trig expressions
      data (xsplit(l),l=0,6) /1.0d-04,5.0d-03,5.0d-02,1.0d-01,
     +     7.5d-01,1.0d0,1.0d0/
c
      if (lcur .eq. 0) then
        if (qr .lt. xsplit(lcur)) then
          fac1 = + 1.0d0
          fac2 = - 1.0d0 / 6.0d0
          fac3 = + 1.0d0 / 120.0d0
          bessel = (fac1+qr**2*(fac2+qr**2*fac3))
        else
          bessel = sin(qr)/qr
          endif
      elseif (lcur .eq. 1) then
        if (qr .lt. xsplit(lcur)) then
          fac1 = + 1.0d0 / 3.0d0
          fac2 = - 1.0d0 / 30.0d0
          fac3 = + 1.0d0 / 840.0d0
          bessel = qr*(fac1+qr**2*(fac2+qr**2*fac3))
        else
          bessel = (-cos(qr)/qr)+(sin(qr)/(qr**2))
        endif
      elseif (lcur .eq. 2) then
        if (qr .lt. xsplit(lcur)) then
          fac1 = + 1.0d0 / 15.0d0
          fac2 = - 1.0d0 / 210.0d0
          fac3 = + 1.0d0 / 7560.0d0
          bessel = qr**2*(fac1+qr**2*(fac2+qr**2*fac3))
        else
          bessel = (((3.0d0/qr**2)-1.0d0)*sin(qr)/qr)
     +             -3.0d0*cos(qr)/(qr**2)
        endif
      elseif (lcur .eq. 3) then
        if (qr .lt. xsplit(lcur)) then
          fac1 = + 1.0d0 / 105.0d0
          fac2 = - 1.0d0 / 1890.0d0
          fac3 = + 1.0d0 / 83160.0d0
          bessel = qr**3*(fac1+qr**2*(fac2+qr**2*fac3))
        else
          qr2=qr*qr
          bessel=( sin(qr)*(15.0d0/(qr2*qr)-6.0d0/qr)
     +            +cos(qr)*(1.0d0-15.0d0/qr2)           )/qr
        endif
      elseif (lcur .eq. 4) then
        if (qr .lt. xsplit(lcur)) then
          fac1 = + 1.0d0 / 945.0d0
          fac2 = - 1.0d0 / 20790.0d0
          fac3 = + 1.0d0 / 1081080.0d0
          fac4 = - 1.0d0 / 97297200.0d0
          bessel = qr**4*(fac1+qr**2*(fac2+qr**2*(fac3+qr**2*fac4)))
        else
          qr2=qr*qr
          bessel=( sin(qr)*(105.0d0/(qr2*qr2)-45.0d0/qr2+1.0d0)
     +             +cos(qr)*(10.0d0/qr-105.0d0/(qr2*qr)) )/qr
        endif
      elseif (lcur .eq. 5) then
        if (qr .lt. xsplit(lcur)) then
          fac1 = + 1.0d0 / 10395.0d0
          fac2 = - 1.0d0 / 270270.0d0
          fac3 = + 1.0d0 / 16216200.0d0
          fac4 = - 1.0d0 / 1654052400.0d0
          bessel = qr**5*(fac1+qr**2*(fac2+qr**2*(fac3+qr**2*fac4)))
        else
          qr2 = qr * qr
          bessel = ((-945.0d0*qr + 105.0d0*qr**3 - qr**5)*cos(qr) 
     +           + (945.0d0 - 420.0d0*qr**2 + 15.0d0*qr2**2)*sin(qr))
     +           /qr2**3
        endif
      elseif (lcur .eq. 6) then
        if (qr .lt. xsplit(lcur)) then
          fac1 = + 1.0d0 / 135135.0d0
          fac2 = - 1.0d0 / 4054050.0d0
          fac3 = + 1.0d0 / 275675400.0d0
          fac4 = - 1.0d0 / 31426995600.0d0
          bessel = qr**6*(fac1+qr**2*(fac2+qr**2*(fac3+qr**2*fac4)))
        else
          qr2 = qr * qr
          bessel = ((-10395.d0*qr+1260.d0*qr**3-21.d0*qr**5)*cos(qr)+
     +        (10395.d0-4725.d0*qr2+210.d0*qr2**2-qr2**3)*sin(qr))/qr**7
        endif
      else
        write(iout,*) '***error in function bessel'
        write(iout,*) 'l = ',lcur,' is out of programmed range'
        call exit(1)
      endif
c
      return
      end
c
