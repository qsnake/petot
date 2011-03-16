c----------------------------------------------------------------------------
c
      program reform
c
c     converts pseudopotential files between formatted and unformatted
c     forms
c
c     this version is compatible only with pseudopotential program
c     version 3.0.0 or higher.  it assumes logarithmic mesh.
c
c=======================================================================
c
c     directions for reforming potentials:
c
c     a typical set of input parameters for reform might be:
c
c     1                              # itype
c     022-Ti-ca-sp-vgrp.uspp         # file to be read from
c     022-Ti-ca-sp-vgrp.txt          # file to be written to
c    
c     The first number 'itype' is a key:
c
c        itype=1 means convert unformatted to formatted
c        itype=2 means convert formatted to unformatted
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
      parameter(idim1=1000,idim3=10,idim4=5,idim5=4,idim6=2*idim5-1)
      parameter(idim8=20)
c
c.....radial mesh information
      dimension r(idim1),rab(idim1)
c.....charge densities
      dimension rsatom(idim1),rspsco(idim1)
c.....pseudo quantum numbers, energies occupancies and cutoff radii
      dimension nnlzps(idim4),eeps(idim4),wwnlps(idim4),rc(idim5)
      dimension wf(idim1,idim4)
c.....beta and q functions and q pseudization coefficients
      dimension beta(idim1,idim3),qfunc(idim1,idim3,idim3),
     +qfcoef(idim8,idim6,idim3,idim3),vloc(idim1),vloc0(idim1),
     +rinner(idim6)
c.....indexing for the beta functions
      dimension nbl0(idim5),nbl(idim5)
c.....angular momenta, reference energy, qij and dij of vanderbilt scheme
      dimension lll(idim3),eee(idim3),iptype(idim3),qqq(idim3,idim3),
     +ddd0(idim3,idim3),ddd(idim3,idim3)
c.....logic to keep track of logarithmic vs herman skillman meshes
      logical tlog
c.....version numbers and date
      dimension iver(3),idmy(3)
c
c.....title
      character*20 title,xctype
c.....file name array
      character*36 infil,outfil
c
      data input,iout,iops /5,6,14/
c
c-----------------------------------------------------------------------
c
c        r e a d  i n  t h e  p a r a m e t e r s  f o r  r u n
c
c        itype             (operation type)
c        infil             (input file name)
c        outfil            (output file name)
c
c
      write(iout,10)
   10 format(' input parameter itype',/,
     +' ( itype = 1 : unformatted --- >   formatted',/,
     +'   itype = 2 :   formatted --- > unformatted',/)
      read (input,*) itype
      write(iout,*) 'selected itype =',itype
      write(iout,*) 'name of input file?'
      read (input,20) infil
      write(iout,*) 'reading from ',infil
      write(iout,*) 'name of output file?'
      read (input,20) outfil
      write(iout,*) 'writing to   ',outfil
   20 format (a36)
c
      write(iout,*) 'working on basis that file has version prefix'
c
c
c-----------------------------------------------------------------------
c
c              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c              u n f o r m a t t e d  t o  f o r m a t t e d
c              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      if ( itype .eq. + 1 ) then
c
c
c           r e a d  p s e u d o p o t e n t i a l  f r o m  f i l e
c
        open( unit = iops , file = infil , status = 'old' ,
     +    form = 'unformatted' )
c
        read (iops) (iver(i),i=1,3),(idmy(i),i=1,3)
c
c       specify logarithmic mesh
        tlog = .true.
c
        read (iops) title,zps,zvps,exftps,nvalps,mesh,etot
c
c       check that the array dimensions are sufficient
        if (mesh.gt.idim1) call error(' reform',' mesh.gt.idim1',mesh)
        if (nvalps.gt.idim4) call error(' reform',' nvalps.gt.idim4',
     +    nvalps)
c
        read (iops) (nnlzps(i),wwnlps(i),eeps(i),i=1,nvalps)
        read (iops) keyps,ifpcor,rinner(1)
c
        if ( iver(1) .ge. 3 ) then
          read(iops) nang,lloc,eloc,ifqopt,nqf,qtryc
        endif

        if (10*iver(1)+iver(2).ge.51) then
          read (iops) (rinner(i),i=1,nang*2-1)
       endif
c
        if ( iver(1) .ge. 4 ) then
          read(iops) irel
        endif
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
        if (nang.gt.idim5) call error(' reform',' nang.gt.idim5',nang)
        if (nqf.gt.idim8) call error(' reform',' nqf.gt.idim8',nqf)
c
        read (iops) (rc(i),i=1,nang)
        read (iops) nbeta,kkbeta
c
        if (nbeta.gt.idim3) call error(' reform',' nbeta.gt.idim3',
     +    nbeta)
c
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
           read (iops) rpcor
           read (iops) (rspsco(i),i=1,mesh)
        endif
c
        read (iops) (vloc(i),i=1,mesh)
        read (iops) (rsatom(i),i=1,mesh)
c
c       possible further read in of information in log mesh case
        if ( tlog ) then
          read (iops) (r(i),i=1,mesh)
          read (iops) (rab(i),i=1,mesh)
        endif
c
        if (iver(1) .ge. 6) then
c           nchi = nvales
           nchi = nvalps
           if (iver(1) .ge. 7)  read (iops) nchi
           if (nchi.gt.idim4)
     +          call error(' reform',' nchi.gt.idim4',nchi)
           read (iops) ((wf(i,j), i=1,mesh),j=1,nchi)
        endif
c
        close (iops)
c
        write (iout,190)
  190   format ('1pseudopotential has been read in')
c
c       read in complete
c
c
c           w r i t e  p s e u d o p o t e n t i a l  t o  f i l e
c
        open( unit = iops , file = outfil , status = 'new' ,
     +    form = 'formatted')
c
        write (iops,890) (iver(i),i=1,3),(idmy(i),i=1,3)
        write (iops,900) title,zps,zvps,exftps,nvalps,mesh,etot
        write (iops,910) (nnlzps(i),wwnlps(i),eeps(i),i=1,nvalps)
        write (iops,920) keyps,ifpcor,rinner(1)
        if ( iver(1) .ge. 3 ) then
          write(iops,925) nang,lloc,eloc,ifqopt,nqf,qtryc
        endif
        if (10*iver(1)+iver(2).ge.51) then
          write (iops,930) (rinner(i),i=1,nang*2-1)
       endif

        if ( iver(1) .ge. 4 ) then
          write(iops,926) irel
        endif
        write (iops,930) (rc(i),i=1,nang)
        write (iops,920) nbeta,kkbeta
        do 300 j=1,nbeta
          write (iops,940) lll(j),eee(j),(beta(i,j),i=1,kkbeta)
          do 300 k=j,nbeta
            write (iops,930) ddd0(j,k),ddd(j,k),qqq(j,k),
     +      (qfunc(i,j,k),i=1,kkbeta),
     +      ((qfcoef(i,lp,j,k),i=1,nqf),lp=1,nlc)
  300   continue
c
        if (10*iver(1)+iver(2).ge.72) then
          write (iops,890) (iptype(j),j=1,nbeta)
          write (iops,910) npf,ptryc
        endif
c
        write (iops,930) rcloc,(vloc0(i),i=1,mesh)
        if (ifpcor.gt.0) then
           write (iops,930) rpcor
           write (iops,930) (rspsco(i),i=1,mesh)
        endif
        write (iops,930) (vloc(i),i=1,mesh)
        write (iops,930) (rsatom(i),i=1,mesh)
c
        if ( tlog ) then
          write (iops,930) (r(i),i=1,mesh)
          write (iops,930) (rab(i),i=1,mesh)
        endif

        if (iver(1) .ge. 6) then
           if (iver(1) .ge. 7)  write (iops,*) nchi
           write (iops,930) ((wf(i,j), i=1,mesh),j=1,nchi)
        endif

c
        close (iops)
c
        write (iout,310)
  310   format ('1pseudopotential has been written out')
c
  890   format (6i5)
  900   format (a20,3f15.9/2i5,1pe19.11)
  910   format (i5,2f15.9)
  920   format (2i5,f15.9)
  925   format (2i5,f9.5,2i5,f9.5)
  926   format (i5)
  930   format (1p4e19.11)
  940   format (i5/(1p4e19.11))
c
      endif
c
c-----------------------------------------------------------------------
c
c              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c              f o r m a t t e d  t o  u n f o r m a t t e d
c              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      if ( itype .eq. + 2 ) then
c
c
c           r e a d  p s e u d o p o t e n t i a l  f r o m  f i l e
c
        open( unit = iops , file = infil , status = 'old' ,
     +    form = 'formatted' )
c
        read (iops,945) (iver(i),i=1,3),(idmy(i),i=1,3)
c
c       specify logarithmic mesh
        tlog = .true.
c
        read (iops,950) title,zps,zvps,exftps,nvalps,mesh,etot
c
c       check that the array dimensions are sufficient
        if (mesh.gt.idim1) call error(' reform',' mesh.gt.idim1',mesh)
        if (nvalps.gt.idim4) call error(' reform',' nvalps.gt.idim4',
     +    nvalps)
c
        read (iops,960) (nnlzps(i),wwnlps(i),eeps(i),i=1,nvalps)
        read (iops,970) keyps,ifpcor,rinner(1)
c
        if ( iver(1) .ge. 3 ) then
          read(iops,975) nang,lloc,eloc,ifqopt,nqf,qtryc
        endif
c
        if (10*iver(1)+iver(2).ge.51) then
          read (iops,980) (rinner(i),i=1,nang*2-1)
       endif

        if ( iver(1) .ge. 4 ) then
          read(iops,976) irel
        endif
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
        if (nang.gt.idim5) call error(' reform',' nang.gt.idim5',nang)
        if (nqf.gt.idim8) call error(' reform',' nqf.gt.idim8',nqf)
c
        read (iops,980) (rc(i),i=1,nang)
        read (iops,970) nbeta,kkbeta
c
        if (nbeta.gt.idim3) call error(' reform',' nbeta.gt.idim3',
     +    nbeta)
c
        do 400 j=1,nbeta
          read (iops,990) lll(j),eee(j),(beta(i,j),i=1,kkbeta)
          do 400 k=j,nbeta
            read (iops,980) ddd0(j,k),ddd(j,k),qqq(j,k),
     +      (qfunc(i,j,k),i=1,kkbeta),
     +      ((qfcoef(i,lp,j,k),i=1,nqf),lp=1,nlc)
            ddd0(k,j)=ddd0(j,k)
            ddd (k,j)=ddd (j,k)
            qqq (k,j)=qqq (j,k)
            do 400 i=1,kkbeta
              qfunc(i,k,j)=qfunc(i,j,k)
  400   continue
        do 410 i=1,nqf
        do 410 lp = 1,nlc
          qfcoef(i,lp,k,j)=qfcoef(i,lp,j,k)
  410   continue
c
        if (10*iver(1)+iver(2).ge.72) then
          read (iops,945) (iptype(j),j=1,nbeta)
          read (iops,960) npf,ptryc
        endif
c
        read (iops,980) rcloc,(vloc0(i),i=1,mesh)
c
c       set index arrays nbl and nbl0
        do 450 lp = 1,nang
          nbl(lp) = 0
  450   continue
        do 460 ib = 1,nbeta
          lp = lll(ib) + 1
          nbl(lp) = nbl(lp) + 1
  460   continue
        lmaxps = lll(nbeta)
        nbl0(1) = 0
        do 470 i = 1,lmaxps
          nbl0(i+1) = nbl0(i) + nbl(i)
  470   continue
c
c       possible readin of the frozen core correction
        if (ifpcor.gt.0) then
           read (iops,980) rpcor
           read (iops,980) (rspsco(i),i=1,mesh)
        endif
c
        read (iops,980) (vloc(i),i=1,mesh)
        read (iops,980) (rsatom(i),i=1,mesh)
c
c       possible further read in of information in log mesh case
        if ( tlog ) then
          read (iops,980) (r(i),i=1,mesh)
          read (iops,980) (rab(i),i=1,mesh)
        endif
        if (iver(1) .ge. 6) then
c     nchi = nvales
           nchi = nvalps
           if (iver(1) .ge. 7)  read (iops,*) nchi
           read (iops,980) ((wf(i,j), i=1,mesh),j=1,nchi)
        endif
c
c
        close (iops)
c
        write (iout,490)
  490   format ('1pseudopotential has been read in')
c
  945   format (6i5)
  950   format (a20,3f15.9/2i5,1pe19.11)
  960   format (i5,2f15.9)
  970   format (2i5,f15.9)
  975   format (2i5,f9.5,2i5,f9.5)
  976   format (i5)
  980   format (1p4e19.11)
  990   format (i5/(1p4e19.11))
c
c       read in complete
c
c
c           w r i t e  p s e u d o p o t e n t i a l  t o  f i l e
c
        open( unit = iops , file = outfil , status = 'new' ,
     +    form = 'unformatted')
c
        write (iops) (iver(i),i=1,3),(idmy(i),i=1,3)
        write (iops) title,zps,zvps,exftps,nvalps,mesh,etot
        write (iops) (nnlzps(i),wwnlps(i),eeps(i),i=1,nvalps)
        write (iops) keyps,ifpcor,rinner(1)
        if ( iver(1) .ge. 3 ) then
          write(iops) nang,lloc,eloc,ifqopt,nqf,qtryc
        endif
        if (10*iver(1)+iver(2).ge.51) then
          write (iops) (rinner(i),i=1,nang*2-1)
       endif

        if ( iver(1) .ge. 4 ) then
          write(iops) irel
        endif
        write (iops) (rc(i),i=1,nang)
        write (iops) nbeta,kkbeta
        do 500 j=1,nbeta
          write (iops) lll(j),eee(j),(beta(i,j),i=1,kkbeta)
          do 500 k=j,nbeta
            write (iops) ddd0(j,k),ddd(j,k),qqq(j,k),
     +      (qfunc(i,j,k),i=1,kkbeta),
     +      ((qfcoef(i,lp,j,k),i=1,nqf),lp=1,nlc)
  500   continue
c
        if (10*iver(1)+iver(2).ge.72) then
          write (iops) (iptype(j),j=1,nbeta),npf,ptryc
        endif
c
        write (iops) rcloc,(vloc0(i),i=1,mesh)
        if (ifpcor.gt.0) then
           write (iops) rpcor
           write (iops) (rspsco(i),i=1,mesh)
        endif
        write (iops) (vloc(i),i=1,mesh)
        write (iops) (rsatom(i),i=1,mesh)
c
        if ( tlog ) then
          write (iops) (r(i),i=1,mesh)
          write (iops) (rab(i),i=1,mesh)
        endif
        if (iver(1) .ge. 6) then
           if (iver(1) .ge. 7)  write (iops) nchi
           write (iops) ((wf(i,j), i=1,mesh),j=1,nchi)
        endif
c
c
        close (iops)
c
        write (iout,510)
  510   format ('1pseudopotential has been written out')
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
        do 1140 ib=1,nbeta
          lp=lll(ib)+1
          if (10*iver(1)+iver(2).lt.72) iptype(ib)=-1
c           -1 means iptype was not recorded and is unknown
          write (iout,1150) ib,lll(ib),eee(ib),rc(lp),iptype(ib)
 1140   continue
 1150   format(4x,'|',6x,i2,5x,i2,4x,2f7.2,6x,i1,18x,'|')
c
        if (10*iver(1)+iver(2).ge.72) then
          ifip2=0
          do ib=1,nbeta
            if (iptype(ib).eq.2) ifip2=1
          end do
          if (ifip2.eq.1) then
            write (iout,1155) npf,ptryc
 1155       format(4x,'|  npf    =',i2,'  ptryc =',f8.3,29x,'|')
          endif
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
      end
c
c-----------------------------------------------------------------------
      subroutine error(a,b,n)
c-----------------------------------------------------------------------
      character*(*) a,b
      write(6,1) a,b,n
    1 format(//' program ',a,':',a,'.',8x,i8,8x,'stop')
      stop
      end
c 
