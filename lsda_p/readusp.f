       program readusp

c      subroutine rwps(iread,title,zps,zvps,exftps,nvalps,irel,
c     +  z,zv,exfact,nvales,mesh,etot,nnlzps,wwnlps,eeps,
c     +  nnlz,wwnl,ee,ncores,keyps,ifpcor,rpcor,rinner,rc,nbeta,
c     +  kkbeta,lll,eee,iptype,npf,ptryc,
c     +  beta,ddd0,ddd,qqq,qfunc,qfcoef,
c     +  rcloc,vloc0,vloc,snl,lloc,eloc,rsatom,rspsco,flname,
c     +  ifqopt,nqf,qtryc,nang,nbl0,nbl,r,rab,iver,idmy,
c     +  idim1,idim2,idim3,idim4,idim5,idim6,idim8)

c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)

      parameter( idim1 = 1000      , idim2 = 26         )
      parameter( idim3 = 10        , idim4 = 5          )
      parameter( idim5 = 4         , idim6 = 2*idim5-1  )
      parameter( idim7 = 4         , idim8 = 20         )

c     idim1  .ge.  no. of points in the radial mesh
c     idim2  .ge.  no. of shells in the atom
c     idim3  .ge.  no. of beta functions in vanderbilt potential
c     idim4  .ge.  no. of valence states for pseudization
c     idim5  .ge.  lmax + 1, lmax is maximum l value of potential
c     idim7  .ge.  no. of reference states for each angular momentum value
c     idim8  .ge.  no. of terms in taylor expansion of pseudized q_i,j



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
c      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c
c--------------------------------------------------------------------------
c
c       r e a d  p s e u d o p o t e n t i a l  f r o m  f i l e
c

        iops=10
        iout=6
        stderr=6
        write(6,*) "input the name of unformatted ultrasoft pseud file"
        read(5,*) flname(6)
      

c       save data so compatibility checks can be run later
c        mold  = mesh
c        r2    = r(2)
c        rmesh = r(mesh)
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
c           nchi = nvales        ! warning, used before defined ! use nvalps ?
           nchi = nvalps        ! warning, used before defined ! use nvalps ?
           if (iver(1) .ge. 7)  read (iops) nchi
           read (iops) ((snl(i,j), i=1,mesh),j=ncores+1,ncores+nchi)
        endif
c
        close (iops)
c
        write (iout,'(/a)') (' pseudopotential has been read in')
c
c
c
c       nvales can be smaller, e.g. psp has s,p,d but testing
c       neutral config so only want to test s,p
        nvalesmin= nvales
c        if ( nvales .gt. nvalps ) then
c           nvalesmin= nvalps
c          write(iout,*) '***error in subroutine rwps'
c          write(iout,*) 'nvales=',nvales,'.gt. nvalps=',nvalps
c          call exit(1)
c        endif
c
c       i=1,nvales must label s, p, d respectively
c        do 200 i=1,nvalesmin
c          if ( nnlz(ncores+i) .ne. nnlzps(i) ) then
c            write(iout,*) '***error in subroutine rwps'
c            write(iout,*) 'nnlz',nnlz(ncores+i),'.ne. nnlzps',nnlzps(i)
c            call exit(1)
c          endif
c  200   continue
c
c       check that relativity switches compatible
c        if ( iver(1) .lt. 4 ) then
c         potentials generated prior to version 4 are never relativistic
c          if ( irel .ne. 0 ) then
c            write(iout,*) '***error in subroutine rwps'
c            write(iout,*) 'irel=',irel,' .and. iver=',iver(1),' illegal'
c            call exit(1)
c          endif
c        else
c         only proceed if relativistic switches of ae and ps agree
c          if ( irelps .ne. irel ) then
c            write(iout,*) '***error in subroutine rwps'
c            write(iout,*) 'irel =',irel,' .ne. irelps =',irelps
c            call exit(1)
c          endif
c        endif
c
c       read in complete
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
cc ----------------------------------------------------------------------
cc ----------------------------------------------------------------------
cc ----------------------------------------------------------------------
      write(6,*) "kkbeta=",kkbeta
      write(6,*) "mesh=", mesh
      write(6,*) "nbeta=",nbeta

      write(6,*) "****** ddd ******"
      do ii=1,nbeta
      write(6,2011) (ddd(jj,ii),jj=1,nbeta)
      enddo
      write(6,*) "***************"
      write(6,*) "****** ddd0 ******"
      do ii=1,nbeta
      write(6,2011) (ddd0(jj,ii),jj=1,nbeta)
      enddo
      write(6,*) "***************"
      write(6,*) "****** qqq ******"
      do ii=1,nbeta
      write(6,2011) (qqq(jj,ii),jj=1,nbeta)
      enddo
      write(6,*) "***************"
2011  format(6(E12.4,1x))

      open(12,file="graph")
      rewind(12)
c      do i=2,mesh
c      write(12,2010) r(i),rsatom(i)/r(i)**2,
c     &     vloc(i)/r(i),vloc0(i)/r(i)
c      enddo
cccccccccccc
c      do i=2,kkbeta
c      write(12,2010) r(i),(beta(i,j)/r(i),j=1,nbeta)
c      enddo
cccccccccccc
c      do i=2,mesh
c      write(12,2010) r(i),(snl(i,j)/r(i),j=1,nvalps)
c      enddo
c      close(12)
      do i=2,kkbeta
      write(12,2010) r(i),qfunc(i,1,1),qfunc(i,2,2),
     &  qfunc(i,1,2),qfunc(i,3,3),qfunc(i,1,3),qfunc(i,i,4)
      enddo
      close(12)
2010  format(7(E12.5,1x))

cccccccccccccccccccccccccccccccccccccccccccccc
ccc  definitions of the read-in functions:
ccc    do 100 i=1,mesh
ccc    r(i)=a*(exp(b*(i-1))-1)
ccc    rab(i)=b*(r(i)+a)
ccc    sqr(i)=exp(b*i/2)
ccc    enddo
ccc
ccc    local ionic potential:     vloc0(i)/r(i)          ! in Ryd
ccc    total screened potential:  vloc(i)/r(i)           ! in Ryd
ccc    atom charge:               rsatom(i)/r(i)**2
ccc    valence el. wave func:     snl(i,j)/r(i), j=1,nvalps
ccc    projector wave func:       beta(i,j)/r(i), j=1,nbeta, i=1,kkbeta (real space)
ccc    The nolocal potential is:
ccc    V_NL=sum_L,m sum_i,j  ddd(i,j) |beta(:,i)/r Y_Lm> <beta(:,j)/r Y_Lm|
ccc    note the overlap between L and i,j. Also, ddd(i,j) is zero if L(i).ne.L(j)
ccc    The overlap operator is:
ccc    S=1+sum_L,m sum_i,j qqq(i,j) |beta(:,i)/r Y_Lm> <beta(:.j)/r Y_Lm|
ccc    The orthnormal condition is:
ccc    <psi_i1| S | psi_i2> = delta_i1,i2
ccc    The Schroedinger's equationi s:
ccc    H psi_i1 = E_i1 S psi_i1
ccc    The charge density is:
ccc    rho(r) = sum_i1 |psi_i1(r)|^2+ sum_i,j  qfunc(r,i,j) X 
ccc        sum_i1 <beta(:,i)/r Y_Lm| psi_i1>  < psi_i1|beta(:,j)/r Y_Lm>  
cccccccccccccccccccccccccccccccccccc







      return
      end
c
c----------------------------------------------------------------------------
