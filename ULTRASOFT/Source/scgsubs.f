c
c Copyright (C) 2002 Vanderbilt group
c This file is distributed under the terms of the GNU General Public
c License as described in the file 'License' in the current directory.
c
c-------------------------------------------------------------------------
c
      subroutine scpgef(r,rab,sqr,a,b,rlogd,klogd,irel,ruae,
     +  zz,vzero,mesh,ifprt,lmaxps,nbeta,rcloc,rc,rinner,
     +  lll,eee,iptype,npf,ptryc,iploctype,nbl,nbl0,
     +  nang,ddd,ddd0,kkbeta,beta,qqq,qfunc,qfcoef,ifqopt,nqf,qtryc,
     +  lloc,eloc,sss,ttt,bbb,bbbi,iwork,ibfix,rfix,qfix,nfix,
     +  psi,phi,chi,vloc,yy,bndmat,npivot,dum1,dum2,dum3,
     +  ikeyee,ikeyeeloc,refwf,
     +  idim1,idim3,idim5,idim6,idim7,idim8)
c
c----------------------------------------------------------------------------
c     routine for generating soft core pseudopotentials
c     using the generalised eigenvalue formulation
c
c     see d.h. vanderbilt, phys. rev. b 41 p.7892 (1990)
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1),sqr(idim1)
c.....potentials for the all electron calculation
      dimension ruae(idim1)
c.,...work arrays for the tcm routine call
      dimension bndmat(idim1,5),npivot(idim1)
c.....cutoff radii for the chi and ql functions
      dimension rc(idim5),rinner(idim6)
c.....beta and q functions and q pseudization coefficients
      dimension beta(idim1,idim3),qfunc(idim1,idim3,idim3),
     +qfcoef(idim8,idim6,idim3,idim3)
c.....extra fixed points to avoid negative densities when pseudizing q_ij
      dimension ibfix(idim3),rfix(idim3),qfix(idim3)
c.....angular momenta, reference energy, qij and dij of vanderbilt scheme
      dimension lll(idim3),eee(idim3),iptype(idim3),qqq(idim3,idim3),
     +ddd0(idim3,idim3),ddd(idim3,idim3),sss(idim3,idim3),
     +ttt(idim3,idim3),bbb(idim3,idim3),bbbi(idim3,idim3),iwork(idim3)
c.....indexing for the beta functions
      dimension nbl0(idim5),nbl(idim5)
c.....psi,phi and chi functions for generation of vanderbilt potential
      dimension psi(idim1,idim3),phi(idim1,idim3),chi(idim1,idim3)
c.....work space for relativistic wavefunctions
      dimension yy(idim1,2)
c.....type of wave function
      dimension ikeyee(idim3)
      dimension refwf(idim1,idim3)
c.....scratch arrays
      dimension vloc(idim1),dum1(idim1),dum2(idim1),dum3(idim1),con(10)
c
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c     input data for the cubic spline call
      data as1,bs1,asn,bsn,isx / 0.0d0,0.0d0,0.0d0,0.0d0,0 /
c
c---------------------------------------------------------------------------
c
c                     i n i t i a l i s a t i o n
c
      write(iout,'(/1x,47(1h-))')
      write(iout,'(1x,a)') 'begin subroutine scpgef'
      write(iout,'(1x,47(1h-)/)')
c
c     check that rlogd is greater than all the pseudization radii 
      rcmax = rcloc
      do 10 i = 1,nang
        rcmax = dmax1(rcmax,rc(i))
   10 continue
c
      if ( rcmax .gt. rlogd ) then
        write(iout,*) '***error in subroutine scpgef'
        write(iout,*) 'rcmax =',rcmax,' .gt. rlogd =',rlogd
        call exit(1)
      endif
c
c     set kkbeta to be a bit larger than klogd
      kkbeta = klogd + 5
c     ensure that kkbeta is odd for future calls to radin
      if ( kkbeta - 2 * ( kkbeta / 2 ) .ne. 1 ) kkbeta = kkbeta + 1
      rbeta = r(kkbeta)
c
      write(iout,9999) rcmax,rbeta,rlogd,kkbeta,klogd
 9999 format('   rcmax, rbeta, rlogd = ',3f10.5,/,
     +       '          kkbeta,klogd = ',10x,2i10)
c
c     check that kkbeta has not exceeded the mesh size
      if ( kkbeta .gt. mesh ) then
        write(iout,*) '***error in subroutine scpgef'
        write(iout,*) 'kkbeta =',kkbeta,' .gt. mesh =',mesh
        call exit(1)
      endif
c
c---------------------------------------------------------------------------
c
c        p s e u d i z e  t h e  l o c a l  p o t e n t i a l
c
c
c     fill vloc with the all electron potential
      do 100 i = 2,mesh
        vloc(i) = ruae(i)
  100 continue
c
      if ( lloc .eq. -1 ) then
c
c       ============================
c       original pseudization method
c       ============================
c
        write(iout,'(3x,a)')
     +    'lloc .eq. -1 so original pseudization for vloc'
c
        call psrv(vloc,rcloc, +1 ,r,mesh,ifprt,idim1)
c
c
      elseif ( lloc .ge. 0 .or. lloc .le. 4 ) then
c
c       =======================
c       new pseudization scheme
c       =======================
c
        write(iout,'(3x,a,i2,a,f8.3)') 'lloc =',lloc,
     +    ' so new pseudization scheme for vloc with eloc = ',eloc
c
        call psv(ruae,vloc,rcloc,rcmax,r,rab,sqr,a,b,mesh,
     +   irel,yy,bndmat,
     +   kkbeta,klogd,dum1,dum2,dum3,lloc,eloc,ifprt,idim1,iploctype,
     +   ikeyeeloc,refwf,
     +   nang,npf,ptryc,idim8,idim3)
c
      else
c
        write(iout,*) '***error in subroutine scpgef'
        write(iout,*) 'lloc is unreasonable, lloc =',lloc
        call exit(1)
c
      endif
c
c     print out of the local potential
      call prvloc(ifprt,vloc,rcloc,r,a,b,idim1)
c
c---------------------------------------------------------------------------
c
c           c o m p u t e  t h e  c h i  f u n c t i o n s
c
c     loop over the reference states
      do 200 ibeta = 1,nbeta
c
c       initialise angular momenta, energy and cutoff radius
        lcur = lll(ibeta)
        ecur = eee(ibeta)
        ipsu = iptype(ibeta)
        lp   = lcur + 1
        rcut = rc(lp)
        llp  = lcur * ( lcur + 1 )
c
        write(stderr,*) 'constructing ibeta =',ibeta
        if ( ifprt .ge. 0) then
          write(iout,'(/3x,43(1h-))')
          write(iout,'(3x,a,i2)') 'constructing ibeta =',ibeta
          write(iout,'(3x,43(1h-))')
          write(iout,'(/3x,a,2i5,f8.4,i5)')
     +     'ibeta,lcur,ecur,ipsu',
     +     ibeta,lcur,ecur,ipsu
        endif
c
c       obtain all electron wavefunction at energy ecur
        if (ikeyee(ibeta) .lt. 0) then
           if ( ifprt .ge. 0) write(iout,'(5x,a)')
     +        'taking wavefunction from input'
           do i=1,kkbeta
              psi(i,ibeta) = refwf(i,3-ikeyee(ibeta))
           enddo
           do i=1,kkbeta
            if (r(i) .ge. rcut) goto 203
         enddo
 203     if (psi(i,ibeta) .lt. 0.0) then
            do i=1,kkbeta
               psi(i,ibeta) = -1.0*psi(i,ibeta)
            enddo
         endif

        elseif ( irel .eq. 0 ) then
c
          write(iout,'(3x,a)') 'calling phase'
          zz = 0.0d0
          vzero = ruae(2) / r(2)
          call phase(zz,vzero,ecur,lcur,b,r,rab,sqr,
     +      ruae,dum1,psi(1,ibeta),mesh,kkbeta,dlwf,idim1)
c
c         renormalise psi and remove square root type factor
          factor = 1.0d0 / ( sqr(klogd) * psi(klogd,ibeta) )
          do 202 i = 1,kkbeta
            psi(i,ibeta) = psi(i,ibeta) * sqr(i) * factor
  202     continue
c
        elseif ( irel .eq. 2 ) then
c
          write(iout,'(3x,a)') 'calling phakoe'
          call phakoe(ruae,r,rab,bndmat,ecur,lcur,
     +      yy,dum1,kkbeta,dlwf,idim1)
c
c         copy major component of koelling equation into psi
c         and renormalise
          factor = 1.0d0 / yy(klogd,2)
          do 205 ir = 1,kkbeta
            psi(ir,ibeta) = yy(ir,2) * factor
  205     continue
c
        else
c
          write(iout,*) '***error in scpgef'
          write(iout,*) 'irel =',irel,' is not allowed'
          call exit(1)
c
        endif
c
c       set phi equal to psi in preparation for pseudization
        do 210 i = 1,kkbeta
          phi(i,ibeta) = psi(i,ibeta)
  210   continue
c
c       pseudise phi
        if ( ipsu .eq. 0 ) then
c
c         polynomial pseudisation
          write(iout,'(3x,a)') 'calling pswf - polynomial pseudization'
          call pswf(phi(1,ibeta),dum1,rcut,lcur,r,rab,con,kkbeta,
     +      ifprt,idim1)
c
        elseif ( ipsu .eq. 1 ) then
c
          write(iout,'(3x,a)')
     +      'calling pswf1 - exponential pseudization'
          call pswf1(phi(1,ibeta),dum1,rcut,lcur,r,rab,con,kkbeta,
     +      ifprt,idim1)
          npfl=4
c
        elseif ( ipsu .eq. 2 ) then
c
c         polynomial optimization
          write(iout,'(3x,a)')
     +      'calling pswf2 - polynomial optimization'
          call pswf2(phi(1,ibeta),dum1,rcut,lcur,r,rab,con,kkbeta,
     +      ifprt,nang,npf,ptryc,idim1,idim8)
          npfl=npf
c     
        elseif ( (ipsu .eq. 3 ).or.( ipsu .eq. 4)) then
c
c         normconserving pseudozing
          write(iout,'(3x,a)')
     +      'calling pswfnormc - normconserving pseudization'
          if (ikeyee(ibeta) .lt. 0) then
             call pswfnormc(phi(1,ibeta),dum1,
     +            refwf(1,1),rcut,lcur,ecur,r,rab,
     +            b,con,kkbeta,ifprt,idim1,10)
          else
             call pswfnormc(phi(1,ibeta),dum1,ruae,rcut,lcur,ecur,r,rab,
     +            b,con,kkbeta,ifprt,idim1,idim8)
          endif
          rc(lcur+1)=rcut
          npfl=7
c
        else
c
          write(iout,*) '***error in subroutine scpgef'
          write(iout,*) 'ipsu =',ipsu,' is out of range'
          call exit(1)
c
        endif
c
c       generate analytic value for chi(2,ibeta)
        call polyn(r(2),d2p,r(1),phi(1,ibeta),-1)
        call polyn(r(2),d2p,r(1),phi(1,ibeta),+2)
c        write(*,*) 'polyn2',d2p
c        do i=1,12
c           write(*,*) r(i),phi(i,ibeta)
c        enddo
c
        if (ikeyee(ibeta) .lt. 0) then
           do i=1,kkbeta
              dum2(i) = vloc(i)+refwf(i,2)-refwf(i,3)
           enddo
        else
           do i=1,kkbeta
              dum2(i) = vloc(i)
           enddo
        endif

        chi(2,ibeta) = + d2p - 
     +  ( dum2(2) / r(2) + dble(llp) / r(2)**2.0 - ecur ) * phi(2,ibeta)
c
        dum1(2) = rab(2) * rab(2) * chi(2,ibeta) / sqr(2)
        dum1(kkbeta) = 0.0d0
c
c
c       set the sum of potentials and centrifugal terms
        call setv(r,rab,dum2,b,ecur,llp,dum3,kkbeta,idim1)

        do 220 i = 1,kkbeta
          dum2(i) = phi(i,ibeta) / sqr(i)
  220   continue

c
c       invert the differential equation
        call difinv(dum3,dum2,dum1,3,kkbeta-1,bndmat,npivot,idim1)
CCCCCCCCCCCCCCCCCCCCCCCCCC
C     write (19,'(4e18.9)') (r(j),dum3(j),dum2(j),dum1(j),j=1,20)
CCCCCCCCCCCCCCCCCCCCCCCCCC
c
c       compute the chi function for this ibeta
        chi(1,ibeta) = 0.0d0
        do 230 i = 2,kkbeta
          chi(i,ibeta) = sqr(i) * dum1(i) / ( rab(i) * rab(i) )
  230   continue

c       zero chi leaving one point beyond max(rcut,rcloc)
        chi(1,ibeta) = 0.0
c        write(*,*) r(1),chi(1,ibeta)
        do 240 i = 1,kkbeta-1
           if (( r(i) .ge. rcut ) .and. (r(i) .ge. rcloc))
     $          chi(i+1,ibeta) = 0.0d0
  240   continue

        if((ipsu .eq. 3) .or.( ipsu .eq. 4)) then
c for numerical stability calculate the projectors analytical 
c this might not be nessesary
        do i = 2,kkbeta
           if ( r(i) .le. rcut ) then
               t1 = 0.d0
               t2 = 0.d0
               do j =2,npfl
                  t1 = t1 + dble(2*(j-1))*con(j)*r(i)**dble(2*j-3)
                  t2 = t2 + dble((2*j-2)*(2*j-3))*con(j)*r(i)**(2*j-4)
               enddo
               vtemp = ecur + dble(2*lcur+2)*t1/r(i) + t1**2.0 + t2
c     remember that vloc is potential * r
            else
c     no pseudization necessary beyond rcloc
               vtemp = ruae(i)/r(i)
               if (ikeyee(ibeta) .lt. 0) then
                  vtemp = refwf(i,1)/r(i)
               endif
            endif
c        write(*,*) r(i),chi(i,ibeta),(vtemp-vloc(i))*phi(i,ibeta)/r(i)
            chi(i,ibeta) = (vtemp-vloc(i)/r(i))*phi(i,ibeta)
            if (ikeyee(ibeta) .lt. 0) then
               chi(i,ibeta) = (vtemp-(vloc(i)+refwf(i,2)-refwf(i,3))
     $              /r(i))*phi(i,ibeta)
            endif
         enddo
      endif
c
c *** here is a mystery: in a5.1.3 this was chi**2, but it
c *** was changed by Kurt evidently to phi**2.  I don't
c *** understand why; I guess it shouldn't matter much, but
c *** as I want to be able to compare runs in detail with a5.1.3,
c *** I have changed it back to chi**2.  -- DV 7/30/97
c
c c       balance by normalizing to unit charge
c         do 250 i = 1,kkbeta
c           dum1(i) = phi(i,ibeta) ** 2.0
c   250   continue
c c
c       balance by normalizing to unit charge
        do 250 i = 1,kkbeta
          dum1(i) = chi(i,ibeta) ** 2
  250   continue
c
        call radin(kkbeta,dum1,rab,asum,idim1)
c
        factor = 1.0d0 / sqrt(asum)
        do 260 i = 1,kkbeta
          psi(i,ibeta) = psi(i,ibeta) * factor
          phi(i,ibeta) = phi(i,ibeta) * factor
          chi(i,ibeta) = chi(i,ibeta) * factor
  260   continue
c
c       possible print out of results
        if ( ifprt .ge. 1 ) then
c
          write(iout,9992) ibeta
 9992     format(/,3x,'results for ibeta =',i2,/,
     +            '    psi - unpseudised wavefunction',/,
     +            '    phi - pseudised wavefunction',/,
     +            '    chi - projector function',//,
     +      7x,'r',9x,'psi',10x,'phi',10x,'chi')
          nind = 70
          rinc = rcmax / dble(nind)
          rcur = -0.75d0 * rinc
          do 270 i = 1,nind
            rcur = rcur + rinc
            in = int( dlog( rcur / a + 1.0d0 ) / b ) + 1
            write(iout,9991) r(in),psi(in,ibeta),
     +        phi(in,ibeta),chi(in,ibeta)
  270     continue
c
        endif
 9991   format(2x,f9.4,3f13.8)
c
        
  200 continue
c
c---------------------------------------------------------------------------
c
c         s e t  b b b  q q q  q f u n c  a n d  d d d
c
      if ( ifprt .ge. 0) then
        write(iout,'(/3x,43(1h-))')
        write(iout,'(3x,a)') 
     +   'construct beta functions and matrices'
        write(iout,'(3x,43(1h-))')
      endif
c
c no qfunction for normconserving pseudos
      do ib=1,nbeta
         if (iptype(ib) .eq. 3) then
            do i=1,kkbeta
               psi(i,ib) = phi(i,ib)
            enddo
         endif
      enddo

      do 300 ib = 1,nbeta
c
c       compute the kinetic operator on chi using spline routine
        call splift(r,chi(1,ib),dum1,dum2,kkbeta,bndmat,ierr,
     +    isx,as1,bs1,asn,bsn)
CCCCCCCCCCCCCCCCCCCCCCCCCC
C     write (iout,'(a,i5)') 'dum2',ib
C     write (iout,'(3e16.7)') (r(j),chi(j,ib),dum2(j),j=1,10)
C     write (iout,'(3e16.7)') (r(j),chi(j,ib),dum2(j),j=20,kkbeta,40)
C     write (17,*) kkbeta
C     write (17,'(3e18.10)') (r(j),chi(j,ib),dum2(j),j=1,kkbeta)
CCCCCCCCCCCCCCCCCCCCCCCCCC
c
        if ( ierr .ne. 1 ) then
          write(iout,*) '***error in subroutine scpgef'
          write(iout,*) 'abnormal exit from splift, ierr=',ierr
          call exit(1)
        endif
c
c       dum2 contains second derivative of chi so multiply by -1 for ke
        do 310 i = 1,kkbeta
          dum2(i) = - dum2(i)
  310   continue
        if((ipsu .eq. 3) .or. (ipsu .eq. 4)) then
          do i = 2,kkbeta
            if ( r(i) .le. rcut ) then
                t1 = 0.d0
                t2 = 0.d0
                do j =2,npfl
                  t1 = t1 + dble(2*(j-1))*con(j)*r(i)**(2*j-3)
                  t2 = t2 + dble((2*j-2)*(2*j-3))*con(j)*r(i)**(2*j-4)
                enddo
                vtemp = ecur + dble(2*lcur+2)*t1/r(i) + t1**2 + t2
c               remember that vloc is potential * r
             else
c              no pseudization necessary beyond rcloc
               vtemp = ruae(i)/r(i)
             endif
             dum2(i) = (ecur-vtemp)*phi(i,ibeta)
           enddo
         endif
c
c
        do 300 jb = 1,nbeta
c
          li = lll(ib)
          lj = lll(jb)
c
c         set qfunc
          do 320 i = 1,kkbeta
            qfunc(i,ib,jb) = psi(i,ib)*psi(i,jb) - phi(i,ib)*phi(i,jb)
  320     continue
c
c         set bbb, sss, ttt and qqq
c
          if ( li .ne. lj ) then
            qqq(ib,jb) = 0.0d0
            sss(ib,jb) = 0.0d0
            ttt(ib,jb) = 0.0d0
            bbb(ib,jb) = 0.0d0
          else
c
            do 330 i = 1,kkbeta
              dum1(i) = qfunc(i,ib,jb)
  330       continue
            call radin(kkbeta,dum1,rab,qqq(ib,jb),idim1)
c
            do 340 i = 1,kkbeta
              dum1(i) = chi(i,ib) * chi(i,jb)
  340       continue
            call radin(kkbeta,dum1,rab,sss(ib,jb),idim1)
c
            do 350 i = 1,kkbeta
              dum1(i) = dum2(i) * chi(i,jb)
  350       continue
            call radin(kkbeta,dum1,rab,ttt(ib,jb),idim1)
c
            do 360 i = 1,kkbeta
              dum1(i) = phi(i,ib) * chi(i,jb)
  360       continue
            call radin(kkbeta,dum1,rab,bbb(ib,jb),idim1)
c
          endif
c
c         set ddd
c
          ddd(ib,jb) = bbb(ib,jb) + eee(jb) * qqq(ib,jb)
c
  300 continue
c
c---------------------------------------------------------------------------
c
c     p r i n t  /  p l o t  t h e  c h i  f u n c t i o n s
c
      call prbeta(ifprt,'chi ',chi,kkbeta,nbeta,lll,r,a,b,idim1,idim3)
c     print and plot fourier transform of chi functions
      call fanalb(ifprt,'chi ',chi,kkbeta,nbeta,lll,r,rab,dum1,
     +   idim1,idim3)
c
c---------------------------------------------------------------------------
c
      call prmat(bbb,nbeta,'b',idim3)
      call prmat(qqq,nbeta,'q',idim3)
      call prmat(ddd,nbeta,'d',idim3)
CCCCCCCCCCCCCCCCCCCCCCCCCC
C     call prmat(sss,nbeta,'s',idim3)
C     call prmat(ttt,nbeta,'t',idim3)
CCCCCCCCCCCCCCCCCCCCCCCCCC
c
c     symmetrize ddd
c
      do 380 ib = 1,nbeta
        do 390 jb = ib+1,nbeta
          xx = 0.5d0 * ( ddd(ib,jb) + ddd(jb,ib) )
          ddd(ib,jb) = xx
          ddd(jb,ib) = xx
  390   continue
  380 continue
c
      call prmat(ddd,nbeta,'d',idim3)
c
c     -------
c     num rec - with my driver nrinv (see nrsubs)
c     -------
      call nrinv(bbb,idim3,nbeta,bbbi,idim3,iwork)
c     -------
c
      call prmat(bbbi,nbeta,'i',idim3)
c
      write (stderr,*) 'construction of qqq and ddd complete'
c
c---------------------------------------------------------------------------
c
c              f o r m  t h e  b e t a  f u n c t i o n s
c
      do 400 ibeta=1,nbeta
c
c       zeroing
        do 410 i=1,mesh
          beta(i,ibeta) = 0.0d0
  410   continue
c
        do 400 jbeta=1,nbeta
          binv = bbbi(jbeta,ibeta)
          do 400 i=1,kkbeta
            beta(i,ibeta)=beta(i,ibeta)+binv*chi(i,jbeta)
  400 continue
c
c---------------------------------------------------------------------------
c
c     p r i n t  /  p l o t  t h e  b e t a  f u n c t i o n s
c
      call prbeta(ifprt,'beta',beta,kkbeta,nbeta,lll,r,a,b,idim1,idim3)
c     print and plot fourier transform of beta functions
      call fanalb(ifprt,'beta',beta,kkbeta,nbeta,lll,r,rab,dum1,
     +   idim1,idim3)
c
c-----------------------------------------------------------------------
c
c     p o s i t i v e  p s e u d i z a t i o n  o f  q f u n c
c
      if ( (ifqopt.eq.2 .or. ifqopt.eq.3) .and. ifprt .ge. 0) then
        write(iout,'(/3x,43(1h-))')
        write(iout,'(3x,a)')
     +   'pseudize q-functions'
        write(iout,'(3x,43(1h-))')
        write(iout,*)
      endif
c
c     check that overlaps are as expected
      do 620 lp = 1,nang
        istart = nbl0(lp) + 1
        istop = istart + nbl(lp) - 1
        do 630 ibeta = istart,istop
        do 630 jbeta = istart,istop
          do 640 ir = 1,mesh
            dum1(ir) = phi(ir,ibeta) * beta(ir,jbeta)
  640     continue
          call radin(mesh,dum1,rab,asum,idim1)
c          write(iout,*) 'ibeta,jbeta',ibeta,jbeta,' asum',asum
  630   continue
  620 continue
c
c     check that rinner is suitable for pseudization
      if ( rinner(1) .gt. r(3) ) then
        do 600 ib = 1,nbeta
        write(stderr,'(a,i2)') ' pseudizing qfuncs for beta=',ib
        do 600 jb = 1,nbeta
           if (iptype(ib) .eq. 3 .and. iptype(jb) .eq. 3 ) then
              do i=1,idim8
                 do j=1,idim6
                    qfcoef(i,j,ib,jb) = 0.0
                 enddo
              enddo
           else
          lmin = iabs ( lll(ib) - lll(jb) )
          lmax = iabs ( lll(ib) + lll(jb) )
          if ( ifqopt.eq.2 .or. ifqopt.eq.3 ) then
            call psqf2(qfunc(1,ib,jb),qfcoef(1,1,ib,jb),dum1,lll,
     +        rinner,mesh,r,rab,a,b,kkbeta,ifprt,nang,
     +        nqf,qtryc,ibfix,rfix,qfix,nfix,ifqopt,
     +        phi(1,ib),phi(1,jb),ib,jb,idim1,idim3,idim6,idim8)
c
c           
c           ======================================================
c           test whether there might be a negative density problem
c           ======================================================
c
            if ( ib .eq. jb ) then
              ineg = 0
              idom = 0
              dneg = 0.0d0
              rhocur = 1.0d0
              do 610 ir = 2,kkbeta
                 if ( r(ir) .gt. rinner(2*nang-1) ) goto 610
                 rhoold = rhocur
                 rhocur = qfunc(ir,ib,jb) + phi(ir,ib)*phi(ir,jb)
                 if ( rhocur .lt. 0.0d0 ) then
c
c                  another point of negative density has been found
                   ineg = ineg + 1
c
c                  is this a new domain of negative density?
                   if ( rhocur * rhoold .lt. 0.0d0 ) idom = idom + 1
c
c                  is this the most negative desnsity so far encountered?
                   if ( rhocur .lt. dneg ) then
                     iworst = ir
                     dneg = rhocur
                   endif
                 endif
  610         continue
c
              if ( ifprt .ge. 2) write (iout,'(5x,a,i5,a,2i3/5x,a,i5)')
     +          'number of negative densities = ',ineg,
     +          ' for ib,jb = ',ib,jb,
     +          'number of domains of negative density = ',idom
              if ( ineg .gt. 0 ) then
                write (iout,'(5x,a,i5,a,f9.5,a,f9.5)')
     +            'most negative density at mesh point',iworst,
     +            ' at ',r(iworst),' au with value ',dneg
                write(iout,'(5x,a,f9.5)') 
     +            'phi(iworst,ib)*phi(iworst,jb) =',
     +            phi(iworst,ib)*phi(iworst,jb)
              endif
            endif
          endif
       endif
  600   continue
      endif
c
      if ( ifqopt.eq.2 .or. ifqopt.eq.3 )
     +    write(stderr,*) 'qfunc pseudized'
c
c----------------------------------------------------------------------
c
c     t r a n s f o r m  b e t a  s o  o r t h o g o n a l
c
      if ( ifprt .ge. 0) then
        write(iout,'(/3x,43(1h-))')
        write(iout,'(3x,a)')
     +   'transform so betas are orthogonal'
        write(iout,'(3x,43(1h-))')
        write(iout,*)
      endif
c
c     do it block by block so that blocks do not get rearranged
c
      do 470 lp=1,nang
c
        nn=nbl(lp)
        ns=nbl0(lp)+1
        if (nn.le.0) go to 470
c
c       transform so that beta are orthonormal and kinetic energy
c       operator is diagonal
c       on input,  ttt carries info about mat el of kin energy op
c       on output, ttt carries info about linear transformation
        call qdiag(nn,bbbi(ns,ns),qqq(ns,ns),sss(ns,ns),ttt(ns,ns),
     +    ddd(ns,ns),beta(1,ns),kkbeta,idim1,idim3)
c
  470 continue
c
      write (iout,*) ' after transformation :'
      call prmat(qqq,nbeta,'q',idim3)
      call prmat(ddd,nbeta,'d',idim3)
c
c---------------------------------------------------------------------------
c
c     p r i n t  /  p l o t  t h e  b e t a  f u n c t i o n s
c
      call prbeta(ifprt,'beta',beta,kkbeta,nbeta,lll,r,a,b,idim1,idim3)
c     print and plot fourier transform of beta functions
      call fanalb(ifprt,'beta',beta,kkbeta,nbeta,lll,r,rab,dum1,
     +   idim1,idim3)
c
c---------------------------------------------------------------------------
c
      write (stderr,*) 'transformation of q and d complete'
c
c     call prmat(ttt,nbeta,'i',idim3)
      do 490 i=1,kkbeta
        call abat(nbeta,idim1,1,ttt,qfunc(i,1,1),idim3)
  490 continue
c
      write (stderr,*) 'transformation of qfunc complete'
c
c     if we have already pseudized q_ij then we must transform qfcoef
      if ( ifqopt.eq.2 .or. ifqopt.eq.3 ) then
c
        do 495 nq = 1,nqf
        do 495 ltp = 1,2*nang-1
          call abat(nbeta,idim8*idim6,1,ttt,qfcoef(nq,ltp,1,1),idim3)
  495   continue
c
        write(stderr,*) 'transformation of qfcoef complete'
c
      endif
c
c---------------------------------------------------------------------------
c
c     p o s s i b l e  p s e u d i z a t i o n  o f  q f u n c
c
      if ( (ifqopt.eq.0 .or. ifqopt.eq.1) .and. ifprt .ge. 0) then
        write(iout,'(/3x,43(1h-))')
        write(iout,'(3x,a)')
     +   'pseudize q-functions'
        write(iout,'(3x,43(1h-))')
        write(iout,*)
      endif
c
c     check that rinner is suitable for pseudization
      if ( rinner(1) .gt. r(3) ) then
        do 500 ib = 1,nbeta
        do 500 jb = 1,nbeta
          lmin = iabs ( lll(ib) - lll(jb) )
          if ( ifqopt .eq. 0 ) then
            call psqf(qfunc(1,ib,jb),qfcoef(1,1,ib,jb),dum1,lmin,
     +        rinner(1),r,rab,a,b,kkbeta,ifprt,nang,idim1,idim6,idim8)
          elseif ( ifqopt .eq. 1 ) then
            call psqf1(qfunc(1,ib,jb),qfcoef(1,1,ib,jb),dum1,lmin,
     +        rinner(1),mesh,r,rab,a,b,kkbeta,ifprt,nang,
     +        nqf,qtryc,idim1,idim6,idim8)
          endif
  500   continue
      endif
c
      if ( ifqopt .eq. 0 .or. ifqopt .eq. 1 ) write(stderr,*) 
     +  'qfunc pseudized'
c
      write(iout,'(/1x,47(1h-))')
      write(iout,'(1x,a)') 'end subroutine scpgef'
      write(iout,'(1x,47(1h-)/)')
c
      return
      end
c
c---------------------------------------------------------------------------
      subroutine difinv(g,p,h,jj1,jj2,bndmat,npivot,idim1)
c---------------------------------------------------------------------------
c     routine solves for h(x) in differential equation of form
c
c         2
c        d p(x)
c        ------  = g(x) p(x) + h(x)
c           2
c         dx
c
c     by inverting the method due to noumerov.
c     method taken from g.w. pratt, phys. rev. 88 p.1217 (1952).
c     see also w.w. piper phys. rev. 123 p.1281 (1961).
c
c     routine integrates from jj1 to jj2 and can cope with both cases
c     jj1 < jj2 and jj1 > jj2.
c
c     *******************************************************
c     n.b.  values of h(min(jj1,jj2)-1) and h(max(jj1,jj2)+1)
c     and   values og g(min(jj1.jj2)-1) and g(max(jj1,jj2)+1)
c     must be set in calling program.
c     *******************************************************
c---------------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
c
c.....differential equation arrays
      dimension g(idim1),p(idim1),h(idim1)
c.....matrix inversion arrays
      dimension bndmat(idim1,5),npivot(idim1)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c---------------------------------------------------------------------------
c
c     i n i t i a l i s a t i o n
c
      r12 = 1.0d0 / 12.0d0
c
c     inversion limits
      kk1 = min(jj1,jj2)
      kk2 = max(jj1,jj2)
c
c     run some test just to be conservative
      if ( kk1 .lt. 2 .or. kk2 .gt. idim1 - 1 ) then
        write(iout,10) kk1,kk2,idim1
        call exit(1)
      endif
c
  10  format(' ***error in subroutine difinv',/,
     +' kk1 =',i5,' kk2 =',i5,' idim1 =',i5,
     +' are not allowed')
c
c     set up for the inversion of the tridiagonal matrix
      emach = 1.0d-12
      diag = 10.0d0 / 12.0d0
      npts = kk2 - kk1 + 1
      do 20 i = 1,npts
        bndmat(i,1) = r12
        bndmat(i,2) = diag
        bndmat(i,3) = r12
        bndmat(i,4) = 0.0d0
        bndmat(i,5) = 0.0d0
   20 continue
      bndmat(1,1)    = 0.0d0
      bndmat(npts,3) = 0.0d0
c
c     factorise the bndmat matrix using tcmlib routine
      call bndglm(bndmat,idim1,npts,5,npivot,emach)
c
c-------------------------------------------------------------------------
c
c     b e g i n  i n v e r s i o n  l o o p
c
c     initialise the starting values of z
      i = kk1 - 1
      zold = p(i) - r12 * p(i) * g(i)
      i = kk1
      zcur = p(i) - r12 * p(i) * g(i)
c
      do 100 i = kk1,kk2
        znew = p(i+1) - r12 * p(i+1) * g(i+1)
        h(i) = znew - 2.0d0 * zcur + zold - p(i) * g(i)
        zold = zcur
        zcur = znew
  100 continue
c
c-------------------------------------------------------------------------
c
c     s o l u t i o n  o f  s i m u l t a n e o u s  e q u a t i o n s
c
      h(kk1) = h(kk1) - r12 * h(kk1-1)
      h(kk2) = h(kk2) - r12 * h(kk2+1)
c
c     tcm routine solves equation
      call bndsub(bndmat,idim1,npts,5,npivot,h(kk1))
c
c
      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine schgef(ecur,lcur,b,r,rab,sqr,vloc,v,
     +  snlo,mesh,imatch,dlwf,beta,nbl,nbl0,gamma,alpha,
     +  eta,kkbeta,ddd,qqq,idim1,idim3,idim5,idim7)
c
c     analogue of subroutine phase for the generalised
c     eigenvalue parts of the formulation
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1),sqr(idim1)
c.....wavefunctions
      dimension snlo(idim1)
c.....beta and q functions and q pseudization coefficients
      dimension beta(idim1,idim3)
c.....angular momenta, reference energy, qij and dij of vanderbilt scheme
      dimension qqq(idim3,idim3),ddd(idim3,idim3)
c.....work space for solution to inhomogeneous equation and overlaps
      dimension eta(idim1,idim7),alpha(idim3)
c.....indexing for the beta functions
      dimension nbl0(idim5),nbl(idim5)
c.....local component of the pseudopotential
      dimension vloc(idim1),v(idim1)
c.....scratch arrays
      dimension gamma(idim1)
      dimension gg(4,4),ff(4),coef(4)
c
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c--------------------------------------------------------------------------
c
c           i n i t i a l i s a t i o n
c
      imloc = imatch
      if ( imloc - 2 * ( imloc / 2 ) .ne. 1 ) imloc = imloc + 1
      if ( imloc .gt. mesh ) then
        write(iout,*) '***error in subroutine schgef'
        write(iout,*) 'imloc =',imloc,' .gt. mesh =',mesh
        call exit(1)
      endif
c
      lp  = lcur + 1
      llp = lcur * ( lcur + 1 )
c
c     routine hard wired to deal only with nbl(lp) .le. 4
      if (lp.le.idim5) then
      if ( nbl(lp) .gt. 4 ) then
        write(iout,*) '***error in subroutine schgef'
        write(iout,*) 'nbl(',lp,') =',nbl(lp),' is too large'
        call exit(1)
      endif
      endif
c
      ifprt = 0
c
c---------------------------------------------------------------------------
c
c     s o l v e  t h e  h o m o g e n e o u s  e q u a t i o n
c
c     compute the potential for this angular momentum and energy
      call setv(r,rab,vloc,b,ecur,llp,v,mesh,idim1)
c
c     analytic starting value for snlo using taylor expansion
      zz = 0.0d0
      vzero = vloc(2) / r(2)
      call orgsnl(r,snlo,1,3,lcur,zz,vzero,ecur,idim1)
c
c     convert to the phi function
      snlo(1) = snlo(1) / sqr(1)
      snlo(2) = snlo(2) / sqr(2)
      snlo(3) = snlo(3) / sqr(3)
c
c     integrate the schroedinger equation outwards
      call difsol(v,snlo,4,imloc,idim1)
c      
c     check whether generalised eigenvalue parts required
      if ( lp .gt. idim5 ) go to 400
      if ( nbl(lp) .le. 0 ) goto 400
c
c---------------------------------------------------------------------------
c
c     s o l v e  t h e  i m h o m o g e n e o u s  e q u a t i o n
c
      do 100 ib = 1,nbl(lp)
c
        ibeta = nbl0(lp) + ib
c
c       construct gamma
        do 110 i = 1,imloc
          gamma(i) = 0.0d0
  110   continue
c
        do 120 jb = 1,nbl(lp)
          jbeta = nbl0(lp) + jb
          con = ddd(jbeta,ibeta) - ecur * qqq(jbeta,ibeta)
          do 130 i = 1,imloc
            gamma(i) = gamma(i) - con * beta(i,jbeta)
  130     continue
  120   continue
c
c       generate a starting guess for the eta function
c       zero value and slope at the origin seems to work fine (dv 7/98)
        do 140 i = 1,3
          eta(i,ib) = 0.d0
  140   continue
c
c       transform gamma into form appropriate for logarithmic mesh
        do 150 i = 1,imloc
          gamma(i) = - gamma(i) * rab(i) * rab(i) / sqr(i)
  150   continue
c
c       solve the inhomogeneous equation
        call difsli(v,eta(1,ib),gamma,4,imloc,idim1)
c
c       form the gg and ff arrays
        do 160 jb = 1,nbl(lp)
          jbeta = nbl0(lp) + jb
c
          do 170 i = 1,imloc
            gamma(i) = beta(i,jbeta) * eta(i,ib) * sqr(i)
  170     continue
c
          call radin(imloc,gamma,rab,gg(jb,ib),idim1)
c
  160   continue
c
        do 180 i = 1,imloc
          gamma(i) = beta(i,ibeta) * snlo(i) * sqr(i)
  180   continue
c
        call radin(imloc,gamma,rab,ff(ib),idim1)
c
  100 continue
c
c     write (iout,190) ((ib,jb,gg(ib,jb),jb=1,nbl(lp)),ib=1,nbl(lp))
c 190 format ('gg ='/(2i3,e20.12))
c
c-------------------------------------------------------------------------
c
c     i n v e r t  t h e  g g  m a t r i x  a n d  f o r m  s n l o       
c
c     solve system of linear equations
      do 200 ib = 1,nbl(lp)
      do 200 jb = 1,nbl(lp)
        gg(ib,jb) = - gg(ib,jb)
  200 continue
      do 210 ib = 1,nbl(lp)
        gg(ib,ib) = 1.0d0 + gg(ib,ib)
  210 continue
c
c     numerical recipies to invert the equation
      call simeq(4,nbl(lp),gg,ff,coef)
c
c     write (iout,*) ' coef=',(coef(ib),ib=1,nbl(lp))
c
c     reconstruct wavefunction and load the overlap vector
      do 250 ib = 1,nbl(lp)
        ibeta = nbl0(lp) + ib
        alpha(ibeta) = coef(ib)
        do 250 i  = 1,imloc
          snlo(i) = snlo(i) + coef(ib) * eta(i,ib)
 250  continue
c
c-------------------------------------------------------------------------
c
c     c o m p u t e  t h e  l o g  d e r i v a t i v e s
c
  400 continue
c
      snlom = snlo(imatch-1) * sqr(imatch-1)
      snlop = snlo(imatch)   * sqr(imatch)
c
      dwf =   snlop - snlom
      wf  = ( snlop + snlom ) * 0.5d0
      dr  =   r(imatch) - r(imatch-1)
ctest
c     write(6,*)'dwf,dr,wf',dwf,dr,wf
      dlwf = dwf / ( dr * wf )
c
      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine difsli(g,p,h,jj1,jj2,idim1)
c
c---------------------------------------------------------------------------
c     routine for solving inhomogeneous second order differential equation
c
c         2
c        d p(x)
c        ------  = g(x) p(x) + h(x)
c           2
c         dx
c
c     using the method due to noumerov.
c     method taken from g.w. pratt, phys. rev. 88 p.1217 (1952).
c     see also w.w. piper phys. rev. 123 1281 (1961).
c
c     routine integrates from jj1 to jj2 and can cope with both cases
c     jj1 < jj2 and jj1 > jj2.  first two starting values of p must be
c     provided by the calling program.
c---------------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
c
      dimension g(idim1),p(idim1),h(idim1)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c     i n i t i a l i s a t i o n
c
      r12 = 1.0d0 / 12.0d0
c
c     decide whether integrating from:
c            left to right ---> isgn = + 1
c        or  right to left ---> isgn = - 1
c
      isgn = ( jj2 - jj1 ) / iabs( jj2 - jj1 )
c
c     run some test just to be conservative
      if ( isgn .eq. + 1 ) then
        if ( jj1 .le. 2 .or. jj2 .gt. idim1 ) then
          write(iout,10) isgn,jj1,jj2,idim1
          call exit(1)
        endif
      elseif ( isgn .eq. - 1 ) then
        if ( jj1 .ge. ( idim1 - 1 ) .or. jj2 .lt. 1 ) then
          write(iout,10) isgn,jj1,jj2,idim1
          call exit(1)
        endif
      else
        write(iout,10) isgn,jj1,jj2,idim1
      endif
c
  10  format(' ***error in subroutine difsli',/,
     +' isgn =',i2,' jj1 =',i5,' jj2 =',i5,' idim1 =',i5,
     +' are not allowed')
c
c     initialise current second derivative of p (d2p), ycur and yold
c
      i = jj1 - isgn
      d2p = g(i) * p(i) + h(i)
      ycur = p(i) - r12 * ( g(i) * p(i) + h(i) )
      i = i - isgn
      yold = p(i) - r12 * ( g(i) * p(i) + h(i) )
c
c
c     b e g i n  i n t e g r a t i o n  l o o p
c
      do 100 i = jj1,jj2,isgn
        ynew = 2.0d0 * ycur - yold + d2p
        p(i) = ( ynew + r12 * h(i) ) / ( 1.0d0 - r12 * g(i) )
        d2p  = g(i) * p(i) + h(i)
        yold = ycur
        ycur = ynew
100   continue
c
      return
      end
c
c---------------------------------------------------------------------------
c
c     ***************************************************************
      subroutine qdiag(nbeta,bbbi,qqq,sss,ttt,ddd,beta,kkbeta,
     +  idim1,idim3)
c     ***************************************************************
c     transform beta's to be orthogonal and to diagonalize the
c     kinetic energy matrix ttt. This tends to order beta functions
c     according to the number of nodes.
c ---------------------------------------------------------------
      implicit double precision (a-h,o-z)
c ---------------------------------------------------------------
c
      parameter ( idim3l = 10 )
c
      dimension bbbi(idim3,idim3),qqq(idim3,idim3),sss(idim3,idim3),
     +  ddd(idim3,idim3),ttt(idim3,idim3)
      dimension beta(idim1,idim3)
      dimension uuu(idim3l,idim3l),uuui(idim3l,idim3l)
      dimension r(idim3l),w1(idim3l),w2(idim3l),iwork(idim3l)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c-----------------------------------------------------------------------
c
c                 i n i t i a l i z a t i o n
c
c     check that the local idim3l and global idim3 are identical
      if ( idim3l .ne. idim3) then
        write(iout,*) '***error in subroutine qdiag'
        write(iout,*) 'idim3l=',idim3l,' .ne. idim3 =',idim3
        call exit(1)
      endif
c
c     check that real action is required
      if ( nbeta .le. 0 ) return
c
c-----------------------------------------------------------------------
c
c         t r a n s f o r m  t h e  b e t a  f u n c t i o n s
c
c     recall sss is < chi_i | chi_j >
      call abat(nbeta,1,-1,bbbi,sss,idim3)
c
c     puts bbb_inv_tr dot sss dot bbb_inv onto sss
c     sss is < beta_i | beta_j >
c
c     recall ttt is < chi_i | t | chi_j >
      call abat(nbeta,1,-1,bbbi,ttt,idim3)
c     puts bbb_inv_tr dot ttt dot bbb_inv onto ttt
c     ttt is < beta_i | t | beta_j >
c
c     call prmat(sss,nbeta,'s',idim3)
c
c     solve gen eigenvalue problem (ttt - lambda_k sss) uuu_k = 0
c     -------
c     nag lib
c     -------
c     ifail=0
c     call f02aef(ttt,idim3,sss,idim3,nbeta,r,uuu,idim3,w1,w2,ifail)
c     if ( ifail .ne. 0 ) then
c       write(iout,*) '***error in subroutine qdiag'
c       write(iout,*) 'naglib ifail =',ifail,' .ne. 0'
c       call exit(1)
c     endif
c     -------
c     eispack
c     -------
      ifvec=1
      call rsg(idim3,nbeta,ttt,sss,r,ifvec,uuu,w1,w2,ierr)
      if ( ierr .ne. 0 ) then
        write(iout,*) '***error in subroutine qdiag'
        write(iout,*) 'eispack ierr =',ierr,' .ne. 0'
        call exit(1)
      endif
c     -------
c     note sss and ttt are corrupted
c
c     transform beta's
      do 150 i=1,kkbeta
      do 140 ib=1,nbeta
      w1(ib)=0.
      do 140 jb=1,nbeta
  140 w1(ib)=w1(ib)+uuu(jb,ib)*beta(i,jb)
      do 150 ib=1,nbeta
  150 beta(i,ib)=w1(ib)
c
c     call prmat(uuu,nbeta,'u',idim3)
c
c     invert uuu
c     -------
c     num rec - with my driver nrinv (see nrsubs)
c     -------
      call nrinv(uuu,idim3,nbeta,uuui,idim3,iwork)
c     -------
c
c     call prmat(uuui,nbeta,'i',idim3)
c
c     transform matrices
      call abat(nbeta,1,1,uuui,qqq,idim3)
      call abat(nbeta,1,1,uuui,ddd,idim3)
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c skip ghost analysis here, since now done at very end
c dv 7/29/97
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c c     get "ghost eigenvalues," defined to be eigenvalues of
c c       non-local piece alone.
c c     if one 'ghost eigenvalues' is very large and negative, this
c c       may indicates the existence of a spurious solution (ghost).
c c     if one 'ghost eigenvalues' is very large and positive,
c c       iterative methods may converge poorly due to inadequate
c c       preconditioning
c       do 162 i=1,nbeta
c       do 160 j=1,nbeta
c       ttt(i,j)=ddd(i,j)
c   160 sss(i,j)=qqq(i,j)
c   162 sss(i,i)=sss(i,i)+1.
c c     solve gen eigenvalue problem (ttt - lambda_k sss) uuu_k = 0
c c     -------
c c     nag lib
c c     -------
c c     ifail=0
c c     call f02aef(ttt,idim3,sss,idim3,nbeta,r,uuu,idim3,w1,w2,ifail)
c c     if ( ifail .ne. 0 ) then
c c       write(iout,*) '***error in subroutine qdiag'
c c       write(iout,*) 'naglib ifail =',ifail,' .ne. 0'
c c       call exit(1)
c c     endif
c c     -------
c c     eispack
c c     -------
c       ifvec=0
c       call rsg(idim3,nbeta,ttt,sss,r,ifvec,uuu,w1,w2,ierr)
c       if ( ierr .ne. 0 ) then
c         write(iout,*) '***error in subroutine qdiag'
c         write(iout,*) 'eispack ierr =',ierr,' .ne. 0'
c         call exit(1)
c       endif
c c     -------
c c
c       write (iout,910) (r(i),i=1,nbeta)
c   910 format (' ghost eigenvalues : ',4f12.5)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c     save uuui into ttt for use in transforming qfunc
      do 180 i=1,nbeta
      do 180 j=1,nbeta
  180 ttt(i,j)=uuui(i,j)
c
      return
      end
c
c
c     ***************************************************************
      subroutine qdiag1(lp,nbeta,bbbi,qqq,sss,ttt,ddd,beta,kkbeta,
     +  idim1,idim3)
c     ***************************************************************
c     find the energy eigenvalues of the nonlocal part alone
c ---------------------------------------------------------------
      implicit double precision (a-h,o-z)
c ---------------------------------------------------------------
c
      parameter ( idim3l = 10 )
c
      dimension bbbi(idim3,idim3),qqq(idim3,idim3),sss(idim3,idim3),
     +  ddd(idim3,idim3),ttt(idim3,idim3)
      dimension beta(idim1,idim3)
      dimension uuu(idim3l,idim3l)
      dimension r(idim3l),w1(idim3l),w2(idim3l)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c-----------------------------------------------------------------------
c
c                 i n i t i a l i z a t i o n
c
c     check that the local idim3l and global idim3 are identical
      if ( idim3l .ne. idim3) then
        write(iout,*) '***error in subroutine qdiag1'
        write(iout,*) 'idim3l=',idim3l,' .ne. idim3 =',idim3
        call exit(1)
      endif
c
c     check that real action is required
      if ( nbeta .le. 0 ) return
c
c
c     get "ghost eigenvalues," defined to be eigenvalues of
c       non-local piece alone.
c     if one 'ghost eigenvalues' is very large and negative, this
c       may indicates the existence of a spurious solution (ghost).
c     if one 'ghost eigenvalues' is very large and positive,
c       iterative methods may converge poorly due to inadequate
c       preconditioning
      do 162 i=1,nbeta
      do 160 j=1,nbeta
      ttt(i,j)=ddd(i,j)
  160 sss(i,j)=qqq(i,j)
  162 sss(i,i)=sss(i,i)+1.
c     solve gen eigenvalue problem (ttt - lambda_k sss) uuu_k = 0
c     -------
c     nag lib
c     -------
c      ifail=0
c      call f02aef(ttt,idim3,sss,idim3,nbeta,r,uuu,idim3,w1,w2,ifail)
c      if ( ifail .ne. 0 ) then
c        write(iout,*) '***error in subroutine qdiag1'
c        write(iout,*) 'naglib ifail =',ifail,' .ne. 0'
c        call exit(1)
c      endif
c     -------
c     eispack
c     -------
      ifvec=0
      call rsg(idim3,nbeta,ttt,sss,r,ifvec,uuu,w1,w2,ierr)
      if ( ierr .ne. 0 ) then
       write(iout,*) '***error in subroutine qdiag1'
       write(iout,*) 'eispack ierr =',ierr,' .ne. 0'
      call exit(1)
      endif
c     -------
c
      write(iout,910) lp-1,(r(i),i=1,nbeta)
  910     format(3x,'l=',i2,4x,4f11.6)
c
      do 180 i=1,nbeta
  180 ttt(i,1)=r(i)
c
      return
      end
c
c     ***************************************************************
      subroutine abat(n,idim,key,a,b,idim3)
c     ***************************************************************
c
c     if key.gt.0 set b equal to a . b . a-transpose
c     if key.lt.0 set b equal to a-transpose . b . a
c
c ---------------------------------------------------------------
      implicit double precision (a-h,o-z)
c ---------------------------------------------------------------
c
      parameter ( idim3l = 10 )
c
c.....routine inputs
      dimension a(idim3,idim3),b(idim,idim3,idim3)
c.....scratch arrays
      dimension t(idim3l,idim3l)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c     check that the local idim3l and global idim3 are identical
      if ( idim3l .ne. idim3) then
        write(iout,*) '***error in subroutine qdiag'
        write(iout,*) 'idim3l=',idim3l,' .ne. idim3 =',idim3
        call exit(1)
      endif
c
      do 40 i=1,n
      do 40 k=1,n
      sum=0.
      if (key.gt.0) then
        do 20 j=1,n
   20   sum=sum+a(i,j)*b(1,j,k)
      else
        do 30 j=1,n
   30   sum=sum+a(j,i)*b(1,j,k)
      endif
   40 t(i,k)=sum
      do 80 i=1,n
      do 80 k=1,n
      sum=0.
      if (key.gt.0) then
        do 60 j=1,n
   60   sum=sum+t(i,j)*a(k,j)
      else
        do 70 j=1,n
   70   sum=sum+t(i,j)*a(j,k)
      endif
   80 b(1,i,k)=sum
      return
      end
c
c-------------------------------------------------------------------------
c
      subroutine scheqg(vloc,r,rab,sqr,a,b,v,z,ee,rcloc,rc,
     +  nnlz,wwnl,nkkk,snl,snlo,ncores,ncspvs,mesh,
     +  thresh,beta,nbl,nbl0,gamma,alpha,eta,klogd,kkbeta,
     +  ddd,qqq,ifprt,lsweep,idim1,idim2,idim3,idim5,idim7,icflag)
c
c     routine to compute non relativistic solutions to
c     schroedinger equation
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1),sqr(idim1)
c.....potentials in the calculation
      dimension vloc(idim1),v(idim1)
c.....wavefunctions
      dimension snl(idim1,idim2),snlo(idim1)
c.....shell labels, energies, occupancies and real space cutoff
      dimension nnlz(idim2),ee(idim2),wwnl(idim2),nkkk(idim2)
c.....beta functions
      dimension beta(idim1,idim3)
c.....core cut of radii
      dimension rc(idim5)
c.....work space for solution to inhomogeneous equation and overlaps
      dimension eta(idim1,idim7),alpha(idim3,idim2)
c.....qij and dij of vanderbilt scheme
      dimension qqq(idim3,idim3),ddd(idim3,idim3)
c.....indexing for the beta functions
      dimension nbl0(idim5),nbl(idim5)
c.....logic to determine whether to scan energies
      logical lsweep
c.....scratch
      dimension yval(-1:1),gamma(idim1)
      logical lscan
c
c.....variable for file = 0
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c--------------------------------------------------------------------------
c
c               r o u t i n e  i n i t i a l i s a t i o n
c
      if ( ifprt .ge. 1) write (iout,'(/1x,a)')
     +     'subroutine scheqg: solve schroedinger equation'
c
c     set the maximum number of iterations for improving wavefunctions
      itmax = 100
c
      do 10 j = ncores+1,ncspvs
      do 10 i = 1,mesh
        snl(i,j) = 0.0d0
  10  continue
c
c--------------------------------------------------------------------------
c
c               m a i n  l o o p  o v e r  o r b i t a l s
c
      do 100 iorb = ncores+1,ncspvs
c
c       check that this orbital is required
        if ( icflag .gt. 0 .and. wwnl(iorb) .le. 0.0d0 ) then
          goto 100
        endif
c
c       compute the quantum numbers
        ncur = nnlz(iorb) / 100
        lcur = nnlz(iorb) / 10 - ncur * 10
        lp   = lcur + 1
c
c       set the starting energies for the iterative loop
        ecur = ee(iorb)
        decurold = -1.0d12
c
c       possible scan of energy range to prevent instabities on
c       initial iterations in ionic tests
        if ( lsweep ) then
          lscan = .true.
          ecur  = 2.0d0 * ecur
        else
          lscan = .false.
        endif
c
c       effective infinity
        ninf = nkkk(iorb)
c       set nctp to a point beyond the range of beta
        rcmax = dmax1(rcloc,rc(lp))
        nctp = 4 + int( dlog ( rcmax / a + 1.0d0 ) / b )
c
        if ( ifprt .ge. 1) write(iout,'(4x,a,3i5)')
     +    'nnlz,nctp,ninf =',nnlz(iorb),nctp,ninf
c
        if ( nctp .gt. ninf - 2 ) then
          write(iout,*) '***error in subroutine scheqg'
          write(iout,*) 'nctp,ninf=',nctp,ninf,' unreasonable'
          call exit(1)
        endif
c
c-------------------------------------------------------------------------
c
c            i t e r a t i o n  o n  e i g e n f u n c t i o n s
c
        do 200 iter = 1,itmax
c
c               o u t w a r d  i n t e g r a t i o n 
c
          call schgef(ecur,lcur,b,r,rab,sqr,vloc,v,
     +      snlo,mesh,nctp,dlwf,beta,nbl,nbl0,gamma,alpha(1,iorb),
     +      eta,kkbeta,ddd,qqq,idim1,idim3,idim5,idim7)
c
          snlctp = snlo(nctp)
c
c               i n w a r d  i n t e g r a t i o n
c
c         exponential start up of the wavefunction
          dum = sqrt(v(ninf))
          snlo(ninf) = exp( - dum * dble(ninf - nctp ) )
          dum = sqrt(v(ninf-1))
          snlo(ninf-1) = exp ( - dum * dble(ninf - nctp - 1) )
c
          call difsol(v,snlo,ninf-2,nctp,idim1)
c
c--------------------------------------------------------------------------
c
c         c h e c k  t h a t  t h e  n u m b e r  o f  n o d e s  o k
c
c         renormalise snlo so continuous at nctp
          factor = snlctp / snlo(nctp)
          do 350 i = nctp,ninf
            snlo(i) = snlo(i) * factor
  350     continue
c
c         count the number of nodes
          call nodeno(snlo,2,ninf,nodes,idim1)
c
c--------------------
c
c         print out of results
c         write(iout,500) lcur
c 500     format(//,'wavefunction and potential for l =',i2,/,
c    +      5x,'r',9x,'psi',10x,'v')
c         nind = 55
c         rinc = r(ninf) / dble(nind)
c         rcur = -0.75d0 * rinc
c         do 510 i = 1,nind
c           rcur = rcur + rinc
c           in = int( dlog( rcur / a + 1.0d0 ) / b ) + 1
c           write(iout,520) r(in),snlo(in),v(in)
c 510     continue
c 520     format(f9.4,3f13.8)
c
c-------------------
c
c         pseudo case so the number of nodes should be zero
c         if ( nodes .ne. 0 ) then
c           if ( iter .lt. itmax ) then
c             try sweeping down to lower energy
c             decur = - ecur * 0.04d0
c             write(iout,*) 'nodes=',nodes,' ecur,decur=',ecur,decur
c             ecur = ecur + decur
c             goto 200
c           else
c             maximum number of iterations and not nodeless so stop
c             write(iout,*) '***error in subroutine scheqg'
c             write(iout,*) 'nodes =',nodes,' .ne. 0'
c             call exit(1)
c           endif
c         endif
c
c---------------------------------------------------------------------------
c
c         v a r i a t i o n a l  i m p ro v e m e n t  o f  e c u r
c
c         find normalisation of psi using simpson's rule
          factor = 0.0d0
          xw = 4.0d0
          do 360 i = 2,ninf
            factor = factor + xw * rab(i) * (snlo(i)*sqr(i))**2.0
            xw = 6.0d0 - xw
  360     continue
          factor = factor / 3.0d0
c
c         do the generalised part of the normalisation
          do 365 ib = 1,nbl(lp)
          do 365 jb = 1,nbl(lp)
            ibeta = ib + nbl0(lp)
            jbeta = jb + nbl0(lp)
            factor = factor + qqq(ibeta,jbeta) * 
     +               alpha(ibeta,iorb) * alpha(jbeta,iorb)
  365    continue
c
c         compute the noumerov discontinuity in the phi function
          do 370 i = -1,1
            yval(i) = snlo(nctp+i) * ( 1.0d0 - v(nctp+i) / 12.0d0 )
  370     continue
          d2p = v(nctp) * snlctp
          disc = - yval(-1) + 2.0d0*yval(0) - yval(1) + d2p
c
c         variational estimate for the change in ecur
          decur = snlctp*disc*sqr(nctp)**2.0 / (factor*rab(nctp))
c change by Kurt Stokbro
c          write(*,*) decur,decurold
          if (decur/decurold .gt. 2.0) then
             write(iout,*) 'WARNING correcting decur'
             decur = -1.5*decurold
          endif
          decurold = decur
c
c
          if ( ifprt .ge. 1) write (iout,'(6x,a,f14.8,e12.4)')
     +       'ecur,decur',ecur,decur
          if ( abs(decur) .lt. thresh ) goto 400
c
c         possible modification of decur to prevent instabilities
          if (5.0d0*abs(decur).gt.-ecur) decur = sign(ecur,decur)/5.0d0
c
          if ( lscan ) then
            if ( decur .lt. 0.0d0 ) then
c             root below this energy so stop scanning
              lscan = .false.
              ecur  = ecur + decur
            else
c             no root found below current energy so sweep upwards
              ecur = ecur * 0.98d0
            endif
          else
            ecur = ecur + decur
          endif
c
          if ( iter .eq. itmax ) then
c           eigenfunction has not converged in allowed number of iterations
            write(iout,*) '***error in subroutine scheqg'
            write(iout,530) nnlz(iorb),ee(iorb),ecur
  530       format(' state',i4,' could not be converged.',/,
     +      ' starting energy for calculation was',f10.5,
     +      ' and end value =',f10.5)
            call exit(1)
          endif
c
c       close iteration loop
  200   continue
c
c-----------------------------------------------------------------------------
c
  400   continue
c
c       update the energy
        ee(iorb) = ecur
c
c       update the wavefunction array
        factor = 1.0d0 / sqrt( factor )
        do 410 i = 1,ninf
          snl(i,iorb) = snlo(i) * sqr(i) * factor
  410   continue
c
c       renormalise the alpha array
        do 420 ib = 1,nbl(lp)
          ibeta = nbl0(lp) + ib
          alpha(ibeta,iorb) = alpha(ibeta,iorb) * factor
  420   continue
c
c     close loop over bands
  100 continue
c
      if ( ifprt .ge. 1) write (iout,'(1x,a)') 'leaving scheqg'
c
      return
      end
c
c--------------------------------------------------------------------------
c
      subroutine descrn(vloc,vloc0,ddd,ddd0,ifpcor,mesh,
     +  idim5,ikeyee,ikeyeeloc,refwf,ruae,
     +  ruexch,ruhar,reexch,rscore,rsvale,rsatom,rspsco,
     +  exfact,r,rab,xwa,etot,wwnl,ee,ncores,ncspvs,zv,
     +  keyps,nbeta,kkbeta,lll,qfunc,ifprt,idim1,idim2,idim3)
c
c     routine for descreening the pseudopotential
c
      implicit double precision (a-h,o-z)
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1)
c.....hartree, exchange and correlation potentials
      dimension ruhar(idim1),reexch(idim1),ruexch(idim1)
c.....charge densities
      dimension rscore(idim1),rsvale(idim1),rsatom(idim1),rspsco(idim1)
c.....shell labels, energies, occupancies and real-space cutoff index
      dimension ee(idim2),wwnl(idim2)
c.....beta and q functions and q pseudization coefficients
      dimension qfunc(idim1,idim3,idim3),vloc(idim1),vloc0(idim1)
c.....angular momenta, reference energy, qij and dij of vanderbilt scheme
      dimension lll(idim3),ddd0(idim3,idim3),ddd(idim3,idim3)
c.....type of wave function
      dimension ikeyee(idim3)
      dimension refwf(idim1,idim3),ruae(idim1)
c.....scratch arrays
      dimension xwa(idim1,5)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c--------------------------------------------------------------------------
c
c                   i n i t i a l i s a t i o n
c
      write (iout,'(/a)') ' subroutine descrn: descreen the psp'
c
      if ( keyps .ne. 3 ) then
        write(iout,*) '***error in subroutine descrn'
        write(iout,*) 'keyps =',keyps,' is out of programed range'
        call exit(1)
      endif
c
c--------------------------------------------------------------------------
c
c           d e s c r e e n  l o c a l  p o t e n t i a l
c
c     compute the valence hartree and exchange correlation potentials
      ipass = 2
      call vhxc(ipass,ifpcor,ruexch,ruhar,reexch,rscore,
     +  rsvale,rsatom,rspsco,xwa(1,1),exfact,r,rab,
     +  xwa(1,2),xwa(1,3),xwa(1,4),ehar,emvxc,mesh,idim1)
c
c     subtract from the local potential
      do i = 1,mesh
         vloc0(i) = vloc(i) - ruhar(i) - ruexch(i)
      end do
c
c     compute alpha of ihm, zunger and cohen j. phys. c 12 4409 (1979)
      do i = 1,mesh
        xwa(i,1) = r(i) * ( vloc0(i) + 2.0d0 * zv )
      end do
      call radin(mesh,xwa(1,1),rab,alpha,idim1)
c
      if ( ifprt .ge. 1 ) then
        write (iout,950) 
        isk=20
        if ( ifprt .ge. 4 ) isk=1
        do i = isk,mesh,isk
          write(iout,960) i,r(i),vloc0(i),vloc(i),ruhar(i),ruexch(i),
     +       rsvale(i),xwa(i,1)
        end do
      endif
  950 format(/4x,'i',13x,'r',9x,'vloc0',10x,'vloc',9x,'ruhar',
     +          8x,'ruexch',8x,'rsvale',9x,'alpha')
  960 format(i5,5f14.8,d14.6,f14.8)
c
c---------------------------------------------------------------------------
c     This test is no longer needed; once the calculation
c     of the xc potential was cleaned up so that it truly
c     goes to zero at small density, the problem went away
c c
c c              c h e c k   c o n v e r g e n c e    o f   a l p h a
c c
c c     find last mesh point at which rsvale exceeds 1.d-18
c       do 160 i=1,mesh
c         i15=mesh+1-i
c         if (rsvale(i15).gt.1.d-18) goto 161
c   160 continue
c   161 alconv=(xwa(mesh,1)-xwa(i15,1))/xwa(mesh,1)
c c  161 alconv=(xwa(mesh,1)-xwa(i15,1))
c c
c       if (abs(alconv).ge.1.d-06) then
c          write(iout,*) '***error subroutine descrn'
c          write(iout,*) '-----inadequate convergence of alpha--------'
c          write(iout,*) '-----try decreasing tolinf in aesubs----------'
c          write(iout,*) '-----or reconstructing your radial mesh-------'
c          write(iout,*) '-----or using irel=2 for ae calculation-------'
c          write(iout,'(a,i5,f12.3,e12.3,f13.8)') 
c      +    'r,rsvale,xwa  at i=',i15 ,r(i15 ),rsvale(i15 ),xwa(i15 ,1)
c          write(iout,'(a,i5,f12.3,e12.3,f13.8)') 
c      +    'r,rsvale,xwa  at i=',mesh,r(mesh),rsvale(mesh),xwa(mesh,1)
c          write(iout,*) 'diff betw xwa values required to be < 1.d-06'
c c         call exit(1)
c       endif
c c
c---------------------------------------------------------------------------
c
c              d e s c r e e n  t h e  d d d  m a t r i x
c
c     take different screening of the angular momenta into account
c
c     note: use of work arrays 'xwa' a little bit tricky:
c             xwa(i,1) is just a dummy array that is used by dddscrl
c             xwa(i,lpp) is the local potential
c                   lpp runs from 2 to l_max+2
c
      do lpp=2,5
         do i=1,kkbeta
            xwa(i,lpp) = vloc(i)
         enddo
      enddo
c
c     allow for specification keyee<0 which means a second reference
c     state was used to construct some of the pseudo wavefunctions
c
      do ibeta = 1,nbeta
        if (ikeyee(ibeta) .lt. 0) then
           lpp=lll(ibeta)+2
           indwf=-2*ikeyee(ibeta)-1
           do i=1,kkbeta
              xwa(i,lpp)=vloc(i)-ruae(i)+refwf(i,indwf)
           enddo
        endif
      enddo

      keydir = - 1
      call dddscrl(keydir,r,rab,xwa(1,2),qfunc,xwa(1,1),kkbeta,nbeta,
     +  lll,ddd,ddd0,mesh,ifprt,idim1,idim3,idim5)
c
c---------------------------------------------------------------------------
c
c             c o m p u t e  t h e  t o t a l  e n e r g y 
c
      etot = 0.0d0
      do 200 i = ncores+1,ncspvs
        etot = etot + wwnl(i) * ee(i)
  200 continue
      etot = etot - ehar + emvxc
c
c     get the value of the screened pseudopotential at origin
      ruat0 = vloc(2) / r(2)
c
      write(iout,999) ruat0,etot
  999 format(/' value of screened s potential at origin =',f14.7//
     +       ' total energy =',f12.7)
c
      write (iout,'(a)') ' leaving descrn'
c
      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine dddscr(keydir,r,rab,vloc,qfunc,dum,kkbeta,nbeta,
     +  lll,ddd,ddd0,mesh,ifprt,idim1,idim3)
c
c     routine for descreening the ddd matrix ( keydir < 0 ) or
c     computing ionic, hartree and exch-corr parts of ddd ( keydir > 0 )
c
      implicit double precision (a-h,o-z)
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1)
c.....beta and q functions and q pseudization coefficients
      dimension qfunc(idim1,idim3,idim3),vloc(idim1)
c.....angular momenta, reference energy, qij and dij of vanderbilt scheme
      dimension lll(idim3),ddd0(idim3,idim3),ddd(idim3,idim3)
c.....scratch arrays
      dimension dum(idim1)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c
      do 100 ib = 1,nbeta
      do 100 jb = 1,nbeta
        dscr = 0.0d0
        if ( lll(ib) .eq. lll(jb) ) then
          do 110 i = 1,mesh
            dum(i) = 0.0d0
  110     continue
          do 120 i = 2,kkbeta
            dum(i) = qfunc(i,ib,jb) * vloc(i) / r(i)
  120     continue
          call radin(mesh,dum,rab,dscr,idim1)
        endif
        if ( keydir .ge. 0 ) then
c         add dscr to ddd0 to get ddd
          ddd(ib,jb) = ddd0(ib,jb) + dscr
        else
c         subtract dscr from ddd to get ddd0
          ddd0(ib,jb) = ddd(ib,jb) - dscr
        endif
  100 continue
c
      if ( ifprt .ge. 1 ) then
        if ( keydir .ge. 0 ) then
          call prmat(ddd,nbeta,'d',idim3)
        else
          call prmat(ddd0,nbeta,'0',idim3)
        endif
      endif
c
      return
      end
c
c---------------------------------------------------------------------------
c

      subroutine dddscrl(keydir,r,rab,vloc,qfunc,dum,kkbeta,nbeta,
     +  lll,ddd,ddd0,mesh,ifprt,idim1,idim3,idim5)
c
c     routine for descreening the ddd matrix ( keydir < 0 ) or
c     computing ionic, hartree and exch-corr parts of ddd ( keydir > 0 )
c
      implicit double precision (a-h,o-z)
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1)
c.....beta and q functions and q pseudization coefficients
      dimension qfunc(idim1,idim3,idim3),vloc(idim1,idim5)
c.....angular momenta, reference energy, qij and dij of vanderbilt scheme
      dimension lll(idim3),ddd0(idim3,idim3),ddd(idim3,idim3)
c.....scratch arrays
      dimension dum(idim1)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c
      do 100 ib = 1,nbeta
      do 100 jb = 1,nbeta
        dscr = 0.0d0
        if ( lll(ib) .eq. lll(jb) ) then
          do 110 i = 1,mesh
            dum(i) = 0.0d0
  110     continue
          do 120 i = 2,kkbeta
            dum(i) = qfunc(i,ib,jb) * vloc(i,lll(ib)+1) / r(i)
  120     continue
          call radin(mesh,dum,rab,dscr,idim1)
        endif
        if ( keydir .ge. 0 ) then
c         add dscr to ddd0 to get ddd
          ddd(ib,jb) = ddd0(ib,jb) + dscr
        else
c         subtract dscr from ddd to get ddd0
          ddd0(ib,jb) = ddd(ib,jb) - dscr
        endif
  100 continue
c
      if ( ifprt .ge. 1 ) then
        if ( keydir .ge. 0 ) then
          call prmat(ddd,nbeta,'d',idim3)
        else
          call prmat(ddd0,nbeta,'0',idim3)
        endif
      endif
c
      return
      end
c
c---------------------------------------------------------------------------
c
