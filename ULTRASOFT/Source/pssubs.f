c
c Copyright (C) 2002 Vanderbilt group
c This file is distributed under the terms of the GNU General Public
c License as described in the file 'License' in the current directory.
c
c----------------------------------------------------------------------------
c
c     ******************************************************
      subroutine psrv(rho,rcut,ifv,r,mesh,ifprt,idim1)
c     ******************************************************
c
c     this subroutine pseudises the function rho to a taylor 
c     series in r**2 between origin and rcut.  the function
c     and its derivative are made continuous at rcut.
c
c     designed for use either on partial core density (ifv=0)
c     or local potential (ifv=1)
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c
c.....logarithmic radial mesh information
      dimension r(idim1)
c.....function to be interpolated
      dimension rho(idim1)
c.....scratch
      dimension targ(2),con(2),bb(2,2)
c     index on targ and 2nd index of bb is (value,deriv)
c     index on con, and 1st index of bb taylor coeff of ( 1 , r^2 )
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c----------------------------------------------------------------------------
c
c                     i n i t i a l i s a t i o n
c
      if (ifprt.gt.-2) write (iout,'(/a)') ' entering subroutine psrv'
c
c     ifv = 1 if this is really a potential, not a charge density
c
c     rho is really 4*pi*r**2*rho
c     v   is really r*v
c     body of subroutine assumes proportional to v or rho
      do 10 i=2,mesh
        rho(i) = rho(i) / r(i)**(2-ifv)
   10 continue
c
c----------------------------------------------------------------------------
c
c           f i n d  i n t e r p o l a t i o n  p o i n t s
c
      do 100 i2 = mesh,1,-1
        if ( r(i2) .lt. rcut ) goto 110
  100 continue
c
  110 continue
c
c     check that the value of i2 is reasonable
      if ( i2 .lt. 3 .or. i2 .gt. mesh-2 ) then
        write(iout,*) '***error in subroutine psrv'
        write(iout,*) 'rcut =',rcut,' generated i2 =',i2,
     +  ' which is illegal with mesh =',mesh
        call exit(0)
      endif
c
      i1 = i2 - 1
      i3 = i2 + 1
      i4 = i2 + 2
c
      if ( ifprt .ge. 5) then
        write (iout,*) ' i1, i2, i3, i4 =',i1,i2,i3,i4
        write (iout,120) ' r1, r2, r3, r4 =',(r(i),i=i1,i4)
        write (iout,120) 'rh1,rh2,rh3,rh4 =',(rho(i),i=i1,i4)
  120   format(a18,4f8.4)
        write (iout,130) rcut
  130   format(' rcut =',f8.4)
      endif
c
c----------------------------------------------------------------------------
c
c     s e t  u p  t h e  p o l y n o m i a l  i n t e r p o l a t i o n
c
      call polyn(rcut,val,r(i1),rho(i1),-1)
      call polyn(rcut,val,r(i1),rho(i1), 0)
      call polyn(rcut, d1,r(i1),rho(i1),+1)
c
c----------------------------------------------------------------------------
c
c         s o l v e  s i m u l t a n e o u s  e q u a t i o n s
c
      targ(1) = val
      targ(2) = d1
c
      bb(1,1) = 1.0d0
      bb(2,1) = 0.0d0
      bb(1,2) = rcut**2
      bb(2,2) = 2.0d0 * rcut
c
c
      if (ifprt .ge. 5) then
        write (iout,150) val,d1
  150   format (' val,d1=',2f15.5)
        write (iout,155) 'val ',(bb(1,j),j=1,2),targ(1)
        write (iout,155) 'd1  ',(bb(2,j),j=1,2),targ(2)
  155   format(1x,a4,3e15.5)
      endif
c
      call simeq(2,2,bb,targ,con)
c
      if ( ifprt .ge. 5) write (iout,160) (con(i),i=1,2)
  160 format (' constants from simeq are ',2e15.5)
c
c----------------------------------------------------------------------------
c
c          f i l l  r h o  w i t h  p s e u d i s e d  f o r m
c
      do 200 i=2,i2
        rr=r(i)**2
        rho(i)=con(1)+rr*con(2)
  200 continue
c
      if ( ifprt .gt. 5) then
        write (iout,*) ' i1, i2, i3, i4 =',i1,i2,i3,i4
        write (iout,120) ' r1, r2, r3, r4 =',(r(i),i=i1,i4)
        write (iout,120) 'rh1,rh2,rh3,rh4 =',(rho(i),i=i1,i4)
      endif
c
      do 210 i=2,mesh
        rho(i) = rho(i) * r(i)**(2-ifv)
  210 continue
c
      if (ifprt.gt.-2) write (iout,'(a)') ' exiting subroutine psrv'
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine psv(ruae,vloc,rcloc,rcmax,r,rab,sqr,a,b,mesh,
     + irel,yy,bndmat,
     + kkbeta,klogd,dum1,dum2,dum3,lloc,eloc,ifprt,idim1,iploctype,
     +  ikeyeeloc,refwf,
     + nang,npf,ptryc,idim8,idim3)
c
c     pseudization of the local potential to obtain correct
c     logarithmic derivative in channel lloc at energy eloc
c
      implicit double precision(a-h,o-z)
c
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1),sqr(idim1)
c.....potentials for the all electron calculation
      dimension ruae(idim1)
c.....local potential
      dimension vloc(idim1)
c.....for the phakoe call
      dimension bndmat(idim1,5),yy(idim1,2)
c.....type of wave function
      dimension refwf(idim1,idim3)
c.....scratch arrays
      dimension dum1(idim1),dum2(idim1),dum3(idim1),con(10)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c-----------------------------------------------------------------------
c
c                 i n i t i a l i z a t i o n
c
      if ( ifprt .ge. 1 ) write (iout,'(/3x,a)')
     +    'subroutine psv: pseudize the local potential'
      write (iout,'(/3x,a,i2,a,f9.5)') 'generating vloc with lloc ',
     +   lloc,' eloc =',eloc
c
      if ( lloc .lt. 0 .or. lloc .gt. 3 ) then
        write(iout,*) '***error subroutine psv'
        write(iout,*) 'lloc is not reasonable, lloc =',lloc
        call exit(1)
      endif
c
c-----------------------------------------------------------------------
c
c          g e n e r a t e  t h e  l o c a l  p o t e n t i a l
c
c     obtain all electron wavefunction at energy eloc in dum2
      if (ikeyeeloc .lt. 0) then
         do i=1,mesh
            dum2(i) = refwf(i,3-ikeyeeloc)
         enddo
         do i=1,mesh
            if (r(i) .gt. rcloc) goto 203
         enddo
 203     if (dum2(i) .lt. 0.0) then
            do i=1,mesh
               dum2(i) = -1.0*dum2(i)
            enddo
         endif

      elseif ( irel .eq. 0 ) then
c
          write(iout,'(3x,a)') 'calling phase'
          zz = 0.0d0
          vzero = ruae(2) / r(2)
          call phase(zz,vzero,eloc,lloc,b,r,rab,sqr,
     +      ruae,dum1,dum2,mesh,kkbeta,dlwf,idim1)
c
c         renormalise psi and remove square root type factor
          factor = 1.0d0 / ( sqr(klogd) * dum2(klogd) )
          do  i = 1,kkbeta
            dum2(i) = dum2(i) * sqr(i) * factor
         enddo
c
        elseif ( irel .eq. 2 ) then
c
          write(iout,'(3x,a)') 'calling phakoe'
          call phakoe(ruae,r,rab,bndmat,eloc,lloc,
     +      yy,dum1,kkbeta,dlwf,idim1)
c
c         copy major component of koelling equation into psi
c         and renormalise
          factor = 1.0d0 / yy(klogd,2)
          do 205 ir = 1,kkbeta
            dum2(ir) = yy(ir,2) * factor
  205     continue
c
        else
c
          write(iout,*) '***error in psw'
          write(iout,*) 'irel =',irel,' is not allowed'
          call exit(1)
c
        endif

c
c     nodeless pseudization of the wavefunction
      if(iploctype .eq. 1) then
         write(iout,'(3x,a)') 'calling pswf1'
         call pswf1(dum2,dum1,rcloc,lloc,r,rab,con,kkbeta,
     +        ifprt,idim1)
c
c        load vloc using analytic inversion of schroedinger equation
         vloc(1) = 0.0d0
         do 120 i = 2,mesh
            if ( r(i) .lt. rcloc ) then
               t1 = 2.0d0*con(2)*r(i) + 4.0d0*con(3)*r(i)**3 + 
     +              6.0d0*con(4)*r(i)**5
               t2 = 2.0d0*con(2) + 12.0d0*con(3)*r(i)**2 +
     +              30.0d0*con(4)*r(i)**4
               vtemp = eloc + dble(2*lloc+2)*t1/r(i) + t1**2 + t2
c              remember that vloc is potential * r
               vloc(i) = r(i) * vtemp
            else
c              no pseudization necessary beyond rcloc
               vloc(i) = ruae(i)
            endif
 120     continue
c
      elseif(iploctype .eq. 2) then
         write(iout,'(3x,a)') 'Not implemented'
         call exit(1)
c
c     this is a little mysterious.  it seems kurt stokbro put in
c     the capability for iploctype=2 in version a6.0 (it was not
c     there in a5.1.3, where pswf1 was always used to generate the
c     pseudo wf as the intermediate step to pseudizing the local
c     potential).  but for some reason kurt put this 'Not
c     implemented' abort into the code.  maybe he just never tested
c     it, or there was a bug he never bothered to track down.
c                          - dv 8/17/97
c
c!         write(iout,*) ' calling pswf2'
c!         call pswf2(dum2,dum1,rcloc,lloc,r,rab,con,kkbeta,
c!     +        ifprt,nang,npf,ptryc,idim1,10)
c!
c!         load vloc using analytic inversion of schroedinger equation
c!         vloc(1) = 0.0d0
c!         do i = 2,mesh
c!            if ( r(i) .lt. rcloc ) then
c!               t1 = 0.d0
c!               t2 = 0.d0
c!               do j =2,npf
c!                  t1 = t1 + dble(2*(j-1))*con(j)*r(i)**(2*j-3)
c!                  t2 = t2 + dble((2*j-2)*(2*j-3))*con(j)*r(i)**(2*j-4)
c!               enddo
c!               vtemp = eloc + dble(2*lloc+2)*t1/r(i) + t1**2 + t2
c!               remember that vloc is potential * r
c!               vloc(i) = r(i) * vtemp
c!            else
c!               no pseudization necessary beyond rcloc
c!               vloc(i) = ruae(i)
c!            endif
c!       enddo
c
      elseif (iploctype .eq. 3) then
c
         do i = 1,mesh
               if (ikeyeeloc .lt. 0) then
                  vloc(i) = refwf(i,1)
               else
                  vloc(i) = ruae(i)
               endif
         enddo
         write(iout,'(3x,a)') 'calling pswfnormc'
         call pswfnormc(dum2,dum1,vloc,rcloc,lloc,eloc,r,rab,b,
     +           con,kkbeta,ifprt,idim1,10)
         nqtr = 7
c     
c     load vloc using analytic inversion of schroedinger equation
         write(iout,'(3x,a,f8.3)') 'rcloc is now',rcloc
         vloc(1) = 0.0d0
         do i = 2,mesh
            if ( r(i) .le. rcloc ) then
                  t1 = 0.d0
                  t2 = 0.d0
                  do j =2,nqtr
                   t1 = t1 + dble(2*(j-1))*con(j)*r(i)**(2*j-3)
                   t2 = t2 + dble((2*j-2)*(2*j-3))*con(j)*r(i)**(2*j-4)
                  enddo
                  vtemp = eloc + dble(2*lloc+2)*t1/r(i) + t1**2 + t2
c                 y = r(i)*r(i)
c                 pr = (((((con(7)*y+con(6))*y+con(5))*y+con(4))*y+
c     +                con(3))*y+con(2))*y + con(1)
c                 t1=exp(pr)*r(i)**(lloc+1)
c                 write(*,*) r(i),vloc(i),r(i) * vtemp,dum2(i),t1
c                 remember that vloc is potential * r
                  vloc(i) = r(i) * vtemp
                  if (ikeyeeloc .lt. 0) then
                     vloc(i) = vloc(i)-refwf(i,2)+refwf(i,3)
                  endif
            else
c              no pseudization necessary beyond rcloc
               vloc(i) = ruae(i)
c              write(*,*) r(i),vloc(i),vloc(i),dum2(i),dum2(i)
            endif
         enddo
c
      else
            write(iout,*) '***error in subroutine psv'
            write(iout,*) 'iploctype=',iploctype,'not implemented'
            call exit(1)
c
      endif
c
c     check the results by generating wavefunction with vloc
      write(iout,'(3x,a)') 'calling phase'
      zz = 0.0d0
      vzero = vloc(2) / r(2)
      call phase(zz,vzero,eloc,lloc,b,r,rab,sqr,
     +  vloc,dum1,dum3,mesh,kkbeta,dlwf,idim1)
c
c     renormalise psi to remove square root type factor
      factor = 1.0d0 / ( sqr(klogd) * dum3(klogd) )
      do 130 i = 1,kkbeta
        dum3(i) = dum3(i) * sqr(i) * factor
  130 continue
c
c     possible print out of results
      if ( ifprt .ge. 1) then
        write(iout,9980)
 9980   format(/3x,'results for the local potential vloc',/,
     +    8x,'r',9x,'vloc',9x,'dum2',9x,'dum3')
        nind = 55
        rinc = rcmax / dble(nind)
        rcur = -0.75d0 * rinc
        do 190 i = 1,nind
          rcur = rcur + rinc
          in = int( dlog( rcur / a + 1.0d0 ) / b ) + 1
          write(iout,9979) r(in),vloc(in),dum2(in),dum3(in)
  190   continue
      endif
 9979 format(3x,f9.4,3f13.8)
c
      if ( ifprt .ge. 1 ) write (iout,'(3x,a)') 'leaving psv'
      return
      end
c
c----------------------------------------------------------------------------
c
c     *************************************************************
      subroutine pswf(psi,rho,rcut,lcur,r,rab,con,mesh,ifprt,idim1)
c     *************************************************************
c     this subroutine pseudises the  radial wave function to
c     as a taylor series
c
c     snl(r) =  r**(l+1) * ( a + b*r**2 + c*r**4 + d*r**6 )
c
c     matches value and 1st through third derivatives at rcut
c
c     snl(r) = snl(r) + r**(l+1) * tor * ( r**2 - rcut**2 )**4
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1)
c.....wavefunction to be pseudised (psi is really r*psi on input)
      dimension psi(idim1),rho(idim1)
c.....scratch
      dimension targ(4),con(4),bb(4,4)
c     index on targ and 2nd index of bb is (val,d1,d2,d3)
c     index on con, and 1st index of bb taylor coeff of 1, r^2, ...
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c----------------------------------------------------------------------------
c
c                    i n i t i a l i s a t i o n
c
      if ( ifprt .ge. 0) write (iout,'(/3x,a)') 'subroutine pswf'
c
c     compute the weight under psi
      rho(1) = 0.0d0
      do 10 i = 2,mesh
        rho(i) = psi(i)**2
   10 continue
      call radin(mesh,rho,rab,asum,idim1)
c
      if (ifprt .ge. 5) write (iout,20) asum
   20 format (' integral under unpseudised psi is ',f15.10)
c
      lp = lcur + 1
      do 30 i=2,mesh
        psi(i) = psi(i) / r(i)**lp
   30 continue
c     body of subroutine assumes proportional to psi/r**l
c
c----------------------------------------------------------------------------
c
c           f i n d  i n t e r p o l a t i o n  p o i n t s
c
      do 100 i2 = mesh,1,-1
        if ( r(i2) .lt. rcut ) goto 110
  100 continue
c
  110 continue
c
c     check that the value of i2 is reasonable
      if ( i2 .lt. 3 .or. i2 .gt. mesh-2 ) then
        write(iout,*) '***error in subroutine pswf'
        write(iout,*) 'rcut =',rcut,' generated i2 =',i2,
     +  ' which is illegal with mesh =',mesh
        call exit(0)
      endif
c
      i1 = i2 - 1
      i3 = i2 + 1
      i4 = i2 + 2
c
      if ( ifprt .ge.5) then
        write (iout,*) ' i1, i2, i3, i4 =',i1,i2,i3,i4
        write (iout,120) ' r1, r2, r3, r4 =',(r(i),i=i1,i4)
        write (iout,120) 'ps1,ps2,ps3,ps4 =',(psi(i),i=i1,i4)
  120   format(a18,4f8.4)
        write (iout,130) rcut
  130   format(' rcut =',f8.4)
      endif
c
c----------------------------------------------------------------------------
c
c     s e t  u p  t h e  p o l y n o m i a l  i n t e r p o l a t i o n
c
      call polyn(rcut,val,r(i1),psi(i1),-1)
      call polyn(rcut,val,r(i1),psi(i1), 0)
      call polyn(rcut, d1,r(i1),psi(i1),+1)
      call polyn(rcut, d2,r(i1),psi(i1),+2)
      call polyn(rcut, d3,r(i1),psi(i1),+3)
c
      if ( ifprt .ge.5) write (iout,140) val,d1,d2,d3
  140 format (' val,d1,d2,d3=',4f15.5)
c
c----------------------------------------------------------------------------
c
c         s o l v e  s i m u l t a n e o u s  e q u a t i o n s
c
      targ(1) = val
      targ(2) = d1
      targ(3) = d2
      targ(4) = d3
c
      bb(1,1) = 1.0d0
      bb(2,1) = 0.0d0
      bb(3,1) = 0.0d0
      bb(4,1) = 0.0d0
      bb(1,2) = rcut**2
      bb(2,2) = 2.0d0 * rcut
      bb(3,2) = 2.0d0
      bb(4,2) = 0.0d0
      bb(1,3) = rcut**4
      bb(2,3) = 4.0d0 * rcut**3
      bb(3,3) = 12.0d0 * rcut**2
      bb(4,3) = 24.0d0 * rcut
      bb(1,4) = rcut**6
      bb(2,4) = 6.0d0 * rcut**5
      bb(3,4) = 30.0d0 * rcut**4
      bb(4,4) = 120.0d0 * rcut**3
c
      if ( ifprt .ge.5) then
        write (iout,150) 'val ',(bb(1,j),j=1,4),targ(1)
        write (iout,150) 'd1  ',(bb(2,j),j=1,4),targ(2)
        write (iout,150) 'd2  ',(bb(3,j),j=1,4),targ(3)
        write (iout,150) 'd3  ',(bb(4,j),j=1,4),targ(4)
      endif
  150 format (1x,a4,5e15.5)
c
      call simeq(4,4,bb,targ,con)
c
      if ( ifprt .ge.5) write (iout,160) 'con ',con
  160 format (1x,a4,4e15.5)
c
c----------------------------------------------------------------------------
c
c          f i l l  p s i  w i t h  p s e u d i s e d  f o r m
c
      do 200 i=2,i2
        rr = r(i)**2
        psi(i) = con(1) + con(2)*rr + con(3)*rr**2 + con(4)*rr**3
  200 continue
c
      if ( ifprt .ge.5) then
        write (iout,*) ' i1, i2, i3, i4 =',i1,i2,i3,i4
        write (iout,120) ' r1, r2, r3, r4 =',(r(i),i=i1,i4)
        write (iout,120) 'ps1,ps2,ps3,ps4 =',(psi(i),i=i1,i4)
      endif
c
c     run a self-consistency check
      call polyn(rcut,val,r(i1),psi(i1),-1)
      call polyn(rcut,val,r(i1),psi(i1), 0)
      call polyn(rcut, d1,r(i1),psi(i1),+1)
      call polyn(rcut, d2,r(i1),psi(i1),+2)
      call polyn(rcut, d3,r(i1),psi(i1),+3)
c
      if ( ifprt .ge.5) write (iout,140) val,d1,d2,d3
c
c-----------------------------------------------------------------------
c
c     r e f o r m  r a d i a l  f u n c t i o n  a n d  w e i g h t
c
      psi(1) = 0.0d0
      do 210 i=2,mesh
        psi(i) = psi(i) * r(i)**lp
  210 continue
c
      rho(1) = 0.0d0
      do 220 i = 2,mesh
        rho(i) = psi(i)**2
  220 continue
      call radin(mesh,rho,rab,asum,idim1)
c
      if ( ifprt .ge.5) write (iout,230) asum
  230 format (' integral under pseudised psi is ',f15.10)
c
      if ( ifprt .ge. 0) write (iout,'(3x,a)') 'leaving pswf'
c
      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine psqf(qfunc,qfcoef,rho,lmin,
     +  rinner,r,rab,a,b,kkbeta,ifprt,nang,idim1,idim6,idim8)
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c     index on targ and 2nd index of bb is (value,deriv,norm)
c     1st indices of qcoef and bb are taylor coeff of 1, r^2, r^4
c
c     ------------------------------------------------------
c     notes on qfunc and qfcoef:
c     ------------------------------------------------------
c     since q_ij(r) is the product of two orbitals like
c     psi_{l1,m1}^star * psi_{l2,m2}, it can be decomposed by
c     total angular momentum l, where l runs over | l1-l2 | ,
c     | l1-l2 | +2 , ... , l1+l2.  (l=0 is the only component
c     needed by the atomic program, which assumes spherical
c     charge symmetry.)
c     
c     recall  qfunc(r) = y1(r) * y2(r)  where y1 and y2 are the
c     radial parts of the wave functions defined according to
c     
c       psi(r-vec) = (1/r) * y(r) * y_lm(r-hat)  .
c     
c     for each total angular momentum l, we pseudize qfunc(r)
c     inside rc as:
c     
c       qfunc(r) = r**(l+2) * [ a_1 + a_2*r**2 + a_3*r**4 ]
c     
c     in such a way as to match qfunc and its 1'st derivative at
c     rc, and to preserve
c     
c       integral dr r**l * qfunc(r)   ,
c     
c     i.e., to preserve the l'th moment of the charge.  the array
c     qfunc has been set inside rc to correspond to this pseudized
c     version using the minimal l, namely l = | l1-l2 | (e.g., l=0
c     for diagonal elements).  the coefficients a_i (i=1,2,3)
c     are stored in the array qfcoef(i,l+1,j,k) for each l so that
c     the correctly pseudized versions of qfunc can be reconstructed
c     for each l.  (note that for given l1 and l2, only the values
c     l = | l1-l2 | , | l1-l2 | +2 , ... , l1+l2 are ever used.)
c     ------------------------------------------------------
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1)
c.....q functions and q pseudization coefficients
      dimension qfunc(idim1),qfcoef(idim8,idim6)
c.....scratch
      dimension rho(idim1)
      dimension targ(3),bb(3,3)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c
      if ( ifprt .ge. 0) write (iout,'(/3x,a)') 'subroutine psqf'
c
c----------------------------------------------------------------------------
c
c           f i n d  i n t e r p o l a t i o n  p o i n t s
c
      do 100 i2 = kkbeta,1,-1
        if ( r(i2) .lt. rinner ) goto 110
  100 continue
c
  110 continue
c
c     check that the value of i2 is reasonable
      if ( i2 .lt. 3 .or. i2 .gt. kkbeta-2 ) then
        write(iout,*) '***error in subroutine psqf'
        write(iout,*) 'rinner =',rinner,' generated i2 =',i2,
     +  ' which is illegal with kkbeta =',kkbeta
        call exit(0)
      endif
c
      i1 = i2 - 1
      i3 = i2 + 1
      i4 = i2 + 2
c
      if ( ifprt .ge.5) then
        write (iout,*) ' i1, i2, i3, i4 =',i1,i2,i3,i4
        write (iout,120) ' r1, r2, r3, r4 =',(r(i),i=i1,i4)
        write (iout,120) 'qf1,qf2,qf3,qf4 =',(qfunc(i),i=i1,i4)
  120   format(a18,4f8.4)
        write (iout,130) rinner
  130   format(' rinner =',f8.4)
      endif
c
c----------------------------------------------------------------------------
c
c     b e g i n  l o o p  o v e r  l t o t
c
      do 200 ltot = 0,2*(nang-1)
c
        ltp = ltot + 1
        ll2 = 2 * ltot + 2
        ll3 = 2 * ltot + 3
        ll5 = 2 * ltot + 5
        ll7 = 2 * ltot + 7
c
        do 210 i = 2,kkbeta
          rho(i) = qfunc(i) / r(i)**(ltot+2)
  210   continue
c
        call polyn(rinner,val,r(i1),rho(i1),-1)
        call polyn(rinner,val,r(i1),rho(i1), 0)
        call polyn(rinner, d1,r(i1),rho(i1),+1)
c
        targ(1) = val
        targ(2) = d1
c
        bb(1,1) = 1.0d0
        bb(2,1) = 0.0d0
        bb(1,2) = rinner**2
        bb(2,2) = 2.0d0 * rinner
        bb(1,3) = rinner**4
        bb(2,3) = 4.0d0*rinner**3
c
c       integral of r**(2*l+2)*rho should be preserved
c
c       bb(3,j) = integral of r**(2*l+2)*rho 
        bb(3,1) = rinner**ll3/ll3
        bb(3,2) = rinner**ll5/ll5
        bb(3,3) = rinner**ll7/ll7
c
c       simpsons rule integration up to rtrue
c
        rho(1)=0.d0
        do 220 i = 2,kkbeta
          rho(i) = qfunc(i) * r(i) ** ltot
  220   continue
c
        call simp(rho,r,rab,a,b,rinner,kkbeta,asum,idim1)
c
        targ(3)=asum
c
        if ( ifprt .ge.5) then
          write (iout,230) val,d1,asum
  230     format ('val,d1,asum=',3f15.5)
          write (iout,240) 'val ',(bb(1,j),j=1,3),targ(1)
          write (iout,240) 'd1  ',(bb(2,j),j=1,3),targ(2)
          write (iout,240) 'asum',(bb(3,j),j=1,3),targ(3)
  240     format (1x,a4,4e15.5)
        endif
c
        call simeq(3,3,bb,targ,qfcoef(1,ltp))
c
        if ( ifprt .ge.5) write (iout,250) 'con ',(qfcoef(i,ltp),i=1,3)
  250   format (1x,a4,3e15.5)
c
c     end loop over tot ang momentum components
  200 continue
c
c-----------------------------------------------------------------------------
c
c     s e t  q f u n c  f o r  l t o t  =  l m i n
c
      lp = lmin+1
      do 300 i = 1,i2
        rr = r(i)**2
        rho(i) = qfcoef(1,lp) + rr*(qfcoef(2,lp) + rr*qfcoef(3,lp))
        qfunc(i) = rho(i)*r(i)**(lmin+2)
  300 continue
c
      if ( ifprt .ge.5) then
        write (iout,*) ' i1, i2, i3, i4 =',i1,i2,i3,i4
        write (iout,120) ' r1, r2, r3, r4 =',(r(i),i=i1,i4)
        write (iout,120) 'qf1,qf2,qf3,qf4 =',(qfunc(i),i=i1,i4)
      endif
c
      if (ifprt.gt.0) write (iout,*) ' exiting subroutine psqf'
c
      return
      end
c
c-----------------------------------------------------------------------------
c
      subroutine simp(rho,r,rab,a,b,rinner,kkbeta,asum,idim1)
c
c     routine for simpson's rule integration up to rinner
c
      implicit double precision (a-h,o-z)
c
      dimension rho(idim1),r(idim1),rab(idim1)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
      x = dlog( rinner / a + 1.0d0 ) / b + 1.0d0
c
c     check that rinner will not require rho beyond kkbeta
      if ( int(x) + 2 .gt. kkbeta ) then
        write(iout,*) '***error in subroutine simp'
        write(iout,*) 'int(x)+2=',int(x)+2,' .gt. kkbeta =',kkbeta
        call exit(1)
      endif
c
      f3 = rho(1) * rab(1)
      asum = 0.0d0
c
c     commence the integration loop
      do 100 i = 1,kkbeta,2
        f1 = f3
        f2 = rho(i+1) * rab(i+1)
        f3 = rho(i+2) * rab(i+2)
        if ( i + 2 .gt. int(x) ) goto 110
        asum = asum + ( f1 + 4.0d0 * f2 + f3 ) / 3.0d0
  100 continue
c
  110 continue
c
      x = x - dble(i)
      aa3 = ( f1 - 2.0d0 * f2 + f3 ) / 6.0d0
      aa2 = ( -3.0d0 * f1 + 4.0d0 * f2 - f3 ) / 4.0d0
      aa1 = f1
      asum = asum + x * ( aa1 + x * ( aa2 + x * aa3 ) )
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine psqf1(qfunc,qfcoef,rho,lmin,
     +  rinner,mesh,r,rab,a,b,kkbeta,ifprt,nang,
     +  nqf,qtryc,idim1,idim6,idim8)
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c     modified by X.-P. Li, in Oct. 1991
c
c     targ are dimensioned nqf+3, the first nqf of them
c     are for taylor coeff of 1, r**2, r**4 etc, the last are
c     (value,deriv,norm) corresponding to 3 lagrange multipliers
c
c     ------------------------------------------------------
c     notes on qfunc and qfcoef:
c     ------------------------------------------------------
c     since q_ij(r) is the product of two orbitals like
c     psi_{l1,m1}^star * psi_{l2,m2}, it can be decomposed by
c     total angular momentum l, where l runs over | l1-l2 | ,
c     | l1-l2 | +2 , ... , l1+l2.  (l=0 is the only component
c     needed by the atomic program, which assumes spherical
c     charge symmetry.)
c     
c     recall  qfunc(r) = y1(r) * y2(r)  where y1 and y2 are the
c     radial parts of the wave functions defined according to
c     
c       psi(r-vec) = (1/r) * y(r) * y_lm(r-hat)  .
c     
c     for each total angular momentum l, we pseudize qfunc(r)
c     inside rc as:
c     
c       qfunc(r) = r**(l+2) * [ a_1 + a_2*r**2 + a_3*r**4 + ...]
c     
c     in such a way as to match qfunc and its 1'st derivative at
c     rc, and to preserve
c     
c       integral dr r**l * qfunc(r)   ,
c     
c     i.e., to preserve the l'th moment of the charge.  the rest of
c     the coeffs are determined by minimizing 
c
c       integral (from qtryc to infinite) dq q**2 * qfunc(q)**2
c     where
c       qfunc(q)=intgral dr qfunc(r) * j_l(qr)
c
c     where qtryc and the number of coeffs nqf are input.
c     
c     qfunc has been set inside rc to correspond to this pseudized
c     version using the minimal l, namely l = | l1-l2 | (e.g., l=0
c     for diagonal elements).  the coefficients a_i (i=1,2,3)
c     are stored in the array qfcoef(i,l+1,j,k) for each l so that
c     the correctly pseudized versions of qfunc can be reconstructed
c     for each l.  (note that for given l1 and l2, only the values
c     l = | l1-l2 | , | l1-l2 | +2 , ... , l1+l2 are ever used.)
c
c     ------------------------------------------------------
c
c.....in order to avoid passing many unnecessary arrays, use following
      parameter(idm8=13,idm1=101,idm9=3)
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1)
c.....q functions and q pseudization coefficients
      dimension qfunc(idim1),qfcoef(idim8,idim6)
c.....scratch
      dimension rho(idim1)
      dimension targ(idm8),bb(idm8,idm8,idm9),targinv(idm8)
      dimension qn(idm1,idm8,idm9),qlf(idm1),g(idm1)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c     to speed up a little
c
      save qn,bb,g
      data ifinit/0/
c
      if ( ifprt .ge. 0) write (iout,'(/3x,a,l)')
     +   'subroutine psqf1, ifinit = ',ifinit
c
      if(idm8.lt.nqf+3) then
        write(iout,*) '***error: in subroutine psqf1, idm8 too small'
        write(iout,*) 'idm8, nqf+3 = ',idm8,nqf+3
        call exit(1)
      endif
c
      if(idm9.lt.nang) then
        write(iout,*) '***error: in subroutine psqf1, idm9 too small'
        write(iout,*) 'idm9, nang = ',idm9,nang
        call exit(1)
      endif
c
c----------------------------------------------------------------------------
c
c           f i n d  i n t e r p o l a t i o n  p o i n t s
c
      do 100 i2 = kkbeta,1,-1
        if ( r(i2) .lt. rinner ) goto 110
  100 continue
c
  110 continue
c
c     check that the value of i2 is reasonable
      if ( i2 .lt. 3 .or. i2 .gt. kkbeta-2 ) then
        write(iout,*) '***error in subroutine psqf1'
        write(iout,*) 'rinner =',rinner,' generated i2 =',i2,
     +  ' which is illegal with kkbeta =',kkbeta
        call exit(0)
      endif
c
      i1 = i2 - 1
      i3 = i2 + 1
      i4 = i2 + 2
c
      if ( ifprt .ge.5) then
        write (iout,*) ' i1, i2, i3, i4 =',i1,i2,i3,i4
        write (iout,120) ' r1, r2, r3, r4 =',(r(i),i=i1,i4)
        write (iout,120) 'qf1,qf2,qf3,qf4 =',(qfunc(i),i=i1,i4)
  120   format(a18,4f8.4)
        write (iout,130) rinner
  130   format(' rinner =',f8.4)
      endif
c
      pi=3.14159265358979
      tpi=2.0*pi
      tpi3=tpi**3
      fpi=4.0*pi
c
c     minimization if nqf gt 3
c
      if(nqf.gt.3) then
c
c----------------------------------------------------------------------------
c
c     prepare things beforehand
c
      if(ifinit.eq.0) then
c
c     prepare g's
c
        dg=qtryc/dfloat(idm1-1)
        do 2010 iq=1,idm1
 2010   g(iq)=dble(iq-1)*dg
c
c     prepare qn's 
c
        facl=1.0
        do 2050 il=1,nang
         ll=il-1
         fll=dble(ll)
         facl=facl*(2.0*fll+1.0)
        do 2050 in=1,nqf
         lp2n=ll+2*in+1
         tlnp=dble(2*(ll+in)+1)
         if(ll.eq.0) then
           qn(1,in,il)=fpi*rinner**lp2n/facl/tlnp
         else
           qn(1,in,il)=0.0d0
         endif
c
        do 2050 iq=2,idm1
          qdr=g(iq)*rinner
          hqdr2=0.5*qdr*qdr
c
          termth=1.0d0
          ith=0
 2020     continue
          ith=ith+1
          di=dble(ith)
          termth=termth*hqdr2/di/(2.0*(fll+di)+1.0)
          if(abs(termth).gt.1.0d-20) goto 2020
          nterm=ith
c
          termth=0.0d0
          do 2030 ith=nterm,1,-1
            di=dble(ith)
            termth=hqdr2/di/(2.0*(fll+di)+1.0)
     *             *(1.0/(tlnp+2.0*di)-termth)
 2030     continue
          qn(iq,in,il)=fpi*rinner**lp2n*qdr**ll/facl
     *                 *(1.0/tlnp-termth)
 2050   continue
c
        if(ifprt.gt.0) then
          write(iout,'(5x,a,i5)') 'notice!! idm1 = ',idm1
          write(iout,'(5x,a,i5)') 'qns for l = 2, n = ',nqf
          write(iout,250) (qn(iq,nqf,3),iq=1,idm1)
        endif
c
c       now construct matrices bb
c
        do 2100 ltp=1,nang
c
        do 2060 iqf=1,nqf+3
        do 2060 jqf=1,nqf+3
 2060   bb(iqf,jqf,ltp)=0.0
c
        do 2070 iqf=1,nqf
          iqft=2*(iqf-1)
          bb(nqf+1,iqf,ltp)=rinner**iqft
          bb(nqf+2,iqf,ltp)=dfloat(iqft)*rinner**(iqft-1)
c       integral of r**(2*l+2)*rho should be preserved
          lmtp=2*(ltp-1+iqf)+1
          bb(nqf+3,iqf,ltp)=rinner**lmtp/dfloat(lmtp)
c       symmetrical part
          bb(iqf,nqf+1,ltp)=bb(nqf+1,iqf,ltp)
          bb(iqf,nqf+2,ltp)=bb(nqf+2,iqf,ltp)
          bb(iqf,nqf+3,ltp)=bb(nqf+3,iqf,ltp)
 2070   continue
c
        do 2090 iqf=1,nqf
        do 2090 jqf=1,nqf
          do 2080 iq=1,idm1
 2080     rho(iq)=g(iq)*g(iq)*qn(iq,iqf,ltp)*qn(iq,jqf,ltp)
          call simpagn(idm1,g,rho,asum,0)
          lijtm=2*(ltp-1+iqf+jqf)-1
          bb(iqf,jqf,ltp)=tpi3*rinner**lijtm/dfloat(lijtm)-asum
 2090   continue
c
 2100   continue
c
      endif
c
c
c     b e g i n  l o o p  o v e r  l t o t
c
      do 200 ltot = 0,2*(nang-1)
c
        ltp = ltot + 1
c
c     since only ltot=lmin is used, set the rest to 0
c
       if(ltot.ne.lmin) then
c
        do 2110 iqf=1,nqf
 2110   qfcoef(iqf,ltp)=0.0d0
c
       else
c
        if(ltp.gt.nang) then
          write(iout,*) '***error! you must have skipped an l.'
          write(iout,*) 'not programed. it will work with some change'
          call exit(1)
        endif
c
c     construct qlf 
c
        do 140 iq=1,idm1
          do 150 i=1,mesh
 150      rho(i)=qfunc(i)*bessel(g(iq)*r(i),ltot)
          call simp(rho,r,rab,a,b,rinner,kkbeta,asum2,idim1)
          call radin(mesh,rho,rab,asum1,idim1)
          qlf(iq)=fpi*(asum1-asum2)
 140    continue
c
        do 152 iq=1,idm1
 152    rho(iq)=(g(iq)*qlf(iq))**2
        call simpagn(idm1,g,rho,tqlf,0)
c
        rho(1)=0.0
        do 154 i=2,mesh
 154    rho(i)=(qfunc(i)/r(i))**2
        call simp(rho,r,rab,a,b,rinner,kkbeta,asum2,idim1)
        call radin(mesh,rho,rab,asum1,idim1)
        tqr=tpi3*(asum1-asum2)
c
c     now construct the vector targ which depends on qfunc
c
        do 210 i = 2,kkbeta
          rho(i) = qfunc(i) / r(i)**(ltot+2)
  210   continue
c
        call polyn(rinner,val,r(i1),rho(i1),-1)
        call polyn(rinner,val,r(i1),rho(i1), 0)
        call polyn(rinner, d1,r(i1),rho(i1),+1)
c
        targ(nqf+1) = val
        targ(nqf+2) = d1
c
c       simpsons rule integration up to rtrue
c
        rho(1)=0.d0
        do 220 i = 2,kkbeta
          rho(i) = qfunc(i) * r(i) ** ltot
  220   continue
c
        call simp(rho,r,rab,a,b,rinner,kkbeta,asum,idim1)
c
        targ(nqf+3)=asum
c
        do 240 iqf=1,nqf
          do 235 iq=1,idm1
 235      rho(iq)=g(iq)*g(iq)*qn(iq,iqf,ltp)*qlf(iq)
          call simpagn(idm1,g,rho,asum,0)
          targ(iqf)=asum
 240    continue
c
c     now invert bb for qfcoef
c
        call simeq(idm8,nqf+3,bb(1,1,ltp),targ,targinv)
        do 245 iqf=1,nqf
 245    qfcoef(iqf,ltp)=targinv(iqf)
c
c     print out information if necessary
c
        if ( ifprt .ge.4) then
          write(iout,*) ' the bb matrix = '
          do 246 i=1,nqf
 246      write(iout,250) (bb(i,j,ltp),j=1,nqf)
          write(iout,*) ' the targ vector is '
          write(iout,250) (targ(i),i=1,nqf)
          write (iout,*) ' qfcoef for ltp = ',ltp
          write (iout,250) (qfcoef(i,ltp),i=1,nqf)
          write (iout,*) ' the lagrangian mutipliers = '
          write (iout,250) (targinv(i),i=nqf+1,nqf+3)
c     calculate and print residual chisq
          chisq=tqr-tqlf
          do 247 i=1,nqf
            chisq=chisq-2.0*targ(i)*targinv(i)
          do 247 j=1,nqf
            chisq=chisq+targinv(i)*bb(i,j,ltp)*targinv(j)
 247      continue
          write (iout,*) ' tqr, tqlf = ',tqr,tqlf
          write (iout,*) ' residual chi-square = ',chisq
        endif
  250   format (5x,8(1pd13.5))
c
       endif
c
c     end loop over tot ang momentum components
  200 continue
c
c     for nqf=3, use the original algorithm
c
      else if(nqf.eq.3) then
c
c
      do 500 ltot = 0,2*(nang-1)
c
        ltp = ltot + 1
        ll2 = 2 * ltot + 2
        ll3 = 2 * ltot + 3
        ll5 = 2 * ltot + 5
        ll7 = 2 * ltot + 7
c
        do 510 i = 2,kkbeta
          rho(i) = qfunc(i) / r(i)**(ltot+2)
  510   continue
c
        call polyn(rinner,val,r(i1),rho(i1),-1)
        call polyn(rinner,val,r(i1),rho(i1), 0)
        call polyn(rinner, d1,r(i1),rho(i1),+1)
c
        targ(1) = val
        targ(2) = d1
c
        bb(1,1,1) = 1.0d0
        bb(2,1,1) = 0.0d0
        bb(1,2,1) = rinner**2
        bb(2,2,1) = 2.0d0 * rinner
        bb(1,3,1) = rinner**4
        bb(2,3,1) = 4.0d0*rinner**3
c
c       integral of r**(2*l+2)*rho should be preserved
c
c       bb(3,j,1) = integral of r**(2*l+2)*rho 
        bb(3,1,1) = rinner**ll3/ll3
        bb(3,2,1) = rinner**ll5/ll5
        bb(3,3,1) = rinner**ll7/ll7
c
c       simpsons rule integration up to rtrue
c
        rho(1)=0.d0
        do 520 i = 2,kkbeta
          rho(i) = qfunc(i) * r(i) ** ltot
  520   continue
c
        call simp(rho,r,rab,a,b,rinner,kkbeta,asum,idim1)
c
        targ(3)=asum
c
        if ( ifprt .ge.4) then
          write (iout,530) val,d1,asum
  530     format ('val,d1,asum=',3f15.5)
          write (iout,540) 'val ',(bb(1,j,1),j=1,3),targ(1)
          write (iout,540) 'd1  ',(bb(2,j,1),j=1,3),targ(2)
          write (iout,540) 'asum',(bb(3,j,1),j=1,3),targ(3)
  540     format (1x,a4,4e15.5)
        endif
c
        call simeq(idm8,3,bb(1,1,1),targ,qfcoef(1,ltp))
c
        if ( ifprt .ge.4) write (iout,550) 'con ',(qfcoef(i,ltp),i=1,3)
  550   format (1x,a4,3e15.5)
c
c     end loop over tot ang momentum components
  500 continue
c
c-----------------------------------------------------------------------------
c
c     otherwise, error
c
      else
c
      write(iout,*) '***error!!!'
      write(iout,*) 'nqf not appropriate, nqf = ',nqf
      call exit(1)
c
      endif
c
c     s e t  q f u n c  f o r  l t o t  =  l m i n
c
      lp = lmin+1
      do 300 i = 1,i2
        rr = r(i)**2
        rho(i)=qfcoef(1,lp)
        do 295 iqf=2,nqf
 295    rho(i)=rho(i)+qfcoef(iqf,lp)*rr**(iqf-1)
        qfunc(i) = rho(i)*r(i)**(lmin+2)
  300 continue
c
      if ( ifprt .ge.5) then
        write (iout,*) ' i1, i2, i3, i4 =',i1,i2,i3,i4
        write (iout,120) ' r1, r2, r3, r4 =',r(i1),r(i2),r(i3),r(i4)
        write (iout,120) 'qf1,qf2,qf3,qf4 =',qfunc(i1),qfunc(i2),
     +       qfunc(i3),qfunc(i4)
      endif
c
      if (ifprt.gt.0) write (iout,*) ' exiting subroutine psqf1'
      ifinit=1
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine psqf2(qfunc,qfcoef,rho,lll,
     +  rinner,mesh,r,rab,a,b,kkbeta,ifprt,nang,
     +  nqf,qtryc,ibfix,rfix,qfix,nfix,ifqopt,
     +  phi1,phi2,ib,jb,idim1,idim3,idim6,idim8)
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c     written from psqf1 by dominic king-smith in november 1992
c     modified by kurt stokbro
c     modified further by d vanderbilt 8/97
c
c     arrays like 'targ' that are dimensioned idm8 have the
c     following structure:
c       1 to nqf             taylor functions 1, r**2, r**4 etc
c       nqf+1 to nqf+ncon-1  normal constraints (see below)
c       nqf+ncon             optional extra constraint set by nfix
c     where ncon=4 and ncon=6 for ifqopt=2 and 3 respectively:
c       ifqopt=2:  3 normal constraints are (value,d1,norm)
c       ifqopt=3:  5 normal constraints are (value,d1,d2,d3,norm)
c     where d1,d2,d3 are 1st, 2nd, 3rd derivative
c
c     ------------------------------------------------------
c     notes on qfunc and qfcoef:
c     ------------------------------------------------------
c     since q_ij(r) is the product of two orbitals like
c     psi_{l1,m1}^star * psi_{l2,m2}, it can be decomposed by
c     total angular momentum l, where l runs over | l1-l2 | ,
c     | l1-l2 | +2 , ... , l1+l2.  (l=0 is the only component
c     needed by the atomic program, which assumes spherical
c     charge symmetry.)
c     
c     recall  qfunc(r) = y1(r) * y2(r)  where y1 and y2 are the
c     radial parts of the wave functions defined according to
c     
c       psi(r-vec) = (1/r) * y(r) * y_lm(r-hat)  .
c     
c     for each total angular momentum l, we pseudize qfunc(r)
c     inside rc as:
c     
c       qfunc(r) = r**(l+2) * [ a_1 + a_2*r**2 + a_3*r**4 + ...]
c     
c     in such a way as to match qfunc and its
c             1'st derivative at rc   (ifqopt=2)
c             1'st-3'rd derivs at rc  (ifqopt=3)
c     and to preserve
c     
c       int_{0}^{r_c} dr r**l * qfunc(r)   ,
c     
c     i.e., to preserve the l'th moment of the charge.  also,
c     depending on nfix,ibfix,etc, there may be one more constraint
c     to obtain a particular value at a particular radius (as a trick
c     for avoiding negative charge densities).  the rest of
c     the coeffs are determined by minimizing 
c
c       integral (from qtryc to infinity) dq q**2 * qfunc(q)**2
c     where
c       qfunc(q)=int_{0}^{r_c} dr qfunc(r) * j_l(qr)
c
c     where qtryc and the number of coeffs nqf are input.
c     
c     qfunc has been set inside rc to correspond to this pseudized
c     version using the minimal l, namely l = | l1-l2 | (e.g., l=0
c     for diagonal elements).  the coefficients a_i (i=1,2,3)
c     are stored in the array qfcoef(i,l+1,j,k) for each l so that
c     the correctly pseudized versions of qfunc can be reconstructed
c     for each l.  (note that for given l1 and l2, only the values
c     l = | l1-l2 | , | l1-l2 | +2 , ... , l1+l2 are ever used.)
c
c     ------------------------------------------------------
c
c.....in order to avoid passing many unnecessary arrays, use following
      parameter(idm6=7,nqfx=10,nconx=6,ngx=101,idm8=nqfx+nconx)
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1)
c.....q functions and q pseudization coefficients
      dimension qfunc(idim1),qfcoef(idim8,idim6)
c.....q functions and q pseudization cutoffs
      dimension rinner(idim6),icut(idm6)
c.....q functions and l-dependent pseudization arrays for plot
      dimension qfplot(1000,idm6)
c.....extra fixed points to avoid negative densities when pseudizing q_ij
      dimension ibfix(idim3),rfix(idim3),qfix(idim3)
c.....the phi functions
      dimension phi1(idim1),phi2(idim1)
c.....angular momentum information
      dimension lll(idim3)
c.....scratch
      dimension rho(idim1)
      dimension targ(idm8),bb(idm8,idm8,idm6),targinv(idm8)
      dimension qn(ngx,nqfx,idm6),qlf(ngx),g(ngx)
      logical   tinit,tcon
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c     data to be saved to avoid reinitialization
c
      save qn,bb,g,icut
      data tinit/.false./
c
c-----------------------------------------------------------------------
c
c                      i n i t i a l i z a t i o n
c
      if ( ifprt .ge. 2) write (iout,'(/3x,a,l)')
     +  'entering subroutine psqf2, tinit = ',tinit
c
c     set up qfplot for original unpseudized qfs
      if (ib .eq. jb) then
        nline = 1
        do 90 i=1,idim1
          qfplot(i,nline)=qfunc(i)
 90     continue
      endif
c
c     compute maximum and minimum angular momenta
c
      lmin = iabs ( lll(ib) - lll(jb) )
      lmax = iabs ( lll(ib) + lll(jb) )
c
      if ( ifprt .ge.1) write(iout,'(3x,100a)') ('=',i=1,100)
      if ( ifprt .ge.1) write(iout,'(3x,6(a,i2))')
     +  ' ib =',ib,' , lll(ib) =',lll(ib),' , jb = ',jb,
     +  ' , lll(jb) =',lll(jb),' , lmin =',lmin,' , lmax =',lmax
      if ( ifprt .ge.2) write(iout,'(3x,100a)') ('-',i=1,100)
c
c     =========================================================
c     set variables that control array index idm8
c     =========================================================
c
      if (ifqopt.eq.2) ncon=4
      if (ifqopt.eq.3) ncon=6
c
      nnorm=nqf+ncon-1
      nxtra=nqf+ncon
c
c     ================
c     ckeck dimensions
c     ================
c
      if ( nqf .le. ncon-1 ) then
        write(iout,*) '***error in subroutine psqf2'
        write(iout,*) 'nqf must be .gt.',ncon-1,' but nqf=',nqf
        call exit(1)
      endif
c
      if(nqfx.lt.nqf) then
        write(iout,*) '***error: in subroutine psqf2, nqfx too small'
        write(iout,*) 'nqfx, nqf = ',nqfx,nqf
        call exit(1)
      endif
c
      if(nconx.lt.ncon) then
        write(iout,*) '***error: in subroutine psqf2, nconx too small'
        write(iout,*) 'nconx, ncon = ',nconx,ncon
        call exit(1)
      endif
c
      if(idm6.lt.2*nang-1) then
        write(iout,*) '***error: in subroutine psqf2, idm6 too small'
        write(iout,*) 'idm6, 2*nang-1 = ',idm6,2*nang-1
        call exit(1)
      endif
c
c     set the value of pi
c
      pi=4.0d0*atan(1.0d0)
      tpi=2.0*pi
      tpi3=tpi**3
      fpi=4.0*pi
c
c     =========================
c     possible extra constraint
c     =========================
c
      tcon = .false.
      if ( ib .eq. jb ) then
c
        do 10 ifix = 1,nfix
c
          if ( ib .eq. ibfix(ifix) ) then
c
            if ( ifprt .ge.2)
     +        write (iout,'(5x,a/5x,a,i2/5x,a,f12.7,a,f12.7)')
     +          'extra constraint q_ij pseudization',
     +          'constraint for ibeta =',ib,
     +          'will ensure q_ij passes through',qfix(ifix),
     +          ' at r =',rfix(ifix)
c
            ifcur = ifix
c
            tcon = .true.
            if ( nqf .le. ncon ) then
              write(iout,*) '***error subroutine psqf2'
              write(iout,*) 'nqf .le.',ncon,' will produce singular bb'
              call exit(1)
            endif
c
c           set up the last row and column of the bb matrix
c           corresponding to the extra constraint
c
            do 20 iqf = 1,nqf
c
              iqft = 2*(iqf-1)
c
              bb(nxtra,iqf,1) = rfix(ifix) ** iqft
c
   20       continue
c
            do 30 iqf = nqf+1,nxtra
              bb(nxtra,iqf,1) = 0.0d0
   30       continue
c
            do 40 iqf = 1,nxtra
              bb(iqf,nxtra,1) = bb(nxtra,iqf,1)
   40       continue
c
            targ(nxtra) = qfix(ifix) / rfix(ifix) ** 2
c
          endif
c
   10   continue
c
      endif
c
c     =============================
c     one time initialization tasks
c     =============================
c
      if (.not. tinit) then
c
c     ==========================================
c     prepare the uniformly spaced g-vector list
c     ==========================================
c
      dg=qtryc/dfloat(ngx-1)
      do 200 iq=1,ngx
        g(iq)=dble(iq-1)*dg
  200 continue
c
c     ====================
c     compute the qn array
c     ====================
c
c     write(iout,*)'nqf,nqfx,idm6',nqf,nqfx,idm6
      call setqn(qn,g,rinner,nqf,2*nang-1,ngx,nqfx,idm6)
c     write(iout,*)'exit setqn psqf2'
c
      do 250 ltp=1,2*nang-1
c
c     =========================
c     find interpolation points
c     =========================
c
      do 100 i2 = kkbeta,1,-1
        if ( r(i2) .lt. rinner(ltp) ) goto 110
  100 continue
c
  110 icut(ltp) = i2
c
c     =====================
c     construct matrices bb
c     =====================
c
        do 260 iqf=1,nnorm
        do 260 jqf=1,nnorm
          bb(iqf,jqf,ltp)=0.0d0
  260   continue
c
c       ========================================================
c       set parts of bb arising from imposing normal constraints
c       ========================================================
c
        do 270 iqf=1,nqf
c
          iqft=2*(iqf-1)
c
c         (1) value at rinner is preserved
c
          bb(nqf+1,iqf,ltp)=rinner(ltp)**iqft
c
c         (2) derivative at rinner is preserved
c
          bb(nqf+2,iqf,ltp)=iqft*rinner(ltp)**(iqft-1)
c
          if (ifqopt.eq.3) then
c
c           (3) 2nd derivative at rinner is preserved
c
            bb(nqf+3,iqf,ltp)=(iqft*(iqft-1))*rinner(ltp)**(iqft-2)
c
c           (4) 3rd derivative at rinner is preserved
c
            bb(nqf+4,iqf,ltp)=(iqft*(iqft-1)*(iqft-2))*
     +           rinner(ltp)**(iqft-3)
c
          endif

c         (5) integral of r**(2*l+2)*rho should be preserved
c
          lmtp=2*(ltp-1+iqf)+1
          bb(nnorm,iqf,ltp)=rinner(ltp)**lmtp/dfloat(lmtp)
c
c         symmetrical part
          do jqf=nqf+1,nnorm
            bb(iqf,jqf,ltp)=bb(jqf,iqf,ltp)
          end do
c
  270   continue
c
c       set parts of bb to do with smoothing of qfunc
c
        do 280 iqf=1,nqf
        do 280 jqf=1,nqf
          do 290 iq=1,ngx
            rho(iq)=g(iq)*g(iq)*qn(iq,iqf,ltp)*qn(iq,jqf,ltp)
  290     continue
          call simpagn(ngx,g,rho,asum,0)
          lijtm=2*(ltp-1+iqf+jqf)-1
          bb(iqf,jqf,ltp)=tpi3*rinner(ltp)**lijtm/dfloat(lijtm)-asum
  280   continue
c
  250 continue
c
c     end one time initialization events
c
      tinit = .true.
      endif
c
c-----------------------------------------------------------------------
c
c     b e g i n  l o o p  o v e r  l t o t
c
      do 300 ltot = 0,2*(nang-1)
c
        ltp = ltot + 1
c
c       only need to pseudize for lmin,lmin+2 ...lmax-2,lmax
c
        if ( ltot .lt. lmin .or. ltot .gt. lmax .or. 
     +    mod(ltot-lmin,2) .eq. 1 ) then
c
          do 310 iqf=1,nqf
            qfcoef(iqf,ltp)=0.0d0
  310     continue
c
        else
c
c         =============
c         construct qlf 
c         =============
c
c         note that 
c
c           qlf(g) = 4 * pi * int_{rinner}^{r_c} qfunc(r) j_l(gr) dr
c
          do 320 iq=1,ngx
            do 330 i=1,mesh
              rho(i)=qfunc(i)*bessel(g(iq)*r(i),ltot)
  330       continue
            call simp(rho,r,rab,a,b,rinner(ltp),kkbeta,asum2,idim1)
            call radin(mesh,rho,rab,asum1,idim1)
            qlf(iq)=fpi*(asum1-asum2)
  320     continue
c
c         ============
c         compute tqlf
c         ============
c
c         note that
c
c           tqlf = int_{0}^{qtryc} g**2 qlf(g)**2 dg
c
          do 340 iq=1,ngx
            rho(iq)=(g(iq)*qlf(iq))**2
  340     continue
          call simpagn(ngx,g,rho,tqlf,0)
c
c         ===========
c         compute tqr
c         ===========
c
c         note that
c
c           tqr = 8 * pi**3 int_{rinner}^{r_c} qfunc(r)**2 / r**2 dr
c
          rho(1)=0.0
          do 350 i=2,mesh
            rho(i)=(qfunc(i)/r(i))**2
  350     continue
          call simp(rho,r,rab,a,b,rinner(ltp),kkbeta,asum2,idim1)
          call radin(mesh,rho,rab,asum1,idim1)
          tqr=tpi3*(asum1-asum2)
c
c         ================================================
c         construct the vector targ which depends on qfunc
c         ================================================
c
          do 400 i = 2,kkbeta
            rho(i) = qfunc(i) / r(i)**(ltot+2)
  400     continue
c
          if (ifqopt.eq.2) then
c
            i2 = icut(ltp)
            i1 = i2 - 1
            i3 = i2 + 1
            i4 = i2 + 2
c
            call polyn(rinner(ltp),val,r(i1),rho(i1),-1)
            call polyn(rinner(ltp),val,r(i1),rho(i1), 0)
            call polyn(rinner(ltp), d1,r(i1),rho(i1),+1)
c
            targ(nqf+1) = val
            targ(nqf+2) = d1
c
c           possibly save val and d1 for later print out
            if ( ltot .eq. lmin ) then
              valmin = val
              d1min = d1
            endif
c
          elseif (ifqopt.eq.3) then
c
            i2 = icut(ltp)
            i1 = i2 - 5
            i3 = i2 + 1
            i4 = i2 + 6
c
            call polyn2(rinner(ltp),val,r(i1),rho(i1),-1)
            call polyn2(rinner(ltp),val,r(i1),rho(i1), 0)
            targ(nqf+1)=val
            do i = 1,ncon-3
               call polyn2(rinner(ltp), d1,r(i1),rho(i1),i)
               targ(nqf+1+i) = d1
            enddo
c
c           possibly save val and d1 for later print out
            if ( ltot .eq. lmin ) then
              valmin = val
              d1min = targ(nqf+2)
              d2min = targ(nqf+3)
              d3min = targ(nqf+4)
            endif
c
          else
            call exit(1)
          endif
c
          rho(1) = 0.d0
          do 410 i = 2,kkbeta
            rho(i) = qfunc(i) * r(i) ** ltot
  410     continue
c
          call simp(rho,r,rab,a,b,rinner(ltp),kkbeta,conae,idim1)
c
          targ(nnorm)=conae
c
          do 420 iqf=1,nqf
            do 430 iq=1,ngx
              rho(iq)=g(iq)*g(iq)*qn(iq,iqf,ltp)*qlf(iq)
  430       continue
            call simpagn(ngx,g,rho,asum,0)
            targ(iqf)=asum
  420     continue
c
c         ====================
c         invert bb for qfcoef
c         ====================
c
          if ( tcon .and. ltot .eq. 0 ) then
            ncons = ncon
          else
            ncons = ncon-1
          endif
c
          call simeq(idm8,nqf+ncons,bb(1,1,ltp),targ,targinv)
          do 435 iqf=1,nqf
            qfcoef(iqf,ltp)=targinv(iqf)
  435     continue
c
c         =================================
c         possiple print out of information
c         =================================
c
          if ( ifprt .ge.5) then
            write(iout,'(//,a,i3,a,2i3,//)') 'psqf2: results ltot =',
     +        ltot,' with ibeta,jbeta =',ib,jb
            write(iout,*) ' the bb matrix = '
            do 440 i=1,nqf
              write(iout,490) (bb(i,j,ltp),j=1,nqf)
  440       continue
            write (iout,*) ' the targ vector is '
            write (iout,490) (targ(i),i=1,nqf)
            write (iout,*) ' qfcoef for ltp = ',ltp
            write (iout,490) (qfcoef(i,ltp),i=1,nqf)
            write (iout,*) ' the lagrange mutipliers = '
            write (iout,490) (targinv(i),i=nqf+1,nqf+ncons)
  490       format (8(1pd13.5))
c
          endif
c
c         ==================================
c         calculate and print residual chisq
c         ==================================
c
c         note: chisq = fpi**2 int_{qtryc}^{inf} g**2 qfunc(g)**2 dg
c
          chisq=tqr-tqlf
          do 450 i=1,nqf
            chisq=chisq-2.0*targ(i)*targinv(i)
            do 460 j=1,nqf
              chisq=chisq+targinv(i)*bb(i,j,ltp)*targinv(j)
  460       continue
  450     continue
c
c         =====================================
c         compute percentage weight above qtryc
c         =====================================
c
c         note: total = fpi**2 int_{0}^{inf}  g**2 qfunc(g)**2 dg
c               total = tpi**3 int_{0}^{r_c}  qfunc(r)**2 / r**2 dr
c
c         load rho for r <= r(i2) with qfunc
          call setqf(qfcoef,rho,r,nqf,ltot,i2,idim1,idim6,idim8)
c
          do 600 ir = i3,mesh
            rho(ir) = qfunc(ir)
  600     continue
c
          rho(1) = 0.0d0
          do 610 ir = 2,mesh
            rho(ir) = ( rho(ir) / r(ir) ) ** 2
  610     continue
          call radin(mesh,rho,rab,total,idim1)
          total = total * tpi**3
c
c         ==================================================
c         check constraint int_{0}^{r_c} r**ltot qfunc(r) dr
c         ==================================================
c
c         load rho with pseudized form for r <= r(i2)
          call setqf(qfcoef,rho,r,nqf,ltot,i2,idim1,idim6,idim8)
c
          rho(1)=0.d0
          do 630 ir = 2,i2
            rho(ir) = rho(ir) * r(ir) ** ltot
  630     continue
c
          do 620 ir = i2,kkbeta
            rho(ir) = qfunc(ir) * r(ir) ** ltot
  620     continue
c
          call simp(rho,r,rab,a,b,rinner(ltp),kkbeta,conps,idim1)
c
c         summary of results
          if ( ifprt .ge.2) then
            write(iout,'(5x,a,8x,a,13x,a,12x,a,11x,a,13x,a)') 
     +        'ltot','chisq','total','residual','conae','conps'
            write(iout,'(5x,i3,3x,5(g15.7,3x))') 
     +         ltot,chisq,total,chisq/total,conae,conps
          endif
c
c--------------------------------------------------------------
c     set up qfplot for l-dependent pseudized q functions
c
          if (ib .eq. jb) then
            nline = nline + 1
            call setqf(qfcoef,rho,r,nqf,ltot,i2,
     +               idim1,idim6,idim8)
            do 91 i=1,i2
              qfplot(i,nline)=rho(i)
 91         continue
c           write(iout,*) 'psqf2',nline,qfplot(1,nline)
            do 92 i=i3,idim1
              qfplot(i,nline)=qfunc(i)
 92         continue
          endif
c
        endif
c
c     end loop over tot ang momentum components
  300 continue
c
c       print the real-space qfunctions and Fourier transforms
c
      if (ib .eq. jb) then
        call prqf1(ifprt,qfplot,kkbeta,nline,r,a,b,idim1,idim3)
c
c       fourier transform and print reciprocal-space q functions
        call fanalq1(r,rab,qfplot,nline,mesh,ifprt,
     +    idim1,idim3)
c
      endif
c
c-----------------------------------------------------------------------
c
c     s e t  q f u n c  f o r  l t o t  =  l m i n
c
      lminp = lmin +1
      i2 = icut(lminp)

      if (ifqopt.eq.2) then
c
        i1 = i2 - 1
        i3 = i2 + 1
        i4 = i2 + 2
c
c       copy elements up to i2 into qfunc
        call setqf(qfcoef,qfunc,r,nqf,lmin,i2,idim1,idim6,idim8)
c
c       copy elements up to i4 into rho
        call setqf(qfcoef,  rho,r,nqf,lmin,i4,idim1,idim6,idim8)
c
c       run some checks
        do ir = i1,i4
          rho(ir) = rho(ir) / r(ir)**(lmin+2)
        end do
c
        lminp = lmin +1
        call polyn(rinner(lminp),val,r(i1),rho(i1),-1)
        call polyn(rinner(lminp),val,r(i1),rho(i1), 0)
        call polyn(rinner(lminp), d1,r(i1),rho(i1),+1)
c
        if ( ifprt .ge. 2) write (iout,'(5x,a,2f9.5)')
     +    'target and achieved values of qfunc ',valmin,val,
     +    'target and achieved values of dqfunc',d1min,d1
c
      elseif (ifqopt.eq.3) then
c
        i1 = i2 - 5
        i3 = i2 + 1
        i4 = i2 + 6
c
c       copy elements up to i2 into qfunc
        call setqf(qfcoef,qfunc,r,nqf,lmin,i2,idim1,idim6,idim8)
c
c       copy elements up to i4 into rho
        call setqf(qfcoef,  rho,r,nqf,lmin,i4,idim1,idim6,idim8)
c
c       run some checks
        do ir = i1,i4
          rho(ir) = rho(ir) / r(ir)**(lmin+2)
        end do
c
        lminp = lmin +1
        call polyn2(rinner(lminp),val,r(i1),rho(i1),-1)
        call polyn2(rinner(lminp),val,r(i1),rho(i1), 0)
        call polyn2(rinner(lminp), d1,r(i1),rho(i1),+1)
        call polyn2(rinner(lminp), d2,r(i1),rho(i1),+2)
        call polyn2(rinner(lminp), d3,r(i1),rho(i1),+3)
c
        if ( ifprt .ge. 5) then
          write(iout,*) 'for lp=',lminp
          write(iout,'(a,2f11.5)')
     +      'target and achieved values of qfunc',valmin,val
          write(iout,'(a,2f11.5)')
     +      'target and achieved values of dqfunc',d1min,d1
          write(iout,'(a,2f11.5)')
     +      'target and achieved values of d2qfunc',d2min,d2
          write(iout,'(a,2f11.5)')
     +      'target and achieved values of d3qfunc',d3min,d3
        endif
c
      else
        call exit(1)
      endif
c
c     possible test of the constraint agreement
      if ( ifprt .ge. 2 .and. tcon ) then
c
        do 560 j3 = 3,i2
          if ( rfix(ifcur) .lt. r(j3) ) goto 570
  560   continue
c
  570   continue
c
        call polyn(rfix(ifcur),val,r(j3-2),qfunc(j3-2),-1)
        call polyn(rfix(ifcur),val,r(j3-2),qfunc(j3-2), 0)
c
        write(iout,'(a,4f9.5)') 'values of r about rfix    ',
     +    (r(j),j=j3-2,j3+1)
        write(iout,'(a,4f9.5)') 'values of q_ij about qfix ',
     +    (qfunc(j),j=j3-2,j3+1)
        write(iout,'(2a,2f12.7)') 'extra constraint:',
     +    ' target and achieved',qfix(ifcur),val
c
      endif
c
      if ( ifprt .ge. 5) then
        write (iout,*) ' i1, i2, i3, i4 =',i1,i2,i3,i4
        write (iout,120) ' r1, r2, r3, r4 =',r(i1),r(i2),r(i3),r(i4)
        write (iout,120) 'qf1,qf2,qf3,qf4 =',qfunc(i1),qfunc(i2),
     +       qfunc(i3),qfunc(i4)
  120   format(a18,4f8.4)
      endif
c
      return
      end

c-----------------------------------------------------------------------
c
      subroutine pspcor(rspsco,qfcoef,rho,rpcor,mesh,r,rab,a,b,
     +          ifprt,nqf,qtryc,idim1,idim8)

c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c     written from psqf2 by Kurt Stokbro February 1996
c
c     targ are dimensioned nqf+nconx-1, the first nqf of them
c     are for taylor coeff of 1, r**2, r**4 etc, the last are
c     (value,deriv,norm) corresponding to 5 lagrange multipliers
c
c     ------------------------------------------------------
c     notes on rspsco and qfcoef:
c     ------------------------------------------------------
c     since q_ij(r) is the product of two orbitals like
c     psi_{l1,m1}^star * psi_{l2,m2}, it can be decomposed by
c     total angular momentum l, where l runs over | l1-l2 | ,
c     | l1-l2 | +2 , ... , l1+l2.  (l=0 is the only component
c     needed by the atomic program, which assumes spherical
c     charge symmetry.)
c     
c     recall  rspsco(r) = y1(r) * y2(r)  where y1 and y2 are the
c     radial parts of the wave functions defined according to
c     
c       psi(r-vec) = (1/r) * y(r) * y_lm(r-hat)  .
c     
c     for each total angular momentum l, we pseudize rspsco(r)
c     inside rc as:
c     
c       rspsco(r) = r**(l+2) * [ a_1 + a_2*r**2 + a_3*r**4 + ...]
c     
c     in such a way as to match rspsco and its 1'st derivative at
c     rc, and to preserve
c     
c       int_{0}^{r_c} dr r**l * rspsco(r)   ,
c     
c     i.e., to preserve the l'th moment of the charge.  the rest of
c     the coeffs are determined by minimizing 
c
c       integral (from qtryc to infinity) dq q**2 * rspsco(q)**2
c     where
c       rspsco(q)=int_{0}^{r_c} dr rspsco(r) * j_l(qr)
c
c     where qtryc and the number of coeffs nqf are input.
c     
c     rspsco has been set inside rc to correspond to this pseudized
c     version using the minimal l, namely l = | l1-l2 | (e.g., l=0
c     for diagonal elements).  the coefficients a_i (i=1,2,3)
c     are stored in the array qfcoef(i,l+1,j,k) for each l so that
c     the correctly pseudized versions of rspsco can be reconstructed
c     for each l.  (note that for given l1 and l2, only the values
c     l = | l1-l2 | , | l1-l2 | +2 , ... , l1+l2 are ever used.)
c
c     ------------------------------------------------------
c
c.....in order to avoid passing many unnecessary arrays, use following
      parameter(idm6=1,nqfx=10,nconx=4,ngx=101,idm8=nqfx+nconx)
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1)
c.....q functions and q pseudization coefficients
      dimension rspsco(idim1),qfcoef(idim8,idm6)
c.....scratch
      dimension rho(idim1),rinner(1),icut(1)
      dimension targ(idm8),bb(idm8,idm8,idm6),targinv(idm8)
      dimension qn(ngx,nqfx,idm6),qlf(ngx),g(ngx)
      logical   tinit
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c     data to be saved to avoid reinitialization
c
      save qn,bb,g,icut
      data tinit/.false./
c
c-----------------------------------------------------------------------
c
c                      i n i t i a l i z a t i o n
c
      if (ifprt .ge. 2 .or. (.not. tinit))
     +   write (iout,'(/3x,a,l)')
     +   'subroutine pspcor: tinit = ',tinit
c
ctest ?
c      write(6,*)'nqf,idim1,idim8',nqf,idm1,idm8
c
c     =========================================================
c     this routine is not appropriate for case where nqf .le. nconx-1
c     =========================================================
c
      if ( nqf .le. nconx ) then
        write(iout,*) '***error in subroutine pspcor'
        write(iout,*) 'nqf must be .gt. nconx, but nqf=',nqf
        call exit(1)
      endif
c
c     ================
c     ckeck dimensions
c     ================
c
      if(nqfx.lt.nqf) then
        write(iout,*) '***error: in subroutine pspcor, nqfx too small'
        write(iout,*) 'nqfx, nqf = ',nqfx,nqf
        call exit(1)
      endif
c
c     set the value of pi
c
      pi=4.0d0*atan(1.0d0)
      tpi=2.0*pi
      tpi3=tpi**3
      fpi=4.0*pi
c
c
c     =============================
c     one time initialization tasks
c     =============================
c
      if (.not. tinit) then
c
c     ==========================================
c     prepare the uniformly spaced g-vector list
c     ==========================================
c
      dg=qtryc/dfloat(ngx-1)
      do 200 iq=1,ngx
        g(iq)=dble(iq-1)*dg
  200 continue
c
c     ====================
c     compute the qn array
c     ====================
c
      rinner(1) = rpcor
ctest?
c      write(6,*)'nqf,nqfx,idm6',nqf,nqfx,idm6
      call setqn(qn,g,rinner,nqf,1,ngx,nqfx,1)
c      write(6,*)'exit setqn in pspcor'
c
      ltp=1
c
c     =========================
c     find interpolation points
c     =========================
c
      do 100 i2 = mesh,1,-1
        if ( r(i2) .lt. rinner(ltp) ) goto 110
  100 continue
c
  110 icut(ltp) = i2
c
c     =====================
c     construct matrices bb
c     =====================
c
        do 260 iqf=1,nqf+nconx
        do 260 jqf=1,nqf+nconx
          bb(iqf,jqf,ltp)=0.0d0
  260   continue
c
c       ===================================================
c       set parts of bb arising from imposing nconx constraints
c       ===================================================
c
        do 270 iqf=1,nqf
c
          iqft=2*(iqf-1)
c
c         (1) value at rinner is preserved
c
ctest?
c          write(6,*)'rinner(ltp)**dble(iqft)',rinner(ltp),dble(iqft)
          bb(nqf+1,iqf,ltp)=rinner(ltp)**dble(iqft)
c
c         (2) derivative at rinner is preserved
c
          bb(nqf+2,iqf,ltp)=dble(iqft)*rinner(ltp)**(iqft-1)
c
c         (3) 2. derivative at rinner is preserved
c
          bb(nqf+3,iqf,ltp)=dble(iqft*(iqft-1))*rinner(ltp)**(iqft-2)
c
c         (4) 3. derivative at rinner is preserved
c
          bb(nqf+4,iqf,ltp)=dble(iqft*(iqft-1)*(iqft-2))*
     +         rinner(ltp)**(iqft-3)
c
c
c         symmetrical part
          do  jqf=nqf+1,nqf+nconx
             bb(iqf,jqf,ltp)=bb(jqf,iqf,ltp)
          enddo
c
  270   continue
c
c       set parts of bb to do with smoothing of rspsco
c
        do 280 iqf=1,nqf
        do 280 jqf=1,nqf
          do 290 iq=1,ngx
            rho(iq)=g(iq)*g(iq)*qn(iq,iqf,ltp)*qn(iq,jqf,ltp)
  290     continue
          call simpagn(ngx,g,rho,asum,0)
          lijtm=2*(ltp-1+iqf+jqf)-1
          if ( ifprt .ge. 5) write(iout,*) 'lijtm=',lijtm
          bb(iqf,jqf,ltp)=
     1      tpi3*rinner(ltp)**dble(lijtm)/dfloat(lijtm)-asum
  280   continue
c
c
c     end one time initialization events
c
      tinit = .true.
      endif
c
c-----------------------------------------------------------------------
c
c     b e g i n  l o o p  o v e r  l t o t
c
      ltot=0
c
        ltp = ltot + 1
c
        i2 = icut(ltp)
        i1 = i2 - 5
        i3 = i2 + 1
        i4 = i2 + 6
c
c
c         =============
c         construct qlf 
c         =============
c
c         note that 
c
c           qlf(g) = 4 * pi * int_{rinner}^{r_c} rspsco(r) j_l(gr) dr
c
          do 320 iq=1,ngx
            do 330 i=1,mesh
              rho(i)=rspsco(i)*bessel(g(iq)*r(i),ltot)
  330       continue
            call simp(rho,r,rab,a,b,rinner(ltp),mesh,asum2,idim1)
            call radin(mesh,rho,rab,asum1,idim1)
            qlf(iq)=fpi*(asum1-asum2)
  320     continue
c
c         ============
c         compute tqlf
c         ============
c
c         note that
c
c           tqlf = int_{0}^{qtryc} g**2 qlf(g)**2 dg
c
          do 340 iq=1,ngx
            rho(iq)=(g(iq)*qlf(iq))**2.0
  340     continue
          call simpagn(ngx,g,rho,tqlf,0)
c
c         ===========
c         compute tqr
c         ===========
c
c         note that
c
c           tqr = 8 * pi**3 int_{rinner}^{r_c} rspsco(r)**2 / r**2 dr
c
          rho(1)=0.0
          do 350 i=2,mesh
            rho(i)=(rspsco(i)/r(i))**2.0
  350     continue
          call simp(rho,r,rab,a,b,rinner(ltp),mesh,asum2,idim1)
          call radin(mesh,rho,rab,asum1,idim1)
          tqr=tpi3*(asum1-asum2)
c
c         ================================================
c         construct the vector targ which depends on rspsco
c         ================================================
c
          do 400 i = 2,mesh
            rho(i) = rspsco(i) / r(i)**dble(ltot+2)
  400     continue
c
          call polyn2(rinner(ltp),val,r(i1),rho(i1),-1)
          call polyn2(rinner(ltp),val,r(i1),rho(i1), 0)
          targ(nqf+1)=val
          do i = 1,nconx-1
             call polyn2(rinner(ltp), d1,r(i1),rho(i1),i)
             targ(nqf+1+i) = d1
          enddo
c
c
c         possibly save val and d1 for later print out
          valmin = val
          d1min = targ(nqf+2)
          d2min = targ(nqf+3)
          d3min = targ(nqf+4)
c
          do 410 i = 1,mesh
            rho(i) = rspsco(i) * r(i) ** dble(ltot)
  410     continue
c
c
          do 420 iqf=1,nqf
            do 430 iq=1,ngx
              rho(iq)=g(iq)*g(iq)*qn(iq,iqf,ltp)*qlf(iq)
  430       continue
            call simpagn(ngx,g,rho,asum,0)
            targ(iqf)=asum
  420     continue
c
c         ====================
c         invert bb for qfcoef
c         ====================
c
          ncon = nconx

c
          call simeq(idm8,nqf+ncon,bb(1,1,ltp),targ,targinv)
          do 435 iqf=1,nqf
            qfcoef(iqf,ltp)=targinv(iqf)
  435     continue
c
c         =================================
c         possiple print out of information
c         =================================
c
          if ( ifprt .ge. 3) then
c
            write(iout,*) ' the bb matrix = '
            do 440 i=1,nqf
              write(iout,490) (bb(i,j,ltp),j=1,nqf)
  440       continue
            write (iout,*) ' the targ vector is '
            write (iout,490) (targ(i),i=1,nqf)
            write (iout,*) ' rspscoef for ltp = ',ltp
            write (iout,490) (qfcoef(i,ltp),i=1,nqf)
            write (iout,*) ' the lagrange mutipliers = '
            write (iout,490) (targinv(i),i=nqf+1,nqf+ncon)
  490       format (8(1pd13.5))
c
          endif
c
c         ==================================
c         calculate and print residual chisq
c         ==================================
c
c         note: chisq = fpi**2 int_{qtryc}^{inf} g**2 rspsco(g)**2 dg
c
          chisq=tqr-tqlf
          do 450 i=1,nqf
            chisq=chisq-2.0*targ(i)*targinv(i)
            do 460 j=1,nqf
              chisq=chisq+targinv(i)*bb(i,j,ltp)*targinv(j)
  460       continue
  450     continue
c
c         =====================================
c         compute percentage weight above qtryc
c         =====================================
c
c         note: total = fpi**2 int_{0}^{inf}  g**2 rspsco(g)**2 dg
c               total = tpi**3 int_{0}^{r_c}  rspsco(r)**2 / r**2 dr
c
c         load rho for r <= r(i2) with rspsco
          call setqf(qfcoef,rho,r,nqf,ltot,i2,idim1,1,idim8)
c
          do 600 ir = i3,mesh
            rho(ir) = rspsco(ir)
  600     continue
c
          rho(1) = 0.0d0
          do 610 ir = 2,mesh
            rho(ir) = ( rho(ir) / r(ir) ) ** 2.0
  610     continue
          call radin(mesh,rho,rab,total,idim1)
          total = total * tpi**3
c
c         ==================================================
c         check constraint int_{0}^{r_c} r**ltot rspsco(r) dr
c         ==================================================
c
c         load rho with pseudized form for r <= r(i2)
          call setqf(qfcoef,rho,r,nqf,ltot,i2,idim1,1,idim8)
c
          do 630 ir = 1,i2
            rho(ir) = rho(ir) * r(ir) ** dble(ltot)
  630     continue
c
          do 620 ir = i2,mesh
            rho(ir) = rspsco(ir) * r(ir) ** dble(ltot)
  620     continue
c     end loop over tot ang momentum components

c
c-----------------------------------------------------------------------
c
c     s e t  q f u n c  f o r  l t o t  =  l m i n
c
      lmin=0
      lminp = lmin +1
      i2 = icut(lminp)
      i1 = i2 - 5
      i3 = i2 + 1
      i4 = i2 + 6
      call setqf(qfcoef,rspsco,r,nqf,lmin,i2,idim1,1,idim8)
c
c     run some checks
      do 520 ir = i1,i4
        rho(ir) = rspsco(ir) / r(ir)**dble(lmin+2)
  520 continue
c
      lminp = lmin +1
      call polyn2(rinner(lminp),val,r(i1),rho(i1),-1)
      call polyn2(rinner(lminp),val,r(i1),rho(i1), 0)
      call polyn2(rinner(lminp), d1,r(i1),rho(i1),+1)
      call polyn2(rinner(lminp), d2,r(i1),rho(i1),+2)
      call polyn2(rinner(lminp), d3,r(i1),rho(i1),+3)
      
c
      if ( ifprt .ge. 4) then
c
        write(iout,*) 'For lp=',lminp
        write(iout,'(a,2f11.5)')
     +    'target and achieved values of rspsco',valmin,val
        write(iout,'(a,2f11.5)')
     +    'target and achieved values of drspsco',d1min,d1
        write(iout,'(a,2f11.5)')
     +    'target and achieved values of d2rspsco',d2min,d2
        write(iout,'(a,2f11.5)')
     +    'target and achieved values of d3rspsco',d3min,d3
c
      endif
c
      if ( ifprt .ge. 5) then
c
        write (iout,*) ' i1, i2, i3, i4 =',i1,i2,i3,i4
        write (iout,120) ' r1, r2, r3, r4 =',r(i1),r(i2),r(i3),r(i4)
        write (iout,120) 'qf1,qf2,qf3,qf4 =',rspsco(i1),rspsco(i2),
     +       rspsco(i3),rspsco(i4)
  120   format(a18,4f8.4)
c
      endif
c
      if ( ifprt .ge. 2) write (iout,'(3x,a)') 'leaving pspcor'
c
      return
      end


c
c
c-----------------------------------------------------------------------
c
      subroutine setqn(qn,g,rinner,nqf,lmaxp,ngx,nqfx,idm6)
c
c-----------------------------------------------------------------------
c
c     ====================
c     compute the qn array
c     ====================
c
c     qn(g,n,l) is a moment of a spherical bessel
c     function between 0 and rinner
c
c     qn(g,n,l) = 4 * pi * int_{0}^{rinner} r**(2n+l) j_l(gr) dr
c
c     for n = 1 to nqf and l = 0 to lmaxp-1
c
c     programmed as a recurrence relation using power series for j_l
c
c-----------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
      dimension qn(ngx,nqfx,idm6),g(ngx),rinner(idm6)
c
ctest?
c      write(6,*)'qn,g,rinner,nqf',qn,g(1),rinner(1),nqf
c      write(6,*)'lmaxp,ngx,nqfx,idm6',lmaxp,ngx,nqfx,idm6
c
      fpi = 16.0d0 * atan(1.0d0)
c
      facl=1.0d0
c
c     loop over angular momenta
      do 200 lp=1,lmaxp
c
        ll=lp-1
        fll=dble(ll)
c       note: facl is double factorial
        facl=facl*(2.0d0*fll+1.0d0)
c
c       loop over the terms in the taylor expansion
        do 210 in=1,nqf
c
          lp2n=ll+2*in+1
          tlnp=dble(2*(ll+in)+1)
          if(ll.eq.0) then
ctest?
c             write(6,*)'rinner(lp),lp2n',rinner(lp),lp2n
c             write(6,*)'facl/tlnp',facl,tlnp
c
            qn(1,in,lp)=fpi*rinner(lp)**lp2n/facl/tlnp
          else
            qn(1,in,lp)=0.0d0
          endif
c
c         loop over points in wavevector grid
          do 220 ig=2,ngx
c
            qdr=g(ig)*rinner(lp)
            hqdr2=0.5d0*qdr*qdr
c
            termth=1.0d0
            ith=0
c
c           establish the number of terms in power series
  250       continue
c
            ith=ith+1
            di=dble(ith)
            termth=termth*hqdr2/di/(2.0*(fll+di)+1.0)
ctest?

            if(abs(termth).gt.1.0d-20) goto 250
            nterm=ith
c
            termth=0.0d0
            do 240 ith=nterm,1,-1
              di=dble(ith)
c
              termth=hqdr2/di/(2.0*(fll+di)+1.0)
     +           *(1.0/(tlnp+2.0*di)-termth)
  240       continue
c
c           finally set the value of qn
            qn(ig,in,lp)=fpi*rinner(lp)**lp2n*qdr**ll/facl
     +                 *(1.0/tlnp-termth)
c

c         close loop over ig
  220     continue
c       close loop over in
ctest?
c          write(6,*)'end loop in'
  210   continue
c     close loop over lp
ctest?
c          write(6,*)'end loop lp'
  200 continue
c
      
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine setqf(qfcoef,rho,r,nqf,ltot,i2,idim1,idim6,idim8)
c
c-----------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c.....logarithmic radial mesh information
      dimension r(idim1)
c.....q pseudization coefficients
      dimension qfcoef(idim8,idim6)
c.....scratch
      dimension rho(idim1)
c
c-----------------------------------------------------------------------
c
c     s e t  q f u n c  f o r  t h i s  l t o t  =  l m i n
c
      ltotp = ltot + 1
c
      do 500 ir = 1,i2
        rr = r(ir)**2.0
        rho(ir)=qfcoef(1,ltotp)
        do 510 iqf=2,nqf
          rho(ir)=rho(ir)+qfcoef(iqf,ltotp)*rr**dble(iqf-1)
  510   continue
        rho(ir) = rho(ir)*r(ir)**(ltot+2)
  500 continue
c
c
      return
      end
c
c-----------------------------------------------------------------------------
c
      subroutine polyn(x,val,r,f,imode)
c
c     routine fits a cubic polynomial through function f
c     returns 0 to 3rd derivative of f at point x
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
      dimension r(4),f(4),c(4)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
      save c
c
c     polynomial definitions:
      p0(x,y1,y2,y3) = (x-y1) * (x-y2) * (x-y3)
      p1(x,y1,y2,y3) = (x-y2)*(x-y3) + (x-y1)*(x-y3) + (x-y1)*(x-y2)
      p2(x,y1,y2,y3) = 2.0d0 * ( (x-y1) + (x-y2) + (x-y3) )
      p3(x,y1,y2,y3) = 6.0d0
c
c----------------------------------------------------------------------------
c
      if ( imode .eq. - 1 ) then
c       initialisation
        c(1) = f(1) / p0(r(1),r(2),r(3),r(4))
        c(2) = f(2) / p0(r(2),r(3),r(4),r(1))
        c(3) = f(3) / p0(r(3),r(4),r(1),r(2))
        c(4) = f(4) / p0(r(4),r(1),r(2),r(3))
      elseif ( imode .eq. 0 ) then
c       return the zeroth derivative
        val = c(1)*p0(x,r(2),r(3),r(4)) + c(2)*p0(x,r(3),r(4),r(1)) +
     +        c(3)*p0(x,r(4),r(1),r(2)) + c(4)*p0(x,r(1),r(2),r(3))
      elseif ( imode .eq. 1 ) then
c       return the first derivative
        val = c(1)*p1(x,r(2),r(3),r(4)) + c(2)*p1(x,r(3),r(4),r(1)) +
     +        c(3)*p1(x,r(4),r(1),r(2)) + c(4)*p1(x,r(1),r(2),r(3))
      elseif ( imode .eq. 2 ) then
c       return the second derivative
        val = c(1)*p2(x,r(2),r(3),r(4)) + c(2)*p2(x,r(3),r(4),r(1)) +
     +        c(3)*p2(x,r(4),r(1),r(2)) + c(4)*p2(x,r(1),r(2),r(3))
      elseif ( imode .eq. 3 ) then
c       return the third derivative
        val = c(1)*p3(x,r(2),r(3),r(4)) + c(2)*p3(x,r(3),r(4),r(1)) +
     +        c(3)*p3(x,r(4),r(1),r(2)) + c(4)*p3(x,r(1),r(2),r(3))
      endif
c
      return
      end
c----------------------------------------------------------------------------
c
c     **************************************************************
      subroutine pswf1(psi,rho,rcut,lcur,r,rab,con,mesh,ifprt,idim1)
c     **************************************************************
c     this subroutine pseudises the  radial wave function to
c     as a taylor series
c
c     snl(r) =  r**(l+1) * exp ( a + b*r**2 + c*r**4 + d*r**6 )
c
c     matches value and 1st through third derivatives at rcut
c
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1)
c.....wavefunction to be pseudised (psi is really r*psi on input)
      dimension psi(idim1),rho(idim1)
c.....scratch
      dimension targ(4),con(4),bb(4,4)
c     index on targ and 2nd index of bb is (val,d1,d2,d3)
c     index on con, and 1st index of bb taylor coeff of 1, r^2, ...
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c----------------------------------------------------------------------------
c
c                    i n i t i a l i s a t i o n
c
      if ( ifprt .ge.0) write (iout,'(/3x,a)') 'subroutine pswf1'
c
c     compute the weight under psi
      rho(1) = 0.0d0
      do 10 i = 2,mesh
        rho(i) = psi(i)**2.0
   10 continue
      call radin(mesh,rho,rab,asum,idim1)
c
      if ( ifprt .ge.5) write (iout,20) asum
   20 format (' integral under unpseudised psi is ',f15.10)
c
      lp = lcur + 1
      do 30 i=2,mesh
        psi(i) = psi(i) / r(i)**dble(lp)
   30 continue
c     body of subroutine assumes proportional to psi/r**l
c
c----------------------------------------------------------------------------
c
c           f i n d  i n t e r p o l a t i o n  p o i n t s
c
      do 100 i2 = mesh,1,-1
        if ( r(i2) .lt. rcut ) goto 110
  100 continue
c
  110 continue
c
c     check that the value of i2 is reasonable
      if ( i2 .lt. 3 .or. i2 .gt. mesh-2 ) then
        write(iout,*) '***error in subroutine pswf1'
        write(iout,*) 'rcut =',rcut,' generated i2 =',i2,
     +  ' which is illegal with mesh =',mesh
        call exit(0)
      endif
c
      i1 = i2 - 1
      i3 = i2 + 1
      i4 = i2 + 2
c
      if ( ifprt .ge.5) then
        write (iout,*) ' i1, i2, i3, i4 =',i1,i2,i3,i4
        write (iout,120) ' r1, r2, r3, r4 =',(r(i),i=i1,i4)
        write (iout,120) 'ps1,ps2,ps3,ps4 =',(psi(i),i=i1,i4)
  120   format(a18,4f8.4)
        write (iout,130) rcut
  130   format(' rcut =',f8.4)
      endif
c
c----------------------------------------------------------------------------
c
c     s e t  u p  t h e  p o l y n o m i a l  i n t e r p o l a t i o n
c
      call polyn(rcut,val,r(i1),psi(i1),-1)
      call polyn(rcut,val,r(i1),psi(i1), 0)
      call polyn(rcut, d1,r(i1),psi(i1),+1)
      call polyn(rcut, d2,r(i1),psi(i1),+2)
      call polyn(rcut, d3,r(i1),psi(i1),+3)
c
      if ( ifprt. ge.5) write (iout,140) val,d1,d2,d3
  140 format (' val,d1,d2,d3=',4f15.5)
c
c----------------------------------------------------------------------------
c
c     s o l v e  s i m u l t a n e o u s  e q u a t i o n s
c
c     the basis for the following is that if f=exp(g) then
c       ln f = g
c       f'   = g' f
c       f''  = g'' f + g' f'
c       f''' = g''' f + 2 g'' f' + g' f''
c
c     NOTE BUG FIX BELOW !  see README and top of runatom.f
c     for discussion of effect of bug
c
c     val must be positive to prevent error when taking log
      sig = val / abs(val)
c
c     possible flip of signs
      val = sig * val
      d1  = sig * d1
      d2  = sig * d2
      d3  = sig * d3
c
      targ(1) = log(val)
      targ(2) = d1
      targ(3) = d2
      targ(4) = d3
c
      bb(1,1) = 1.0d0
      bb(2,1) = 0.0d0
      bb(3,1) = 0.0d0
      bb(4,1) = 0.0d0
      bb(1,2) = rcut**2
      bb(2,2) = 2.0d0*val*rcut
      bb(3,2) = 2.0d0*val + 2.0d0*d1*rcut
c bug!    bb(4,2) = 4.0d0*val + 2.0d0*d2*rcut
      bb(4,2) = 4.0d0*d1 + 2.0d0*d2*rcut
      bb(1,3) = rcut**4
      bb(2,3) = 4.0d0*val*rcut**3
      bb(3,3) = 12.0d0*val*rcut**2 + 4.0d0*d1*rcut**3
      bb(4,3) = 24.0d0*val*rcut + 24.0d0*d1*rcut**2 + 4.0d0*d2*rcut**3
      bb(1,4) = rcut**6
      bb(2,4) = 6.0d0*val*rcut**5
      bb(3,4) = 30.0d0*val*rcut**4 + 6.0d0*d1*rcut**5
      bb(4,4) = 120.0d0*val*rcut**3 + 60.0d0*d1*rcut**4 + 
     +                                          6.0d0*d2*rcut**5
c
      if ( ifprt. ge.5) then
        write (iout,150) 'val ',(bb(1,j),j=1,4),targ(1)
        write (iout,150) 'd1  ',(bb(2,j),j=1,4),targ(2)
        write (iout,150) 'd2  ',(bb(3,j),j=1,4),targ(3)
        write (iout,150) 'd3  ',(bb(4,j),j=1,4),targ(4)
      endif
  150 format (1x,a4,5e15.5)
c
      call simeq(4,4,bb,targ,con)
c
      if ( ifprt. ge.5) write (iout,160) 'con ',con
  160 format (1x,a4,4e15.5)
c
c----------------------------------------------------------------------------
c
c          f i l l  p s i  w i t h  p s e u d i s e d  f o r m
c
c     save two values
      psii3=psi(i3)
      psii4=psi(i4)
c
c     do 200 i=2,i2
      do 200 i=2,i4
        rr = r(i)**2
        psi(i) = con(1) + con(2)*rr + con(3)*rr**2 + con(4)*rr**3
        if (psi(i).gt.1.0d2) then
          write(iout,*) '***error in subroutine pswf1'
          write(iout,*) 'i,r(i),psi(i)=',i,r(i),psi(i)
          call exit(1)
        endif
        psi(i) = sig * exp( psi(i) )
  200 continue
c
      if ( ifprt. ge.5) then
        write (iout,*) ' i1, i2, i3, i4 =',i1,i2,i3,i4
        write (iout,120) ' r1, r2, r3, r4 =',(r(i),i=i1,i4)
        write (iout,120) 'ps1,ps2,ps3,ps4 =',(psi(i),i=i1,i4)
      endif
c
c     run a self-consistency check
      call polyn(rcut,val,r(i1),psi(i1),-1)
      call polyn(rcut,val,r(i1),psi(i1), 0)
      call polyn(rcut, d1,r(i1),psi(i1),+1)
      call polyn(rcut, d2,r(i1),psi(i1),+2)
      call polyn(rcut, d3,r(i1),psi(i1),+3)
c
      if ( ifprt. ge.5) write (iout,140) val,d1,d2,d3
c
c     restore two values
      psi(i3)=psii3
      psi(i4)=psii4
c
c-----------------------------------------------------------------------
c
c     r e f o r m  r a d i a l  f u n c t i o n  a n d  w e i g h t
c
      psi(1) = 0.0d0
      do 210 i=2,mesh
        psi(i) = psi(i) * r(i)**dble(lp)
  210 continue
c
      rho(1) = 0.0d0
      do 220 i = 2,mesh
        rho(i) = psi(i)**2
  220 continue
      call radin(mesh,rho,rab,asum,idim1)
c
      if ( ifprt .ge.5) write (iout,230) asum
  230 format (' integral under pseudised psi is ',f15.10)
c
      if ( ifprt .ge. 0) write (iout,'(3x,a)') 'leaving pswf1'
c
      return
      end
c
c 
c----------------------------------------------------------------------------
c
c     *************************************************************
      subroutine pswf2(psi,rho,rcut,lcur,r,rab,dummy,mesh,ifprt
     +  ,nang,npf,ptryc,idim1,idim8)
c     *************************************************************
c
c     modified from pswf by X.-P. Li
c
c     this subroutine pseudises the  radial wave function to
c     as a taylor series
c
c     snl(r) =  r**(l+1) * ( a + b*r**2 + c*r**4 + d*r**6 + ... )
c
c     matches value and 1st through third derivatives at rcut
c     and the rest of the coeffients are determined by optimizing
c     the smoothness of the pseudowavefunction
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
      parameter(idm8=24,idm1=201,idm9=5)
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1)
c.....wavefunction to be pseudised (psi is really r*psi on input)
      dimension psi(idim1),rho(idim1)
c.....scratch
      dimension targ(idm8),bb(idm8,idm8,idm9),qn(idm1,idm8,idm9)
     +         ,qlf(idm1),g(idm1),targin(idm8),ifin(idm9),con(idm8)
     +         ,dvs(4),dummy(4)
c
c     bug fix by d.v. 5/96: argument 'dummy' of subroutine pswf2
c       used to be 'con' but really only allocated as con(4)
c       now actally allocate scratch array in this subroutine,
c       as for other arrays like targin
c
c     index on con, and 1st index of bb taylor coeff of 1, r^2, ...
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
      save qn,bb,g
      data ifinit,ifin/0,idm9*0/
c
c
      if ( ifprt .ge. 0) write (iout,'(/3x,a)') 'subroutine pswf2'
c
      if(idm8.lt.npf+4) then
        write(iout,*) '***error: in subroutine pswf2, idm8 too small'
        write(iout,*) 'idm8, npf+4 = ',idm8,npf+4
        call exit(1)
      endif
c
      if(idm9.lt.nang) then
        write(iout,*) '***error: in subroutine psqf1, idm9 too small'
        write(iout,*) 'idm9, nang = ',idm9,nang
        call exit(1)
      endif
c
c
c     compute the weight under psi
      rho(1) = 0.0d0
      do 10 i = 2,mesh
        rho(i) = psi(i)**2
   10 continue
      call radin(mesh,rho,rab,asum,idim1)
c
      if ( ifprt .ge.5) write (iout,20) asum
   20 format (' integral under unpseudised psi is ',f15.10)
c
c----------------------------------------------------------------------------
c
c           f i n d  i n t e r p o l a t i o n  p o i n t s
c
      do 100 i2 = mesh,1,-1
        if ( r(i2) .lt. rcut ) goto 110
  100 continue
c
  110 continue
c
c     check that the value of i2 is reasonable
      if ( i2 .lt. 6 .or. i2 .gt. mesh-6 ) then
        write(iout,*) '***error in subroutine pswf2'
        write(iout,*) 'rcut =',rcut,' generated i2 =',i2,
     +  ' which is illegal with mesh =',mesh
        call exit(0)
      endif
c
      i1 = i2 - 5
      i3 = i2 + 1
      i4 = i2 + 6
c
      if ( ifprt. ge.5) then
        write (iout,*) ' i1, i2, i3, i4 =',i1,i2,i3,i4
        write (iout,120) ' r1, r2, r3, r4 =',r(i1),r(i2),r(i3),r(i4)
        write (iout,120) 'ps1,ps2,ps3,ps4 ='
     +            ,psi(i1),psi(i2),psi(i3),psi(i4)
  120   format(a18,4f8.4)
        write (iout,130) rcut
  130   format(' rcut =',f8.4)
      endif
c
      lp = lcur + 1
c     do 30 i=1,12
c       psitt(i) = psi(i+i1-1) / r(i+i1-1)**lp
c  30 continue
c
      pi=3.14159265358979
      tpi=2.0*pi
      fpi=4.0*pi
      qwc=ptryc/2.0
      npfw=npf
      if(npfw.lt.4) npfw=4
c
c     optimization if npfw gt 4
c
      if(npfw.ge.4) then
c
c     prepare things beforehand
c
      if(ifinit.eq.0) then
c
c     prepare g's
c
        dg=ptryc/dfloat(idm1-1)
        do 2010 iq=1,idm1
 2010   g(iq)=qwc+dfloat(iq-1)*dg
c
        ifinit=1
      endif
c
      if(ifin(lp).eq.0) then
c
c     prepare qn's
c
        ll=lp-1
        fll=dble(ll)
        if(ll.eq.0) then
          facl=1.0
        else if(ll.eq.1) then
          facl=3.0
        else if(ll.eq.2) then
          facl=3.0*5.0
        else if(ll.eq.3) then
          facl=3.0*5.0*7.0
        endif
        do 2050 in=1,npfw
         lp2n=ll+2*in+1
         tlnp=dble(2*(ll+in)+1)
c
        do 2050 iq=1,idm1
          qdr=g(iq)*rcut
          hqdr2=0.5*qdr*qdr
c
          termth=1.0d0
          ith=0
 2020     continue
          ith=ith+1
          di=dble(ith)
          termth=termth*hqdr2/di/(2.0*(fll+di)+1.0)
          if(abs(termth).gt.1.0d-20) goto 2020
          nterm=ith
c
          termth=0.0d0
          do 2030 ith=nterm,1,-1
            di=dble(ith)
            termth=hqdr2/di/(2.0*(fll+di)+1.0)
     *             *(1.0/(tlnp+2.0*di)-termth)
 2030     continue
          qn(iq,in,lp)=fpi*rcut**dble(lp2n)*qdr**dble(ll)/facl
     *                 *(1.0/tlnp-termth)
 2050   continue
c
        if((ifprt .ge.1).and.(lp.eq.3)) then
          write(iout,'(5x,a,i5)') 'notice!! idm1 = ',idm1
          write(iout,'(5x,a,i5)') 'qns for l = 2, n = ',npfw
          write(iout,'(5x,8(1pd13.5))') (qn(iq,npfw,3),iq=1,idm1)
        endif
c
c       now construct matrices bb
c
        do 2060 iqf=1,npfw+4
        do 2060 jqf=1,npfw+4
 2060   bb(iqf,jqf,lp)=0.0
c
        do 2070 iqf=1,npfw
          iqft=2*(iqf-1)+lp
          bb(npfw+1,iqf,lp)=rcut**dble(iqft)
          bb(npfw+2,iqf,lp)=dble(iqft)*rcut**(iqft-1)
          bb(npfw+3,iqf,lp)=dble(iqft*(iqft-1))*rcut**(iqft-2)
          bb(npfw+4,iqf,lp)=dble(iqft*(iqft-1)*(iqft-2))*rcut**(iqft-3)
c       symmetrical part
          bb(iqf,npfw+1,lp)=bb(npfw+1,iqf,lp)
          bb(iqf,npfw+2,lp)=bb(npfw+2,iqf,lp)
          bb(iqf,npfw+3,lp)=bb(npfw+3,iqf,lp)
          bb(iqf,npfw+4,lp)=bb(npfw+4,iqf,lp)
 2070   continue
c
        do 2090 iqf=1,npfw
        do 2090 jqf=1,npfw
          do 2080 iq=1,idm1
 2080     rho(iq)=g(iq)**4*qn(iq,iqf,lp)*qn(iq,jqf,lp)
c2080     rho(iq)=g(iq)*g(iq)*qn(iq,iqf,lp)*qn(iq,jqf,lp)
          call simpagn(idm1,g,rho,asum,0)
          bb(iqf,jqf,lp)=asum
 2090   continue
c
        ifin(lp)=1
      endif
c
c     construct the vector targ
c
      call polyn2(rcut,val,r(i1),psi(i1),-1)
      call polyn2(rcut,val,r(i1),psi(i1), 0)
      call polyn2(rcut, d1,r(i1),psi(i1),+1)
      call polyn2(rcut, d2,r(i1),psi(i1),+2)
      call polyn2(rcut, d3,r(i1),psi(i1),+3)
c
      if ( ifprt. ge.5) write (iout,140) val,d1,d2,d3
  140 format (' val,d1,d2,d3=',4f15.5)
c
      targ(npfw+1) = val
      targ(npfw+2) = d1
      targ(npfw+3) = d2
      targ(npfw+4) = d3
c
      dvs(1)=val
      dvs(2)=d1
      dvs(3)=d2
      dvs(4)=d3
c
      do 2100 iq=1,idm1
 2100 qlf(iq)=fpi*expjl(rcut,g(iq),dvs,lp)
      do 2120 in=1,npfw
        do 2110 iq=1,idm1
 2110   rho(iq)=g(iq)**4*qlf(iq)*qn(iq,in,lp)
c2110   rho(iq)=g(iq)*g(iq)*qlf(iq)*qn(iq,in,lp)
        call simpagn(idm1,g,rho,asum,0)
        targ(in)=-asum
 2120 continue
c
c----------------------------------------------------------------------------
c
c         s o l v e  s i m u l t a n e o u s  e q u a t i o n s
c
      call simeq(idm8,npfw+4,bb(1,1,lp),targ,targin)
c
      do 2130 in=1,npfw
 2130 con(in)=targin(in)
      if ( ifprt. ge.5) then
        write (iout,*) 'con: '
        write (iout,160) (con(in),in=1,npfw)
        write (iout,*) 'Lagrange multipliers: '
        write (iout,160) (targin(in),in=npfw+1,npfw+4)
  160   format (8(1pd13.5))
c
c       residual
c
        chisq=0.0
        do 2170 in1=1,npfw
        do 2170 in2=1,npfw
 2170   chisq=chisq+con(in1)*bb(in1,in2,lp)*con(in2)
        do 2180 in=1,npfw
 2180   chisq=chisq-2.0*con(in)*targ(in)
        do 2190 iq=1,idm1
 2190   rho(iq)=(qlf(iq)*g(iq)**2)**2
c2190   rho(iq)=(qlf(iq)*g(iq))**2
        call simpagn(idm1,g,rho,asum,0)
        chisq=chisq+asum
        write(iout,*) 'residual chi-sq = ',chisq
      endif
c
      else
c
      write(iout,*) 'npfw not correct, stop'
      call exit(1)
c
      endif
c
c----------------------------------------------------------------------------
c
c          f i l l  p s i  w i t h  p s e u d i s e d  f o r m
c

      do 210 i=2,i2
        rr = r(i)**2
        psi(i)=con(1)
        do 200 in=2,npfw
          psi(i) = psi(i)+con(in)*rr**(in-1)
  200   continue
        psi(i)=psi(i)*r(i)**dble(lp)
  210 continue

c
      if ( ifprt. ge.5) then
        write (iout,*) ' i1, i2, i3, i4 =',i1,i2,i3,i4
        write (iout,120) ' r1, r2, r3, r4 =',r(i1),r(i2),r(i3),r(i4)
        write (iout,120) 'ps1,ps2,ps3,ps4 ='
     +                   ,psi(i1),psi(i2),psi(i3),psi(i4)
      endif
c
c     run a self-consistency check
c
      if ( ifprt. ge.5) then 
c
c       do 2140 i=1,12
c         psitt(i) = psi(i+i1-1) / r(i+i1-1)**lp
c2140   continue
        call polyn2(rcut,val,r(i1),psi(i1),-1)
        call polyn2(rcut,val,r(i1),psi(i1), 0)
        call polyn2(rcut, d1,r(i1),psi(i1),+1)
        call polyn2(rcut, d2,r(i1),psi(i1),+2)
        call polyn2(rcut, d3,r(i1),psi(i1),+3)
        write (iout,140) val,d1,d2,d3
c
      endif
c
      rho(1) = 0.0d0
      do 220 i = 2,mesh
        rho(i) = psi(i)**2
  220 continue
      call radin(mesh,rho,rab,asum,idim1)
c
      if ( ifprt .ge.5) write (iout,230) asum
  230 format (' integral under pseudised psi is ',f15.10)
c
      if ( ifprt .ge. 0) write (iout,'(3x,a)') 'leaving pswf2'
c
      return
      end
c
c-----------------------------------------------------------------------------
c
      subroutine polyn2(x,val,r,f,imode)
c
c     routine fits a 5th order polynomial through function f using 12 points
c     returns 0 to 3rd derivative of f at point x
c
c     by X.-P. Li
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
      dimension r(12),f(12),bb(6,6),c(6),cc(6)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
      save c
c
c----------------------------------------------------------------------------
c
      if ( imode .eq. - 1 ) then
c
c       initialize bb
        do n1=1,6
        do n2=1,6
          npow=n1+n2-2
          if (npow.eq.0) then
            bb(n1,n2)=12.0d0
          else
            bb(n1,n2)=0.0
            do i=1,12
              bb(n1,n2)=bb(n1,n2)+(r(i)-x)**npow
            end do
          endif
        end do
        end do
c
c       initialize cc
        do n1=1,6
          if (n1.eq.1) then
            cc(n1)=0.0
            do i=1,12
              cc(n1)=cc(n1)+f(i)
            end do
          else
            cc(n1)=0.0
            do i=1,12
              cc(n1)=cc(n1)+f(i)*(r(i)-x)**(n1-1)
            end do
          endif
        end do
c
c       solve set of equations
        call simeq(6,6,bb,cc,c)
c
      elseif ( imode .eq. 0 ) then
c       return the zeroth derivative
        val=c(1)
c
      elseif ( imode .eq. 1 ) then
c       return the first derivative
        val=c(2)
c
      elseif ( imode .eq. 2 ) then
c       return the second derivative
        val=2.0*c(3)
c
      elseif ( imode .eq. 3 ) then
c       return the third derivative
        val=6.0*c(4)
c
      endif
c
      return
      end
c
      double precision function expjl(rc,q,dvs,lp)
c
c     the integral from rc to inf. r*psi(r)*jl(q*r) with large q
c     and assuming that psi(inf.) is not important
c     expand to q**(-5)
c
c     by X.-P. Li
c
c     l=3 added by Chris J. Pickard February 1998
c
      implicit double precision (a-h,o-z)
      dimension dvs(4)
c
      qdr=q*rc
      cqdr=cos(qdr)
      sqdr=sin(qdr)
      qinv=1.0/q
      qinv2=qinv*qinv
      if(lp.eq.1) then
        expjl=cqdr*qinv2*(dvs(1)-dvs(3)*qinv2)
     +       -sqdr*qinv2*qinv*(dvs(2)-dvs(4)*qinv2)
      else if(lp.eq.2) then
        expjl=sqdr*qinv2*(dvs(1)-qinv2*(dvs(3)+dvs(2)/rc
     +                                 -dvs(1)/(rc*rc)))
     +       +cqdr*qinv2*qinv*(dvs(2)+dvs(1)/rc-qinv2*(dvs(4)
     +          +dvs(3)/rc-2.0*(dvs(2)-dvs(1)/rc)/(rc*rc)))
      else if(lp.eq.3) then
        expjl=cqdr*qinv2*(-dvs(1)+qinv2*(dvs(3)+3.0*dvs(2)/rc))
     +       +sqdr*qinv2*qinv*(dvs(2)+3.0*dvs(1)/rc
     +          -qinv2*(dvs(4)+3.0*(dvs(3)-dvs(2)/rc)/rc))
      else if(lp.eq.4) then
         expjl=dvs(1)*(cqdr*(-6.d0*qinv**3/rc-3.d0*qinv**5/rc**3)+
     +        sqdr*(9.d0*qinv**4/rc**2-qinv2)) +
     +        dvs(2)*(cqdr*(3.d0*qinv**5/rc**2-qinv**3)+
     +        sqdr*6.d0*qinv**4/rc) +
     +        dvs(3)*(6.d0*cqdr*qinv**5/rc+sqdr*qinv**4)+
     +        dvs(4)*cqdr*qinv**5
      else
        stop 99
      endif
c
      return
      end
