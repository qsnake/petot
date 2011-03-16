      subroutine intpnt(wsum,
     +     ninp,nout,
     +     ainp,binp,b,
     +     ntrans,mtrans,icof,mtrx,
     +     nrk,nband,w,rk,indk,fvec,ppflag,akshift,
     +     iflag_kgrid)
c     
c     subroutine reads or computes k points for
c     integration over the brillouin zone.
c     $Id: intpnt.x,v 1.2 89/04/14 20:57:21 sverre Exp $
c
c     $Log:	intpnt.x,v $
c     Revision 1.2  89/04/14  20:57:21  sverre
c     checked in with -k by jimb at 89.05.02.18.03.28.
c     
c     Revision 1.2  89/04/14  20:57:21  sverre
c     Added symmetry check for grid generation vectors.
c     
c     Revision 1.1  89/03/11  19:09:16  sverre
c     Initial revision
c     
      implicit DOUBLE PRECISION(a-h,o-z)
#ifdef PARAM_H
#include PARAM_H
#else
      parameter ( nrkp = nrkp )
#endif
      parameter ( nxkp = 32768 )
      dimension ainp(3,3),binp(3,3),b(3,3),
     +     icof(48),mtrx(3,3,48),
     +     nband(nrkp),w(nrkp),rk(3,nrkp),indk(6,nrkp),fvec(3,4)
c     
c O   wsum                the sum of the weight factors w(i)
c     
      dimension xk(3,nxkp),xw(nxkp),imap(nxkp),
     +     gvec(3,3),rkmod(3,4),rkcar(3),bm(3,3),rkx(3)
      dimension rkpp(3,nxkp)	
      common /work3/ xk,xw,imap
c
      parameter ( zero = 0.0D0 )
      parameter ( eps  = 1.0D-5 )
      parameter ( one  = 1.0D0 )
      parameter ( huge = 2000.D0 )
c
c     rcs id string - allows use of ident to identify binaries
c
      integer ppflag
      real*8 akshift(3)
      character rcsid*50
      rcsid = '$RCSfile: intpnt.x,v $$Revision: 1.2 $'
c     
 
      iflag_kgrid=0

      irand = 0
      icalc = 0
      indk(1,1) = 0
      do 110 i=1,4
         do 100 j=1,3
            fvec(j,i) = 0
  100    continue
  110 continue

c     
c     ****************************************************************
c     read the number of k points to be used in the integration over
c     the brillouin zone and the number of bands to be computed
c     
c     nrk > 0 nrk points will be read from input
c     nrk = 0 points will be computed from parameters on the next line
c     nrk < 0 (-nrk) points will be chosen at random
c     ****************************************************************
c     
ccc      read(ninp,*) nrk,nbandi
      nrk=0
      nbandi=6
cccccccccc LWW

      if (iabs(nrk) .gt. nxkp) call warn(191,xdum,nxkp)
      if (nrk .gt. 0) then
c     
c     ****************************
c     read nrk k points from input
c     ****************************
c     
         wsum = zero
         n = nrk
         do 140 i=1,n
c     
c     ******************************************************
c     read in number of bands to be computed (if zero, the
c     overall nband read above will be used), an optional
c     integer integration weight, a real integration weight,
c     and the reduced k vector (in ibc units) --
c     if nonzero, the integer weight is normalized
c     and used in place of the real weight
c     ******************************************************
c     
ccccc            read(ninp,*) nband(i),iwi,xw(i),(rkx(j),j=1,3)
ccccc LWW
            iwi=1       ! not really used
ccccc LWW
  120       format(2i5,4f10.3)
            call cnvrt(3,0,binp,rkx,xk(1,i))
            if (nband(i) .eq. 0) nband(i) = nbandi
c     
c     kludge to get nonzero nband for expansion points
c     
            if (nbandi .lt. nband(i)) nbandi = nband(i)
c     
c     translate to interval (0-1) (<1)
c     
            do 130 j=1,3
               kdel = xk(j,i) + eps
               if (xk(j,i) .lt. -eps) kdel = kdel - 1
               xk(j,i) = xk(j,i) - DBLE(kdel)
  130       continue
c     
c     initialize imap and sum up integer weights (if used)
c     
            imap(i) = 0
            if (iwi .ne. 0) then
               xw(i) = DBLE(iwi)
               wsum = wsum + xw(i)
            end if
  140    continue
c     
c     normalize if integer weights used
c     
         if (wsum .ne. zero) then
            do 150 i=1,n
               xw(i) = xw(i) / wsum
  150       continue
         end if
c     
c     reduce to irreducible zone using the symmetry operations
c     
         call intsub(n,xk,xw,imap,nred,nexp,
     +        ntrans,mtrans,icof,mtrx,
     +        nrk,w,rk)
c     
      else if (nrk .eq. 0) then
c     
c     **********************************
c     calculate k points from parameters
c     **********************************
c     
         icalc = 1
c     
c     ****************************************************************
c     read four vectors (in ibc units) which will be used to compute
c     the integration k points
c     
c     xk(*,i) = fvec(*,1) i1 + fvec(*,2) i2 + fvec(*,2) i3 + fvec(*,4)
c     
c     where i1-i3 are any set of integers such that xk is within the
c     interval (0-1) (<1)
c     ****************************************************************
c     
      do 170 i=1,4
ccc         read(ninp,*) (rkx(j),j=1,3)
       if(i.eq.1) then
       rkx(1)=1.d0
       rkx(2)=0.d0
       rkx(3)=0.d0
       endif
       if(i.eq.2) then
       rkx(1)=0.d0
       rkx(2)=1.d0
       rkx(3)=0.d0
       endif
       if(i.eq.3) then
       rkx(1)=0.d0
       rkx(2)=0.d0
       rkx(3)=1.d0
       endif
       if(i.eq.4) then
       rkx(1)=0.5*akshift(1)
       rkx(2)=0.5*akshift(2)
       rkx(3)=0.5*akshift(3)
       endif
cccccccccc LWW

  160    format(3f10.6)
c
c     convert to units of wave number space basis vectors
c     and make each component a rational number
c
         call cnvrt(3,0,binp,rkx,fvec(1,i))
         maxd = 50
         call ration(3,fvec(1,i),maxd)
         if (maxd .eq. 0) call warn(192,xdum,i)
  170 continue
c     
c     find reciprocal (gvec) and check that fvec(1-3) are
c     linearly independent
c     
      call recip(one, fvec, gvec, xvol)
      if (abs(xvol) .lt. eps) call warn(193,xvol,idum)
c
c     check that each gvec is equal to some lattice vector
c
      do 190 i=1,3
         do 180 j=1,3
            if (abs(anint(gvec(j,i)) - gvec(j,i)) .gt. eps)
     +           call warn(194,gvec(j,i),i)
  180    continue
  190 continue
c     
c     symmetry check:
c     trans(gvec) * op * fvec(1-3) = int matrix
c     trans(gvec) * op * fvec(4) - trans(gvec) * fvec(4) = int matrix
c     (rkmod is used as scratch array here)
c
      do 240 i = 1,mtrans
         do 210 j = 1,3
            do 200 k = 1,4
               rkmod(j,k) = DBLE(mtrx(j,1,i))*fvec(1,k)
     +              + DBLE(mtrx(j,2,i))*fvec(2,k)
     +              + DBLE(mtrx(j,3,i))*fvec(3,k)
  200       continue
  210    continue
         rkmod(1,4) = rkmod(1,4) - fvec(1,4)
         rkmod(2,4) = rkmod(2,4) - fvec(2,4)
         rkmod(3,4) = rkmod(3,4) - fvec(3,4)
         do 230 k = 1,4
            do 220 j = 1,3
               xm = gvec(1,j)*rkmod(1,k)
     +              + gvec(2,j)*rkmod(2,k)
     +              + gvec(3,j)*rkmod(3,k)
c
c     if result has non-integral element(s), call warn
c
cccccc L.W. Wang               if (abs(anint(xm) - xm) .gt. eps) call warn(195,xm,k)
cccccccc Let it passes the test, just print out a warning
         if (abs(anint(xm) - xm) .gt. eps)  then
        iflag_kgrid=1
        write(6,*) 
     &"**WARNING, k-grid in kpgen.input has lower sym. than atom.conf"
        endif
cccccccc L.W. Wang, July 3, 2001
   
  220       continue
  230    continue
  240 continue
c     
c     generate all k points within a sphere with radius 1
c     and with origin at 0.5, 0.5, 0.5
c     select k points within the interval (0-1)
c     
         ni = int(sqrt(gvec(1,1)**2 + gvec(2,1)**2 + gvec(3,1)**2)) + 1
         i0 = int(gvec(1,1) * (0.5D0 - fvec(1,4))
     +          + gvec(2,1) * (0.5D0 - fvec(2,4))
     +          + gvec(3,1) * (0.5D0 - fvec(3,4)) + 0.5D0)
         nj = int(sqrt(gvec(1,2)**2 + gvec(2,2)**2 + gvec(3,2)**2)) + 1
         j0 = int(gvec(1,2) * (0.5D0 - fvec(1,4))
     +          + gvec(2,2) * (0.5D0 - fvec(2,4))
     +          + gvec(3,2) * (0.5D0 - fvec(3,4)) + 0.5D0)
         nk = int(sqrt(gvec(1,3)**2 + gvec(2,3)**2 + gvec(3,3)**2)) + 1
         k0 = int(gvec(1,3) * (0.5D0 - fvec(1,4))
     +          + gvec(2,3) * (0.5D0 - fvec(2,4))
     +          + gvec(3,3) * (0.5D0 - fvec(3,4)) + 0.5D0)
c     
c     set up array of k points
c     
         n = 1
         imax = 2*ni + 1
         jmax = 2*nj + 1
         kmax = 2*nk + 1
         do 280 ii=1,imax
            i = ii - ni - 1 + i0
            do 270 jj=1,jmax
               j = jj - nj - 1 + j0
               do 260 kk=1,kmax
                  k = kk - nk - 1 + k0
                  if (n .gt. nxkp) call warn(191,xdum,nxkp)
                  do 250 l=1,3
                     xk(l,n) = fvec(l,1) * DBLE(i)
     +                    + fvec(l,2) * DBLE(j)
     +                    + fvec(l,3) * DBLE(k)
     +                    + fvec(l,4)
                     if (xk(l,n) .lt. -eps) goto 260
                     if (xk(l,n) .gt. one-eps) goto 260
  250             continue
c--------------	Check of unreduced k point set
		  write(nout, 913) (xk(lll,n),lll=1,3)
  913    	  format(3f10.7)
                  n = n + 1
  260          continue
  270       continue
  280    continue
         n = n - 1
c     
c     initialize imap and set weights
c     imap is used to map reducible to irreducible k points
c     
         do 290 i=1,n
            imap(i) = 0
            xw(i)   = one / DBLE(n)
  290    continue
c     
c     reduce to irreducible zone using the symmetry operations
c     
         call intsub(n,xk,xw,imap,nred,nexp,
     +        ntrans,mtrans,icof,mtrx,
     +        nrk,w,rk)
c     
c     **************************************************
c     set up array indk
c     indk(1-6,i) gives the index for irreducible points
c     that are neighbors of point i
c     **************************************************
c     
c     cannot be done with current grid generation
c     
      else
c     
c     *****************************
c     generate (-nrk) random points
c     *****************************
c     
         irand = 1
ccccc         read(ninp,*) iseed
ccccccc LWW
         iseed=0
ccccccc LWW

  300    format(2i5)
#ifdef MFECRAY
         if (iseed .ne. 0) call ranset(iseed)
#else
#ifdef CYBER
         if (iseed .ne. 0) call ranset(iseed)
#else
#ifdef SUN
         if (iseed .ne. 0) xdum = drand(iseed)
#else
#ifdef SYSV
         if (iseed .ne. 0) call srand(iseed)
#endif
#endif
#endif
#endif
         nrk = -nrk
         dw = one / DBLE(nrk)
         do 320 i=1,nrk
            w(i) = dw
            do 310 j=1,3
#ifdef MFECRAY
               rk(j,i) = ranf()
#else
#ifdef CYBER
               rk(j,i) = ranf()
#else
#ifdef SUN
               rk(j,i) = drand(0)
#else
#ifdef SYSV
               rk(j,i) = rand()
#endif
#endif
#endif
#endif
  310       continue
  320    continue
         indk(1,1) = 0
      end if
c     
c     ******************************************
c     initialize nband, find rk**2, and printout
c     ******************************************
c     
  330 write(nout,340)
  340 format('integration k-points',/,1x,20('*'),//,
     +     '   i',5x,'nband',7x,'w',16x,'rk (ibc units)',
     +     22x,'rk (a.u.)',17x,'ekin',/)
c     
c     need b metric
c     
      call metric(b,bm)
      wsum = zero
      do 420 i=1,nrk
         if (irand+icalc .ne. 0) nband(i) = nbandi
         if (nband(i) .le. 0) nband(i) = nbandi
         wsum = wsum + w(i)
c     
c     find rk**2
c     
         do 360 j=1,4
            do 350 k=1,3
               rkmod(k,j) = rk(k,i)
  350       continue
  360    continue
         rkmod(1,2) = rkmod(1,2)-one
         rkmod(2,3) = rkmod(2,3)-one
         rkmod(3,4) = rkmod(3,4)-one
         ek = huge
         do 390 j=1,4
            ek1 = zero
            ek2 = zero
            do 380 k=1,3
               do 370 l=1,3
                  ek1 = ek1 + rkmod(k,j)*bm(k,l)*rkmod(l,j)
                  ek2 = ek2 + (rkmod(k,j)-one)
     +                 *bm(k,l)*(rkmod(l,j)-one)
  370          continue
  380       continue
            if (ek1 .lt. ek) then
               ek = ek1
               jcar = j
            end if
            if (ek2 .lt. ek) then
               ek = ek2
               jcar = -j
            end if
  390    continue
c     
c     rk in cartesian coordinates
c     
         neg = 0
         do 400 k=1,3
            jabs = iabs(jcar)
            rkcar(k) = b(k,1)*rkmod(1,jabs)
     +               + b(k,2)*rkmod(2,jabs)
     +               + b(k,3)*rkmod(3,jabs)
            if (jcar .lt. 0) rkcar(k) = rkcar(k)-b(k,1)-b(k,2)-b(k,3)
            if (rkcar(k) .lt. zero) neg = neg + 1
  400    continue
         if (neg .le. 1) neg = 1
         if (neg .gt. 1) neg = -1
c     
c     printout
c     
         call cnvrt(3,1,ainp,rk(1,i),rkx)

         write(nout,410) i,nband(i),w(i),(rkx(j),j=1,3),
     +        (DBLE(neg)*rkcar(j),j=1,3),ek
  410    format(1x,i3,i8,3x,f10.6,3x,3f10.6,3x,3f10.6,3x,f10.6)

C----HELP-LOOP for saving rk's for output in PP format (SM,25 Nov.99,NREL)
 	 Do 9211 ipp=1,3
	   rkpp(ipp,i) = rkx(ipp)
9211	 Continue
C-----------------------------------------------------
	 
  420 continue
c     
c    printout in LAPW format
c
	write(nout,500)
  500   format(//'integration k-points in LAPW format',/,1x,35('*')//,
     +   8x,'rk (rec units)',8x,'     w     ','   i'/) 
	do 510 i=1,nrk
	  write(nout,520) (rk(j,i),j=1,3),w(i)*nred,i
  510   continue

C-NEW PART: ALL WEIGHTS, SUM NORMED TO 1.0 
C-(FOR PSEUDOPOT-PLANEWAVE-CODE)
C WPSUM: Control-Parameter: wpsum must be 1.0!!! (SM,NREL,14 Oct98)

	wpsum = zero


        open(17,file="kpt.file")
        rewind(17)

        write(17,*) nrk

	If (ppflag.eq.1) then
cccccccc the result is wrong !!!
        
        a0=1.d0
        write(17,517) ppflag,a0

 517    format(2x,i3,2x,f12.6,1x, 
     &       "       !iflag,a0,   Cart. Coord !!WRONG RESULT!!")
C
	write(nout,516)
  516   format(/,/'k-points in PP-PW-format(Cart. Coord.)',
     +  /,1x,35('*')/,/,8x,'rk (rec units)',8x,'     w     ','   i'/) 

	do 515 ii=1,nrk
	  wpsum = wpsum + w(ii)
	  write(nout,521) (rkpp(j,ii),j=1,3),w(ii)
          write(17,521) (rkpp(j,ii),j=1,3),w(ii)
  515	continue

	else if (ppflag.eq.2) then 

        a0=1.d0
        write(17,519) ppflag,a0
 519    format(2x,i3,2x,f12.6,1x, 
     &       "    !iflag,a0,     Latt. Coord")

	write(nout,526)
  526   format(/,/'k-points in PP-PW-format(Latt. Coord.)',
     +  /,1x,35('*')/,/,8x,'rk (rec units)',8x,'     w     ','   i'/) 

	do 525 ii=1,nrk
	  wpsum = wpsum + w(ii)
	  write(nout,521) (rk(j,ii),j=1,3),w(ii)
	  write(17,521) (rk(j,ii),j=1,3),w(ii)
  525	continue

	endif
        close(17)

  	write(nout,518) wpsum
  518	format(2X,'CONTROL: WEIGHT-SUM: WPSUM = ',f10.3)
	
  520   format(3f10.7,f10.4,i5)
  521   format(3f10.5,f10.5)

      if (irand .ne. 1) write(nout,430) nred,nexp
  430 format(/,' corresponding to',i5,' distinct reducible vectors',
     +     /,' number of expansion points =',i4)
      if (icalc .eq. 1) then
         write(nout,440)
  440    format(' k points generated by program from parameters')
         do 460 i=1,4
            call cnvrt(3,1,ainp,fvec(1,i),rkx)
            write(nout,450) (rkx(j),j=1,3)
  450       format(5x,3f10.6)
  460    continue
      end if
      if (irand .eq. 1) write(nout,470) iseed
  470 format(/,' k points generated at random from seed =',i5)
      if (abs(wsum) .lt. eps) wsum = zero
      if (wsum .eq. zero) return
      dw = wsum - one
      if (abs(dw) .gt. eps) call warn(196,dw,idum)
      return
      end
