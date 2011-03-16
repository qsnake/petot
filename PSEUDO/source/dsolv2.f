      subroutine dsolv2(iter,iconv,ispp,ifcore,lmax,nr,a,b,r,
     1 rab,norb,ncore,no,lo,so,zo,znuc,cdd,cdu,cdc,viod,
     2 viou,vid,viu,ev,ek,ep,wk1,wk2,wk3,wk4,v,ar,br,evi)
c
c  dsolv2 finds the (non) relativistic wave function using
c  difnrl to intgrate the Scroedinger equation or
c  difrel to intgrate the Dirac equation.
c  The energy level from the previous iteration is used
c  as initial guess, and it must therefore be reasonable
c  accurate.
c
c  njtj
c  ###  Cray conversions  
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      character*1 ispp
c
      parameter (zero=0.D0,smev=1.D-4)
Cray      parameter (zero=0.0,smev=1.E-4)
c
      dimension r(nr),rab(nr),no(norb),lo(norb),so(norb),zo(norb),
     1 cdd(nr),cdu(nr),cdc(nr),viod(lmax,nr),viou(lmax,nr),
     2 vid(nr),viu(nr),ev(norb),ek(norb),ep(norb),evi(norb),
     3 wk1(nr),wk2(nr),wk3(nr),wk4(nr),v(nr),ar(nr),br(nr)
c
c  Initialize arrays for charge density.
c
      do 5 i=1,nr
        cdd(i) = zero
 5    continue
      do 10 i=1,nr
        cdu(i) = zero
 10   continue
      if (ifcore .eq. 0) then
        do 15 i=1,nr
          cdc(i)= zero
 15     continue
      endif
c
c  Start the loop over orbitals.
c  Note that spin zero is treated as down.
c
      do 50 i=1,norb
        if (no(i) .le. 0) goto 50
        if (zo(i) .eq. 0.0 .and. iconv .eq. 0) goto 50
        if (ev(i) .ge. 0.0) ev(i)=-smev
c
c  Set up the potential, set the wave functionc array to zero-ar.
c
        lp  = lo(i)+1
        llp = lo(i)*lp
        do 17 j=1,nr
          ar(j)=zero
 17     continue
        if (so(i) .lt. 0.1) then
          do 18 j=2,nr
            v(j) = viod(lp,j)/r(j) + vid(j)
 18       continue
        else
          do 19 j=2,nr
            v(j) = viou(lp,j)/r(j) + viu(j)
 19       continue
        endif
        if (ispp .ne. 'r') then
          do 20 j=2,nr
            v(j) = v(j) + llp/r(j)**2
 20       continue
        endif
c
c  Call the integration routine.
c
        if (ispp .ne. 'r') then
          call difnrl(iter,i,v,ar,br,lmax,nr,a,b,r,
     1     rab,norb,no,lo,so,znuc,viod,viou,vid,viu,
     2     ev,iflag,wk1,wk2,wk3,evi)
        else
          call difrel(iter,i,v,ar,br,lmax,nr,a,b,r,
     1     rab,norb,no,lo,so,znuc,viod,viou,vid,viu,
     2     ev,wk1,wk2,wk3,wk4,evi)
        endif
c
c  Add to the charge density.
c
       if (ispp .eq. 'r') then
         if (so(i) .lt. 0.1) then
           do 30 j=1,nr
             denr = zo(i) *(br(j) * br(j) + ar(j) * ar(j))
             cdd(j) = cdd(j) + denr
 30        continue
         else
           do 31 j=1,nr
             denr = zo(i) *(br(j) * br(j) + ar(j) * ar(j))
             cdu(j) = cdu(j) + denr
 31        continue
         endif
       else
         if (so(i) .lt. 0.1) then
           do 32 j=1,nr
             denr = zo(i) * ar(j) * ar(j)
             cdd(j) = cdd(j) + denr
 32        continue
         else
           do 33 j=1,nr
             denr = zo(i) * ar(j) * ar(j)
             cdu(j) = cdu(j) + denr
 33        continue
         endif
       endif
       if (ifcore .eq. 0 .and. i .le. ncore) then
         do 34 j=1,nr
           denr = zo(i) * ar(j) * ar(j)
           cdc(j)=cdc(j)+denr
 34      continue
       endif
c
c  Compute various quantitities if last iteration.
c
        if (iconv .eq. 1) call orban(ispp,i,ar,br,
     1   lmax,nr,a,b,r,rab,norb,no,lo,zo,so,viod,viou,
     2   vid,viu,ev,ek,ep)
 50   continue
c
c  End loop over orbitals.
c
      return
      end
