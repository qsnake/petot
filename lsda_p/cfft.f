      subroutine cfft(n1,n2,n3,chdr,chdi,wrk,lwrk,mode)
c     
c     computes complex 3d fft of chd
c     $Id: cfft.f,v 1.1.1.1 1992/05/08 19:44:47 jimb Exp $
c
c     $Log: cfft.f,v $
c Revision 1.1.1.1  1992/05/08  19:44:47  jimb
c
c     Revision 1.1  89/03/11  19:09:31  sverre
c     Initial revision
c     
c     
      implicit DOUBLE PRECISION(a-h,o-z)
      dimension chdr(n1,n2,n3),chdi(n1,n2,n3),wrk(1)
c     level 2, chdr,chdi
c     
c I   n1,n2,n3            grid size
c IO  chd(r,i)(i,j,k)     i=1,n1,j=1,n2,k=1,n3.
c     .                   array to be fourier transformed
c     wrk(i)              i=1,6*max(n1,n2,n3)+15
c     .                   work array
c     lwrk                length of wrk
c I   mode                determines transform mode
c     .                   mode > 0   forward
c     .                   mode < 0   backward
c
c     rcs id string - allows use of ident to identify binaries
c
      character rcsid*50
      rcsid = '$RCSfile: cfft.f,v $$Revision: 1.1.1.1 $'
c     
      np = 2 * max0(n1,n2,n3) + 1
      if (lwrk .lt. 3 * np + 12) then
	write (*,*) 'work array in cfft is too short'
	stop
      end if

      call dcffti(n3,wrk(np))
      do 130 i=1,n1
         do 120 j=1,n2
cdir$       ivdep
            do 100 k=1,n3
               wrk(2*k-1) = chdr(i,j,k)
               wrk(2*k)   = chdi(i,j,k)
  100       continue
            if (mode .gt. 0) call dcfftf(n3,wrk,wrk(np))
            if (mode .lt. 0) call dcfftb(n3,wrk,wrk(np))
            do 110 k=1,n3
               chdr(i,j,k) = wrk(2*k-1)
               chdi(i,j,k) = wrk(2*k)
  110       continue
  120    continue
  130 continue
      if (n2 .ne. n3) call dcffti(n2,wrk(np))
      do 170 i=1,n3
         do 160 j=1,n1
cdir$       ivdep
            do 140 k=1,n2
               wrk(2*k-1) = chdr(j,k,i)
               wrk(2*k)   = chdi(j,k,i)
  140       continue
            if (mode .gt. 0) call dcfftf(n2,wrk,wrk(np))
            if (mode .lt. 0) call dcfftb(n2,wrk,wrk(np))
            do 150 k=1,n2
               chdr(j,k,i) = wrk(2*k-1)
               chdi(j,k,i) = wrk(2*k)
  150       continue
  160    continue
  170 continue
      if (n1 .ne. n2) call dcffti(n1,wrk(np))
      do 210 i=1,n2
         do 200 j=1,n3
cdir$       ivdep
            do 180 k=1,n1
               wrk(2*k-1) = chdr(k,i,j)
               wrk(2*k)   = chdi(k,i,j)
  180       continue
            if (mode .gt. 0) call dcfftf(n1,wrk,wrk(np))
            if (mode .lt. 0) call dcfftb(n1,wrk,wrk(np))
            do 190 k=1,n1
               chdr(k,i,j) = wrk(2*k-1)
               chdi(k,i,j) = wrk(2*k)
  190       continue
  200    continue
  210 continue
c     
c     normalize
c     
ccc      if (mode .lt. 0) return
      if (mode .gt. 0) return
      factor = 1.D0 / DBLE(n1*n2*n3)
      do 240 i=1,n1
         do 230 j=1,n2
            do 220 k=1,n3
               chdr(i,j,k) = factor * chdr(i,j,k)
               chdi(i,j,k) = factor * chdi(i,j,k)
  220       continue
  230    continue
  240 continue
      return
      end
