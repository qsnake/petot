      subroutine ration (n, v, maxd)
c     
c     replace each component of the vector v
c     with the closest rational number which has
c     a denominator less than or equal to maxd
c     $Id: ration.x,v 1.2 89/04/14 19:16:20 sverre Exp $
c
c     $Log:	ration.x,v $
c     Revision 1.2  89/04/14  19:16:20  sverre
c     checked in with -k by jimb at 89.05.02.18.03.32.
c     
c     Revision 1.2  89/04/14  19:16:20  sverre
c     Corrected value for delta.
c     
c     Revision 1.1  89/03/11  19:09:20  sverre
c     Initial revision
c     
      integer n, maxd
      DOUBLE PRECISION v(n)
c     
c     input:
c     
c     n        length of vector v
c     v(n)     vector to adjust
c     maxd     is the largest denominator to check
c     
c     output:
c     
c     v(n)     adjusted vector
c     maxd     is largest denominator needed, zero if none found
c
      parameter ( zero = 0.D0, one = 1.D0 )
c     
c     local variables
c     
      integer maxdt, i, j
      DOUBLE PRECISION delta, vrat
c
c     rcs id string - allows use of ident to identify binaries
c
      character rcsid*50
      rcsid = '$RCSfile: ration.x,v $$Revision: 1.2 $'
c
      delta = one / DBLE (maxd * maxd)
      maxdt = 0
c     
c     for each element in v
c     
      do 110 i = 1,n
c     
c     choose lowest possible denominator
c     
         do 100 j = 1,maxd
            vrat = DBLE (j)
            vrat = anint (vrat * v(i)) / vrat
            if (abs (vrat - v(i)) .lt. delta) then
               v(i) = vrat
               if (maxdt .lt. j) maxdt = j
               goto 110
            end if
  100    continue
c
c     failed
c
         maxdt = maxd + 1
  110 continue
      if (maxdt .eq. maxd + 1) then
         maxd = 0
      else
         maxd = maxdt
      end if
      return
      end
