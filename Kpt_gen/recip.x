      subroutine recip (fac, a, b, vol)
c     
c     generate reciprocal space vectors
c     $Id: recip.x,v 1.1 89/03/11 19:09:03 sverre Exp $
c
c     $Log:	recip.x,v $
c     Revision 1.1  89/03/11  19:09:03  sverre
c     Initial revision
c     
c     Revision 1.1  89/03/11  19:09:03  sverre
c     Initial revision
c     
c     
      DOUBLE PRECISION fac, a(3,3), b(3,3), vol
c     
c     input:
c     
c     fac      scaling factor (see below)
c     a(3,3)   three linearly independent column vectors
c     
c     output:
c     
c     b(3,3)   three column vectors forming the reciprocal of a
c     .        i.e. satisfying a(*,i) dot b(*,j) = delta(i,j) * fac
c     vol      the volume of the unit cell formed by a.
c     
c     local variables
c     
      integer p(3), i, j, i1, i2, j1, j2
      DOUBLE PRECISION  delta, factor
c
c     rcs id string - allows use of ident to identify binaries
c
      character rcsid*50
      rcsid = '$RCSfile: recip.x,v $$Revision: 1.1 $'
c
      data p / 2, 3, 1 /
      data delta / 1.0D-8 /
c     
c     do vector products
c     
      do 110 i = 1,3
         i1 = p(i)
         i2 = p(i1)
         do 100 j = 1,3
            j1 = p(j)
            j2 = p(j1)
            b(j,i) = a(j1,i1)*a(j2,i2) - a(j2,i1)*a(j1,i2)
  100    continue
  110 continue
c     
c     find volume of unit cell a
c     
      vol = a(1,1)*b(1,1) + a(2,1)*b(2,1) + a(3,1)*b(3,1)
c
      if (abs(vol) .gt. delta) then
         factor = fac / vol
      else
         factor = fac
      end if
c     
c     rescale b
c     
      do 130 i = 1,3
         do 120 j = 1,3
            b(j,i) = factor * b(j,i)
  120    continue
  130 continue
c     
c     volume can be negative for left handed coordinates
c     
      vol = abs(vol)
      return
      end
