      subroutine metric (a, am)
c     
c     compute metric am from unit cell vectors a
c     $Id: metric.x,v 1.1 89/03/11 19:09:10 sverre Exp $
c
c     $Log:	metric.x,v $
c     Revision 1.1  89/03/11  19:09:10  sverre
c     Initial revision
c     
c     Revision 1.1  89/03/11  19:09:10  sverre
c     Initial revision
c     
c     
      DOUBLE PRECISION a(3,3), am(3,3)
c     
c     input:
c     
c     a(3,3)   unit cell vectors stored as column vectors
c     
c     output:
c     
c     am(n,n)  metric for computing lengths
c     am(i,j) = a(*,i) dot a(*,j)
c     
c     local variables
c     
      integer i, j, k
c
c     rcs id string - allows use of ident to identify binaries
c
      character rcsid*50
      rcsid = '$RCSfile: metric.x,v $$Revision: 1.1 $'
c     
c     find metric
c     
      do 120 i = 1,3
         do 110 j = 1,3
            am(i,j) = 0.D0
            do 100 k = 1,3
               am(i,j) = am(i,j) + a(k,i) * a(k,j)
  100       continue
  110    continue
  120 continue
      return
      end
