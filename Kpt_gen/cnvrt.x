      subroutine cnvrt (n, tflag, a, vin, vout)
c     
c     vector vout = matrix a times vector vin
c     $Id: cnvrt.x,v 1.1 89/03/11 19:09:04 sverre Exp $
c
c     $Log:	cnvrt.x,v $
c     Revision 1.1  89/03/11  19:09:04  sverre
c     Initial revision
c     
c     Revision 1.1  89/03/11  19:09:04  sverre
c     Initial revision
c     
c     
      integer n, tflag
      DOUBLE PRECISION a(n,n), vin(n), vout(n)
c     
c     input:
c     
c     n        dimension of matrix and vector
c     tflag    transpose flag (transpose matrix a if tflag positive)
c     a(n,n)   matrix
c     vin(n)   input vector 
c     
c     output:
c     
c     vout(n)  result vector
c     
c     local variables
c     
      integer i, j
c
c     rcs id string - allows use of ident to identify binaries
c
      character rcsid*50
      rcsid = '$RCSfile: cnvrt.x,v $$Revision: 1.1 $'
c     
c     initialize
c     
      do 100 i = 1,n
         vout(i) = 0.D0
  100 continue
c     
c     multiply with a
c     
      if (tflag .le. 0) then
         do 120 i = 1,n
            do 110 j = 1,n
               vout(j) = vout(j) + a(j,i) * vin(i)
  110       continue
  120    continue
      else
c     
c     or a transpose
c     
         do 140 i = 1,n
            do 130 j = 1,n
               vout(j) = vout(j) + a(i,j) * vin(i)
  130       continue
  140    continue
      endif
      return
      end
