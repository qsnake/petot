      subroutine dottvv(a,b,c,n)
c
      implicit real*8 (a-h,o-z)
      parameter (zero=0.0d0)
c
      dimension a(n),b(n)
c
      c=zero
c
      do 10 i=1,n
        c=c+a(i)*b(i)
 10   continue
      return
      end
