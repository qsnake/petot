      subroutine trnsvv(a,b,c,n) 
c
      implicit real*8 (a-h,o-z)
c
      dimension a(n),b(n) 
c
      do 10 i=1,n
        a(i)=a(i)+c*b(i)
 10   continue
      return
      end
