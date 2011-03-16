      subroutine prdiff(nconf,econf)
c   
c   Prints out the energy differences between 
c   different atomic configurations.
c
c   njtj  ***  modifications  ***
c     econf is able to handle larger numbers 
c     of configurations.
c   njtj  ***  modifications  ***                              
c
c  njtj
c  ###  Cray conversions  
c  ###    1)Comment out implicit double precision.
c  ###  Cray conversions
c  njtj
c  
      implicit double precision (a-h,o-z)
c
      dimension econf(100)
c
      write(6,10) (i,i=1,nconf)
      do 30 i=1,nconf
        write(6,20) i,(econf(i)-econf(j),j=1,i)
 30   continue
 10   format(/,' total energy difference',//,2x,9i9)
 20   format(1x,i2,1x,9f9.4)
      return
      end
