      subroutine system_csfft(init,isign,n,scale,x,y,table,work,
     &              isys,ntable,nwork)
cccc input x (complex)
cccc output y (real)
      implicit none
      integer init,isign,n,isys,ntable,nwork
      real*8 scale,table(1),work(1)
      real*8 x(1),y(1)
ccccccc  complex to real 1D FFT

cccccc for T3E
      if(init.eq.1) then
      call csfft(0,n,scale,0,0,table,0,0)
      else
      call csfft(isign,n,scale,x,y,table,work,isys)
      endif
cccccccccccccccccccccc
cccccc for IBM SP2  
c      if(init.eq.1) then
c      call dcrft(1,0,0,0,0,
c     &    n,1,-isign,scale,table,ntable,0,0)
c      else
c      call dcrft(0,x,0,y,0,n,1,-isign,scale,
c     &    table,ntable,work,nwork) 
c      endif
cccccccccccccccccccccccccccccccccccc
      return
      end

 
      
     

