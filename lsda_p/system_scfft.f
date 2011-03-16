      subroutine system_scfft(init,isign,n,scale,x,y,table,work,
     &              isys,ntable,nwork)
cccc input: x  (real)
cccc output: y  (complex)
      implicit none
      integer init,isign,n,isys,ntable,nwork
      real*8 scale,table(1),work(1)
      real*8 x(1),y(1)
ccccccc real to complex 1D fft

cccccc for T3E
      if(init.eq.1) then
      call scfft(0,n,scale,0,0,table,0,0)
      else
      call scfft(isign,n,scale,x,y,table,work,isys)
      endif
cccccccccccccccccccccc
cccccc for IBM SP2  !!!  NEED to be changed, dcft to ???
c      if(init.eq.1) then
c      call drcft(1,0,0,0,0,
c     &    n,1,-isign,scale,table,ntable,0,0)
c      else
c      call drcft(0,x,0,y,0,n,1,-isign,scale,
c     &    table,ntable,work,nwork) 
c      endif
cccccccccccccccccccccccccccccccccccc
      return
      end

 
      
     

