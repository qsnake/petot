      subroutine system_scfft(init,isign,n,scale,x,y,table,work,
     &              isys,ntable,nwork)
*************************************************************************
*** Written by Lin-Wang Wang, 2001
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

cccc input: x  (real)
cccc output: y  (complex)
      implicit none
      integer init,isign,n,isys,ntable,nwork
      real*8 scale,table(1),work(1)
      real*8 x(1),y(1)
ccccccc real to complex 1D fft

cccccc for T3E
c      if(init.eq.1) then
c      call scfft(0,n,scale,0,0,table,0,0)
c      else
c      call scfft(isign,n,scale,x,y,table,work,isys)
c      endif
cccccccccccccccccccccc
cccccc for IBM SP2, essl  
      if(init.eq.1) then
      call drcft(1,0,0,0,0,
     &    n,1,-isign,scale,table,ntable,0,0)
      else
      call drcft(0,x,0,y,0,n,1,-isign,scale,
     &    table,ntable,work,nwork) 
      endif
cccccccccccccccccccccccccccccccccccc
cccccc for Intel KML
c       if(init.eq.1) then
c       call dzfft1d(y,n,0,table)
c       else
c       do i=1,n
c       y(i)=x(i)*scale
c       enddo
c       call dzfft1d(y,n,-1,table)    ! or -isign
ccccc dzfft1d(y,n,1,table) and dzfft1d(y,n,-1,table) are the same !
c       if(isign.eq.1) then
c       do i=2,n+2,2
c       y(i)=-y(i)
c       enddo
c       endif
c       endif
cccccccccccccccccccccccccccccccccccccccccccc


      return
      end

 
      
     

