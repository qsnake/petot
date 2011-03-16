      subroutine system_csfft(init,isign,n,scale,x,y,table,work,
     &              isys,ntable,nwork)
*************************************************************************
*** Written by Lin-Wang Wang, 2001
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

cccc input x (complex)
cccc output y (real)
      implicit none
      integer init,isign,n,isys,ntable,nwork
      real*8 scale,table(1),work(1)
      real*8 x(1),y(1)
ccccccc  complex to real 1D FFT

cccccc for T3E
c      if(init.eq.1) then
c      call csfft(0,n,scale,0,0,table,0,0)
c      else
c      call csfft(isign,n,scale,x,y,table,work,isys)
c      endif
cccccccccccccccccccccc
cccccc for IBM SP2  
      if(init.eq.1) then
      call dcrft(1,0,0,0,0,
     &    n,1,-isign,scale,table,ntable,0,0)
      else
      call dcrft(0,x,0,y,0,n,1,-isign,scale,
     &    table,ntable,work,nwork) 
      endif
cccccccccccccccccccccccccccccccccccc
ccccc    for Intel KML
c       if(init.eq.1) then
c       call zdfft1d(y,n,0,table)
c       else
c       call zdfft1d(x,n,1,table)    ! isign, or -isign
ccccc call zdfft1d(x,n,-1,table) is wrong ! The same result as zdfft1d(x,n,1,table)
c        if(isign.eq.-1) then
c        y(1)=x(1)*n*scale
c        do i=2,n
c        y(i)=x(n+2-i)*n*scale
c        enddo
c        else
c        do i=1,n
c        y(i)=x(i)*n*scale
c        enddo
c        endif
c       endif
ccccccccccccccccccccccccccccccccccccccccc

      return
      end

 
      
     

