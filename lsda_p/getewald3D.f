      subroutine getewald3D(fatom,xatom,AL,
     &    ityatom,ewald)
***************************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************

ccccc 
ccccc The final result is fatom=+dE/dR

       use fft_data
       use load_data
       use data

      implicit double precision (a-h,o-z)
      include 'param.escan_real'
      include "mpif.h"

      real*8 AL(3,3)
      real*8 xatom(3,matom),fatom(3,matom)
      real*8 zatom(mtype)
      integer ityatom(matom)
      integer, allocatable, dimension(:) :: ncount_tmp
      integer icoul
      real*8 xcoul(3)

      common /comzatom/zatom
      common /comcoul/icoul,xcoul

cccccccccccccccccccccccccccccccccccccccccccccc
ccccc assigne atoms to each node
     
      ewald=0.d0

      do 2000 ia1=1,natom

      fatom(1,ia1)=0.d0
      fatom(2,ia1)=0.d0
      fatom(3,ia1)=0.d0

      x1=xatom(1,ia1)
      if(x1.ge.1.d0) x1=x1-1.d0
      if(x1.lt.0.d0) x1=x1+1.d0
      if(x1.gt.xcoul(1)) x1=x1-1.d0
      y1=xatom(2,ia1)
      if(y1.ge.1.d0) y1=y1-1.d0
      if(y1.lt.0.d0) y1=y1+1.d0
      if(y1.gt.xcoul(2)) y1=y1-1.d0
      z1=xatom(3,ia1)
      if(z1.ge.1.d0) z1=z1-1.d0
      if(z1.lt.0.d0) z1=z1+1.d0
      if(z1.gt.xcoul(3)) z1=z1-1.d0

      do 1000 ia2=1,natom

      if(ia2.eq.ia1) goto 1000

      x2=xatom(1,ia2)
      if(x2.ge.1.d0) x2=x2-1.d0
      if(x2.lt.0.d0) x2=x2+1.d0
      if(x2.gt.xcoul(1)) x2=x2-1.d0
      y2=xatom(2,ia2)
      if(y2.ge.1.d0) y2=y2-1.d0
      if(y2.lt.0.d0) y2=y2+1.d0
      if(y2.gt.xcoul(2)) y2=y2-1.d0
      z2=xatom(3,ia2)
      if(z2.ge.1.d0) z2=z2-1.d0
      if(z2.lt.0.d0) z2=z2+1.d0
      if(z2.gt.xcoul(3)) z2=z2-1.d0

      dx=x1-x2     ! the definition is different from getewald
      dy=y1-y2
      dz=z1-z2

      x0=AL(1,1)*dx+AL(1,2)*dy+AL(1,3)*dz
      y0=AL(2,1)*dx+AL(2,2)*dy+AL(2,3)*dz
      z0=AL(3,1)*dx+AL(3,2)*dy+AL(3,3)*dz

      r=dsqrt(x0**2+y0**2+z0**2)

      ch1=zatom(ityatom(ia1))
      ch2=zatom(ityatom(ia2))

      ewald=ewald+0.5*ch1*ch2/r
      fatom(1,ia1)=fatom(1,ia1)-ch1*ch2*x0/r**3
      fatom(2,ia1)=fatom(2,ia1)-ch1*ch2*y0/r**3
      fatom(3,ia1)=fatom(3,ia1)-ch1*ch2*z0/r**3

1000  continue
2000  continue


      return
      end



  
