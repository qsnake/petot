      subroutine system_czheev(w1,w2,n,z,m,EE,workx,
     &  lwork,workrx,info)
*************************************************************************
*** Written by Lin-Wang Wang, 2001
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

ccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      character*1 w1,w2
      integer n,m,lwork,info
      complex*16 workx(1),z(1)
      real*8 workrx(1),EE

ccccc for T3E, lapack routine (cheev is for complex*16 in T3E)
ccccc in Cray T3E, the single precision is defined as real*8 and complex*16
c      call cheev(w1,w2,n,z,m,EE,workx,
c     &  lwork,workrx,info)

ccccc for all the other machines, lapack routine
      call zheev(w1,w2,n,z,m,EE,workx,
     &  lwork,workrx,info)

      return
      end


      
       

