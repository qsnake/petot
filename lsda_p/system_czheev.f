      subroutine system_czheev(w1,w2,n,z,m,EE,workx,
     &  lwork,workrx,info)

      implicit none
      character*1 w1,w2
      integer n,m,lwork,info
      complex*16 workx(1),z(1)
      real*8 workrx(1),EE

ccccc for T3E, lapack routine
      call cheev(w1,w2,n,z,m,EE,workx,
     &  lwork,workrx,info)

ccccc for IBM SP2, lapack routine
c      call zheev(w1,w2,n,z,m,EE,workx,
c     &  lwork,workrx,info)

      return
      end


      
       

