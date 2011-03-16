      subroutine mch_pulay(w_in,w_out,nint0,AA,nreset,islda)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************
ccccc  Now, this is in n1L,n2L,n3L


      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)
      include 'param.escan_real'

      real*8 w_in(mr_nL,islda),w_out(mr_nL,islda)

      real*8 AA(npulay_max,npulay_max)
      real*8, allocatable, dimension(:,:)  :: AA1
      real*8, allocatable, dimension(:)  :: B
      integer nreset

      allocate(AA1(npulay_max,npulay_max))
      allocate(B(npulay_max))

      alpha2=1.d0
cccc alpha2 controls how many recent charge densities
cccc to be used. If alpha2 > 1, then, the very old
cccc charge density is not used effectivly. 
cccc We find that alpha2=1 is O.K.
******************************************************
      if(nint0.eq.1) nreset=0
      nint=nint0-nreset

      if(nint.gt.npulay_max) then
      write(6,*) "restart pulay, nint0,npulay_max",
     &    nint0,npulay_max
      nreset=nint0-1
      nint=1
      endif

      if(nint.eq.1) then

      do iislda=1,islda
      do i=1,nr_nL
      R0(i,iislda)=w_out(i,iislda)-w_in(i,iislda)
      w_in0(i,iislda)=w_in(i,iislda)
      enddo
      enddo

      endif
******************************************************
      if(nint.gt.1) then

      do iislda=1,islda
      do i=1,nr_nL
      dw(i,nint-1,iislda)=w_in(i,iislda)-w_in0(i,iislda)
      dR(i,nint-1,iislda)=w_out(i,iislda)-w_in(i,iislda)-
     &                 R0(i,iislda)
      R0(i,iislda)=w_out(i,iislda)-w_in(i,iislda)
      w_in0(i,iislda)=w_in(i,iislda)
      enddo
      enddo


      do m=1,nint-1
      s=0.d0

      do iislda=1,islda
      do i=1,nr_nL
      s=s+dR(i,m,iislda)*R0(i,iislda)
      enddo
      enddo

      call global_sumr(s)

      s=s*vol/nrL
      B(m)=-s
      enddo

******************************************************
      do m=1,nint-1
      s1=0.d0

      do iislda=1,islda
      do i=1,nr_nL
      s1=s1+dR(i,m,iislda)*dR(i,nint-1,iislda)
      enddo
      enddo

      call global_sumr(s1)

      s1=s1*vol/nrL
      AA(m,nint-1)=s1
      AA(nint-1,m)=s1
      enddo

cccccccccc pulay optimization
**********************************************************
      do m1=1,nint-1
      do m2=1,nint-1
      AA1(m1,m2)=AA(m1,m2)
      enddo
      enddo

      w=1.0d0
      do m=nint-1,1,-1
      AA1(m,m)=AA1(m,m)*w
      w=w*alpha2
      enddo

**********************************************************

      call gaussj(AA1,nint-1,npulay_max,B,1,1)

      do iislda=1,islda
      do i=1,nr_nL
      w_out(i,iislda)=R0(i,iislda)
      enddo
      enddo

      
      do m=1,nint-1
      do iislda=1,islda
      do i=1,nr_nL
      w_in(i,iislda)=w_in(i,iislda)+B(m)*dw(i,m,iislda)
      w_out(i,iislda)=w_out(i,iislda)+B(m)*dR(i,m,iislda)
      enddo
      enddo
      enddo

      do iislda=1,islda
      do i=1,nr_nL
      w_out(i,iislda)=w_out(i,iislda)+w_in(i,iislda)
      enddo
      enddo

      endif

      deallocate(AA1)
      deallocate(B)

      return
      end
      

