      subroutine mch_pulay(w_in,w_out,nint0,AA,nreset,islda)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
******************************************


      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)
      include 'param.escan_real'

      real*8 w_in(mr_n,islda),w_out(mr_n,islda)

      real*8 AA(nmax,nmax),B(nmax)
      real*8 AA1(nmax,nmax)
      integer nreset

      alpha2=1.d0
cccc alpha2 controls how many recent charge densities
cccc to be used. If alpha2 > 1, then, the very old
cccc charge density is not used effectivly. 
cccc We find that alpha2=1 is O.K.
******************************************************
      if(nint0.eq.1) nreset=0
      nint=nint0-nreset

      if(nint.gt.nmax) then
      write(6,*) "restart pulay, nint0,nmax",nint0,nmax
      nreset=nint0-1
      nint=1
      endif

      if(nint.eq.1) then

      do iislda=1,islda
      do i=1,nr_n
      R0(i,iislda)=w_out(i,iislda)-w_in(i,iislda)
      w_in0(i,iislda)=w_in(i,iislda)
      enddo
      enddo

      endif
******************************************************
      if(nint.gt.1) then

      do iislda=1,islda
      do i=1,nr_n
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
      do i=1,nr_n
      s=s+dR(i,m,iislda)*R0(i,iislda)
      enddo
      enddo

      call global_sumr(s)

      s=s*vol/nr
      B(m)=-s
      enddo

******************************************************
      do m=1,nint-1
      s1=0.d0

      do iislda=1,islda
      do i=1,nr_n
      s1=s1+dR(i,m,iislda)*dR(i,nint-1,iislda)
      enddo
      enddo

      call global_sumr(s1)

      s1=s1*vol/nr
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

      call gaussj(AA1,nint-1,nmax,B,1,1)

      do iislda=1,islda
      do i=1,nr_n
      w_out(i,iislda)=R0(i,iislda)
      enddo
      enddo

      
      do m=1,nint-1
      do iislda=1,islda
      do i=1,nr_n
      w_in(i,iislda)=w_in(i,iislda)+B(m)*dw(i,m,iislda)
      w_out(i,iislda)=w_out(i,iislda)+B(m)*dR(i,m,iislda)
      enddo
      enddo
      enddo

      do iislda=1,islda
      do i=1,nr_n
      w_out(i,iislda)=w_out(i,iislda)+w_in(i,iislda)
      enddo
      enddo

      endif

      return
      end
      

