      subroutine mch_kerk(w_in,w_out,workr_n)
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

      real*8 w_in(mr_n),w_out(mr_n)
      real*8 workr_n(mr_n)

      real*8 vi(mr_n)  !work array


      ng2_n=ngtotnod2(inode)
********************************************************
*** Kerker mixing, this is like the precoditioning,
*** and it adds the new component w_out-w_in into w_in
********************************************************
      do i=1,nr_n
      w_out(i)=w_out(i)-w_in(i)
      vi(i)=0.d0
      enddo

      workr_n = w_out

      call d3fft_real2(vi,workr_n,1,0)


      do i=1,ng2_n

********** in SLDA, v(q=0) is not zero, but has the information of
********** relative highs between up and down, so this point should be
********** treated carefully.
       if(inode.eq.iorg2(1).and.i.eq.iorg2(2)) then
       ig = iorg2(2)
       vi(ig*2)=0.3d0*vi(ig*2)       
       vi(ig*2-1)=0.3d0*vi(ig*2-1)
        else
        vi(i*2)= 0.8d0*vi(i*2)*gkk2_n(i)/(gkk2_n(i)+0.5d0)
        vi(i*2-1)= 0.8d0*vi(i*2-1)*gkk2_n(i)/(gkk2_n(i)+0.5d0)
       endif

      enddo


      call d3fft_real2(vi,workr_n,-1,0)

      w_out = workr_n

      do i=1,nr_n
      w_in(i)=w_in(i)+w_out(i)
      enddo
**************************************************
      return
      end
      

