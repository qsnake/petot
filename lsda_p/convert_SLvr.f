      subroutine convert_SLvr(vr_p,vrL_p,iflag)

******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************

ccc    iflag.eq.1, from vr_p to vrL_p
cc     iflag.eq.-1, from vrL_p to vr_p

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'param.escan_real'
      include 'mpif.h'

      real*8 vr_p(mr_n),vrL_p(mr_nL)
      complex*16,allocatable, dimension(:) :: vg_p,vgL_p
      real*8, allocatable, dimension(:)  :: workr_n

*******************************************************
cccccccccccccccccccccccccccccccccccccccccccc
      if(n1.eq.n1L.and.n2.eq.n2L.
     &                 and.n3.eq.n3L) then
       if(iflag.eq.1) vrL_p=vr_p
       if(iflag.eq.-1) vr_p=vrL_p
       return
      endif
cccccccccccccccccccccccccccccccccccccccccccc

       nrL=n1L*n2L*n3L
       nr_nL=nrL/nnodes


      allocate(vg_p(mr_n/2))
      allocate(vgL_p(mr_nL/2))


      if(iflag.eq.1) then     ! from vr_n to vrL_n
   
      allocate(workr_n(mr_n))

      workr_n=vr_p

      call d3fft_real2(vg_p,workr_n,1,0)         ! workr_n will be destroyed

      ng2_n=ngtotnod2(inode)

      vgL_p=dcmplx(0.d0,0.d0)
      do ig=1,ng2_n
      vgL_p(map_StoL(ig))=vg_p(ig)         ! the mapping is on the same node !
      enddo

      call d3fft_real2L(vgL_p,vrL_p,-1,0)

      deallocate(workr_n)

      endif

cccccccccccccccccccccccccccccccccccccccccccccc

      if(iflag.eq.-1) then       ! from vrL_n to vr_n

      allocate(workr_n(mr_nL))
      
      workr_n=vrL_p

      call d3fft_real2L(vgL_p,workr_n,1,0)

      ng2_n=ngtotnod2(inode)

      vg_p=dcmplx(0.d0,0.d0)
      do ig=1,ng2_n
      vg_p(ig)=vgL_p(map_StoL(ig))
      enddo

      call d3fft_real2(vg_p,vr_p,-1,0)
  
      deallocate(workr_n)

      endif
cccccccccccccccccccccccccccccccccccccccccccccc

      deallocate(vg_p)
      deallocate(vgL_p)

      return
      end
 


