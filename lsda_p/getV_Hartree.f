      subroutine getV_Hartree(rho,vtot)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************
***** All the arrays are in nr_nL in this subroutine


      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)
      include 'param.escan_real'

      real*8 vtot(mr_nL),rho(mr_nL)
      real*8, allocatable, dimension(:) :: workr_nL,workr_nL2,vtmp_nL2
      integer icoul
      real*8 xcoul(3)
      
      common /comcoul/icoul,xcoul

**************************************************
ccccc generate the potential vtot from rho(i)

      if(icoul.eq.0) then

      allocate(workr_nL(mr_nL))

      ng2_nL=ngtotnod2L(inode)
      nr_nL=n1L*n2L*n3L/nnodes
      nrL=n1L*n2L*n3L

      vtot = 0.0d0

      do i=1,nr_nL
      workr_nL(i) = rho(i)
      enddo

      call d3fft_real2L(vtot,workr_nL,1,0)


      pi=4*datan(1.d0)

      do 100 i=1,ng2_nL
 
      if(inode.eq.iorg2L(1).and.i.eq.iorg2L(2)) then
      ig = iorg2L(2)
      vtot(ig*2)=0.d0
      vtot(ig*2-1)=0.d0
      else
      vtot(i*2)=vtot(i*2)*2*pi/gkk2_nL(i)
      vtot(i*2-1)=vtot(i*2-1)*2*pi/gkk2_nL(i)
      endif

100   continue

      call d3fft_real2L(vtot,workr_nL,-1,0)

      do i=1,nr_nL
      vtot(i) = workr_nL(i)
      enddo

      deallocate(workr_nL)
      return
      endif

cccccccccccccccccccccccccccccccccccccccccccc

      if(icoul.eq.1.or.icoul.eq.11.or.
     &  icoul.eq.12.or.icoul.eq.13) then

      allocate(workr_nL2(mr_nL2))

      call  convert_2LtoL(rho,workr_nL2,1)

      allocate(vtmp_nL2(mr_nL2))


      call d3fft_real2L2(vtmp_nL2,workr_nL2,1,0)

      ng2_nL2=ngtotnod2L2(inode)

      do i=1,ng2_nL2
      i2=i*2

       vreal= vcoul_nL2(i2-1,0)*vtmp_nL2(i2-1)-
     &        vcoul_nL2(i2,0)*vtmp_nL2(i2)

       vimag= vcoul_nL2(i2-1,0)*vtmp_nL2(i2)+
     &        vcoul_nL2(i2,0)*vtmp_nL2(i2-1)
       vtmp_nL2(i2-1)=vreal
       vtmp_nL2(i2)=vimag

       enddo

      call d3fft_real2L2(vtmp_nL2,workr_nL2,-1,0)

      deallocate(vtmp_nL2)

      call  convert_2LtoL(vtot,workr_nL2,2)

      deallocate(workr_nL2)

      return

      endif      ! icoul.eq.1
ccccccccccccccccccccccccccccccccccccccccccc     
      return

      end

