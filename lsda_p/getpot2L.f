      subroutine getpot2L(rho,vion,rhocr,vtot,v0,E_Hxc,
     & E_coul,E_ion)
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

      real*8 vtot(mr_nL)
      real*8 rho(mr_nL),vion(mr_nL),rhocr(mr_nL)

      real*8 vi(mr_nL)


**************************************************
ccccc generate the potential vtot from rho(i)

      call getV_Hartree(rho,vtot)

      ng2_nL=ngtotnod2L(inode)
      nr_nL=n1L*n2L*n3L/nnodes
      nrL=n1L*n2L*n3L


      s=0.d0
      E_Hxc=0.d0
      E_coul=0.d0
      E_ion=0.d0

      do i=1,nr_nL
      E_Hxc=E_Hxc+0.5d0*vtot(i)*rho(i)
      E_coul=E_coul+0.5d0*vtot(i)*rho(i)
      E_ion=E_ion+vion(i)*rho(i)
      vtot(i)=vtot(i)+vion(i)+UxcCA(rho(i)+rhocr(i),uxc2)

      E_Hxc=E_Hxc+uxc2*(rho(i)+rhocr(i))
      s=s+vtot(i)
      enddo
 

      call global_sumr(s)
      call global_sumr(E_Hxc)
      call global_sumr(E_coul)
      call global_sumr(E_ion)

      E_Hxc=E_Hxc*vol/nrL
      E_coul=E_coul*vol/nrL
      E_ion=E_ion*vol/nrL

      s=s/nrL
      do i=1,nr_nL
      vtot(i)=vtot(i)-s
      enddo

      v0=s
ccccccccccccccccccccc
cccccc ave_vtot=0.d0 is necessary for mch_Broyden and mch_pulay,
cccccc otherwise it doesn't work
ccccccccccccccccccccc


      return
      end
      

