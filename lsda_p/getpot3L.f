      subroutine getpot3L(rho,vion,rhocr,vtot,v0,E_Hxc,
     &  E_coul,E_ion)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************


      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)
      include 'param.escan_real'

      real*8 vtot(mr_nL,2)
      real*8 rho(mr_nL,2),vion(mr_nL),rhocr(mr_nL)

      real*8 vi(mr_nL)

**************************************************
ccccc generate the potential vtot from rho(i)

      ng2_nL=ngtotnod2L(inode)
      nr_nL=n1L*n2L*n3L/nnodes
      nrL=n1L*n2L*n3L




      do i=1,nr_nL
      vtot(i,2) = rho(i,1)+rho(i,2)     ! total charge for Coulomb potential
      enddo

      call getV_Hartree(vtot(1,2),vtot(1,1))

      vtot(:,2)=vtot(:,1)

      
      s=0.d0
      E_Hxc=0.d0
      E_coul=0.d0
      E_ion=0.d0


      do i=1,nr_nL
      E_Hxc=E_Hxc+0.5d0*(vtot(i,1)*rho(i,1)+
     &                   vtot(i,2)*rho(i,2))
      E_coul=E_coul+0.5d0*(vtot(i,1)*rho(i,1)+
     &                   vtot(i,2)*rho(i,2))
      E_ion=E_ion+vion(i)*(rho(i,1)+rho(i,2))
      call UxcCA2(rho(i,1)+rhocr(i)*0.5d0,rho(i,2)+rhocr(i)*0.5d0,
     &    vxc1,vxc2,uxc1,uxc2)
      vtot(i,1)=vtot(i,1)+vion(i)+vxc1
      vtot(i,2)=vtot(i,2)+vion(i)+vxc2
      E_Hxc=E_Hxc+uxc1*(rho(i,1)+rhocr(i)*0.5d0)+
     &         uxc2*(rho(i,2)+rhocr(i)*0.5d0)
      s=s+(vtot(i,1)+vtot(i,2))*0.5d0
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
      vtot(i,1)=vtot(i,1)-s
      vtot(i,2)=vtot(i,2)-s
      enddo

      v0=s
ccccccccccccccccccccc
cccccc ave_vtot=0.d0 is necessary for mch_Broyden and mch_pulay,
cccccc otherwise it doesn't work
ccccccccccccccccccccc


      return
      end
      

