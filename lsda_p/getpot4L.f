      subroutine getpot4L(rho,vion,rhocr,vtot,v0,E_Hxc,
     & E_coul,E_ion)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************

***********************************************************
**** This subroutine calculates the GGA, using PBE-GGA, 
**** by J.P. Perdew, K. Burke and M. Ernzerhof, Phys. Rev. Lett. 77, 3865
**** (1996); ibid. 78, 1396 (1997) (E).
**** It assumes that spin up = spin down, unpolarized calculation
*************************************

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)
      include 'param.escan_real'
      include 'mpif.h'

      real*8 workr_n(mr_nL)
      real*8 vtot(mr_nL)
      real*8 rho(mr_nL),vion(mr_nL),rhocr(mr_nL)


      real*8, allocatable, dimension(:) :: vi,vitmp

      real*8, allocatable, dimension(:) :: drho,ddrho,dotrho

      real*8, allocatable, dimension(:,:) :: vtmp

**************************************************
**************************************************
**************************************************

      call getV_Hartree(rho,vtot)

****** finished the Coulomb potential
**************************************************


ccccc generate the potential vtot from rho(i)
ccccc (mr_n,3): 1, spin up; 2, spin down; 3 spin up+spin down

      allocate(drho(mr_nL))         !  |grad(rho)|
      allocate(dotrho(mr_nL))       !  [grad(rho)]\cdot[grad|grad(rho)|]
      allocate(ddrho(mr_nL))        !  \nablda^2(rho)

ccccccccc drho,dotrho,ddrho are all we need to call the one r point gga subroutines

      drho=0.d0
      dotrho=0.d0
      ddrho=0.d0

      allocate(vi(mr_nL))
      allocate(vitmp(mr_nL))
      allocate(vtmp(mr_nL,3))

****************************************************


      ng2_nL=ngtotnod2L(inode)
      nr_nL=n1L*n2L*n3L/nnodes
      nrL=n1L*n2L*n3L

*************************************************************
*************************************************************
****** now,  calculate drho, ddrho, dotrho
*************************************************************
      rho_max=0.d0
      do i=1,nr_nL
      if(rho(i).gt.rho_max) rho_max=rho(i)
      enddo
  
      call mpi_allreduce(rho_max,rho_maxtmp,1,MPI_REAL8,
     &  MPI_MAX,MPI_COMM_K,ierr)
      rho_max=rho_maxtmp
cccccccccccccccccccccccccccccccccccccccccccccccccccccc


      do i=1,nr_nL
      workr_n(i)=rho(i)+rhocr(i)
      enddo
     
      call d3fft_real2L(vi,workr_n,1,0)

      do 201 ixyz=1,3

      call vitmp_eq_vigkxyz()     ! multiply i*gkx,i*gky,i*gkz to vi

      call d3fft_real2L(vitmp,workr_n,-1,0)

      do i=1,nr_nL
      vtmp(i,ixyz)=workr_n(i)
      enddo

201   continue



      do i=1,nr_nL
      drho(i)=dsqrt(vtmp(i,1)**2+
     &  vtmp(i,2)**2+vtmp(i,3)**2)
      enddo
      
******************************************
      do i=1,ng2_nL
      vitmp(i*2)=2*vi(i*2)*gkk2_nL(i)
      vitmp(i*2-1)=2*vi(i*2-1)*gkk2_nL(i)
      enddo
      call d3fft_real2L(vitmp,workr_n,-1,0)
      do i=1,nr_nL
      ddrho(i)=-workr_n(i)
      enddo
*******************************************
      do i=1,nr_nL
      workr_n(i)=drho(i)
      enddo
     
      call d3fft_real2L(vi,workr_n,1,0)

      do 202 ixyz=1,3

      call vitmp_eq_vigkxyz()     ! multiply i*gkx,i*gky,i*gkz to vi

      call d3fft_real2L(vitmp,workr_n,-1,0)

      do i=1,nr_nL
      dotrho(i)=dotrho(i)+workr_n(i)*vtmp(i,ixyz)
      enddo
202   continue

********************************************
      deallocate(vi)
      deallocate(vitmp)
      deallocate(vtmp)


*********************************************************
**** finished the calculation of drho, ddrho and dotrho,
**** now, we can call GGA subroutine one r point at a time
*********************************************************
ccc      drho(mr_n)        !  |grad(rho)|
ccc      dotrho(mr_n)       !  [grad(rho)]\cdot[grad|grad(rho)|]
ccc      ddrho(mr_n)        !  \nablda^2(rho)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      s=0.d0
      E_Hxc=0.d0
      E_coul=0.d0
      E_ion=0.d0

      do 1000 i=1,nr_nL

cccccccccccccccccccccccccccccccccccccccccccc
cccccccccc Hartree energy

      E_Hxc=E_Hxc+0.5d0*vtot(i)*rho(i)
      E_coul=E_coul+0.5d0*vtot(i)*rho(i)
      E_ion=E_ion+vion(i)*rho(i)


       up=(rho(i)+rhocr(i))*0.5d0
       agrup=drho(i)*0.5d0
       delgrup=dotrho(i)*0.25d0
       uplap=ddrho(i)*0.5d0
       dn=(rho(i)+rhocr(i))*0.5d0
       agrdn=drho(i)*0.5d0
       delgrdn=dotrho(i)*0.25d0
       dnlap=ddrho(i)*0.5d0
       agr=drho(i)
       delgr=dotrho(i)
       lcor=1
       lpot=1

      call easypbe(up,agrup,delgrup,uplap,dn,agrdn,delgrdn,dnlap,
     &           agr,delgr,lcor,lpot,
     &           exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd,
     &           expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91,
     &           expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)

      if(up.lt.rho_max*1.D-4) then        ! special, make it stable
       vxuppbe=0.d0
       vcuppbe=0.d0
       vxdnpbe=0.d0
       vcdnpbe=0.d0
       expbe=0.d0
       ecpbe=0.d0
       endif


ccccccccc dirty trick 
      if(rho(i).gt.rho_max*1.D-2) then
      vtot(i)=vtot(i)+vion(i)+vxuppbe+vcuppbe
      E_Hxc=E_Hxc+(expbe+ecpbe)*(rho(i)+rhocr(i))
      else        ! use LDA formula
      vtot(i)=vtot(i)+vion(i)+UxcCA(rho(i)+rhocr(i),uxc2)
      E_Hxc=E_Hxc+uxc2*(rho(i)+rhocr(i))
      endif

      s=s+vtot(i)
1000  continue

****************************************************
      deallocate(drho)      
      deallocate(dotrho)    
      deallocate(ddrho)     
****************************************************
****************************************************

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

      contains
      
      subroutine vitmp_eq_vigkxyz()

      implicit double precision (a-h,o-z)

        if(ixyz.eq.1) then
        do i=1,ng2_nL
        vitmp(i*2)=-vi(i*2-1)*gkx2_nL(i)
        vitmp(i*2-1)=vi(i*2)*gkx2_nL(i)
        enddo
        endif
        if(ixyz.eq.2) then
        do i=1,ng2_nL
        vitmp(i*2)=-vi(i*2-1)*gky2_nL(i)
        vitmp(i*2-1)=vi(i*2)*gky2_nL(i)
        enddo
        endif
        if(ixyz.eq.3) then
        do i=1,ng2_nL
        vitmp(i*2)=-vi(i*2-1)*gkz2_nL(i)
        vitmp(i*2-1)=vi(i*2)*gkz2_nL(i)
        enddo
        endif

      end subroutine vitmp_eq_vigkxyz

      end
