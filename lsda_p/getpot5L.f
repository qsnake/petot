      subroutine getpot5L(rho,vion,rhocr,vtot,v0,E_Hxc,
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
**** It assumes that there are up and down spins.
*************************************

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)
      include 'param.escan_real'
      include 'mpif.h'


      real*8 workr_n(mr_nL)
      real*8 vtot(mr_nL,2)
      real*8 rho(mr_nL,2),vion(mr_nL),rhocr(mr_nL)


      real*8, allocatable, dimension(:) :: vi,vitmp

      real*8, allocatable, dimension(:,:) :: drho,ddrho,dotrho

      real*8, allocatable, dimension(:,:,:) :: vtmp

**************************************************
**************************************************
      rho_max=0.d0
      do i=1,nr_nL
      vtot(i,2) = rho(i,1)+rho(i,2)     ! total charge for Coulomb potential
      if(vtot(i,2).gt.rho_max) rho_max=vtot(i,2)
      enddo

      call mpi_allreduce(rho_max,rho_maxtmp,1,MPI_REAL8,
     &  MPI_MAX,MPI_COMM_K,ierr)
      rho_max=rho_maxtmp

      call getV_Hartree(vtot(1,2),vtot(1,1))

      vtot(:,2)=vtot(:,1)

****** finished the coulomb potential      
**************************************************
ccccc generate the potential vtot from rho(i)
ccccc (mr_n,3): 1, spin up; 2, spin down; 3 spin up+spin down

      allocate(drho(mr_nL,3))         !  |grad(rho)|
      allocate(dotrho(mr_nL,3))       !  [grad(rho)]\cdot[grad|grad(rho)|]
      allocate(ddrho(mr_nL,2))        !  \nablda^2(rho)

ccccccccc drho,dotrho,ddrho are all we need to call the one r point gga subroutines

      drho=0.d0
      dotrho=0.d0
      ddrho=0.d0

      allocate(vi(mr_nL))
      allocate(vitmp(mr_nL))
      allocate(vtmp(mr_nL,3,2))

****************************************************

      ng2_nL=ngtotnod2L(inode)
      nr_nL=n1L*n2L*n3L/nnodes
      nrL=n1L*n2L*n3L

*************************************************************
****** now,  calculate drho, ddrho, dotrho
*************************************************************

      do 200 iislda=1,2

      do i=1,nr_nL
      workr_n(i)=rho(i,iislda)+rhocr(i)*0.5d0
      enddo
     
      call d3fft_real2L(vi,workr_n,1,0)

      do 201 ixyz=1,3

      call vitmp_eq_vigkxyz()     ! multiply i*gkx,i*gky,i*gkz to vi

      call d3fft_real2L(vitmp,workr_n,-1,0)

      do i=1,nr_nL
      vtmp(i,ixyz,iislda)=workr_n(i)
      enddo

201   continue



      do i=1,nr_nL
      drho(i,iislda)=dsqrt(vtmp(i,1,iislda)**2+
     &  vtmp(i,2,iislda)**2+vtmp(i,3,iislda)**2)
      enddo
      
******************************************
      do i=1,ng2_nL
      vitmp(i*2)=2*vi(i*2)*gkk2_nL(i)
      vitmp(i*2-1)=2*vi(i*2-1)*gkk2_nL(i)
      enddo
      call d3fft_real2L(vitmp,workr_n,-1,0)
      do i=1,nr_nL
      ddrho(i,iislda)=-workr_n(i)
      enddo
*******************************************
      do i=1,nr_nL
      workr_n(i)=drho(i,iislda)
      enddo
     
      call d3fft_real2L(vi,workr_n,1,0)

      do 202 ixyz=1,3

      call vitmp_eq_vigkxyz()     ! multiply i*gkx,i*gky,i*gkz to vi

      call d3fft_real2L(vitmp,workr_n,-1,0)

      do i=1,nr_nL
      dotrho(i,iislda)=dotrho(i,iislda)+workr_n(i)*vtmp(i,ixyz,iislda)
      enddo
202   continue

200   continue
********************************************

      do i=1,nr_nL
      drho(i,3)=dsqrt((vtmp(i,1,1)+vtmp(i,1,2))**2+
     &                (vtmp(i,2,1)+vtmp(i,2,2))**2+
     &                (vtmp(i,3,1)+vtmp(i,3,2))**2)
      enddo
      

      do i=1,nr_nL
      workr_n(i)=drho(i,3)
      enddo
     
      call d3fft_real2L(vi,workr_n,1,0)

      do 203 ixyz=1,3

      call vitmp_eq_vigkxyz()     ! multiply i*gkx,i*gky,i*gkz to vi

      call d3fft_real2L(vitmp,workr_n,-1,0)

      do i=1,nr_nL
      dotrho(i,3)=dotrho(i,iislda)+workr_n(i)*
     &     (vtmp(i,ixyz,1)+vtmp(i,ixyz,2))
      enddo

203   continue

*********************************************************
      deallocate(vi)
      deallocate(vitmp)
      deallocate(vtmp)
*********************************************************
**** finished the calculation of drho, ddrho and dotrho,
**** now, we can call GGA subroutine one r point at a time
*********************************************************
ccccc (mr_n,3): 1, spin up; 2, spin down; 3 spin up+spin down
ccc
ccc      drho(mr_n,3)        !  |grad(rho)|
ccc      dotrho(mr_n,3)       !  [grad(rho)]\cdot[grad|grad(rho)|]
ccc      ddrho(mr_n,2)        !  \nablda^2(rho)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      s=0.d0
      E_Hxc=0.d0
      E_coul=0.d0
      E_ion=0.d0

      do 1000 i=1,nr_nL

cccccccccccccccccccccccccccccccccccccccccccc
cccccccccc Hartree energy

      E_Hxc=E_Hxc+0.5d0*(vtot(i,1)*rho(i,1)+
     &                   vtot(i,2)*rho(i,2))
      E_coul=E_coul+0.5d0*(vtot(i,1)*rho(i,1)+
     &                   vtot(i,2)*rho(i,2))
      E_ion=E_ion+vion(i)*(rho(i,1)+rho(i,2))

       up=rho(i,1)+rhocr(i)*0.5d0
       agrup=drho(i,1)
       delgrup=dotrho(i,1)
       uplap=ddrho(i,1)
       dn=rho(i,2)+rhocr(i)*0.5d0
       agrdn=drho(i,2)
       delgrdn=dotrho(i,2)
       dnlap=ddrho(i,2)
       agr=drho(i,3)
       delgr=dotrho(i,3)
       lcor=1
       lpot=1
      call easypbe(up,agrup,delgrup,uplap,dn,agrdn,delgrdn,dnlap,
     &           agr,delgr,lcor,lpot,
     &           exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd,
     &           expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91,
     &           expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)

       if(up+dn.lt.rho_max*1.D-4) then     
       vxuppbe=0.d0
       vcuppbe=0.d0
       vxdnpbe=0.d0
       vcdnpbe=0.d0
       expbe=0.d0
       ecpbe=0.d0
       endif



      if(up+dn.gt.rho_max*1.D-2) then
      vtot(i,1)=vtot(i,1)+vion(i)+vxuppbe+vcuppbe
      vtot(i,2)=vtot(i,2)+vion(i)+vxdnpbe+vcdnpbe
      E_Hxc=E_Hxc+(expbe+ecpbe)*(rho(i,1)+rho(i,2)+rhocr(i))
      else      ! using LSDA potential, more stable
      call UxcCA2(rho(i,1)+rhocr(i)*0.5d0,rho(i,2)+rhocr(i)*0.5d0,
     &    vxc1,vxc2,uxc1,uxc2)
      vtot(i,1)=vtot(i,1)+vion(i)+vxc1
      vtot(i,2)=vtot(i,2)+vion(i)+vxc2
      E_Hxc=E_Hxc+uxc1*(rho(i,1)+rhocr(i)*0.5d0)+
     &         uxc2*(rho(i,2)+rhocr(i)*0.5d0)
      endif
      
      s=s+(vtot(i,1)+vtot(i,2))*0.5d0
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
      vtot(i,1)=vtot(i,1)-s
      vtot(i,2)=vtot(i,2)-s
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
