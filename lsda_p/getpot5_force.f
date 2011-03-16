      subroutine getpot5_force(rho,vtot,rhocr,workr_n)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
******************************************

***********************************************************
**** This program calculate only the vxc to be used for core correction 
**** force calculation.
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

      real*8 workr_n(mr_n)
      real*8 vtot(mr_n)
      real*8 rho(mr_n,2),rhocr(mr_n)


      real*8, allocatable, dimension(:) :: vi,vitmp

      real*8, allocatable, dimension(:,:) :: drho,ddrho,dotrho

      real*8, allocatable, dimension(:,:,:) :: vtmp

**************************************************
ccccc generate the potential vtot from rho(i)
ccccc (mr_n,3): 1, spin up; 2, spin down; 3 spin up+spin down

      allocate(drho(mr_n,3))         !  |grad(rho)|
      allocate(dotrho(mr_n,3))       !  [grad(rho)]\cdot[grad|grad(rho)|]
      allocate(ddrho(mr_n,2))        !  \nablda^2(rho)

ccccccccc drho,dotrho,ddrho are all we need to call the one r point gga subroutines

      drho=0.d0
      dotrho=0.d0
      ddrho=0.d0

      allocate(vi(mr_n))
      allocate(vitmp(mr_n))
      allocate(vtmp(mr_n,3,2))

****************************************************

      ng2_n=ngtotnod2(inode)

      pi=4*datan(1.d0)

*************************************************************
****** calculate drho, ddrho, dotrho
*************************************************************

      do 200 iislda=1,2

      do i=1,nr_n
      workr_n(i)=rho(i,iislda)+rhocr(i)*0.5d0
      enddo
     
      call d3fft_real2(vi,workr_n,1,0)

      do 201 ixyz=1,3

      call vitmp_eq_vigkxyz()     ! multiply i*gkx,i*gky,i*gkz to vi

      call d3fft_real2(vitmp,workr_n,-1,0)

      do i=1,nr_n
      vtmp(i,ixyz,iislda)=workr_n(i)
      enddo

201   continue



      do i=1,nr_n
      drho(i,iislda)=dsqrt(vtmp(i,1,iislda)**2+vtmp(i,2,iislda)**2+
     &     vtmp(i,3,iislda)**2)
      enddo
      
******************************************
      do i=1,ng2_n
      vitmp(i*2)=2*vi(i*2)*gkk2_n(i)
      vitmp(i*2-1)=2*vi(i*2-1)*gkk2_n(i)
      enddo
      call d3fft_real2(vitmp,workr_n,-1,0)
      do i=1,nr_n
      ddrho(i,iislda)=-workr_n(i)
      enddo
*******************************************
      do i=1,nr_n
      workr_n(i)=drho(i,iislda)
      enddo
     
      call d3fft_real2(vi,workr_n,1,0)

      do 202 ixyz=1,3

      call vitmp_eq_vigkxyz()     ! multiply i*gkx,i*gky,i*gkz to vi

      call d3fft_real2(vitmp,workr_n,-1,0)

      do i=1,nr_n
      dotrho(i,iislda)=dotrho(i,iislda)+workr_n(i)*vtmp(i,ixyz,iislda)
      enddo
202   continue

200   continue
********************************************

      do i=1,nr_n
      drho(i,3)=dsqrt((vtmp(i,1,1)+vtmp(i,1,2))**2+
     &                (vtmp(i,2,1)+vtmp(i,2,2))**2+
     &                (vtmp(i,3,1)+vtmp(i,3,2))**2)
      enddo
      

      do i=1,nr_n
      workr_n(i)=drho(i,3)
      enddo
     
      call d3fft_real2(vi,workr_n,1,0)

      do 203 ixyz=1,3

      call vitmp_eq_vigkxyz()     ! multiply i*gkx,i*gky,i*gkz to vi

      call d3fft_real2(vitmp,workr_n,-1,0)

      do i=1,nr_n
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
      
      do 1000 i=1,nr_n

cccccccccccccccccccccccccccccccccccccccccccc
cccccccccc Hartree energy
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

      vtot(i)=(vxuppbe+vcuppbe+vxdnpbe+vcdnpbe)*0.5d0
1000  continue

****************************************************
      deallocate(drho)      
      deallocate(dotrho)    
      deallocate(ddrho)     
****************************************************
      return

      contains
      
      subroutine vitmp_eq_vigkxyz()

      implicit double precision (a-h,o-z)

        if(ixyz.eq.1) then
        do i=1,ng2_n
        vitmp(i*2)=-vi(i*2-1)*gkx2_n(i)
        vitmp(i*2-1)=vi(i*2)*gkx2_n(i)
        enddo
        endif
        if(ixyz.eq.2) then
        do i=1,ng2_n
        vitmp(i*2)=-vi(i*2-1)*gky2_n(i)
        vitmp(i*2-1)=vi(i*2)*gky2_n(i)
        enddo
        endif
        if(ixyz.eq.3) then
        do i=1,ng2_n
        vitmp(i*2)=-vi(i*2-1)*gkz2_n(i)
        vitmp(i*2-1)=vi(i*2)*gkz2_n(i)
        enddo
        endif

      end subroutine vitmp_eq_vigkxyz

      end
