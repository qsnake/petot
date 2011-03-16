      subroutine getpot(rhoL,vionL,rhocrL,vtot,v0,E_Hxc,
     & E_coul,E_ion,islda,igga)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************
*** input rhoL(nr_nL),vionL(nr_nL),rhocrL(nr_nL), output vtot(nr_n)
*** the nr_nL version of vtot(nr_n): viL(nr_nL) has been deleted
*******************************************

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)
      include 'param.escan_real'

      real*8 vtot(mr_n,islda)
      real*8 rhoL(mr_nL,islda),vionL(mr_nL),rhocrL(mr_nL)

      real*8, allocatable,dimension(:)   :: workr_nL
      real*8, allocatable,dimension(:,:)   :: viL 
     
      real*8, allocatable,dimension(:)   :: vi2L,vitmpL,drhoL,dotrhoL,
     &  ddrhoL.  
      real*8, allocatable,dimension(:,:)   :: vtmpL,drho2L,dotrho2L,
     &  ddrho2L
      real*8, allocatable,dimension(:,:,:)  :: vtmp2L

**************************************************
      allocate(workr_nL(mr_nL))
      allocate(viL(mr_nL,islda))


      ng2_nL=ngtotnod2L(inode)

      viL = 0.0d0

      if(islda.eq.1) then
      do i=1,nr_nL
      workr_nL(i)=rhoL(i,1)
      enddo
      else
      do i=1,nr_nL
      workr_nL(i)=rhoL(i,1)+rhoL(i,2)
      enddo
      endif

      call d3fft_real2L(viL,workr_nL,1,0)

      pi=4*datan(1.d0)

      do 100 i=1,ng2_nL
 
      if(inode.eq.iorg2L(1).and.i.eq.iorg2L(2)) then
      ig = iorg2L(2)
      viL(ig*2,1)=0.d0
      viL(ig*2-1,1)=0.d0
      else
      viL(i*2,1)=viL(i*2,1)*2*pi/gkk2_nL(i)
      viL(i*2-1,1)=viL(i*2-1,1)*2*pi/gkk2_nL(i)
      endif

100   continue

      call d3fft_real2L(viL,workr_nL,-1,0)

      viL(:,1)=workr_nL(:)

ccccc change below according to islda, iGGA
      E_Hxc=0.d0
      E_coul=0.d0
      E_ion=0.d0
ccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccc
ccccc case one, the old getpot2

      if(islda.eq.1.and.igga.eq.0) then
      do i=1,nr_nL
      E_Hxc=E_Hxc+0.5d0*viL(i,1)*rhoL(i,1)
      E_coul=E_coul+0.5d0*viL(i,1)*rhoL(i,1)
      E_ion=E_ion+vionL(i)*rhoL(i,1)
      viL(i,1)=viL(i,1)+vionL(i)+
     &         UxcCA(rhoL(i,1)+rhocrL(i),uxc2)
      E_Hxc=E_Hxc+uxc2*(rhoL(i,1)+rhocr(i))
      enddo

      call convert_SLvr(vtot(1,1),viL(1,1),-1)

      endif      ! (islda.eq.1.and.igga.eq.0)

ccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccc
ccccc case two, the old getpot3

      if(islda.eq.2.and.igga.eq.0) then
      viL(:,2)=viL(:,1)
      do i=1,nr_nL
      E_Hxc=E_Hxc+0.5d0*(viL(i,1)*rhoL(i,1)+
     &                   viL(i,2)*rhoL(i,2))
      E_coul=E_coul+0.5d0*(viL(i,1)*rhoL(i,1)+
     &                   viL(i,2)*rhoL(i,2))
      E_ion=E_ion+vionL(i)*(rhoL(i,1)+rhoL(i,2))
      call UxcCA2(rhoL(i,1)+rhocrL(i)*0.5d0,
     & rhoL(i,2)+rhocrL(i)*0.5d0,vxc1,vxc2,uxc1,uxc2)
      viL(i,1)=viL(i,1)+vionL(i)+vxc1
      viL(i,2)=viL(i,2)+vionL(i)+vxc2
      E_Hxc=E_Hxc+uxc1*(rhoL(i,1)+rhocrL(i)*0.5d0)+
     &         uxc2*(rhoL(i,2)+rhocrL(i)*0.5d0)
      enddo

      call convert_SLvr(vtot(1,1),viL(1,1),-1)
      call convert_SLvr(vtot(1,2),viL(1,2),-1)

      endif    ! (islda.eq.2.and.igga.eq.0)

ccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccc
ccccc case three, the old getpot4

      if(islda.eq.1.and.igga.eq.1) then
      allocate(vi2L(mr_nL))
      allocate(vitmpL(mr_nL))
      allocate(vtmpL(mr_nL,3))
      allocate(drhoL(mr_nL))
      allocate(dotrhoL(mr_nL))
      allocate(ddrhoL(mr_nL))
      drho=0.d0
      dotrho=0.d0
      ddrho=0.d0
*******************

      workr_nL(:)=rhoL(:,1)+rhocrL(:)

      call d3fft_real2L(vi2L,workr_nL,1,0)
      do ixyz=1,3
      call vitmp_eq_vigkxyz()
      call d3fft_real2L(vitmpL,workr_nL,-1,0)
      vtmpL(:,ixyz)=workr_nL(:)
      enddo

      drhoL(:)=dsqrt(vtmpL(:,1)**2+vtmpL(:,2)**2+
     &               vtmpL(:,3)**2)
*********************
      do i=1,ng2_nL
      vitmpL(i*2)=2*vi2L(i*2)*gkk2_nL(i)
      vitmpL(i*2-1)=2*vi2L(i*2-1)*gkk2_nL(i)
      enddo
      call d3fft_real2(vitmpL,workr_nL,-1,0)
      ddrhoL=-workr_nL
**************
      workr_nL=drhoL

      call d3fft_real2L(vi2L,workr_nL,1,0)

      do ixyz=1,3
      call vitmp_eq_vigkxyz()     ! multiply i*gkx,i*gky,i*gkz to vi
      call d3fft_real2L(vitmpL,workr_nL,-1,0)
      dotrhoL(:)=dotrhoL(:)+workr_nL(:)*vtmpL(:,ixyz)
      enddo
*********************************************************
**** finished the calculation of drho, ddrho and dotrho,
**** now, we can call GGA subroutine one r point at a time
*********************************************************
ccc      drhoL(mr_nL)        !  |grad(rho)|
ccc      dotrhoL(mr_nL)       !  [grad(rho)]\cdot[grad|grad(rho)|]
ccc      ddrhoL(mr_nL)        !  \nablda^2(rho)
*********************************************************
       do 1000 i=1,nr_nL
 
      E_Hxc=E_Hxc+0.5d0*viL(i,1)*rhoL(i,1)
      E_coul=E_coul+0.5d0*viL(i,1)*rhoL(i,1)
      E_ion=E_ion+vionL(i)*rhoL(i,1)


       up=(rhoL(i,1)+rhocrL(i))*0.5d0
       agrup=drhoL(i)*0.5d0
       delgrup=dotrhoL(i)*0.25d0
       uplap=ddrhoL(i)*0.5d0
       dn=(rhoL(i)+rhocrL(i))*0.5d0
       agrdn=drhoL(i)*0.5d0
       delgrdn=dotrhoL(i)*0.25d0
       dnlap=ddrhoL(i)*0.5d0
       agr=drhoL(i)
       delgr=dotrhoL(i)
       lcor=1
       lpot=1

      call easypbe(up,agrup,delgrup,uplap,dn,agrdn,delgrdn,dnlap,
     &           agr,delgr,lcor,lpot,
     &           exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd,
     &           expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91,
     &           expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)

ccccccccc dirty trick to avoid the convergency
      if(rhoL(i,1).gt.0.003) then
      viL(i,1)=viL(i,1)+vionL(i)+vxuppbe+vcuppbe
      E_Hxc=E_Hxc+(expbe+ecpbe)*(rhoL(i,1)+rhocrL(i))
      else        ! use LDA formula
      viL(i)=viL(i)+vionL(i)+UxcCA(rhoL(i,1)+rhocrL(i),uxc2)
      E_Hxc=E_Hxc+uxc2*(rhoL(i,1)+rhocrL(i))
      endif
1000  continue
*************************************
      call convert_SLvr(vtot(1,1),viL(1,1),-1)

      deallocate(vi2L)
      deallocate(vitmpL)
      deallocate(vtmpL)
      deallocate(drhoL)
      deallocate(dotrhoL)
      deallocate(ddrhoL)
      endif      ! (islda.eq.1.and.igga.eq.1)

ccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccc
ccccc case four, the old getpot5

      if(islda.eq.2.and.igga.eq.1) then

      allocate(drho2L(mr_nL,3))         !  |grad(rho)|
      allocate(dotrho2L(mr_nL,3))       !  [grad(rho)]\cdot[grad|grad(rho)|]
      allocate(ddrho2L(mr_nL,2))        !  \nablda^2(rho)
      allocate(vi2L(mr_nL))
      allocate(vitmpL(mr_nL))
      allocate(vtmp2L(mr_nL,3,2))
ccccccccc drho,dotrho,ddrho are all we need to call the one r point gga subroutine
      viL(:,2)=viL(:,1)

      drho2L=0.d0
      dotrho2L=0.d0
      ddrho2L=0.d0

******************************************
      do 200 iislda=1,2

      workr_nL(:)=rhoL(:,iislda)+rhocrL(:)*0.5d0

      call d3fft_real2L(vi2L,workr_nL,1,0)

      do ixyz=1,3
      call vitmp_eq_vigkxyz()     ! multiply i*gkx,i*gky,i*gkz to vi
      call d3fft_real2L(vitmpL,workr_nL,-1,0)
      vtmp2L(:,ixyz,iislda)=workr_nL(:)
      enddo

      drho2L(:,iislda)=dsqrt(vtmp2L(:,1,iislda)**2+
     &  vtmp2L(:,2,iislda)**2+vtmp2L(:,3,iislda)**2)

******************************************
      do i=1,ng2_nL
      vitmpL(i*2)=2*vi2L(i*2)*gkk2_nL(i)
      vitmpL(i*2-1)=2*vi2L(i*2-1)*gkk2_nL(i)
      enddo
      call d3fft_real2L(vitmpL,workr_nL,-1,0)
      ddrho2L(:,iislda)=-workr_nL(:)
*******************************************
      workr_nL(:)=drho2L(:,iislda)

      call d3fft_real2L(vi2L,workr_nL,1,0)

      do ixyz=1,3
      call vitmp_eq_vigkxyz()     ! multiply i*gkx,i*gky,i*gkz to vi2L
      call d3fft_real2L(vitmpL,workr_nL,-1,0)
      dotrho2L(:,iislda)=dotrho2L(:,iislda)+
     &            workr_nL(:)*vtmp2L(:,ixyz,iislda)
      enddo

200   continue
********************************************
      drho2L(:,3)=dsqrt((vtmp2L(:,1,1)+vtmp2L(:,1,2))**2+
     &                (vtmp2L(:,2,1)+vtmp2L(:,2,2))**2+
     &                (vtmp2L(:,3,1)+vtmp2L(:,3,2))**2)

      workr_nL(:)=drho2L(:,3)

      call d3fft_real2L(vi2L,workr_nL,1,0)

      do ixyz=1,3
      call vitmp_eq_vigkxyz()     ! multiply i*gkx,i*gky,i*gkz to vi
      call d3fft_real2L(vitmpL,workr_nL,-1,0)
      dotrho2L(:,3)=dotrho2L(:,iislda)+workr_nL(:)*
     &     (vtmp2L(:,ixyz,1)+vtmp2L(:,ixyz,2))
      enddo
********************************************************
**** finished the calculation of drho, ddrho and dotrho,
**** now, we can call GGA subroutine one r point at a time
*********************************************************
ccccc (mr_n,3): 1, spin up; 2, spin down; 3 spin up+spin down
ccc
ccc      drho(mr_n,3)        !  |grad(rho)|
ccc      dotrho(mr_n,3)       !  [grad(rho)]\cdot[grad|grad(rho)|]
ccc      ddrho(mr_n,2)        !  \nablda^2(rho)
*********************************************************

      do 1001 i=1,nr_nL
cccccccccc Hartree energy

      E_Hxc=E_Hxc+0.5d0*(viL(i,1)*rhoL(i,1)+
     &                   viL(i,2)*rhoL(i,2))
      E_coul=E_coul+0.5d0*(viL(i,1)*rhoL(i,1)+
     &                   viL(i,2)*rhoL(i,2))
      E_ion=E_ion+vionL(i)*(rhoL(i,1)+rhoL(i,2))

       up=rhoL(i,1)+rhocrL(i)*0.5d0
       agrup=drho2L(i,1)
       delgrup=dotrho2L(i,1)
       uplap=ddrho2L(i,1)
       dn=rhoL(i,2)+rhocrL(i)*0.5d0
       agrdn=drho2L(i,2)
       delgrdn=dotrho2L(i,2)
       dnlap=ddrho2L(i,2)
       agr=drho2L(i,3)
       delgr=dotrho2L(i,3)
       lcor=1
       lpot=1
      call easypbe(up,agrup,delgrup,uplap,dn,agrdn,delgrdn,dnlap,
     &           agr,delgr,lcor,lpot,
     &           exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd,
     &           expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91,
     &           expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)

      viL(i,1)=viL(i,1)+vionL(i)+vxuppbe+vcuppbe
      viL(i,2)=viL(i,2)+vionL(i)+vxdnpbe+vcdnpbe
      E_Hxc=E_Hxc+(expbe+ecpbe)*(rhoL(i,1)+rhoL(i,2)+rhocrL(i))

1001  continue
      call convert_SLvr(vtot(1,1),viL(1,1),-1)
      call convert_SLvr(vtot(1,2),viL(1,2),-1)

      deallocate(drho2L)         !  |grad(rho)|
      deallocate(dotrho2L)       !  [grad(rho)]\cdot[grad|grad(rho)|]
      deallocate(ddrho2L)        !  \nablda^2(rho)
      deallocate(vi2L)
      deallocate(vitmpL)
      deallocate(vtmp2L)

      endif    !(islda.eq.2.and.igga.eq.1)
ccccccccccccccccccccccccccccccccccccccccc


      call global_sumr(E_Hxc)
      call global_sumr(E_coul)
      call global_sumr(E_ion)

      E_Hxc=E_Hxc*vol/nrL
      E_coul=E_coul*vol/nrL
      E_ion=E_ion*vol/nrL

      s=0.d0
      do iislda=1,islda
      do i=1,nr_nL
      s=s+vtot(i,iislda)
      enddo
      enddo
      s=s/islda
      call global_sumr(s)
      s=s/nrL
      
      do iislda=1,islda
      do i=1,nr_nL
      vtot(i,iislda)=vtot(i,iislda)-s
      enddo
      enddo

      v0=s
ccccccccccccccccccccc
cccccc ave_vtot=0.d0 is necessary for mch_Broyden and mch_pulay,
cccccc otherwise it doesn't work
ccccccccccccccccccccc
      deallocate(workr_nL)
      deallocate(viL)

      return
      contains
cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine vitmp_eq_vigkxyz()

      implicit double precision (a-h,o-z)

        if(ixyz.eq.1) then
        do i=1,ng2_nL
        vitmpL(i*2)=-vi2L(i*2-1)*gkx2_nL(i)
        vitmpL(i*2-1)=vi2L(i*2)*gkx2_nL(i)
        enddo
        endif
        if(ixyz.eq.2) then
        do i=1,ng2_nL
        vitmpL(i*2)=-vi2L(i*2-1)*gky2_nL(i)
        vitmpL(i*2-1)=vi2L(i*2)*gky2_nL(i)
        enddo
        endif
        if(ixyz.eq.3) then
        do i=1,ng2_nL
        vitmpL(i*2)=-vi2L(i*2-1)*gkz2_nL(i)
        vitmpL(i*2-1)=vi2L(i*2)*gkz2_nL(i)
        enddo
        endif

      end subroutine vitmp_eq_vigkxyz

      end

      

