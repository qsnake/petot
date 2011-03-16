      subroutine getpot4(rho,vion,rhocr,vtot,v0,E_Hxc,workr_n,
     & E_coul,E_ion)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
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

      real*8 workr_n(mr_n)
      real*8 vtot(mr_n)
      real*8 rho(mr_n),vion(mr_n),rhocr(mr_n)


      real*8, allocatable, dimension(:) :: vi,vitmp

      real*8, allocatable, dimension(:) :: drho,ddrho,dotrho

      real*8, allocatable, dimension(:,:) :: vtmp

**************************************************
ccccc generate the potential vtot from rho(i)
ccccc (mr_n,3): 1, spin up; 2, spin down; 3 spin up+spin down

      allocate(drho(mr_n))         !  |grad(rho)|
      allocate(dotrho(mr_n))       !  [grad(rho)]\cdot[grad|grad(rho)|]
      allocate(ddrho(mr_n))        !  \nablda^2(rho)

ccccccccc drho,dotrho,ddrho are all we need to call the one r point gga subroutines

      drho=0.d0
      dotrho=0.d0
      ddrho=0.d0

      allocate(vi(mr_n))
      allocate(vitmp(mr_n))
      allocate(vtmp(mr_n,3))

****************************************************

      ng2_n=ngtotnod2(inode)

      vi = 0.0d0

      do i=1,nr_n
      workr_n(i) = rho(i)     ! total charge for Coulomb potential
      enddo

      call d3fft_real2(vi,workr_n,1,0)


      pi=4*datan(1.d0)

      do 100 i=1,ng2_n
 
      if(inode.eq.iorg2(1).and.i.eq.iorg2(2)) then
      ig = iorg2(2)
      vi(ig*2)=0.d0
      vi(ig*2-1)=0.d0
      else
      vi(i*2)=vi(i*2)*2*pi/gkk2_n(i)
      vi(i*2-1)=vi(i*2-1)*2*pi/gkk2_n(i)
      endif

100   continue

      call d3fft_real2(vi,workr_n,-1,0)

      
      do i=1,nr_n
      vtot(i) = workr_n(i)
      enddo

*************************************************************
****** finished the Coulomb potential
*************************************************************
****** now,  calculate drho, ddrho, dotrho
*************************************************************

      do i=1,nr_n
      workr_n(i)=rho(i)+rhocr(i)
      enddo
     
      call d3fft_real2(vi,workr_n,1,0)

      do 201 ixyz=1,3

      call vitmp_eq_vigkxyz()     ! multiply i*gkx,i*gky,i*gkz to vi

      call d3fft_real2(vitmp,workr_n,-1,0)

      do i=1,nr_n
      vtmp(i,ixyz)=workr_n(i)
      enddo

201   continue



      do i=1,nr_n
      drho(i)=dsqrt(vtmp(i,1)**2+
     &  vtmp(i,2)**2+vtmp(i,3)**2)
      enddo
      
******************************************
      do i=1,ng2_n
      vitmp(i*2)=2*vi(i*2)*gkk2_n(i)
      vitmp(i*2-1)=2*vi(i*2-1)*gkk2_n(i)
      enddo
      call d3fft_real2(vitmp,workr_n,-1,0)
      do i=1,nr_n
      ddrho(i)=-workr_n(i)
      enddo
*******************************************
      do i=1,nr_n
      workr_n(i)=drho(i)
      enddo
     
      call d3fft_real2(vi,workr_n,1,0)

      do 202 ixyz=1,3

      call vitmp_eq_vigkxyz()     ! multiply i*gkx,i*gky,i*gkz to vi

      call d3fft_real2(vitmp,workr_n,-1,0)

      do i=1,nr_n
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

      do 1000 i=1,nr_n

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

      vtot(i)=vtot(i)+vion(i)+vxuppbe+vcuppbe
      E_Hxc=E_Hxc+(expbe+ecpbe)*(rho(i)+rhocr(i))

      s=s+vtot(i)
1000  continue

****************************************************
      deallocate(drho)      
      deallocate(dotrho)    
      deallocate(ddrho)     
****************************************************
****************************************************

      do i=1,nr_n
      workr_n(i) = vtot(i)
      enddo

***** only keep the G components within Ecut2
      call d3fft_real2(vtot,workr_n,1,0)
      call d3fft_real2(vtot,workr_n,-1,0)
***** we can save one fft by rearrange the order.

      do i=1,nr_n
      vtot(i) = workr_n(i)
      enddo

      call global_sumr(s)
      call global_sumr(E_Hxc)
      call global_sumr(E_coul)
      call global_sumr(E_ion)

      E_Hxc=E_Hxc*vol/nr
      E_coul=E_coul*vol/nr
      E_ion=E_ion*vol/nr

      s=s/nr
      do i=1,nr_n
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
