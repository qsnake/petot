      subroutine getpot2(rho,vion,rhocr,vtot,v0,E_Hxc,workr_n,
     & E_coul,E_ion)
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

      real*8 workr_n(mr_n)
      real*8 vtot(mr_n)
      real*8 rho(mr_n),vion(mr_n),rhocr(mr_n)

      real*8 vi(mr_n)

**************************************************
ccccc generate the potential vtot from rho(i)

      ng2_n=ngtotnod2(inode)

      vtot = 0.0d0
      vi = 0.0d0

      s = 0.0d0

      do i=1,nr_n
      vtot(i)=rho(i)
      s = s + vtot(i)
      enddo

      call global_sumr(s)

      do i=1,nr_n
      workr_n(i) = vtot(i)
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

      
      s=0.d0
      E_Hxc=0.d0
      E_coul=0.d0
      E_ion=0.d0

      do i=1,nr_n
      E_Hxc=E_Hxc+0.5d0*vtot(i)*rho(i)
      E_coul=E_coul+0.5d0*vtot(i)*rho(i)
      E_ion=E_ion+vion(i)*rho(i)
      vtot(i)=vtot(i)+vion(i)+UxcCA(rho(i)+rhocr(i),uxc2)
      E_Hxc=E_Hxc+uxc2*(rho(i)+rhocr(i))
      s=s+vtot(i)
      enddo
 
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
      end
      

