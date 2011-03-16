      subroutine getVrho(AL,vion_p,vionT_p,xatom,
     &   ntype,iatom,rhocr_p,totNel,ido_rho,
     &   workr_n,islda)

******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
******************************************


ccccc if ido_rho.eq.1, generate rho_n from atom, 
ccccc if ido_rho.eq.0, do not touch rho_n, rho_n is inputted.

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'param.escan_real'
      include 'mpif.h'

      real*8 xatom(3,matom)

      real*8 AL(3,3),ALt(3,3)
      real*8 workr_n(mr_n)    ! use only half of the original working space  
ccc  but workr_n was set up by shpalloc, they still have the same
ccc  address, to be used in shmem_iput

      real*8 vion_p(mr_n)
      real*8 vionT_p(mr_n),rhocr_p(mr_n)

      real*8 qi2(mnq),vq(mnq,mtype),rhoq(mnq,mtype)
      real*8 vqT(mnq,mtype),rhocq(mnq,mtype)
      real*8 occ_t(mtype)
c
      integer iiatom(mtype),iatom(matom),icore(mtype),numref(matom)

      complex*16 cc

***************************************************
****  xatom(1),xatom(2),xatom(3) are the coord in unit of AL(3,3)
****  supercell edges
***************************************************

      common /comNL2/occ_t,iiatom,icore,numref
      common /comVrho/qi2,vq,rhoq,vqT,rhocq

*******************************************************
**** generate the Kleiman-Bylander reference wavefunction
*******************************************************
      icc=0
      do itype=1,ntype
      icc=icc+icore(itype)
      enddo

      vins=1.d0/vol

      nh1=(n1+1)/2+1

      ng2_n=ngtotnod2(inode)

      vion_p = 0.0d0 
      vionT_p = 0.0d0
      rhocr_p = 0.0d0

      if(ido_rho.eq.1) rho_n = 0.0d0

      do 10 i=1,ng2_n

      do 9  itype=1,ntype
        cc=dcmplx(0.d0,0.d0)
        do ia=1,natom

	if(iatom(ia).eq.iiatom(itype)) then
        x1=xatom(1,ia)
        y1=xatom(2,ia)
        z1=xatom(3,ia)
      
        x11=AL(1,1)*x1+AL(1,2)*y1+AL(1,3)*z1
        y11=AL(2,1)*x1+AL(2,2)*y1+AL(2,3)*z1
        z11=AL(3,1)*x1+AL(3,2)*y1+AL(3,3)*z1

        ph=gkx2_n(i)*x11+gky2_n(i)*y11+gkz2_n(i)*z11
        cc=cc+cdexp(dcmplx(0.d0,ph))
	endif

        enddo

      q=dsqrt(gkx2_n(i)**2+gky2_n(i)**2+gkz2_n(i)**2)

      iq=1+q*(mnq-1.d0)/qi2(mnq)

      x=(q-qi2(iq))/(qi2(iq+1)-qi2(iq))    ! assuming equal distance grid

cccccc for even smoother treatment, we should use four points interpolation
cccccc But that might not be necessary if we use mqline=4000

      f1=1-x-0.5d0*x*(1-x)
      f2=x+x*(1-x)
      f3=-0.5d0*x*(1-x)       ! using quadratic interpolation

      y=vq(iq,itype)*f1+vq(iq+1,itype)*f2+vq(iq+2,itype)*f3
      if(q.lt.1.D-6) y=0.d0

      y2=vqT(iq,itype)*f1+vqT(iq+1,itype)*f2+vqT(iq+2,itype)*f3
      y2=y2-y
      if(q.lt.1.D-6) y2=vqT(3,itype)-vq(3,itype)

      y1=rhoq(iq,itype)*f1+rhoq(iq+1,itype)*f2
     &  +rhoq(iq+2,itype)*f3

      i2 = i*2

      vion_p(i2-1)=vion_p(i2-1)+y*dreal(cc)*vins
      vion_p(i2)=vion_p(i2)+y*dimag(cc)*vins

      vionT_p(i2-1)=vionT_p(i2-1)+y2*dreal(cc)*vins
      vionT_p(i2)=vionT_p(i2)+y2*dimag(cc)*vins

      if(ido_rho.eq.1) then
      rho_n(i2-1,1)=rho_n(i2-1,1)+y1*dreal(cc)*vins
      rho_n(i2,1)=rho_n(i2,1)+y1*dimag(cc)*vins
      endif

      if(icore(itype).ne.0) then
      y1=rhocq(iq,itype)*f1+rhocq(iq+1,itype)*f2+
     &      rhocq(iq+2,itype)*f3
      rhocr_p(i2-1)=rhocr_p(i2-1)+y1*dreal(cc)*vins
      rhocr_p(i2)=rhocr_p(i2)+y1*dimag(cc)*vins
      endif

  9   continue
 10   continue

      scale=1.d0

      if(ido_rho.eq.1) then
      call d3fft_real2(rho_n(1,1),workr_n,-1,0)
      if(islda.eq.1) rho_n(:,1) = workr_n
      if(islda.eq.2) then
      rho_n(:,1)=workr_n*0.6d0
      rho_n(:,2)=workr_n*0.4d0
      endif
      endif



      call d3fft_real2(vion_p,workr_n,-1,0)
      vion_p = workr_n


      call d3fft_real2(vionT_p,workr_n,-1,0)
      vionT_p = workr_n


*************************************
      if(icc.gt.0) then
      call d3fft_real2(rhocr_p,workr_n,-1,0)
      s=0.d0
      s0=0.d0
      do i=1,nr_n
      rhocr_p(i) = workr_n(i)
      s0=s0+rhocr_p(i)
      if(rhocr_p(i).lt.0.d0) rhocr_p(i)=0.d0
      s=s+rhocr_p(i)
      enddo
      call global_sumr(s)
      call global_sumr(s0)
      s=s*vol/nr
      s0=s0*vol/nr
      if(inode.eq.1) then
      write(6,*) "The total core charge for core correction",s,s0
      write(22,*) "The total core charge for core correction",s,s0
      endif
      s=s0/s
      do i=1,nr_n
      rhocr_p(i)=rhocr_p(i)*s
      enddo

      else
      rhocr_p=0.d0
      endif
*************************************
      

      s=0.d0
      do iislda=1,islda
      do i=1,nr_n
      s=s+dabs(rho_n(i,iislda))
      enddo
      enddo

      call global_sumr(s)

      s=s*vol/nr

      if(inode.eq.1) then
      write(6,*) "input charge sum=",s
      write(6,*) "correct that to totNel", totNel
      endif

      s=totNel/s
      do iislda=1,islda
      do i=1,nr_n
      rho_n(i,iislda)=dabs(rho_n(i,iislda))*s
      enddo
      enddo

      return
cccccccccccccccccccccccccccccccccccccccc

      end
      
