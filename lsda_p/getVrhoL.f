      subroutine getVrhoL(AL,vion_pL,vionT_pL,xatom,
     &   ntype,iatom,rhocr_pL,totNel,ido_rho,islda)

******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************


ccccc if ido_rho.eq.1, generate rho_n from atom, 
ccccc if ido_rho.eq.0, do not touch rho_n, rho_n is inputted.

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'param.escan_real'
      include 'mpif.h'

      real*8 xatom(3,matom),AL(3,3)

      real*8 vion_pL(mr_nL)
      real*8 vionT_pL(mr_nL),rhocr_pL(mr_nL)

      real*8 qi2(mnq),vq(mnq,mtype),rhoq(mnq,mtype)
      real*8 vqT(mnq,mtype),rhocq(mnq,mtype)
      real*8 occ_t(mtype)
c
      integer iiatom(mtype),iatom(matom),icore(mtype),numref(matom),
     &  ityatom(matom)

      integer icoul
      real*8 xcoul(3)

      character*20 file_tmp

      complex*16 cc

      real*8, allocatable, dimension(:)  :: workr_nL
      real*8, allocatable, dimension(:)  :: workr_nL2,vion_pL2

***************************************************
****  xatom(1),xatom(2),xatom(3) are the coord in unit of AL(3,3)
****  supercell edges
***************************************************

      common /comNL2/occ_t,iiatom,icore,numref,ityatom
      common /comVrho/qi2,vq,rhoq,vqT,rhocq
      common /comcoul/icoul,xcoul

*******************************************************
**** generate the Kleiman-Bylander reference wavefunction
*******************************************************
      allocate(workr_nL(mr_nL))

      icc=0
      do itype=1,ntype
      icc=icc+icore(itype)
      enddo

      vins=1.d0/vol

      nh1L=(n1L+1)/2+1
      nr_nL=n1L*n2L*n3L/nnodes
      nrL=n1L*n2L*n3L

      ng2_nL=ngtotnod2L(inode)

      vion_pL = 0.0d0 
      vionT_pL = 0.0d0
      rhocr_pL = 0.0d0
      

      if(ido_rho.eq.1) rho_nL = 0.0d0

      do 10 i=1,ng2_nL
      y=0.d0

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

        ph=gkx2_nL(i)*x11+gky2_nL(i)*y11+gkz2_nL(i)*z11
        cc=cc+cdexp(dcmplx(0.d0,ph))
	endif

        enddo

      q=dsqrt(gkx2_nL(i)**2+gky2_nL(i)**2+gkz2_nL(i)**2)

      iq=1+q*(mnq-1.d0)/qi2(mnq)

      x=(q-qi2(iq))/(qi2(iq+1)-qi2(iq))    ! assuming equal distance grid

cccccc for even smoother treatment, we should use four points interpolation
cccccc But that might not be necessary if we use mqline=4000

      f1=1-x-0.5d0*x*(1-x)
      f2=x+x*(1-x)
      f3=-0.5d0*x*(1-x)       ! using quadratic interpolation

      if(icoul.eq.0) then
      y=vq(iq,itype)*qi2(iq)**2*f1+vq(iq+1,itype)*qi2(iq+1)**2*f2+
     &   vq(iq+2,itype)*qi2(iq+2)**2*f3
      endif

      y2=vqT(iq,itype)*qi2(iq)**2*f1+vqT(iq+1,itype)*qi2(iq+1)**2*f2
     &   +vqT(iq+2,itype)*qi2(iq+2)**2*f3

      if(q.lt.1.D-6) then
      y=0.d0
      y2=vqT(4,itype)
      else
      y=y/q**2
      y2=y2/q**2
      endif

      y1=rhoq(iq,itype)*f1+rhoq(iq+1,itype)*f2
     &  +rhoq(iq+2,itype)*f3        ! rhoq is the full charge, not just psi^2 

      i2 = i*2

      if(icoul.eq.0) then
      vion_pL(i2-1)=vion_pL(i2-1)+y*dreal(cc)*vins
      vion_pL(i2)=vion_pL(i2)+y*dimag(cc)*vins
      endif

      vionT_pL(i2-1)=vionT_pL(i2-1)+y2*dreal(cc)*vins
      vionT_pL(i2)=vionT_pL(i2)+y2*dimag(cc)*vins

      if(ido_rho.eq.1) then
      rho_nL(i2-1,1)=rho_nL(i2-1,1)+y1*dreal(cc)*vins
      rho_nL(i2,1)=rho_nL(i2,1)+y1*dimag(cc)*vins
      endif

      if(icore(itype).ne.0) then
      y1=rhocq(iq,itype)*f1+rhocq(iq+1,itype)*f2+
     &      rhocq(iq+2,itype)*f3
      rhocr_pL(i2-1)=rhocr_pL(i2-1)+y1*dreal(cc)*vins
      rhocr_pL(i2)=rhocr_pL(i2)+y1*dimag(cc)*vins
      endif

  9   continue
 10   continue

      scale=1.d0


      if(ido_rho.eq.1) then
      call d3fft_real2L(rho_nL(1,1),workr_nL,-1,0)
      if(islda.eq.1) rho_nL(:,1) = workr_nL     ! rho_nL is the full charge
      if(islda.eq.2) then
      rho_nL(:,1)=workr_nL*0.6d0
      rho_nL(:,2)=workr_nL*0.4d0
      endif
      endif


      if(icoul.eq.0) then
      call d3fft_real2L(vion_pL,workr_nL,-1,0)
      vion_pL = workr_nL
      endif

*************************************

      call d3fft_real2L(vionT_pL,workr_nL,-1,0)
      vionT_pL = workr_nL


*************************************
      if(icc.gt.0) then
      call d3fft_real2L(rhocr_pL,workr_nL,-1,0)


      s=0.d0
      s0=0.d0
      do i=1,nr_nL
      rhocr_pL(i) = workr_nL(i)
      s0=s0+rhocr_pL(i)
      if(rhocr_pL(i).lt.0.d0) rhocr_pL(i)=0.d0
      s=s+rhocr_pL(i)
      enddo
      call global_sumr(s)
      call global_sumr(s0)
      s=s*vol/nrL
      s0=s0*vol/nrL
      if(inode_tot.eq.1) then
      write(6,*) "The total core charge for core correction",s,s0
      write(22,*) "The total core charge for core correction",s,s0
      endif
      s=s0/s
      do i=1,nr_nL
      rhocr_pL(i)=rhocr_pL(i)*s
      enddo

      else
      rhocr_pL=0.d0
      endif
*************************************
      s=0.d0
      do iislda=1,islda
      do i=1,nr_nL
      s=s+rho_nL(i,iislda)
      enddo
      enddo

      call global_sumr(s)

      s=s*vol/nrL

      if(inode_tot.eq.1) then
      write(6,*) "input charge sum=",s
      write(6,*) "correct that to totNel", totNel
      endif

      s=totNel/s               ! rho_nL is the full charge, not just \sum psi^2(r)
      do iislda=1,islda
      do i=1,nr_nL
      rho_nL(i,iislda)=rho_nL(i,iislda)*s        ! do not use abs(),
      enddo
      enddo

      deallocate(workr_nL)

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(icoul.eq.0) return
cccc calculate vion_pL for icoul.ne.0
cccc for icoul.ne.0, vion_pL must be calculated based on n1L2,n2L2,n3L2,AL2 

      vins2=1.d0/vol2
      ng2_nL2=ngtotnod2L2(inode)

      allocate(workr_nL2(mr_nL2))
      allocate(vion_pL2(mr_nL2))
    
      vion_pL2=0.d0
      
      do 110 i=1,ng2_nL2

      do 19  itype=1,ntype
        cc=dcmplx(0.d0,0.d0)
        do ia=1,natom

	if(iatom(ia).eq.iiatom(itype)) then   ! even for icoul=11,12,13, shifting x1,y1,z1 is okay
        x1=xatom(1,ia)
         if(x1.gt.1.d0) x1=x1-1.d0
         if(x1.lt.0.d0) x1=x1+1.d0
         if(x1.gt.xcoul(1)) x1=x1-1.d0
        y1=xatom(2,ia)                      ! actually x2, not really y
         if(y1.gt.1.d0) y1=y1-1.d0
         if(y1.lt.0.d0) y1=y1+1.d0
         if(y1.gt.xcoul(2)) y1=y1-1.d0
        z1=xatom(3,ia)                      ! actually x3, not really z
         if(z1.gt.1.d0) z1=z1-1.d0
         if(z1.lt.0.d0) z1=z1+1.d0
         if(z1.gt.xcoul(3)) z1=z1-1.d0
         
         if(icoul.eq.1) then
         x1=x1/2
         y1=y1/2
         z1=z1/2
         endif
         if(icoul.eq.11) then
         x1=x1/2
         endif
         if(icoul.eq.12) then
         y1=y1/2
         endif
         if(icoul.eq.13) then
         z1=z1/2
         endif

        x11=AL2(1,1)*x1+AL2(1,2)*y1+AL2(1,3)*z1
        y11=AL2(2,1)*x1+AL2(2,2)*y1+AL2(2,3)*z1
        z11=AL2(3,1)*x1+AL2(3,2)*y1+AL2(3,3)*z1

        ph=gkx2_nL2(i)*x11+gky2_nL2(i)*y11+gkz2_nL2(i)*z11
        cc=cc+cdexp(dcmplx(0.d0,ph))
	endif

        enddo

      i2 = i*2

      vion_pL2(i2-1)=vion_pL2(i2-1)+
     &        vcoul_nL2(i2-1,itype)*dreal(cc)-
     &        vcoul_nL2(i2,itype)*dimag(cc)
      vion_pL2(i2)=vion_pL2(i2)+
     &        vcoul_nL2(i2-1,itype)*dimag(cc)+
     &        vcoul_nL2(i2,itype)*dreal(cc)

 19   continue


110   continue

      
      call d3fft_real2L2(vion_pL2,workr_nL2,-1,0)
      vion_pL2=workr_nL2

      deallocate(workr_nL2)

      call convert_2LtoL(vion_pL,vion_pL2,2)

      deallocate(vion_pL2)


      return
cccccccccccccccccccccccccccccccccccccccc

      end
      
