      subroutine getrho_only(AL,xatom,
     &   ntype,iatom,totNel,rhonew_nL)

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

      real*8 xatom(3,matom)

      real*8 AL(3,3),ALt(3,3)
      
      real*8 rhonew_nL(mr_nL)

      real*8 qi2(mnq),vq(mnq,mtype),rhoq(mnq,mtype)
      real*8 vqT(mnq,mtype),rhocq(mnq,mtype)
      real*8 occ_t(mtype)
c
      integer iiatom(mtype),iatom(matom),icore(mtype),numref(matom),
     &   ityatom(matom)

      real*8, allocatable, dimension(:)  :: workr_nL

      complex*16 cc

***************************************************
****  xatom(1),xatom(2),xatom(3) are the coord in unit of AL(3,3)
****  supercell edges
***************************************************

      common /comNL2/occ_t,iiatom,icore,numref,ityatom
      common /comVrho/qi2,vq,rhoq,vqT,rhocq

*******************************************************
**** generate the Kleiman-Bylander reference wavefunction
*******************************************************

      vins=1.d0/vol

      nh1L=(n1L+1)/2+1

      ng2_nL=ngtotnod2L(inode)

      allocate(workr_nL(mr_nL))


      rhonew_nL = 0.0d0

      do 10 i=1,ng2_nL

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

      y1=rhoq(iq,itype)*f1+rhoq(iq+1,itype)*f2
     &  +rhoq(iq+2,itype)*f3     ! rhoq is the full charge, not just psi^2

      i2 = i*2

      rhonew_nL(i2-1)=rhonew_nL(i2-1)+y1*dreal(cc)*vins
      rhonew_nL(i2)=rhonew_nL(i2)+y1*dimag(cc)*vins

  9   continue
 10   continue

      scale=1.d0

      call d3fft_real2L(rhonew_nL,workr_nL,-1,0)
      rhonew_nL=workr_nL

      s=0.d0
      do i=1,nr_nL
      s=s+rhonew_nL(i)
      enddo

      call global_sumr(s)

      s=s*vol/nrL

      if(inode.eq.1) then
      write(6,*) "test, inside getrho_only,totNel,sum=",totNel,s
      endif

      deallocate(workr_nL)

      return
cccccccccccccccccccccccccccccccccccccccc

      end
      
