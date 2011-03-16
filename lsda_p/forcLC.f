        subroutine forcLC(AL,rho_nt,vxc_nt,xatom,
     &               ntype,iatom,fatom)
cccccccccccccccccccccccccccccccccccccccccccccc
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc the input rho,dv,vxc, will be destroyed.
cccccccccccccccccccccccccccccccccccco
************************************************************
*** This subroutine calculates
*** \int [\partial Vion[R]/\partial R] rho(r) d^3r + 
*** \int [\partial rho_core[R]/\partial R] v_xc(r) d^3r+
*** \int [\partial rho_atom[R]/\partial R] dv(r) d^3r    ! currently not used
*** They are calculated in q space.
************************************************************
      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'param.escan_real'
      include 'mpif.h'


      real*8 rho_nt(mr_nL),vxc_nt(mr_nL)

      real*8 xatom(3,matom),fatom(3,matom)

      real*8 AL(3,3),ALt(3,3)
      real*8, allocatable, dimension(:) ::  workr_nL

      real*8 qi2(mnq),vq(mnq,mtype),rhoq(mnq,mtype)
      real*8 vqT(mnq,mtype),rhocq(mnq,mtype)
      real*8 occ_t(mtype)

      integer iiatom(mtype),iatom(matom),icore(mtype),numref(matom)
      integer ityatom(matom)

      complex*16 cc,cai

***************************************************
      common /comNL2/occ_t,iiatom,icore,numref,ityatom
      common /comVrho/qi2,vq,rhoq,vqT,rhocq
***************************************************
****  xatom(1),xatom(2),xatom(3) are the coord in unit of AL(3,3)
****  supercell edges
***************************************************

      cai=dcmplx(0.d0,1.d0)
*******************************************************

      allocate(workr_nL(mr_nL))

      workr_nL=rho_nt
      call d3fft_real2L(rho_nt,workr_nL,1,0)

      workr_nL=vxc_nt
      call d3fft_real2L(vxc_nt,workr_nL,1,0)

      deallocate(workr_nL)
*******************************************************
      ng2_nL=ngtotnod2L(inode)

      do 20  ia=1,natom
        x1=xatom(1,ia)
        y1=xatom(2,ia)
        z1=xatom(3,ia)
      
        x11=AL(1,1)*x1+AL(1,2)*y1+AL(1,3)*z1
        y11=AL(2,1)*x1+AL(2,2)*y1+AL(2,3)*z1
        z11=AL(3,1)*x1+AL(3,2)*y1+AL(3,3)*z1
          do j=1,ntype
          if(iatom(ia).eq.iiatom(j)) itype=j
	  enddo

      fac1=(mnq-1)/qi2(mnq)
      fx=0.d0
      fy=0.d0
      fz=0.d0

      do 10 i=1,ng2_nL
        ph=gkx2_nL(i)*x11+gky2_nL(i)*y11+gkz2_nL(i)*z11
        cc=cdexp(dcmplx(0.d0,ph))

      q=dsqrt(gkx2_nL(i)**2+gky2_nL(i)**2+gkz2_nL(i)**2)

      iq=1+q*fac1

      x=(q-qi2(iq))/(qi2(iq+1)-qi2(iq))

      f1=1-x-0.5d0*x*(1-x)
      f2=x+x*(1-x)
      f3=-0.5d0*x*(1-x)

      y=vq(iq,itype)*qi2(iq)**2*f1+vq(iq+1,itype)*qi2(iq+1)**2*f2+
     &    vq(iq+2,itype)*qi2(iq+2)**2*f3

      if(q.lt.1.D-6) then
      y=0.d0
      else
      y=y/q**2
      endif

      y1=rhoq(iq,itype)*f1+rhoq(iq+1,itype)*f2+
     &      rhoq(iq+2,itype)*f3

      y2=rhocq(iq,itype)*f1+rhocq(iq+1,itype)*f2+
     &      rhocq(iq+2,itype)*f3

      s1=y*dreal(cai*cc*dcmplx(rho_nt(2*i-1),-rho_nt(2*i)))
      fx=fx+s1*gkx2_nL(i)
      fy=fy+s1*gky2_nL(i)
      fz=fz+s1*gkz2_nL(i)

      s3=y2*dreal(cai*cc*dcmplx(vxc_nt(2*i-1),-vxc_nt(2*i)))
      fx=fx+s3*gkx2_nL(i)
      fy=fy+s3*gky2_nL(i)
      fz=fz+s3*gkz2_nL(i)
************************************************************
*** Not so sure this will improve the force convergence
************************************************************
c      s2=y1*dreal(cai*cc*dcmplx(dv_n(2*i-1),-dv_n(2*i))
c      fx=fx+s2*gkx2_nL(i)
c      fy=fy+s2*gky2_nL(i)
c      fz=fz+s2*gkz2_nL(i)
 10   continue

       call mpi_allreduce(fx,ftmp,1,MPI_REAL8,
     & MPI_SUM,MPI_COMM_K,ierr)
       fx = ftmp
       call mpi_allreduce(fy,ftmp,1,MPI_REAL8,
     & MPI_SUM,MPI_COMM_K,ierr)
       fy = ftmp
       call mpi_allreduce(fz,ftmp,1,MPI_REAL8,
     & MPI_SUM,MPI_COMM_K,ierr)
       fz = ftmp

      fatom(1,ia)=fatom(1,ia)+fx*2
      fatom(2,ia)=fatom(2,ia)+fy*2
      fatom(3,ia)=fatom(3,ia)+fz*2
 20   continue
****************************************************
*** The phase factor and normalization factor is correct, checked.
ccc It is unclear whether the second term dv will improved the force.
ccc need to be checked later. 
****************************************************

      return
      end
	 

