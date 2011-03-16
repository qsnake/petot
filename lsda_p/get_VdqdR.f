      subroutine get_VdqdR(vr_in_pL,occ,islda,nkpt,nmap_q,
     &   fatom)

******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'param.escan_real'
      include 'mpif.h'

      real*8 vr_in_pL(mr_nL,islda)
      real*8 occ(mst,nkpt,islda)
      real*8 fatom(3,matom),fatom_tmp(3,matom),fatmp(3,matom)
      complex*16,allocatable,dimension(:,:,:)  :: occ_beta,occ_betatmp

      integer iiatom(mtype),iatom(matom),icore(mtype),numref(matom),
     &  ityatom(matom)
      real*8 occ_t(mtype) 
      real*8 q_LMnm(32,32,49,mtype),q_LMnm0(32,32,12,mtype)
      real*8  q_LM(49,matom),q_LM0(12,matom)
      integer nmap_q(matom)

      integer isNLa(9,matom),ipsp_type(mtype),isps_all
      real*8  Dij0(32,32,mtype),Qij(32,32,mtype)

      complex*16 cc

      common /comNL2/occ_t,iiatom,icore,numref,ityatom
      common /comisNLa/isNLa,Dij0,Qij,ipsp_all,ipsp_type

      common /com_qLMnm/q_LMnm,q_LMnm0


      allocate(occ_beta(32,32,natom))
      
       mx_n=mx/nnodes+1

       fatom_tmp=0.d0

      do 200 iislda=1,islda

*****************************
      occ_beta=dcmplx(0.d0,0.d0)

      do 100 kpt=1,nkpt

      if(nkpt.gt.1.or.islda.gt.1) then
       if(icolor.eq.0) then
      call beta_psiIO(beta_psi,kpt,2,0,iislda)    ! let only one group to read beta_psi
       endif
       call mpi_bcast(beta_psi,nref_tot*mx_n,
     &  MPI_DOUBLE_COMPLEX,0,MPI_COMM_N,ierr)       ! bcast from icolor=0 to all the other groups
      endif     ! otherwise, beta_psi is already here

      do 100 im=1,mx_n
      m=(inode-1)*mx_n+im
       if(m.le.mx) then
       iref_tt=0
       do ia=1,natom

       do iref1=1,numref(ia)
       iref_t1=iref_tt+iref1

       do iref2=1,numref(ia)
       iref_t2=iref_tt+iref2

       occ_beta(iref1,iref2,ia)=occ_beta(iref1,iref2,ia)+
     &  beta_psi(iref_t1,im)*dconjg(beta_psi(iref_t2,im))
     &  *occ(m,kpt,iislda)
cccc dconjg should be on iref2, tested.

       enddo
       enddo
       iref_tt=iref_tt+numref(ia)
       enddo
       endif

100    continue

       allocate(occ_betatmp(32,32,natom))
       call mpi_allreduce(occ_beta,occ_betatmp,natom*32*32,
     &   MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_K,ierr)
       occ_beta = occ_betatmp
       deallocate(occ_betatmp)

**************************

       do ia=1,natom
       itype=ityatom(ia)
         if(ipsp_type(itype).eq.2) then
       do LM=2,49       ! LM is a combined index of (L,M), L=0,...,6
       cc=dcmplx(0.d0,0.d0)
       do iref2=1,numref(ia)
       do iref1=1,numref(ia)
       cc=cc+q_LMnm(iref1,iref2,LM,itype)*occ_beta(iref1,iref2,ia)
       enddo
       enddo
ccccc occ_beta(iref1,iref2,ia)=dconjg(occ_beta(iref2,iref1,ia))
       q_LM(LM,ia)=dreal(cc)      ! check, cc should be already real here
       enddo

       do LM0=1,12       ! LM0 is the index for LM=(0,0) and the cases of l(iref1)=l(iref2)
       cc=dcmplx(0.d0,0.d0)
       do iref2=1,numref(ia)
       do iref1=1,numref(ia)
       cc=cc+q_LMnm0(iref1,iref2,LM0,itype)*occ_beta(iref1,iref2,ia)
       enddo
       enddo
       q_LM0(LM0,ia)=dreal(cc)      ! check, cc should be already real here
       enddo


         endif
       enddo

******** d_rho= \sum_LM q_LM(LM,ia)* wmask_q(LM,r-ia)+\sum_LM0 q_LM0(LM0,ia)*wmask_q(LM0,r-ia)

       ico=0
       do ia=1,natom
       do i=1,nmap_q(ia)
       ico=ico+1
       do ixyz=1,3

       sum=0.d0
       do LM=2,49         ! special treatment for LM=1, by LM0
       sum=sum+wmask_dq(LM,ixyz,ico)*q_LM(LM,ia)
       enddo

       do LM0=1,12
       sum=sum+wmask_dq0(LM0,ixyz,ico)*q_LM0(LM0,ia)
       enddo

       fatom_tmp(ixyz,ia)=fatom_tmp(ixyz,ia)+
     &    vr_in_pL(indm_q(ico),iislda)*sum

ccccccccC       test
c       fatom_tmp(ixyz,ia)=fatom_tmp(ixyz,ia)+
c     &    vr_in_pL(indm_q(ico),iislda)*wmask_dq0(1,ixyz,ico)
c     &    vr_in_pL(indm_q(ico),iislda)*wmask_q0(1,ico)
c     &     wmask_q0(1,ico)*wmask_dq0(1,ixyz,ico)

       enddo

       enddo
       enddo

*****************************************8
200    continue

       call mpi_allreduce(fatom_tmp,fatmp,3*matom,
     & MPI_REAL8,MPI_SUM,MPI_COMM_K,ierr)
       fatom_tmp = fatmp

c       fatom=fatom+fatom_tmp*vol/(n1L*n2L*n3L)
       fatom=fatom-fatom_tmp*vol/(n1L*n2L*n3L)

       deallocate(occ_beta)

       return
       end




