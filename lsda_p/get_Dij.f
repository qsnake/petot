      subroutine get_Dij(Dij_t,vtot_nL_t,nmap_q,af)

******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************
cccccccccccccccccccccccccccccccccccc
ccc  possibly this can be improved by dividing the labour (ia=1,natom) among diff icolor. 
ccc  Currently, all the icolor groups are doing the same thing for ia=1,natom



      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'param.escan_real'
      include 'mpif.h'

      real*8 Dij_t(32,32,natom),Dij_t0(32,32,natom)
      real*8 vtot_nL_t(mr_nL)
ccccccccccccccccccc
      integer iiatom(mtype),iatom(matom),icore(mtype),numref(matom),
     &  ityatom(matom)
      real*8 occ_t(mtype)
      real*8 q_LMnm(32,32,49,mtype),q_LMnm0(32,32,12,mtype)
      integer nmap_q(matom)

      integer isNLa(9,matom),ipsp_type(mtype),isps_all
      real*8  Dij0(32,32,mtype),Qij(32,32,mtype)
      real*8, allocatable, dimension (:,:)  :: sum_wv,sum_wv0
      real*8, allocatable, dimension (:,:)  :: sum_wvtmp,sum_wv0tmp

      complex*16 cc
ccccccccccccccccccc

      common /comNL2/occ_t,iiatom,icore,numref,ityatom
      common /comisNLa/isNLa,Dij0,Qij,ipsp_all,ipsp_type

      common /com_qLMnm/q_LMnm,q_LMnm0


      allocate(sum_wv(49,natom))
      allocate(sum_wv0(12,natom))
      allocate(sum_wvtmp(49,natom))
      allocate(sum_wv0tmp(12,natom))
      Dij_t0=Dij_t
      pi=4*datan(1.d0)

c         write(6,*) "inode,wmask",inode,nmap_q(1),wmask_q(1,1),
c     &   wmask_q(1,2),wmask_q(1,nmap_q(1)/2+1)
c         write(6,*) "inode,wmask0",inode,wmask_q0(1,1),
c     &   wmask_q0(1,2),wmask_q0(1,nmap_q(1)/2+1)
c         write(6,*) "inode,vtot",indm_q(1),indm_q(2),
c     &          vtot_nL_t(indm_q(1)),vtot_nL_t(indm_q(2))

       sum_wv=0.d0
       sum_wv0=0.d0
       ico=0
       do ia=1,natom
       do i=1,nmap_q(ia)
       ico=ico+1
       if(indm_q(ico).gt.mr_nL) write(6,*) inode, indm_q(ico),ico
       vtmp=vtot_nL_t(indm_q(ico))
       do LM=1,49
       sum_wv(LM,ia)=sum_wv(LM,ia)+vtmp*wmask_q(LM,ico)
       enddo
       do LM0=1,12
       sum_wv0(LM0,ia)=sum_wv0(LM0,ia)+vtmp*wmask_q0(LM0,ico)
       enddo
       enddo
       enddo


       call mpi_allreduce(sum_wv,sum_wvtmp,natom*49,
     &   MPI_REAL8,MPI_SUM,MPI_COMM_K,ierr)
       sum_wv = sum_wvtmp
       deallocate(sum_wvtmp)

       call mpi_allreduce(sum_wv0,sum_wv0tmp,natom*12,
     &   MPI_REAL8,MPI_SUM,MPI_COMM_K,ierr)
       sum_wv0 = sum_wv0tmp
       deallocate(sum_wv0tmp)

       sum_wv=sum_wv*vol/(n1L*n2L*n3L)
       sum_wv0=sum_wv0*vol/(n1L*n2L*n3L)

ccccccccccccccccccccccccccccccccccc
c       if(inode.eq.1) then
c       write(6,*) "sum_wv",sum_wv(:,1)
c       write(6,*) "**********************************"
c       write(6,*) "sum_wv0",sum_wv0(:,1)
c       write(6,*) "**********************************"
c       endif

       Dij_t=0.d0


       do ia=1,natom
       itype=ityatom(ia)
         if(ipsp_type(itype).eq.2) then
       do LM=2,49       ! LM is a combined index of (L,M), L=0,...,6
       do iref2=1,numref(ia)
       do iref1=1,numref(ia)
       Dij_t(iref1,iref2,ia)=Dij_t(iref1,iref2,ia)+
     &   q_LMnm(iref1,iref2,LM,itype)*sum_wv(LM,ia)
       enddo
       enddo
       enddo

       do LM0=1,12       ! LM0 is a index for (L,M)=(0,0), and (l,m)[iref1]=(l,m)[iref2]
       do iref2=1,numref(ia)
       do iref1=1,numref(ia)
       Dij_t(iref1,iref2,ia)=Dij_t(iref1,iref2,ia)+
     &   q_LMnm0(iref1,iref2,LM0,itype)*sum_wv0(LM0,ia)
       enddo
       enddo
       enddo

         endif
       enddo


       sum1=0.d0
       sum2=0.d0
       diff=0.d0
       do ia=1,natom
       itype=ityatom(ia)
       do iref2=1,numref(ia)
       do iref1=1,numref(ia)
       Dij_t(iref1,iref2,ia)=Dij_t(iref1,iref2,ia)+
     &     Dij0(iref1,iref2,itype)
       sum1=sum1+dabs(Dij_t0(iref1,iref2,ia))
       sum2=sum2+dabs(Dij_t(iref1,iref2,ia))
       diff=diff+dabs(Dij_t(iref1,iref2,ia)-Dij_t0(iref1,iref2,ia))
       enddo
       enddo
       enddo

       Dij_t=Dij_t*af+(1-af)*Dij_t0

       

       call mpi_allreduce(sum1,sumtmp,1,
     &   MPI_REAL8,MPI_SUM,MPI_COMM_K,ierr)
       sum1 = sumtmp

       call mpi_allreduce(sum2,sumtmp,1,
     &   MPI_REAL8,MPI_SUM,MPI_COMM_K,ierr)
       sum2 = sumtmp

       call mpi_allreduce(diff,sumtmp,1,
     &   MPI_REAL8,MPI_SUM,MPI_COMM_K,ierr)
       diff = sumtmp

       if(inode.eq.1) then
        write(6,*) "sum1,sum2,diff",sum1,sum2,diff
       endif

       deallocate(sum_wv)
       deallocate(sum_wv0)

       return
       end

