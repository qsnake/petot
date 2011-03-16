       subroutine forcNLr(workr_n,fatom,
     &   occ,kpt,nkpt,iislda,islda,E_st)
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
       include "mpif.h"

       real*8 fatom(3,matom)

       complex*16 workr_n(mr_n)

       integer nmap(matom)
       integer isNLa(9,matom)
       real*8 Dij0(32,32,mtype),Qij(32,32,mtype)
       real*8 y_tmp1(10000),y_tmp2(10000)
       real*8 sumy1(32),sumy2(32)


       real*8 occ_t(mtype)
       integer iiatom(mtype),icore(mtype),numref(matom),
     &   ityatom(matom),ipsp_type(mtype)
       real*8 occ(mst,nkpt,islda),E_st(mst,nkpt,islda)

       complex*16 s1(32,matom),s2(32,matom,3),stmp(32,matom)
       complex*16 y

       complex*16, allocatable, dimension (:,:,:,:) :: ss,ss2
************************************************
*************************************************
       common /comisNLa/isNLa,Dij0,Qij,ipsp_all,ipsp_type
       common /comnmap/nmap
       common /comNL2/occ_t,iiatom,icore,numref,ityatom

cccccccccccccccccccccccccccccccccccccc
       mref=0
       do ia=1,natom
       if(numref(ia).gt.mref) mref=numref(ia)
       enddo

       if(ipsp_all.eq.2) then
       allocate(ss(mref,mref,natom,3))
       allocate(ss2(mref,mref,natom,3))
       ss=dcmplx(0.d0,0.d0)
       ss2=dcmplx(0.d0,0.d0)
       endif

       do 40 m=1,mx

       call d3fft_comp(ug_n(1,m),workr_n,-1,kpt)

       s1 = dcmplx(0.d0,0.d0)
       s2 = dcmplx(0.d0,0.d0)

       ico1=0
       ico2=0
       do ia=1,natom
       nref=numref(ia)

cccccccccccccccccccccccccccccc
c	do i=1,nmap(ia)
c         ico1=ico1+1
c         y=workr_n(indm(ico1))*cphase(ico1)  
c          do jj = 1,nref
c          ico2=ico2+1
c          s1(jj,ia)=s1(jj,ia)+wmask(ico2)*y
c          s2(jj,ia)=s2(jj,ia)+wmaskX(ico2,ixyz)*y
c          enddo
c        enddo
ccccccccccccccccccccccccccccccccccccccccc
        do i=1,nmap(ia)
        ico1=ico1+1
         y=workr_n(indm(ico1))*cphase(ico1)
         y_tmp1(i)=dreal(y)
         y_tmp2(i)=dimag(y)
        enddo
        
ccccccccccccccccccccccccccccc
cccccccccc wmask is real*8 not complex*16 !
ccccccccccc maybe combine this to a matrix*matrix multiplication ?
        
       call dgemv('N',nref,nmap(ia),1.d0,wmask(ico2+1),nref,
     & y_tmp1,1,0.d0,sumy1,1)
       call dgemv('N',nref,nmap(ia),1.d0,wmask(ico2+1),nref,
     & y_tmp2,1,0.d0,sumy2,1)

         if(nmap(ia).gt.0) then
        do jj=1,nref
        s1(jj,ia)=dcmplx(sumy1(jj),sumy2(jj))
        enddo
          endif

         do ixyz=1,3
       call dgemv('N',nref,nmap(ia),1.d0,wmaskX(ico2+1,ixyz),nref,
     & y_tmp1,1,0.d0,sumy1,1)
       call dgemv('N',nref,nmap(ia),1.d0,wmaskX(ico2+1,ixyz),nref,
     & y_tmp2,1,0.d0,sumy2,1)

         if(nmap(ia).gt.0) then
        do jj=1,nref
        s2(jj,ia,ixyz)=dcmplx(sumy1(jj),sumy2(jj))
        enddo
          endif
         enddo   ! ixyz

         ico2=ico2+nmap(ia)*nref
 
       enddo    ! ia=1,natom
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




       call mpi_barrier(MPI_COMM_K,ierr)

       call mpi_allreduce(s1,stmp,natom*32,
     $     MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_K,ierr)
       s1 = stmp

       do ixyz=1,3
       call mpi_allreduce(s2(1,1,ixyz),stmp,natom*32,
     $     MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_K,ierr)
       s2(:,:,ixyz) = stmp(:,:)
       enddo

       s1=s1*vol/nr
       s2=s2*vol/nr

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       do 21 ixyz=1,3
       do 20 ia=1,natom
       iitype=ityatom(ia)
       nref=numref(ia)

       if(ipsp_type(iitype).eq.2) then
         do j2=1,nref
         do j1=1,nref
         ss(j1,j2,ia,ixyz)=ss(j1,j2,ia,ixyz)+dconjg(s1(j2,ia))*
     &                        s2(j1,ia,ixyz)*occ(m,kpt,iislda)
         ss2(j1,j2,ia,ixyz)=ss2(j1,j2,ia,ixyz)+dconjg(s1(j2,ia))*
     &    s2(j1,ia,ixyz)*occ(m,kpt,iislda)*E_st(m,kpt,iislda)
         enddo
         enddo

        else

         fx=0.d0
         do j1=1,nref
         fx=fx+Dij(j1,j1,ia,iislda)*dconjg(s1(j1,ia))*
     &                         s2(j1,ia,ixyz)*occ(m,kpt,iislda)
         enddo
         fatom(ixyz,ia)=fatom(ixyz,ia)+2*fx
        endif

20      continue
21      continue



40       continue

cccccccccccccccccccccccccccc

      

       do 51 ixyz=1,3
       do 50 ia=1,natom
       iitype=ityatom(ia)
       nref=numref(ia)
         if(ipsp_type(iitype).eq.2) then
         fx=0.d0
         do j1=1,nref
         do j2=1,nref
         fx=fx+Dij(j1,j2,ia,iislda)*ss(j1,j2,ia,ixyz)
     &        -Qij(j1,j2,iitype)*ss2(j1,j2,ia,ixyz)
         enddo
         enddo
         fatom(ixyz,ia)=fatom(ixyz,ia)+2*fx
         endif
50     continue
51     continue

cccccccccccccccccccccccccccccccccccccc
       if(ipsp_all.eq.2) then
       deallocate(ss)
       deallocate(ss2)
       endif

       return
       end
