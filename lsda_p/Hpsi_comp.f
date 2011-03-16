       subroutine Hpsi_comp(wg,wgh,ilocal,vr,workr_n,kpt,
     &  isij,swg,sumdum,iislda)

******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************

********* wgh = H * wg
********* swg= wg+Qij*wg     ! execute only when isij=1 and ipsp_all=2,  used for orth. 
*********                    ! otherwise, swg is a dummy

***** cphase is not installed yet

       use fft_data
       use load_data
       use data

       implicit double precision (a-h,o-z)
       include 'param.escan_real'
       include "mpif.h"


       real*8 vr(mr_n)
       complex*16 workr_n(mr_n)
       complex*16 wg(mg_nx),wgh(mg_nx),swg(1)
       complex*16 wg_dum(mg_nx)
       complex*16, allocatable, dimension(:) :: workr2_n

       integer nmap(matom)

       complex*16 sumdum(32,matom),sumdum2(32),sumdum3(32)
       complex*16 sumdumtmp(32,matom)
       complex*16 y,s,cc,s2
************************************************
****  extra memory for ilocal=2  !!!
*************************************************

c      real*8 xatom(3,matom)
       integer isNLa(9,matom),ityatom(matom),ipsp_type(mtype)
       complex*16 cy

       real*8 occ_t(mtype)
       integer iiatom(mtype),icore(mtype),numref(matom)
       real*8 Dij0(32,32,mtype),Qij(32,32,mtype)

*************************************************
       common /comisNLa/isNLa,Dij0,Qij,ipsp_all,ipsp_type
       common /comnmap/nmap
       common /comNL2/occ_t,iiatom,icore,numref,ityatom


       ng_n=ngtotnod(inode,kpt)


       call d3fft_comp(wg,workr_n,-1,kpt)


       if(ilocal.eq.2) then
       call nonlocal_realsp()
       endif
   
       if(ilocal.ne.2) then
       do i=1,nr_n
       workr_n(i)=dble(vr(i))*workr_n(i)
       enddo
       endif

       call d3fft_comp(wgh,workr_n,1,kpt)

       do i=1,ng_n
       wgh(i)=wgh(i)+gkk_n(i,kpt)*wg(i)
       enddo


       if(ilocal.eq.3) then     
       call nonlocal_qsp()
       endif    

       return

       contains
cccccccccccccccccccccccccccccccc

       subroutine nonlocal_realsp()
       implicit double precision (a-h,o-z)
       complex*16 y,cc1,cc0
       real*8 y_tmp1(10000),y_tmp2(10000)
       real*8 y_tmp31(10000),y_tmp32(10000)
       real*8 sumy1(32),sumy2(32)
       real*8 sumy31(32),sumy32(32)

       cc1=dcmplx(1.d0,0.d0)
       cc0=dcmplx(0.d0,0.d0)

       sumdum = dcmplx(0.d0,0.d0)

       ico1=0
       ico2=0
       do ia=1,natom
        nref=numref(ia)

cccccccccccccccccccccccccc
c	do i=1,nmap(ia)
c         ico1=ico1+1
c         y=dconjg(workr_n(indm(ico1))*cphase(ico1))
c          do jj = 1,nref
c          ico2=ico2+1
c           sumdum(jj,ia)=sumdum(jj,ia)+wmask(ico2)*y
c          enddo
c        enddo
cccccccccccccccccccccccccccccc
        if(nmap(ia).gt.10000) then
        write(6,*) "nmap.gt.10000, stop", nmap(ia),ia
        call  mpi_abort(MPI_COMM_WORLD,ierr)
        endif

	do i=1,nmap(ia)
        ico1=ico1+1
         y=dconjg(workr_n(indm(ico1))*cphase(ico1))
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

          if(nmap(ia).gt.0) then  ! for nmap.eq.0, sumy1 not right, sumdum=0, initialized
        do jj=1,nref
        sumdum(jj,ia)=dcmplx(sumy1(jj),sumy2(jj))
        enddo
        ico2=ico2+nmap(ia)*nref
          endif
cccccccccccccccccccccccccccccccccccccccccccccc
       enddo   ! ia
ccccccccccccccc

       call mpi_barrier(MPI_COMM_K,ierr)

cccccccccccc most often, sumdum is much smaller than natom*32 !
       call mpi_allreduce(sumdum,sumdumtmp,natom*32,
     $     MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_K,ierr)
       sumdum = sumdumtmp
ccccccccccccccccccccccccccccccccccccccccccccccccccc

       sumdum=sumdum*vol/(n1*n2*n3)
cccccccccccccccccccccccccccccccccccccccccc

       if(isij.eq.1.and.ipsp_all.eq.2) then
       allocate(workr2_n(mr_n))
       workr2_n=workr_n
       endif

       do i=1,nr_n
       workr_n(i)=dble(vr(i))*workr_n(i)
       enddo


       ico1=0
       ico2=0
       do ia=1,natom
        nref=numref(ia)
        iitype=ityatom(ia)


       do jj1=1,nref
        if(ipsp_type(iitype).eq.2) then
        cc=dcmplx(0.d0,0.d0)
        do jj2=1,nref
        cc=cc+Dij(jj2,jj1,ia,iislda)*sumdum(jj2,ia)
        enddo
        sumdum2(jj1)=dconjg(cc)
        else
        sumdum2(jj1)=Dij(jj1,jj1,ia,iislda)*
     &                    dconjg(sumdum(jj1,ia))
        endif

        if(isij.eq.1.and.ipsp_type(iitype).eq.2) then
        cc=cmplx(0.d0,0.d0)
        do jj2=1,nref
        cc=cc+Qij(jj2,jj1,iitype)*sumdum(jj2,ia)
        enddo
        sumdum3(jj1)=dconjg(cc)
        endif
       enddo   ! jj1=1,nref

cccccccccccc
c       if(isij.eq.1.and.ipsp_type(iitype).eq.2) then
c         do i=1,nmap(ia)
c         ico1=ico1+1
c         s=dcmplx(0.d0,0.d0)
c         s2=dcmplx(0.d0,0.d0)
c         do jj = 1,nref
c         ico2=ico2+1
c         s=s+sumdum2(jj)*wmask(ico2)
c         s2=s2+sumdum3(jj)*wmask(ico2)
c         enddo

c         workr_n(indm(ico1))=workr_n(indm(ico1))+
c     &         s*dconjg(cphase(ico1))
c         workr2_n(indm(ico1))=workr2_n(indm(ico1))+
c     &         s2*dconjg(cphase(ico1))
c         enddo   
c       else
cccccccccccccccccccccccccccccccccccc
c         do i=1,nmap(ia)
c         ico1=ico1+1
c         s=dcmplx(0.d0,0.d0)
c         do jj = 1,nref
c         ico2=ico2+1
c         s=s+sumdum2(jj)*wmask(ico2)
c         enddo
c         workr_n(indm(ico1))=workr_n(indm(ico1))+
c     &            s*dconjg(cphase(ico1))
c         enddo    ! i=1,nmap
ccccccccccccccccccccccccccccccccccccccccccccccccc
c       endif

       if(isij.eq.1.and.ipsp_type(iitype).eq.2) then

       do jj=1,nref
       sumy1(jj)=dreal(sumdum2(jj))
       sumy2(jj)=dimag(sumdum2(jj))
       sumy31(jj)=dreal(sumdum3(jj))
       sumy32(jj)=dimag(sumdum3(jj))
       enddo
cccccccccccccc maybe we should really use matrix*matrix multiplication
       call dgemv('T',nref,nmap(ia),1.d0,wmask(ico2+1),nref,
     & sumy1,1,0.d0,y_tmp1,1)
       call dgemv('T',nref,nmap(ia),1.d0,wmask(ico2+1),nref,
     & sumy2,1,0.d0,y_tmp2,1)
       call dgemv('T',nref,nmap(ia),1.d0,wmask(ico2+1),nref,
     & sumy31,1,0.d0,y_tmp31,1)
       call dgemv('T',nref,nmap(ia),1.d0,wmask(ico2+1),nref,
     & sumy32,1,0.d0,y_tmp32,1)
        ico2=ico2+nmap(ia)*nref
       do i=1,nmap(ia)
       ico1=ico1+1
       workr_n(indm(ico1))=workr_n(indm(ico1))+
     &       dcmplx(y_tmp1(i),y_tmp2(i))*dconjg(cphase(ico1))
       workr2_n(indm(ico1))=workr2_n(indm(ico1))+
     &       dcmplx(y_tmp31(i),y_tmp32(i))*dconjg(cphase(ico1))
        enddo   

       else

       do jj=1,nref
       sumy1(jj)=dreal(sumdum2(jj))
       sumy2(jj)=dimag(sumdum2(jj))
       enddo
 
       call dgemv('T',nref,nmap(ia),1.d0,wmask(ico2+1),nref,
     & sumy1,1,0.d0,y_tmp1,1)
       call dgemv('T',nref,nmap(ia),1.d0,wmask(ico2+1),nref,
     & sumy2,1,0.d0,y_tmp2,1)
        ico2=ico2+nmap(ia)*nref
       do i=1,nmap(ia)
       ico1=ico1+1
       workr_n(indm(ico1))=workr_n(indm(ico1))+
     &       dcmplx(y_tmp1(i),y_tmp2(i))*dconjg(cphase(ico1))
        enddo   
       endif

cccccccccccccccccccccccccccc

       enddo      ! ia=1,natom


       if(isij.eq.1.and.ipsp_all.eq.2) then
       call d3fft_comp(swg,workr2_n,1,kpt)
       deallocate(workr2_n)
       endif

    
       return
       end subroutine nonlocal_realsp

cccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccc


       subroutine nonlocal_qsp()
       implicit double precision (a-h,o-z)
       complex*16 cc1,cc0,sumy(32)
       cc1=dcmplx(1.d0,0.d0)
       cc0=dcmplx(0.d0,0.d0)

       sumdum = dcmplx(0.d0,0.d0)

       do ia=1,natom

        nref=numref(ia)
        iref_start2=iref_start(ia)

cccccccccccccccccccccccccccccccccccccccc
c        do jj = 1,nref
c        iref=iref_start2+jj
c          do i=1,ng_n
c          sumdum(jj,ia) = sumdum(jj,ia)+wqmask(i,iref)*    
c     &                        dconjg(wg(i))  ! fetch two data. wqmask,wg, instead of one (if we switch (i,iref)). 
c
c          enddo
c        enddo
cccccccccccccccccccccccccccccccccccccccc
       call zgemv('C',ng_n,nref,cc1,wqmask(1,iref_start2+1),
     & mg_nx,wg,1,cc0,sumy,1)

       if(ng_n.gt.0) then
       do jj=1,nref 
       sumdum(jj,ia)=dconjg(sumy(jj))
       enddo
       endif

       enddo      ! ia
cccccccccccccccccccccccccccccccccc

       call mpi_allreduce(sumdum,sumdumtmp,natom*32,
     $     MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_K,ierr)
       sumdum = sumdumtmp

       sumdum=vol*sumdum  

cccccccccccccccccc

       if(isij.eq.1.and.ipsp_all.eq.2) then
       do i=1,ng_n
       swg(i)=wg(i)
       enddo
       endif
    
       do ia=1,natom
        nref=numref(ia)
        iref_start2=iref_start(ia)
        iitype=ityatom(ia)


       do jj1=1,nref
        if(ipsp_type(iitype).eq.2) then
        cc=0.d0
        do jj2=1,nref
        cc=cc+Dij(jj2,jj1,ia,iislda)*sumdum(jj2,ia)
        enddo
        sumdum2(jj1)=dconjg(cc)
        else
        sumdum2(jj1)=Dij(jj1,jj1,ia,iislda)*
     &                    dconjg(sumdum(jj1,ia))
        endif

        if(isij.eq.1.and.ipsp_type(iitype).eq.2) then
        cc=0.d0
        do jj2=1,nref
        cc=cc+Qij(jj2,jj1,iitype)*sumdum(jj2,ia)
        enddo
        sumdum3(jj1)=dconjg(cc)
        endif
       enddo

ccccccccccccccccccccccccccccccccccccc
c        do jj = 1,nref
c        iref=iref_start2+jj
c          if(isij.eq.1.and.ipsp_type(iitype).eq.2) then
c           do i=1,ng_n
c           wgh(i)=wgh(i)+sumdum2(jj)*wqmask(i,iref)
c           swg(i)=swg(i)+sumdum3(jj)*wqmask(i,iref)
c           enddo
c          else
c           do i=1,ng_n
c           wgh(i)=wgh(i)+sumdum2(jj)*wqmask(i,iref)
c           enddo
c          endif
c        enddo
ccccccccccccccccccccccccccccccccccccc

       if(isij.eq.1.and.ipsp_type(iitype).eq.2) then
       call zgemv('N',ng_n,nref,cc1,wqmask(1,iref_start2+1),
     & mg_nx,sumdum2,1,cc1,wgh,1)
       call zgemv('N',ng_n,nref,cc1,wqmask(1,iref_start2+1),
     & mg_nx,sumdum3,1,cc1,swg,1)
       else
       call zgemv('N',ng_n,nref,cc1,wqmask(1,iref_start2+1),
     & mg_nx,sumdum2,1,cc1,wgh,1)
       endif

       enddo     !  do ia

        return

       end subroutine nonlocal_qsp
ccccccccccccccccccccccccc
       end


  

