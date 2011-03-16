       subroutine Hpsi_comp(wg,wgh,ilocal,vr,workr_n,kpt)

******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
******************************************

********* wgh = H * wg

***** cphase is not installed yet

       use fft_data
       use load_data
       use data

       implicit double precision (a-h,o-z)
       include 'param.escan_real'
       include "mpif.h"


       real*8 vr(mr_n)
       complex*16 workr_n(mr_n)
       complex*16 workr2_n(mr_n)
       complex*16 wg(mg_nx),wgh(mg_nx)
       complex*16 wg_dum(mg_nx)

       integer nmap(matom)

       complex*16 sumdum(9,matom)
       complex*16 y,s
************************************************
****  extra memory for ilocal=2  !!!
*************************************************

c      real*8 xatom(3,matom)
       integer isNLa(9,matom)
       complex*16 cy

       real*8 occ_t(mtype)
       integer iiatom(mtype),icore(mtype),numref(matom)

*************************************************
       common /comisNLa/isNLa
       common /comnmap/nmap
       common /comNL2/occ_t,iiatom,icore,numref

       ng_n=ngtotnod(inode,kpt)


       call d3fft_comp(wg,workr_n,-1,kpt)

       if(ilocal.eq.2) then
       do i=1,nr_n
       workr2_n(i)=workr_n(i)
       enddo
       endif

       do i=1,nr_n
       workr_n(i)=dble(vr(i))*workr_n(i)
       enddo

       if(ilocal.eq.2) then

       sumdum = dcmplx(0.d0,0.d0)

       ico1=0
       ico2=0
       do ia=1,natom
        nref=numref(ia)
	do i=1,nmap(ia)
         ico1=ico1+1
         y=workr2_n(indm(ico1))*cphase(ico1)

          do jj = 1,nref
          ico2=ico2+1
           sumdum(jj,ia)=sumdum(jj,ia)+wmask(ico2)*y
          enddo
        enddo


       enddo

       call mpi_barrier(MPI_COMM_WORLD,ierr)

       call mpi_allreduce(sumdum,sumdum,natom*9,
     $     MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_WORLD,ierr)

       ico1=0
       ico2=0
       do ia=1,natom
        nref=numref(ia)
        do jj = 1,nref
          sumdum(jj,ia) = sumdum(jj,ia)*vol/nr*isNLa(jj,ia)
        enddo
        do i=1,nmap(ia)
         ico1=ico1+1
         s=dcmplx(0.d0,0.d0)
         do jj = 1,nref
         ico2=ico2+1
         s=s+sumdum(jj,ia)*wmask(ico2)
         enddo
         workr_n(indm(ico1))=workr_n(indm(ico1))+s*dconjg(cphase(ico1))
        enddo
       enddo

       endif !ilocal.eq.2

       call d3fft_comp(wgh,workr_n,1,kpt)

       do i=1,ng_n
       wgh(i)=wgh(i)+gkk_n(i,kpt)*wg(i)
       enddo

       if(ilocal.eq.3) then

       sumdum = dcmplx(0.d0,0.d0)

       do ia=1,natom

        nref=numref(ia)

	do i=1,ng_n
         cy=dconjg(wg(i))

         do jj = 1,nref
          sumdum(jj,ia) = sumdum(jj,ia)+wqmask(jj,i,ia)*cy
         enddo

	enddo

       enddo

       call mpi_allreduce(sumdum,sumdum,natom*9,
     $     MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_WORLD,ierr)
    
       do ia=1,natom
        nref=numref(ia)
        do jj = 1,nref
          sumdum(jj,ia) = dconjg(sumdum(jj,ia))*vol*isNLa(jj,ia)
        enddo

        do i=1,ng_n
         do jj = 1,nref
          wgh(i)=wgh(i)+sumdum(jj,ia)*wqmask(jj,i,ia)
         enddo
        enddo

       enddo
       endif


       return
       end
