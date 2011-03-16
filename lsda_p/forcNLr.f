       subroutine forcNLr(workr_n,fatom,
     &   occ,kpt,nkpt,ixyz)
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
       include "mpif.h"

       real*8 fatom(3,matom)
       real*8 occ(mst,nkpt)

       complex*16 workr_n(mr_n)

       integer nmap(matom)
       integer isNLa(9,matom)

       real*8 occ_t(mtype)
       integer iiatom(mtype),icore(mtype),numref(matom)

       complex*16 s1(9,matom),s2(9,matom)
       complex*16 y
************************************************
*************************************************
       common /comisNLa/isNLa
       common /comnmap/nmap
       common /comNL2/occ_t,iiatom,icore,numref


       s1=0.d0
       s2=0.d0


       do 20 m=1,mx

       call d3fft_comp(ug_n(1,m),workr_n,-1,kpt)

       s1 = dcmplx(0.d0,0.d0)
       s2 = dcmplx(0.d0,0.d0)

       ico1=0
       ico2=0
       do ia=1,natom
       nref=numref(ia)
	do i=1,nmap(ia)
         ico1=ico1+1
         y=workr_n(indm(ico1))*cphase(ico1)
          do jj = 1,nref
          ico2=ico2+1
          s1(jj,ia)=s1(jj,ia)+wmask(ico2)*y
          s2(jj,ia)=s2(jj,ia)+wmaskX(ico2)*y
          enddo
        enddo
       enddo

       call mpi_barrier(MPI_COMM_WORLD,ierr)

       call mpi_allreduce(s1,s1,natom*9,
     $     MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_WORLD,ierr)

       call mpi_allreduce(s2,s2,natom*9,
     $     MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_WORLD,ierr)

       s1=s1*isNLa*vol/nr
       s2=s2*vol/nr

       do ia=1,natom
       nref=numref(ia)
       Fc=0
       do j=1,nref
       Fc=Fc+2*dreal(s1(j,ia)*dconjg(s2(j,ia)))
       enddo
       fatom(ixyz,ia)=fatom(ixyz,ia)+Fc*occ(m,kpt)
       enddo

20     continue
       return
       end
