      subroutine dens_out(AL,workr_n,kpt_dens,
     & ispin_dens,iw_dens,fdens_out)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************

****************************************
***** this subroutine output the charge density.
******************************************

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'mpif.h'
      include 'param.escan_real'
***********************************************
       complex*16 workr_n(mg_nx)
       real*8,allocatable,dimension(:)  ::  temp_rho
       real*8 AL(3,3)
       character*20 fdens_out
       integer status(MPI_STATUS_SIZE)
       integer kpt_dens(2),ispin_dens(2),iw_dens(2)

**********************************************
       common /com123b/m1,m2,m3,ngb,nghb,ntype,rrcut,msb
      do i=1,nr_n
      rho_n(i,1)=0.d0
      enddo

      do 100 ikpt=kpt_dens(1),kpt_dens(2)        ! this is only called by icolor=0, not by all groups
      call gen_G_comp(ikpt,0)
      call fftprep_comp(n1,n2,n3)

      do 80 iislda=ispin_dens(1),ispin_dens(2)

      call ugIO(ug_n,ikpt,2,0,iislda)

      do 80 m=iw_dens(1),iw_dens(2)

      call d3fft_comp(ug_n(1,m),workr_n,-1,ikpt)

      do i=1,nr_n
      rho_n(i,1)=rho_n(i,1)+cdabs(workr_n(i))**2
      enddo

80    continue
100   continue


      call mpi_barrier(MPI_COMM_K,ierr)

      allocate(temp_rho(mr_n))

       if(inode.eq.1) then
       open(11,file=fdens_out,form="unformatted")
       rewind(11)
       write(11) n1,n2,n3,nnodes
       write(11) AL
       write(11) (rho_n(i,1),i=1,nr_n)
       endif

200    format(10(E14.8,1x))
201    format(3(i5,1x))
202    format(3(f12.6,1x))

       do i=1,nnodes-1
      call mpi_barrier(MPI_COMM_K,ierr)
       if(inode==i+1) then
        do iloop=1,mr_n
        temp_rho(iloop)=rho_n(iloop,1)
        enddo
        call  mpi_send(temp_rho,mr_n,MPI_REAL8,0,
     &   100,MPI_COMM_K,ierr)
       endif
       if(inode.eq.1) then
        call mpi_recv(temp_rho,mr_n,MPI_REAL8,i,
     &   100,MPI_COMM_K,status,ierr)

       do iloop=1,mr_n
       rho_n(iloop,1)=temp_rho(iloop)
       end do
       write(11) (rho_n(j,1),j=1,nr_n)
       endif
      call mpi_barrier(MPI_COMM_K,ierr)

       enddo

       deallocate(temp_rho)

       if(inode.eq.1)  close(11)
*********************************************

      return
      end

