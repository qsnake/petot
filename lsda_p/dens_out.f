      subroutine dens_out(AL,workr_n,m1_out,m2_out,
     &   fdens_out,kpt_dens)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
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

**********************************************
       common /com123b/m1,m2,m3,ngb,nghb,ntype,rrcut,msb

      do i=1,nr_n
      rho_n(i,1)=0.d0
      enddo

      do m=m1_out,m2_out

      call d3fft_comp(ug_n(1,m),workr_n,-1,kpt_dens)

      do i=1,nr_n
      rho_n(i,1)=rho_n(i,1)+cdabs(workr_n(i))**2
      enddo

      enddo


      call mpi_barrier(MPI_COMM_WORLD,ierr)

      allocate(temp_rho(mr_n))

       if(inode.eq.1) then
       open(11,file=fdens_out)
       rewind(11)
       write(11,201) n1,n2,n3
       write(11,202) AL(1,1),AL(2,1),AL(3,1)
       write(11,202) AL(1,2),AL(2,2),AL(3,2)
       write(11,202) AL(1,3),AL(2,3),AL(3,3)
       write(11,*) nnodes
       write(11,200) (rho_n(i,1),i=1,nr_n)
       endif

200    format(10(E10.4,1x))
201    format(3(i5,1x))
202    format(3(f12.6,1x))

       do i=1,nnodes-1
      call mpi_barrier(MPI_COMM_WORLD,ierr)
       if(inode==i+1) then
        do iloop=1,mr_n
        temp_rho(iloop)=rho_n(iloop,1)
        enddo
        call  mpi_send(temp_rho,mr_n,MPI_REAL8,0,
     &   100,MPI_COMM_WORLD,ierr)
       endif
       if(inode.eq.1) then
        call mpi_recv(temp_rho,mr_n,MPI_REAL8,i,
     &   100,MPI_COMM_WORLD,status,ierr)

       do iloop=1,mr_n
       rho_n(iloop,1)=temp_rho(iloop)
       end do
       write(11,200) (rho_n(j,1),j=1,nr_n)
       endif
      call mpi_barrier(MPI_COMM_WORLD,ierr)

       enddo

       deallocate(temp_rho)

       if(inode.eq.1)  close(11)
*********************************************

      return
      end

