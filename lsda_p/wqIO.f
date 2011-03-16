      subroutine wqIO(nkpt,kpt,iflag)

******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************

****************************************
ccccc iflag=1, write, iflag=2, read
****************************************
******************************************
      use fft_data
      use load_data
      use data
      implicit double precision (a-h,o-z)

      include 'mpif.h'
      include 'param.escan_real'

      integer status(MPI_STATUS_SIZE)

      complex*16,allocatable,dimension(:)  :: wqtemp

       character*8 fname
*************************************************
       if(nkpt.eq.1) return

*************************************************
       if(iflag.eq.2) then     ! read the wqmask
       if(inode==nnodes) then

       fname="wqiofile"
       kpt1=mod(kpt,10)
       kpt2=mod((kpt-kpt1)/10,10)
       kpt3=mod((kpt-kpt1-kpt2*10)/100,10)
       kpt4=mod((kpt-kpt1-kpt2*10-kpt3*100)/1000,10)

       open(10,file=
     & fname//char(kpt4+48)//char(kpt3+48)//
     &  char(kpt2+48)//char(kpt1+48),
     &       form="unformatted")

         allocate(wqtemp(mg_nx))

            do i=0,nnodes-2
               do iref=1,nref_tot
                  read(10) wqtemp    ! assuming the mg_nx is the same as before
                  wqmask(:,iref)=wqtemp
               end do

          call mpi_send(wqmask,mg_nx*nref_tot,
     & MPI_DOUBLE_COMPLEX,i,102,MPI_COMM_K,ierr)
            end do
c     Now read in the data for node nnodes
            do iref=1,nref_tot
               read(10) wqtemp
               wqmask(:,iref)=wqtemp
            end do
            deallocate(wqtemp)
            close(10)

         else   ! for other nodes

         call mpi_recv(wqmask,mg_nx*nref_tot,
     &  MPI_DOUBLE_COMPLEX,nnodes-1,102,MPI_COMM_K,status,ierr)

         endif    ! for the nodes
        endif     ! for the iflag, read
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       if(iflag.eq.1) then     ! write the wavefunction

       if(inode==1) then

       fname="wqiofile"
       kpt1=mod(kpt,10)
       kpt2=mod((kpt-kpt1)/10,10)
       kpt3=mod((kpt-kpt1-kpt2*10)/100,10)
       kpt4=mod((kpt-kpt1-kpt2*10-kpt3*100)/1000,10)

       open(10,file=
     & fname//char(kpt4+48)//char(kpt3+48)//
     &  char(kpt2+48)//char(kpt1+48),
     &       form="unformatted")

            allocate(wqtemp(mg_nx))

            do iref=1,nref_tot
               wqtemp=wqmask(:,iref)
               write(10) wqtemp
            end do
        endif    ! inode==1

           call mpi_barrier(MPI_COMM_K,ierr)

        do i=1,nnodes-1
         do iref=1,nref_tot

           call mpi_barrier(MPI_COMM_K,ierr)
           if(inode==i+1) then
           call mpi_send(wqmask(1,iref),mg_nx,
     &  MPI_DOUBLE_COMPLEX,0,102,MPI_COMM_K,ierr)
           endif

           if(inode==1) then
           call mpi_recv(wqtemp,mg_nx,
     & MPI_DOUBLE_COMPLEX,i,102,MPI_COMM_K,status,ierr)
           write(10) wqtemp   
           endif
           call mpi_barrier(MPI_COMM_K,ierr)


          enddo  ! do iref
         enddo    ! do inode

         if(inode==1) then
         close(10)
         deallocate(wqtemp)
         endif

        endif     ! for the iflag, write
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        return
        end

