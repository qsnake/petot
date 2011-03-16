      subroutine ugIO(ug_n_tmp,kpt,iflag,istep,iislda)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
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

       complex*16,allocatable,dimension(:) :: ugtemp
       complex*16 ug_n_tmp(mg_nx,mx)
       character*8 fname
*************************************************

*************************************************
       if(iflag.eq.2) then     ! read the wavefunction
       if(inode==nnodes) then

       fname="ugiofile"
       kpt1=mod(kpt,10)
       kpt2=mod((kpt-kpt1)/10,10)
       kpt3=mod((kpt-kpt1-kpt2*10)/100,10)
       kpt4=mod((kpt-kpt1-kpt2*10-kpt3*100)/1000,10)

       open(10,file=
     & fname//char(iislda+48)//char(istep+48)//char(kpt4+48)//
     &  char(kpt3+48)//char(kpt2+48)//char(kpt1+48),
     &       form="unformatted")

            allocate(ugtemp(mg_nx))

            do i=0,nnodes-2
               do iwavefun=1,mx
                  read(10)ugtemp    ! assuming the mg_nx is the same as before
                  ug_n_tmp(:,iwavefun)=ugtemp
               end do


          call mpi_send(ug_n_tmp,mg_nx*mx,
     & MPI_DOUBLE_COMPLEX,i,102,MPI_COMM_WORLD,ierr)
            end do
c     Now read in the data for node nnodes
            do iwavefun=1,mx
               read(10)ugtemp
               ug_n_tmp(:,iwavefun)=ugtemp
            end do
            deallocate(ugtemp)
            close(10)

         else   ! for other nodes

         call mpi_recv(ug_n_tmp,mg_nx*mx,
     & MPI_DOUBLE_COMPLEX,nnodes-1,102,MPI_COMM_WORLD,status,ierr)

         endif    ! for the nodes
        endif     ! for the iflag, read
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       if(iflag.eq.1) then     ! write the wavefunction

       if(inode==1) then
       fname="ugiofile"
       kpt1=mod(kpt,10)
       kpt2=mod((kpt-kpt1)/10,10)
       kpt3=mod((kpt-kpt1-kpt2*10)/100,10)
       kpt4=mod((kpt-kpt1-kpt2*10-kpt3*100)/1000,10)

       open(10,file=
     & fname//char(iislda+48)//char(istep+48)//char(kpt4+48)//
     &  char(kpt3+48)//char(kpt2+48)//char(kpt1+48),
     &       form="unformatted")

            allocate(ugtemp(mg_nx))

            do iwavefun=1,mx
               ugtemp=ug_n_tmp(:,iwavefun)
               write(10) ugtemp
            end do
        endif    ! inode==1

        call mpi_barrier(MPI_COMM_WORLD,ierr)

        do i=1,nnodes-1
         do iwavefun=1,mx

        call mpi_barrier(MPI_COMM_WORLD,ierr)
           if(inode==i+1) then
           call mpi_send(ug_n_tmp(1,iwavefun),mg_nx,
     &  MPI_DOUBLE_COMPLEX,0,102,MPI_COMM_WORLD,ierr)
           endif

           if(inode==1) then
           call mpi_recv(ugtemp,mg_nx,
     & MPI_DOUBLE_COMPLEX,i,102,MPI_COMM_WORLD,status,ierr)
           write(10) ugtemp   
           endif

        call mpi_barrier(MPI_COMM_WORLD,ierr)

          enddo  ! do iwavefun
         enddo    ! do inode

         if(inode==1) then
         close(10)
         deallocate(ugtemp)
         endif

        endif     ! for the iflag, write
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        return
        end

