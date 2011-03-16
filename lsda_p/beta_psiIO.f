      subroutine beta_psiIO(beta_psi_tmp,kpt,iflag,istep,
     &  iislda)
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

       complex*16 beta_psi_tmp(nref_tot,mx/nnodes+1)
       character*9 fname
*************************************************

       mx_n=mx/nnodes+1
*************************************************
       if(iflag.eq.2) then     ! read the wavefunction
       if(inode==nnodes) then

       fname="bpsiiofil"
       kpt1=mod(kpt,10)
       kpt2=mod((kpt-kpt1)/10,10)
       kpt3=mod((kpt-kpt1-kpt2*10)/100,10)
       kpt4=mod((kpt-kpt1-kpt2*10-kpt3*100)/1000,10)

       open(10,file=
     & fname//char(iislda+48)//char(istep+48)//char(kpt4+48)//
     &  char(kpt3+48)//char(kpt2+48)//char(kpt1+48),
     &       form="unformatted")


            do i=0,nnodes-2
                  read(10) beta_psi_tmp


          call mpi_send(beta_psi_tmp,nref_tot*mx_n,
     & MPI_DOUBLE_COMPLEX,i,102,MPI_COMM_K,ierr)
            end do

c     Now read in the data for node nnodes
               read(10) beta_psi_tmp
            close(10)

         else   ! for other nodes

         call mpi_recv(beta_psi_tmp,nref_tot*mx_n,
     & MPI_DOUBLE_COMPLEX,nnodes-1,102,MPI_COMM_K,status,ierr)

         endif    ! for the nodes
        endif     ! for the iflag, read
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       if(iflag.eq.1) then     ! write the wavefunction

       if(inode==1) then
       fname="bpsiiofil"
       kpt1=mod(kpt,10)
       kpt2=mod((kpt-kpt1)/10,10)
       kpt3=mod((kpt-kpt1-kpt2*10)/100,10)
       kpt4=mod((kpt-kpt1-kpt2*10-kpt3*100)/1000,10)

       open(10,file=
     & fname//char(iislda+48)//char(istep+48)//char(kpt4+48)//
     &  char(kpt3+48)//char(kpt2+48)//char(kpt1+48),
     &       form="unformatted")

               write(10) beta_psi_tmp
        endif    ! inode==1

        call mpi_barrier(MPI_COMM_K,ierr)

        do i=1,nnodes-1

        call mpi_barrier(MPI_COMM_K,ierr)

           if(inode==i+1) then
           call mpi_send(beta_psi_tmp,nref_tot*mx_n,
     &  MPI_DOUBLE_COMPLEX,0,102,MPI_COMM_K,ierr)
           endif

           if(inode==1) then
           call mpi_recv(beta_psi_tmp,nref_tot*mx_n,
     & MPI_DOUBLE_COMPLEX,i,102,MPI_COMM_K,status,ierr)
           write(10) beta_psi_tmp
           endif

         enddo    ! do inode

         if(inode==1) then
         close(10)
         endif

        endif     ! for the iflag, write
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        return
        end

