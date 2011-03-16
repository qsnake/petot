      subroutine init_ug(AL,iwg_in,workr_n,kpt,iranm)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
******************************************

****************************************
cccccc this is just for one kpt. 
****************************************
****  It stores the wavefunction in G space. 
******************************************
      use fft_data
      use load_data
      use data
      implicit double precision (a-h,o-z)

      include 'mpif.h'
      include 'param.escan_real'

      integer status(MPI_STATUS_SIZE)
******************************************
       real*8 AL(3,3),ALt(3,3)
***********************************************
       complex*16 workr_n(mr_n)


*************************************************
****  input ug from file  and check the consistency
*************************************************

       if(iwg_in.eq.1) then
       call read_ug()
       endif
************************************************
*** generate the initial wavefunction ug from random if iwg_in.eq.0
************************************************
      if(iwg_in.eq.0) then
****  this is the one should be used later
	do m=1,mx
	do i=1,nr_n
        x1=ran1(iranm)
        x2=ran1(iranm)
        workr_n(i)=dcmplx(x1-0.5d0,x2-0.5d0)
cc        workr_n(i)=dcmplx(x1-0.5d0,0.d0)
	enddo
************************************************

        call mpi_barrier(MPI_COMM_WORLD,ierr)

        call d3fft_comp(ug_n(1,m),workr_n,1,kpt)

      enddo
      endif
*************************************************
**** end generate the initial wavefunction from random
*************************************************
        do m=1,mx
         call orth_comp(ug_n(1,m),ug_n,m-1,1,kpt)
	enddo
*************************************************
      return
      contains
*************************************************

       subroutine read_ug()

       complex*16,allocatable,dimension(:) :: ugtemp

        call mpi_barrier(MPI_COMM_WORLD,ierr)
cccccc  file (11) should have been opened before the do loop for kpt=1,nkpt

ccccccccccccc check the header.
        if(inode==nnodes.and.kpt.eq.1) then

         read(11)n1t,n2t,n3t,mxt
         read(11)Ecutt
         read(11)ALt
         read(11)nnodes_o

         if(n1t.ne.n1.or.n2t.ne.n2.or.n3t.ne.n3) then
         write(6,*) "n1t,n2t,n3t changed, stop", n1t,n2t,n3t
         call mpi_abort(MPI_COMM_WORLD,ierr)
         endif

         if(Ecut.ne.Ecutt.or.mxt.ne.mx) then
         write(6,*) "Ecutt,mxt changed, stop", Ecutt,mxt
         call mpi_abort(MPI_COMM_WORLD,ierr)
         endif

         diff=0.d0
         do i=1,3
            do j=1,3
               diff=diff+dabs(AL(i,j)-ALt(i,j))
             enddo
         enddo

         if(diff.gt.0.001d0) then
            write(6,*) "the AL.ne.ALt, stop"
         call mpi_abort(MPI_COMM_WORLD,ierr)
         endif

         if(nnodes_o.ne.nnodes) then
           write(6,*) "nnodes_o changed, stop", nnodes_o
         call mpi_abort(MPI_COMM_WORLD,ierr)
         endif

      endif     ! check the header, make sure it is the same
ccccccccccccccccccccccccccccccccccccccccc
         call mpi_barrier(MPI_COMM_WORLD,ierr)

         if(inode==nnodes) then
            allocate(ugtemp(mg_nx))

            do i=0,nnodes-2
               do iwavefun=1,mx
                  read(11)ugtemp    ! assuming the mg_nx is the same as before
                  ug_n(:,iwavefun)=ugtemp
               end do


          call mpi_send(ug_n,mg_nx*mx,MPI_DOUBLE_COMPLEX,i,
     &            102,MPI_COMM_WORLD,ierr)
            end do
c     Now read in the data for node nnodes
            do iwavefun=1,mx
               read(11)ugtemp
               ug_n(:,iwavefun)=ugtemp
            end do
            deallocate(ugtemp)

         else   ! for other nodes

         call mpi_recv(ug_n,mg_nx*mx,MPI_DOUBLE_COMPLEX,nnodes-1,
     &       102,MPI_COMM_WORLD,status,ierr)

         end if

        return
        end subroutine read_ug

        end

