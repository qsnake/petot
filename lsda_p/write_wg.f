      subroutine write_wg(fwg_out,AL,islda,nkpt)

*************************************************************************
*** Written by Lin-Wang Wang, 2001
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************


ccccccccccccccccccccccccccccccccccccccccccc
ccccc This should only be called by one icolor group, it will write all the kpt into one wg file 

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'mpif.h'
                                                                                                          
      include 'param.escan_real'

       real*8 AL(3,3)
       character*20 fwg_out(2)

       complex*16,allocatable,dimension(:) :: ugtemp
       complex*16,allocatable,dimension(:) :: temp_array

       integer status(MPI_STATUS_SIZE)

       do 280 iislda=1,islda
                                                                                                          
        allocate(temp_array(mg_nx))
                                                                                                          
        call mpi_barrier(MPI_COMM_K,ierr)
                                                                                                          
       if(inode==1) then
         allocate(ugtemp(mg_nx))
         open(11,file=fwg_out(iislda),form='unformatted')
         rewind(11)
         write(11)n1,n2,n3,mx
         write(11)Ecut
         write(11)AL
         write(11)nnodes
       end if
                  
       do 180 kpt=1,nkpt
                  
       call ugIO(ug_n,kpt,2,0,iislda)
                  
       if(inode==1) then
         do iwavefun=1,mx
            ugtemp=ug_n(:,iwavefun)
            write(11) ugtemp
         end do
       end if
      call mpi_barrier(MPI_COMM_K,ierr)
                  
      do i=1,nnodes-1
                  
         do iwavefun=1,mx
      call mpi_barrier(MPI_COMM_K,ierr)
            if(inode==i+1) then
               do iloop=1,mg_nx
                  temp_array(iloop)=ug_n(iloop,iwavefun)
               end do
              call mpi_send(temp_array,mg_nx,
     &  MPI_DOUBLE_COMPLEX,0,100,MPI_COMM_K,ierr)
            end if
                  
            if (inode==1) then
              call mpi_recv(temp_array,mg_nx,
     &  MPI_DOUBLE_COMPLEX,i,100,MPI_COMM_K,status,ierr)
                  
               do iloop=1,mg_nx
               ugtemp(iloop)=temp_array(iloop)
               end do
               write(11) ugtemp
            end if
            call mpi_barrier(MPI_COMM_K,ierr)
         end do   ! iwavefun
                  
      end do   ! inodes
                  
180   continue    ! kpt=1,nkpt
              
      if (inode==1) then
      close(11)
      deallocate(ugtemp)
      endif
                  
      deallocate(temp_array)
                  
280   continue   ! iislda=1,islda

      return
      end
              


       

       
      
