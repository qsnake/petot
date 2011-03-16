      subroutine rhoIO(AL,vr_n_tmp,mr_ni,iflag,
     &    n1i,n2i,n3i,fdens_name)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************


cccc when iflag=1, the original vr_n_tmp will not be destroyed.

ccccccc iflag=1, write the charge density vr_n_tmp to fdens_name
ccccccc iflag=2, read the charge density vr_n_tmp from fdens_name

      implicit double precision (a-h,o-z)

      include 'mpif.h'
      include 'param.escan_real'


      real*8,allocatable,dimension(:)   :: tempvr1,tempvr2

      real*8,allocatable,dimension(:)   :: temp_rho
    
      real*8 vr_n_tmp(mr_ni)

      integer status(MPI_STATUS_SIZE)
      real*8 AL(3,3),ALt(3,3)
      character*20 fdens_name

      nri=n1i*n2i*n3i


      if(iflag.eq.2) then   ! read the charge density
      if(icolor.eq.0) then    ! only read from the icolor group
     
      if (inode==1) then
      open(unit=15,file=fdens_name,form='unformatted',
     & status='old',action='read',iostat=ierr)
      if(ierr.ne.0) then
      write(6,*) "file ***",fdens_name,"*** does not exist, stop"
      call mpi_abort(MPI_COMM_WORLD,ierr)
      endif

      read(15)n1t,n2t,n3t,nodes

      if(n1t.ne.n1i.or.n2t.ne.n2i.or.n3t.ne.n3i) then
      write(6,*) "n1,n2,n3 in input dens file not right, stop"
      write(6,*) n1t,n2t,n3t
      stop
      endif

c     Calculate the ratio of the no. of nodes used to generate the
c     potential and the no. used for the escan calculation.
c     We assume iratio=2^n where n can be positive or negative.

         ratio=dble(nodes)/dble(nnodes)

         test=dlog(ratio)/dlog(2.d0)
         itest=test*1.00001d0
         if(dabs(test-itest).gt.0.001d0.and.inode.eq.1) then
       write(6,*) "rhoIO, nodes/nnodes.ne.2**n, stop", nodes,nnodes
         call mpi_abort(MPI_COMM_WORLD,ierr)
         endif
         


         read(15)ALt

      diff=0.d0
      do i=1,3
      do j=1,3
      diff=diff+(ALt(i,j)-AL(i,j))**2
      enddo
      enddo
      if(diff.gt.1.E-5) then
      write(6,*) "AL in the input dens file not right, stop"
      stop
      endif

ccc         write(6,*)'ratio=',ratio
      endif
ccccccccccccccccccccccccccccccc
      call mpi_barrier(MPI_COMM_K,ierr)

      if(inode.eq.1) then

         allocate(tempvr1(nri/nnodes))
         allocate(tempvr2(nri/nodes))
         if (ratio>=1.0d0) then
            iratio=int(ratio+1.0d-06)
            do isend=0,nnodes-1
               do iread=1,iratio
                  read(15) tempvr2
                  tempvr1((iread-1)*nri/nodes+1:iread*nri/nodes)
     &                 =tempvr2
               end do
               if (isend==0) then
                  vr_n_tmp(1:nri/nnodes)=tempvr1
               else
                  call mpi_send(tempvr1,nri/nnodes,MPI_REAL8,isend,
     &                 100,MPI_COMM_K,ierr)
               end if
            end do
         else
            isend=0
            iratio=int(1.0d0/ratio+1.0d-06)
cc            write(6,*)'iratio = ',iratio
            do iloop=1,nodes
               read(15)tempvr2
               do i=1,iratio
                  tempvr1=tempvr2((i-1)*nri/nnodes+1:i*nri/nnodes)
                  if (isend==0) then
                     vr_n_tmp(1:nri/nnodes)=tempvr1
                  else
                     call mpi_send(tempvr1,nri/nnodes,MPI_REAL8,isend,
     &                    100,MPI_COMM_K,ierr)
                  end if
                  isend=isend+1
               end do
            end do
         end if
         deallocate(tempvr1)
         deallocate(tempvr2)
      else
         call mpi_recv(vr_n_tmp,nri/nnodes,MPI_REAL8,0,100,
     &        MPI_COMM_K,status,ierr)
      end if
      call mpi_barrier(MPI_COMM_K,ierr)

234   continue
      if(inode.eq.1) close(15) 

      endif     ! icolor=0 group
ccccccccc  send the data from icolor group to other groups
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call mpi_bcast(vr_n_tmp,nri/nnodes,MPI_REAL8,0,MPI_COMM_N,ierr)

      return
      endif     ! iflag=2, read
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(iflag.eq.1) then     ! write the charge 
      if(icolor.eq.0) then

       allocate(temp_rho(mr_ni))

       if(inode.eq.1) then
       allocate(tempvr1(mr_ni))
       open(11,file=fdens_name,form="unformatted")
       rewind(11)
       write(11) n1i,n2i,n3i,nnodes
       write(11) AL
       write(11) (vr_n_tmp(i),i=1,n1i*n2i*n3i/nnodes)
       endif
       call mpi_barrier(MPI_COMM_K,ierr)

       do i=1,nnodes-1
       call mpi_barrier(MPI_COMM_K,ierr)
       if(inode==i+1) then
        do iloop=1,mr_ni
        temp_rho(iloop)=vr_n_tmp(iloop)
        enddo
        call  mpi_send(temp_rho,mr_ni,MPI_REAL8,0,
     &   100,MPI_COMM_K,ierr)

       endif

       if(inode==1) then
        call mpi_recv(temp_rho,mr_ni,MPI_REAL8,i,
     &   100,MPI_COMM_K,status,ierr)
       do iloop=1,mr_ni
       tempvr1(iloop)=temp_rho(iloop)
       end do
       write(11) tempvr1
       endif
       enddo
       call mpi_barrier(MPI_COMM_K,ierr)

        deallocate(temp_rho)

        if(inode.eq.1) then
        deallocate(tempvr1)
        close(11)
        endif
*****************************************************
        endif    ! icolor.eq.0
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        endif    ! iflag=1, write

        return
        end
      
