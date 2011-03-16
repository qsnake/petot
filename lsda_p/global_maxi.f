c------------------------------------------------------
      subroutine global_maxi(inum)
c------------------------------------------------------
c output=inum 
c this routine calculates the global max of inum 
c answer put on each PE
c
c
      implicit none
c
      include "mpif.h"
      
      integer inum,inum_tmp,ierr
      integer inode,nnodes
      integer inode_tot,nnodes_tot,icolor,num_group,MPI_COMM_K,
     &       MPI_COMM_N

      common /mpi_data/inode,nnodes,inode_tot,nnodes_tot,
     &  icolor,num_group,MPI_COMM_K,MPI_COMM_N

      call mpi_allreduce(inum,inum_tmp,1,
     &  MPI_INTEGER,MPI_MAX,MPI_COMM_K,ierr)

      inum=inum_tmp

      return
      end

