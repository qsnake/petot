       subroutine global_sumr(sum_local_pe)
c--------------------------------------------------------

      implicit none
c
      include "mpif.h"
      real*8    sum_local_pe,res
      integer ierr

      integer inode,nnodes
      integer inode_tot,nnodes_tot,icolor,num_group,MPI_COMM_K,
     &       MPI_COMM_N
                                                                                                          
      common /mpi_data/inode,nnodes,inode_tot,nnodes_tot,
     &  icolor,num_group,MPI_COMM_K,MPI_COMM_N

c
      call mpi_allreduce(sum_local_pe,res,1,
     $     MPI_REAL8,MPI_SUM,MPI_COMM_K,ierr)

      sum_local_pe = res
      return
      end


