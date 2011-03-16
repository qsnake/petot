c------------------------------------------------------
      subroutine global_maxi(inum,nnodes)
c------------------------------------------------------
c output=inum 
c this routine calculates the global max of inum 
c answer put on each PE
c
c
      implicit none
c
      include "mpif.h"
      integer inum,nnodes,inum_tmp,ierr

      call mpi_allreduce(inum,inum_tmp,1,
     &  MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)

      inum=inum_tmp

      return
      end

