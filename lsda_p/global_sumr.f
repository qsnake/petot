       subroutine global_sumr(sum_local_pe)
c--------------------------------------------------------

      implicit none
c
      include "mpif.h"
      real*8    sum_local_pe,res
      integer ierr
c
      call mpi_allreduce(sum_local_pe,res,1,
     $     MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

      sum_local_pe = res
      return
      end


