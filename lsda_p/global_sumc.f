         subroutine global_sumc(sum_local_pe)
c--------------------------------------------------------

      implicit none
c
      include "mpif.h"
      complex*16    sum_local_pe,res
c
      integer while_counter, neighbor,infoa,ierr
      complex*16    sum_other_pe
c
      call mpi_allreduce(sum_local_pe,res,1,
     $     MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_WORLD,ierr)

      sum_local_pe = res
      return
      end


