      subroutine  d3fft_comp(fg_n,fr_n,isign,kpt)

******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
******************************************


***************************************************
****
****    Warning ! fr_n(mr_n) must be statically allocated
****    before the call of this subroutine using:
****    real*8 fr_n(1)
****    pointer(fr_n_p,fr_n)
****    call shpalloc(fr_n_p,mr_n,errcode,-1)
****
****
***************************************************
****    warning !!!
****  for  fr_n --> fg_n, the original fr_n will be destroyed
****  this FFT includes a fg_n(i), i.e the smooth Ecut mask
****  As a result, the wavefunction in real space fr_n(i) are
****  no longer normalized and orthogonal. Use other FFT to
****  get the "real" (without wg_n) real space wavefunction.
***************************************************

***************************************************
****  d3fft_real(fg,fr,isign)
****  isign=1: fr_n --> fg_n   (forward fft)
****  isign=-1: fg_n --> fr_n  (inverse fft)
****  normalization: \sum_i fr_n(i)^2 *vol/(n1*n2*n3)=1
****  \sum_i fg_n(i)^2  * 2 * vol = 1
***************************************************


      use fft_data
      use load_data
      use data

      implicit none

      include "param.escan_real"

      include "mpif.h"
      
      
      complex*16 fg_n(mg_nx)
      complex*16 fr_n(mr_n)   
      complex*16 cdum

      integer i,ig,indepg,n1_inv,n2_inv,indepg_d,isign,
     &     jjnode_dum,kpt,ierr
      

      call mpi_barrier(mpi_comm_world,ierr)

      if(isign.eq.-1) then

c     put into x col format for fft in workr_n

      do i = 1,ncol(inode)*n1
         fr_n(i) = dcmplx(0.d0,0.d0)
      enddo

 
      call mpi_barrier(mpi_comm_world,ierr)

      do ig = 1, ngtotnod(inode,kpt)

         indepg=(jjcol(n2p_n(ig),n3p_n(ig))-1)*n1+n1p_n(ig)
         fr_n(indepg) = fg_n(ig)*wg_n(ig,kpt)

      enddo

      call mpi_barrier(mpi_comm_world,ierr)

      call invcpfft_comp(fr_n,n1,n2,n3)

      call mpi_barrier(mpi_comm_world,ierr)
      return
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc end of inversed FFT 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      if(isign.eq.1) then

      call fwdcpfft_comp(fr_n,n1,n2,n3)

c     put back into load balanced g vector distribution
c     

      do ig = 1, ngtotnod(inode,kpt)
         indepg=(jjcol(n2p_n(ig),n3p_n(ig))-1)*n1+n1p_n(ig)
         fg_n(ig)=wg_n(ig,kpt)*fr_n(indepg)
      enddo


      call mpi_barrier(mpi_comm_world,ierr)
      return
      endif


      return
      end


