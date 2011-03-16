      subroutine  d3fft_real2L2(fg_n,fr_n,isign,ifac2)

******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

***************************************************
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

      include 'param.escan_real'
      include "mpif.h"

      complex*16 fg_n(mr_nL2/2)
      real*8 fr_n(mr_nL2)        
      complex*16 cdum
      real*8 sum

      integer i,j,ig,indepg,n1_inv,n2_inv,indepg_d,isign,
     &     jjnode_dum,ifac2,ierr
      integer ireq(nnodes)
      integer mpistatus(mpi_status_size)
      integer k0npac_max,k0nunp_max

c
c dummy arrays for building k.eq.0 plane for ffts
c
      real*8, allocatable,dimension(:,:) :: k0valL2,frdum



      if(isign.eq.-1) then

      k0npac_max=1
      k0nunp_max=1
      do i=1,nnodes
      if(k0npacL2(i).gt.k0npac_max) k0npac_max=k0npacL2(i)
      if(k0nunpL2(i).gt.k0nunp_max) k0nunp_max=k0nunpL2(i)
      enddo

      allocate(k0valL2(k0npac_max,nnodes))
      allocate(frdum(k0nunp_max,nnodes))

     
c     put into x col format for fft in workr_n

      do i = 1,ncol2L2(inode)*2*n1L2
         fr_n(i) = 0.0d0
      enddo
 
      k0npacL2 = 0

      do ig = 1, ngtotnod2L2(inode)

        indepg=(jjcol2L2(n2p2_nL2(ig),n3p2_nL2(ig))-1)*n1L2+
     &              n1p2_nL2(ig)
         fr_n(2*indepg-1) = dreal(fg_n(ig))
         fr_n(2*indepg) = dimag(fg_n(ig))

c     
c     make up full k.eq.0 plane from half plane for ffts
c     
         if(n3p2_nL2(ig).eq.1) then 
            n1_inv = mod(n1L2-n1p2_nL2(ig)+1,n1L2)+1
            n2_inv = mod(n2L2-n2p2_nL2(ig)+1,n2L2)+1
            jjnode_dum = jjnode2L2(n2_inv,1)
            k0npacL2(jjnode_dum)=k0npacL2(jjnode_dum)+1
            k0valL2(k0npacL2(jjnode_dum),jjnode_dum)=fr_n(2*indepg-1)
            k0npacL2(jjnode_dum)=k0npacL2(jjnode_dum)+1
            k0valL2(k0npacL2(jjnode_dum),jjnode_dum)=-fr_n(2*indepg)
         endif

      enddo
c
c pass out half  k.eq.0 plane
c

      call mpi_barrier(MPI_COMM_K,ierr)

      do i = 1,nnodes
       call mpi_isend(k0valL2(1,i),k0npacL2(i),mpi_real8,i-1,inode,
     &                MPI_COMM_K,ireq(i),ierr)
      enddo

      do i = 1,nnodes
       call mpi_recv(frdum(1,i),k0nunpL2(i),mpi_real8,i-1,i,
     &               MPI_COMM_K,mpistatus,ierr)
      enddo

      do i = 1, nnodes
         call mpi_wait(ireq(i), mpistatus, ierr)
      end do
      call mpi_barrier(MPI_COMM_K,ierr)



      do i = 1,nnodes
        do j = 1,k0nunpL2(i)
          fr_n(k0iunpL2(j,i)) = frdum(j,i)
        enddo
      enddo

c     
c     deal with origin separately (it is its own inverse)

c     
      if(inode.eq.iorg2L2(1)) then
        ig = iorg2L2(2)
        indepg=(jjcol2L2(n2p2_nL2(ig),n3p2_nL2(ig))-1)*n1L2+
     &         n1p2_nL2(ig)
          if(ifac2.eq.1) then
          fr_n(2*indepg-1) = dsqrt(2.d0)*dreal(fg_n(ig))
          else
          fr_n(2*indepg-1) = dreal(fg_n(ig))
          endif
        fr_n(2*indepg) = 0.0
      endif


      call invcpfft2L2(fr_n,n1L2,n2L2,n3L2)


      call mpi_barrier(MPI_COMM_K,ierr)

      deallocate(k0valL2)
      deallocate(frdum)

      return
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc end of inversed FFT 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      if(isign.eq.1) then

      
      call fwdcpfft2L2(fr_n,n1L2,n2L2,n3L2)
c     
c     put back into load balanced g vector distribution
c     
      do ig = 1, ngtotnod2L2(inode)
       indepg=(jjcol2L2(n2p2_nL2(ig),n3p2_nL2(ig))-1)*n1L2+
     &          n1p2_nL2(ig)
       fg_n(ig)=dcmplx(fr_n(2*indepg-1),fr_n(2*indepg))
      enddo

      if(inode.eq.iorg2L2(1)) then
         ig = iorg2L2(2)
        indepg=(jjcol2L2(n2p2_nL2(ig),n3p2_nL2(ig))-1)*n1L2+
     &          n1p2_nL2(ig)
          if(ifac2.eq.1) then
           fg_n(ig)= 1.d0/dsqrt(2.d0)*dcmplx(fr_n(2*indepg-1),0.d0)
           else
           fg_n(ig)= dcmplx(fr_n(2*indepg-1),0.d0)
           endif
      endif

      call mpi_barrier(MPI_COMM_K,ierr)
      return
      endif


      return
      end
