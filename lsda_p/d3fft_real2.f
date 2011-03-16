      subroutine  d3fft_real2(fg_n,fr_n,isign,ifac2)

******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

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

      include 'param.escan_real'
      include "mpif.h"

      complex*16 fg_n(mr_n/2)
      real*8 fr_n(mr_n)        
      complex*16 cdum

      integer i,j,ig,indepg,n1_inv,n2_inv,indepg_d,isign,
     &     jjnode_dum,ifac2,ierr
      integer ireq(nnodes)
      integer mpistatus(mpi_status_size)
      integer k0npac_max,k0nunp_max

c
c dummy arrays for building k.eq.0 plane for ffts
c
      real*8, allocatable,dimension(:,:) :: k0val,frdum




      if(isign.eq.-1) then

      k0npac_max=1
      k0nunp_max=1
      do i=1,nnodes
      if(k0npac(i).gt.k0npac_max) k0npac_max=k0npac(i)
      if(k0nunp(i).gt.k0nunp_max) k0nunp_max=k0nunp(i)
      enddo

      allocate(k0val(k0npac_max,nnodes))
      allocate(frdum(k0nunp_max,nnodes))


c     put into x col format for fft in workr_n

      do i = 1,ncol2(inode)*2*n1
         fr_n(i) = 0.0d0
      enddo
 
      k0npac = 0

      do ig = 1, ngtotnod2(inode)

         indepg=(jjcol2(n2p2_n(ig),n3p2_n(ig))-1)*n1+n1p2_n(ig)
         fr_n(2*indepg-1) = dreal(fg_n(ig))
         fr_n(2*indepg) = dimag(fg_n(ig))

c     
c     make up full k.eq.0 plane from half plane for ffts
c     
         if(n3p2_n(ig).eq.1) then 
            n1_inv = mod(n1-n1p2_n(ig)+1,n1)+1
            n2_inv = mod(n2-n2p2_n(ig)+1,n2)+1
            jjnode_dum = jjnode2(n2_inv,1)
            k0npac(jjnode_dum)=k0npac(jjnode_dum)+1
            k0val(k0npac(jjnode_dum),jjnode_dum)=fr_n(2*indepg-1)
            k0npac(jjnode_dum)=k0npac(jjnode_dum)+1
            k0val(k0npac(jjnode_dum),jjnode_dum)=-fr_n(2*indepg)
         endif

      enddo
c
c pass out half  k.eq.0 plane
c
      call mpi_barrier(MPI_COMM_K,ierr)

      do i = 1,nnodes
       call mpi_isend(k0val(1,i),k0npac(i),mpi_real8,i-1,inode,
     &                MPI_COMM_K,ireq(i),ierr)
      enddo

      do i = 1,nnodes
       call mpi_recv(frdum(1,i),k0nunp(i),mpi_real8,i-1,i,
     &               MPI_COMM_K,mpistatus,ierr)
      enddo

      do i = 1, nnodes
         call mpi_wait(ireq(i), mpistatus, ierr)
      end do

      call mpi_barrier(MPI_COMM_K,ierr)

      do i = 1,nnodes
        do j = 1,k0nunp(i)
          fr_n(k0iunp(j,i)) = frdum(j,i)
        enddo
      enddo

c     
c     deal with origin separately (it is its own inverse)
c     
      if(inode.eq.iorg2(1)) then
        ig = iorg2(2)
        indepg=(jjcol2(n2p2_n(ig),n3p2_n(ig))-1)*n1+n1p2_n(ig)
          if(ifac2.eq.1) then
          fr_n(2*indepg-1) = dsqrt(2.d0)*dreal(fg_n(ig))
          else
          fr_n(2*indepg-1) = dreal(fg_n(ig))
          endif
        fr_n(2*indepg) = 0.0
      endif

      call invcpfft2(fr_n,n1,n2,n3)

      call mpi_barrier(MPI_COMM_K,ierr)

      deallocate(k0val)
      deallocate(frdum)

      return
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc end of inversed FFT 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      if(isign.eq.1) then
      
      call fwdcpfft2(fr_n,n1,n2,n3)
c     
c     put back into load balanced g vector distribution
c     
      do ig = 1, ngtotnod2(inode)
         indepg=(jjcol2(n2p2_n(ig),n3p2_n(ig))-1)*n1+n1p2_n(ig)
         fg_n(ig)=dcmplx(fr_n(2*indepg-1),fr_n(2*indepg))
      enddo

      if(inode.eq.iorg2(1)) then
         ig = iorg2(2)
         indepg=(jjcol2(n2p2_n(ig),n3p2_n(ig))-1)*n1+n1p2_n(ig)
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
