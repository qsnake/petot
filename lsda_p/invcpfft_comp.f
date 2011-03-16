      subroutine invcpfft_comp(psi,nr1,nr2,nr3)

c     Written by A. Canning (CRAY-EPFL) 25th July 94  
c     
c     output = psi(z,y,x)  wavefunction in real space z columns 
c     input  = psi(x,y,z)  x columns of g vectors in load balalanced distr.
c     
c     computes inverse fft specifically for CP algo 
c     ie taking sphere and going to cylinder then two (one for gamma)
c     slices then a complete cube. 
c     see fftprep.f for more details
c     fftprep must be called once before this routine
c     for setups.
c     
c     last revised ADV 7/7/95 to parametrise the no. of 
c     chunks on each PE
c     AMC revised 21/6/97 hardwired for gamma point complex to real fft

      use load_data
      use fft_data
      use data
      
      implicit none

      include "mpif.h"
c     
c     
      real*8,allocatable,dimension(:) :: worknr1,worknr2,worknr3
      complex*16,allocatable,dimension(:) :: psiy,combuf1,combuf2
      
      complex*16 psi(mr_n)

      integer ireq(nnodes)

      integer inode,nnodes
      integer inode_tot,nnodes_tot,icolor,num_group,MPI_COMM_K,
     &       MPI_COMM_N


      common /mpi_data/inode,nnodes,inode_tot,nnodes_tot,
     &  icolor,num_group,MPI_COMM_K,MPI_COMM_N

      integer mpistatus(mpi_status_size)

c     scalars used

      integer i,ib,ic,idum,ii,ilocadd,isc,isign,itar,itaradd,
     c     itnode,iw,ix,iy,izb,j,jcol,ngy,ngyadd,ngz,nr1,nr2,nr3,
     c     nr3u,iloc_dum,ierr
      integer nworknr1,nworknr2,nworknr3


      nworknr1 = 20000+2.28*nr1x
      nworknr2 = 20000+2.28*nr2x
      nworknr3 = 20000+2.28*nr3x

      allocate(worknr1(nworknr1))
      allocate(worknr2(nworknr2))
      allocate(worknr3(nworknr3))
      allocate(psiy(mr_n))
      allocate(combuf1(mr_n))
      allocate(combuf2(mr_n))

c     
      isign = -1
c     
c     do the FFT's in place on the non-zero columns 
c     of psi
c     so this is the x dir FFT
c     


      do jcol = 1,ncol(inode)
         ilocadd = 1+(jcol-1)*nr1

         call system_ccfft(0,isign,nr1,1.0d0,psi(ilocadd),
     &     psi(ilocadd),tabnr1in,worknr1,0,ntabnr1,nworknr1)

c         call ccfft(isign,nr1,1.0,psi(ilocadd),psi(ilocadd),
c     &        tabnr1,worknr1,0)
c        call dcft(0, psi(ilocadd), 1, 0,
c    &                psi(ilocadd), 1, 0,
c    &                nr1, 1, -isign, 1.0,
c    &                tabnr1, ntabnr1, worknr1, nworknr1)
      enddo                     !i
c     
c     transpose to the y slice  y col mode  in psiy
c     
      do ii = 1,ncoly*2*nr2
         psiy(ii) = 0.0
      enddo

      idum = 1
      do i = 1,nnodes
        do j = 1,ivpacn1(i)
          combuf1(idum) = psi(ivpac1(idum))
          idum = idum + 1
        enddo
      enddo

      idum = 1
      do i = 1,nnodes
       call mpi_isend(combuf1(idum),ivpacn1(i),mpi_double_complex,i-1,
     &                inode,mpi_comm_k,ireq(i),ierr)
       idum = idum + ivpacn1(i)
      enddo

      idum = 1
      do i = 1,nnodes
       call mpi_recv(combuf2(idum),ivunpn1(i),mpi_double_complex,i-1,i,
     &               mpi_comm_k,mpistatus,ierr)
       idum = idum + ivunpn1(i)
      enddo

      do i = 1, nnodes
        call mpi_wait(ireq(i), mpistatus, ierr)
      end do

      idum = 1
      do i = 1,nnodes
        do j = 1,ivunpn1(i)
          psiy(ivunp1(idum)) = combuf2(idum)
          idum = idum + 1
        enddo
      enddo
c     
c     do FFT on the y direction on the two (one for gamma) slices 
c     this is the y dir FFT
c     each PE should have ncoly columns for the FFT
c     

      do i = 1,2*ncoly
         ilocadd = 1+(i-1)*nr2
         call system_ccfft(0,isign,nr2,1.0d0,psiy(ilocadd),
     &    psiy(ilocadd),tabnr2in,worknr2,0,ntabnr2,nworknr2)
c         call ccfft(isign,nr2,1.0,psiy(ilocadd),psiy(ilocadd),
c     &        tabnr2,worknr2,0)
c        call dcft(0, psiy(ilocadd), 1, 0,
c    &                psiy(ilocadd), 1, 0,
c    &                nr2, 1, -isign, 1.0,
c    &                tabnr2, ntabnr2, worknr2, nworknr2)
      enddo
c     
c    transpose slice using mpi calls
c     zero values of psi which we do not put into
c
      do i = 1,ncolz
         idum = (i-1)*nr3
         do j = mgz+1,(nr3-mgz)
            psi(idum+j)=0.0
         enddo
      enddo

      idum = 1
      do i = 1,nnodes
        do j = 1,ivpacn2(i)
          combuf1(idum) = psiy(ivpac2(idum))
          idum  = idum + 1
        enddo
      enddo

      idum = 1
      do i = 1,nnodes
       call mpi_isend(combuf1(idum),ivpacn2(i),mpi_double_complex,i-1,
     &                inode,mpi_comm_k,ireq(i),ierr)
       idum = idum + ivpacn2(i)
      enddo

      idum = 1
      do i = 1,nnodes
       call mpi_recv(combuf2(idum),ivunpn2(i),mpi_double_complex,i-1,i,
     &               mpi_comm_k,mpistatus,ierr)
       idum = idum + ivunpn2(i)
      enddo

      do i = 1, nnodes
        call mpi_wait(ireq(i), mpistatus, ierr)
      end do

      idum = 1
      do i = 1,nnodes
        do j = 1,ivunpn2(i)
          psi(ivunp2(idum)) = combuf2(idum)
          idum = idum + 1
        enddo
      enddo

c     
c     now do FFT's on z direction so each proc has nr1*nr2/nnodes 
c     columns of height nr3 complex numbers (real for gamma point)
c     

      do i = 1,ncolz
         ilocadd = 1+(i-1)*nr3
         call system_ccfft(0,isign,nr3,1.0d0,psi(ilocadd),
     &      psi(ilocadd),tabnr3in,worknr3,0,ntabnr3,nworknr3)
c         call ccfft(isign,nr3,1.0,psi(ilocadd),psi(ilocadd),
c     &        tabnr3,worknr3,0)
c        call dcft(0, psi(ilocadd), 1, 0,
c    &                psi(ilocadd), 1, 0,
c    &                nr3, 1, -isign, 1.0,
c    &                tabnr3, ntabnr3, worknr3, nworknr3)
      enddo

      deallocate(worknr1)
      deallocate(worknr2)
      deallocate(worknr3)
      deallocate(psiy)
      deallocate(combuf1)
      deallocate(combuf2)

      end
      


