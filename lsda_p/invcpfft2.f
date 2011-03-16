      subroutine invcpfft2(psi,nr1,nr2,nr3)

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
      
      real*8 combuf1(mr_n),combuf2(mr_n)
      real*8 psi(mr_n)
      real*8 psiy(mr_n)           ! work space
      real*8 psidum(nr3)

      integer inode,nnodes
      integer ireq(nnodes)

      integer jj

      common /mpi_data/inode,nnodes
      integer mpistatus(mpi_status_size)

c     scalars used

      integer i,ib,ic,idum,ii,ilocadd,isc,isign,itar,itaradd,
     c     itnode,iw,ix,iy,izb,j,jcol,ngy,ngyadd,ngz,nr1,nr2,nr3,
     c     nr3u,inode,iloc_dum,ierr
      integer nworknr1,nworknr2,nworknr3

      nworknr1 = 20000+2.28*nr1x
      nworknr2 = 20000+2.28*nr2x
      nworknr3 = 20000+2.28*nr3x

      allocate(worknr1(nworknr1))
      allocate(worknr2(nworknr2))
      allocate(worknr3(nworknr3))


      call mpi_barrier(MPI_COMM_WORLD,ierr)

c     
c     for gamma point fft
c     
      nr3u = nr3/2 + 1
c     
c     memory allocation for local arrays used with SHMEM routines
c     
c     
      isign = -1
c     
c     do the FFT's in place on the non-zero columns 
c     of psi
c     so this is the x dir FFT
c     


      do jcol = 1,ncol2(inode)
         ilocadd = 1+(jcol-1)*(2*nr1)

         call system_ccfft(0,isign,nr1,1.0d0,psi(ilocadd),
     &    psi(ilocadd),tabnr1lin,worknr1,0,ntabnr1,nworknr1)


      enddo                     !i



c     
c     transpose to the slice  y col mode  in psiy
c     
      do ii = 1,ncoly2*2*nr2
         psiy(ii) = 0.0
      enddo

      idum = 1
      do i = 1,nnodes
        do j = 1,ivpacn1l(i)
          combuf1(idum) = psi(ivpac1l(idum))
          idum = idum + 1
        enddo
      enddo

      call mpi_barrier(mpi_comm_world,ierr)

      idum = 1
      do i = 1,nnodes
       call mpi_isend(combuf1(idum),ivpacn1l(i),mpi_real8,i-1,inode,
     &                mpi_comm_world,ireq(i),ierr)
       idum = idum + ivpacn1l(i)
      enddo

      idum = 1
      do i = 1,nnodes
       call mpi_recv(combuf2(idum),ivunpn1l(i),mpi_real8,i-1,i,
     &               mpi_comm_world,mpistatus,ierr)
       idum = idum + ivunpn1l(i)
      enddo

      do i = 1, nnodes
        call mpi_wait(ireq(i), mpistatus, ierr)
      end do

      call mpi_barrier(mpi_comm_world,ierr)

      idum = 1
      do i = 1,nnodes
        do j = 1,ivunpn1l(i)
          psiy(ivunp1l(idum)) = combuf2(idum)
          idum = idum + 1
        enddo
      enddo
c     
c     do FFT on the y direction on the two (one for gamma) slices 
c     this is the y dir FFT
c     each PE should have ncoly columns for the FFT
c     


      do i = 1,ncoly2
         ilocadd = 1+(i-1)*(2*nr2)
         call system_ccfft(0,isign,nr2,1.0d0,psiy(ilocadd),
     c      psiy(ilocadd),tabnr2lin,worknr2,0,ntabnr2,nworknr2)
      enddo

c
c     now use strided put to put the z strips to the correct
c     processors. Put them into psi.
c     zero values of psi which we do not put into

      do i = 1,ncolz2
         idum = (i-1)*2*nr3u
         do j = 2*mgz2+1,2*nr3u
            psi(idum+j)=0.0
         enddo
      enddo

      idum = 1
      do i = 1,nnodes
        do j = 1,ivpacn2l(i)
          combuf1(idum) = psiy(ivpac2l(idum))
          idum  = idum + 1
        enddo
      enddo

      call mpi_barrier(mpi_comm_world,ierr)

      idum = 1
      do i = 1,nnodes
       call mpi_isend(combuf1(idum),ivpacn2l(i),mpi_real8,i-1,inode,
     &                mpi_comm_world,ireq(i),ierr)
       idum = idum + ivpacn2l(i)
      enddo

      idum = 1
      do i = 1,nnodes
       call mpi_recv(combuf2(idum),ivunpn2l(i),mpi_real8,i-1,i,
     &               mpi_comm_world,mpistatus,ierr)
       idum = idum + ivunpn2l(i)
      enddo

      do i = 1, nnodes
        call mpi_wait(ireq(i), mpistatus, ierr)
      end do

      call mpi_barrier(mpi_comm_world,ierr)

      idum = 1
      do i = 1,nnodes
        do j = 1,ivunpn2l(i)
          psi(ivunp2l(idum)) = combuf2(idum)
          idum = idum + 1
        enddo
      enddo
c

c     
c     now do FFT's on z direction so each proc has nr1*nr2/nnodes 
c     columns of height nr3 complex numbers (real for gamma point)
c     


      do i = 1,ncolz
         ilocadd = 1+(i-1)*(2*nr3u)

         call system_csfft(0,isign,nr3,1.0d0,psi(ilocadd),
     &    psidum,tabnr3lcr,worknr3,0,ntabnr3,nworknr3)

         iloc_dum = 1+(i-1)*nr3

         do ii = 1,nr3
            psi(iloc_dum-1+ii)=psidum(ii)
         enddo

      enddo



      call mpi_barrier(mpi_comm_world,ierr)

      deallocate(worknr1)
      deallocate(worknr2)
      deallocate(worknr3)

      end
      


