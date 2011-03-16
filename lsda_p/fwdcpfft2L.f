      subroutine fwdcpfft2L(psi,nr1L,nr2L,nr3L)
c------------------------------------------------------
c
c Written by A. Canning (CRAY-EPFL) 25th July 94  
c 
c
c output = psi(x,y,z)  wavefunction in fourier space load balanced
c                      as small sphere columns in x direction
c input = psi(z,y,x) real space grid with each PE having consecutive 
c                    sets of z columns.
c
c computes forward fft specifically for CP algo
c ie taking cube and going to two slices (one slice for gamma)
c then  to a sphere.  See fftprep.f for more details
c fftprep must be called once before this routine
c for setups.
       
      use fft_data
      use load_data
      use data

      implicit none

      include  'param.escan_real'
      include "mpif.h"

      real*8,allocatable,dimension(:) :: worknr1,worknr2,worknr3
      
      real*8 psi(mr_nL)
      real*8,allocatable,dimension(:) :: psiy,combuf1,combuf2
      real*8 fac,s,sum
      
c scalars used

      integer i,ib,ic,idum,ii,iloc,ilocadd,isc,isign,itar,itaradd,
     c        itnode,iw,ix,iy,izb,j,jcol,ngy,ngyadd,ngz,nr1L,
     c        nr2L,nr3L,nr3uL,iloc_dum,ierr
      integer nworknr1,nworknr2,nworknr3

c mpi version arrays
c
      integer mpistatus(mpi_status_size)
      integer ireq(nnodes)

      call mpi_barrier(mpi_comm_k,ierr)
c
      nworknr1 = 20000+2.28*nr1xL
      nworknr2 = 20000+2.28*nr2xL
      nworknr3 = 20000+2.28*nr3xL

      allocate(worknr1(nworknr1))
      allocate(worknr2(nworknr2))
      allocate(worknr3(nworknr3))
      allocate(psiy(mr_nL))
      allocate(combuf1(mr_nL))
      allocate(combuf2(mr_nL))
      
c
c for gamma point fft
c
      nr3uL = nr3L/2 + 1
c
c memory allocation for local arrays used with SHMEM routines
c with SHMEM routines address on all PE's must be the same for 
c same variable
c
c
      isign = 1
c
c full FFT on the cube psi this is on the z dir
c layout is nz,ny,nx  colulmns of height nr3 complex
c reslt not put directly into psi (problem for scfft)
c due to different sizes of input and output
c


       do i = 1,ncolz2L
        ilocadd = 1+(i-1)*nr3L
       iloc_dum = 1+(i-1)*(2*nr3uL)

        call system_scfft(0,isign,nr3L,1.0d0,psi(ilocadd),
     &    psiy(iloc_dum),tabnr3lrcL,worknr3,0,ntabnr3L,nworknr3)
       enddo
c
c now transpose nz,ny,nx to ny,nz,nx in the two slice mode
c into psiy. Each PE will have ncoly columns in psiy
c
       

      idum = 1
      do i = 1,nnodes
        do j = 1,ivunpn2lL(i)
          combuf1(idum) = psiy(ivunp2lL(idum))
          idum = idum + 1
        enddo
      enddo

      call mpi_barrier(mpi_comm_k,ierr)

      idum = 1
      do i = 1,nnodes
       call mpi_isend(combuf1(idum),ivunpn2lL(i),mpi_real8,i-1,
     &              inode,mpi_comm_k,ireq(i),ierr)
       idum = idum + ivunpn2lL(i)
      enddo

      idum = 1
      do i = 1,nnodes
       call mpi_recv(combuf2(idum),ivpacn2lL(i),mpi_real8,i-1,i,
     &               mpi_comm_k,mpistatus,ierr)
       idum = idum + ivpacn2lL(i)
      enddo

      do i = 1, nnodes
         call mpi_wait(ireq(i), mpistatus, ierr)
      end do

      idum = 1
      do i = 1,nnodes
        do j = 1,ivpacn2lL(i)
          psiy(ivpac2lL(idum)) = combuf2(idum)
          idum = idum + 1
        enddo
      enddo


c
c now do FFT's on the one slice
c
       do i = 1,ncoly2L
        ilocadd = 1+(i-1)*(2*nr2L)
       call system_ccfft(0,isign,nr2L,1.0d0,psiy(ilocadd),
     c     psiy(ilocadd),tabnr2lfwL,worknr2,0,ntabnr2L,nworknr2)
       enddo

c
c now transpose back to format of program into psi
c ie into x columns load balanced
c
      idum = 1
      do i = 1,nnodes
        do j = 1,ivunpn1lL(i)
          combuf1(idum) = psiy(ivunp1lL(idum))
          idum  = idum + 1
        enddo
      enddo

      call mpi_barrier(mpi_comm_k,ierr)

      idum = 1
      do i = 1,nnodes
       call mpi_isend(combuf1(idum),ivunpn1lL(i),mpi_real8,i-1,
     &            inode,mpi_comm_k,ireq(i),ierr)
       idum = idum + ivunpn1lL(i)
      enddo

      idum = 1
      do i = 1,nnodes
       call mpi_recv(combuf2(idum),ivpacn1lL(i),mpi_real8,i-1,i,
     &               mpi_comm_k,mpistatus,ierr)
       idum = idum + ivpacn1lL(i)
      enddo

      do i = 1, nnodes
         call mpi_wait(ireq(i), mpistatus, ierr)
      end do

      idum = 1
      do i = 1,nnodes
        do j = 1,ivpacn1lL(i)
          psi(ivpac1lL(idum)) = combuf2(idum)
          idum = idum + 1
        enddo
      enddo

c
c now do last FFT in place in the format of the program 
c for psi in the x direction 
c
       do jcol = 1,ncol2L(inode)
        ilocadd = 1+(jcol-1)*(2*nr1L)
        call system_ccfft(0,isign,nr1L,1.0d0,psi(ilocadd),
     c    psi(ilocadd),tabnr1lfwL,worknr1,0,ntabnr1L,nworknr1)
       enddo !i
c
c now do the scaling of the FFT
c
      fac=1.d0/dfloat(nr1L*nr2L*nr3L)
 
      do i = 1,ncol2L(inode)*2*nr1L
        psi(i) = psi(i)*fac
      enddo


      call mpi_barrier(MPI_COMM_K,ierr)

      deallocate(worknr1)
      deallocate(worknr2)
      deallocate(worknr3)
      deallocate(psiy)
      deallocate(combuf1)
      deallocate(combuf2)

      end



