      subroutine fwdcpfft_comp(psi,nr1,nr2,nr3)
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

      include "mpif.h"

      real*8,allocatable,dimension(:) :: worknr1,worknr2,worknr3
      complex*16,allocatable,dimension(:) :: psiy,combuf1,combuf2
      
      complex*16 psi(mr_n)

      integer inode,nnodes
      integer inode_tot,nnodes_tot,icolor,num_group,MPI_COMM_K,
     &       MPI_COMM_N

      common /mpi_data/inode,nnodes,inode_tot,nnodes_tot,
     &  icolor,num_group,MPI_COMM_K,MPI_COMM_N

      integer mpistatus(mpi_status_size)
      integer ireq(nnodes)

c scalars used

      integer i,ib,ic,idum,ii,iloc,ilocadd,isc,isign,itar,itaradd,
     c        itnode,iw,ix,iy,izb,j,jcol,ngy,ngyadd,ngz,nr1,nr2,nr3,
     c        nr3u,iloc_dum,ierr
      integer nworknr1,nworknr2,nworknr3

      real*8 fac,s

c
      nworknr1 = 20000+2.28*nr1x
      nworknr2 = 20000+2.28*nr2x
      nworknr3 = 20000+2.28*nr3x

      allocate(worknr1(nworknr1))
      allocate(worknr2(nworknr2))
      allocate(worknr3(nworknr3))
      allocate(psiy(mr_n))
      allocate(combuf1(mr_n))
      allocate(combuf2(mr_n))



      call mpi_barrier(mpi_comm_k,ierr)

c
      isign = 1
c
c full FFT on the cube psi this is on the z dir
c layout is nz,ny,nx  colulmns of height nr3 complex
c reslt not put directly into psi (problem for scfft)
c due to different sizes of input and output
c

       do i = 1,ncolz
        ilocadd = 1+(i-1)*nr3
        call system_ccfft(0,isign,nr3,1.0d0,psi(ilocadd),
     &     psi(ilocadd),tabnr3fw,worknr3,0,ntabnr3,nworknr3)

c        call ccfft(isign,nr3,1.0,psi(ilocadd),psi(ilocadd),
c     &             tabnr3,worknr3,0)
c        call dcft(0, psi(ilocadd), 1, 0,
c    &                psi(ilocadd), 1, 0,
c    &                nr3, 1, -isign, 1.0,
c    &                tabnr3, ntabnr3, worknr3, nworknr3)

       enddo


c
c now transpose nz,ny,nx to ny,nz,nx in the two slice mode
c into psiy. Each PE will have ncoly columns in psiy 
c

      idum = 1
      do i = 1,nnodes
        do j = 1,ivunpn2(i)
          combuf1(idum) = psi(ivunp2(idum))
          idum = idum + 1
        enddo
      enddo

      call mpi_barrier(mpi_comm_k,ierr)

      idum = 1
      do i = 1,nnodes
       call mpi_isend(combuf1(idum),ivunpn2(i),mpi_double_complex,i-1,
     &                inode,mpi_comm_k,ireq(i),ierr)
       idum = idum + ivunpn2(i)
      enddo

      idum = 1
      do i = 1,nnodes
       call mpi_recv(combuf2(idum),ivpacn2(i),mpi_double_complex,i-1,i,
     &               mpi_comm_k,mpistatus,ierr)
       idum = idum + ivpacn2(i)
      enddo

      do i = 1, nnodes
         call mpi_wait(ireq(i), mpistatus, ierr)
      end do

      call mpi_barrier(mpi_comm_k,ierr)

      idum = 1
      do i = 1,nnodes
        do j = 1,ivpacn2(i)
          psiy(ivpac2(idum)) = combuf2(idum)
          idum = idum + 1
        enddo
      enddo
c
c now do FFT's on the two slice
c
       do i = 1,2*ncoly
        ilocadd = 1+(i-1)*nr2
        call system_ccfft(0,isign,nr2,1.0d0,psiy(ilocadd),
     &  psiy(ilocadd),tabnr2fw,worknr2,0,ntabnr2,nworknr2)
c        call ccfft(isign,nr2,1.0,psiy(ilocadd),psiy(ilocadd),
c     c             tabnr2,worknr2,0)
c        call dcft(0, psiy(ilocadd), 1, 0,
c    &                psiy(ilocadd), 1, 0,
c    &                nr2, 1, -isign, 1.0,
c    &                tabnr2, ntabnr2, worknr2, nworknr2)
       enddo



c 
c now transpose back to format of program into psi 
c ie into x columns load balanced 
c 
      idum = 1
      do i = 1,nnodes
        do j = 1,ivunpn1(i)
          combuf1(idum) = psiy(ivunp1(idum))
          idum  = idum + 1
        enddo
      enddo

      call mpi_barrier(mpi_comm_k,ierr)

      idum = 1
      do i = 1,nnodes
       call mpi_isend(combuf1(idum),ivunpn1(i),mpi_double_complex,i-1,
     &                inode,mpi_comm_k,ireq(i),ierr)
       idum = idum + ivunpn1(i)
      enddo

      idum = 1
      do i = 1,nnodes
       call mpi_recv(combuf2(idum),ivpacn1(i),mpi_double_complex,i-1,i,
     &               mpi_comm_k,mpistatus,ierr)
       idum = idum + ivpacn1(i)
      enddo

      do i = 1, nnodes
         call mpi_wait(ireq(i), mpistatus, ierr)
      end do

      call mpi_barrier(mpi_comm_k,ierr)

      idum = 1
      do i = 1,nnodes
        do j = 1,ivpacn1(i)
          psi(ivpac1(idum)) = combuf2(idum)
          idum = idum + 1
        enddo
      enddo

c
c for psi in the x direction 
c
       do jcol = 1,ncol(inode)
        ilocadd = 1+(jcol-1)*nr1
        call system_ccfft(0,isign,nr1,1.0d0,psi(ilocadd),
     &       psi(ilocadd),tabnr1fw,worknr1,0,ntabnr1,nworknr1)
c        call ccfft(isign,nr1,1.0,psi(ilocadd),psi(ilocadd),
c     c             tabnr1,worknr1,0)
c        call dcft(0, psi(ilocadd), 1, 0,
c    &                psi(ilocadd), 1, 0,
c    &                nr1, 1, -isign, 1.0,
c    &                tabnr1, ntabnr1, worknr1, nworknr1)
       enddo !i

c
c now do the scaling of the FFT
c
      fac=1.d0/dfloat(nr1*nr2*nr3)
 
      do i = 1,ncol(inode)*nr1
        psi(i) = psi(i)*fac
      enddo

      call mpi_barrier(mpi_comm_k,ierr)

      deallocate(worknr1)
      deallocate(worknr2)
      deallocate(worknr3)
      deallocate(psiy)
      deallocate(combuf1)
      deallocate(combuf2)

      end



