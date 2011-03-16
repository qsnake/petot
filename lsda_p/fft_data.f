********************************************************************
**   Written by Andrew Canning, 2001
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

c     Data for the Parallel FFT's
      
c     ngcol: number of g vectors per column
c     ngcol0: as ngcol but including empty cols on k=0 plane requied in fft
c     ngycol: y position of column in grid
c     ngzcol: z position of column in grid
c     jjcol: column number on PE of column (y,z)
c     jjnode: node number of column (y,z)
c     ncol: num of columns on each node
c     ngtotnod: tot number of g vecs on each node
c     n1p_n,n2p_n,n3p_n: i,j,k of g vector in n1*n2*n3 grid (local)
c     iycol: y position of column iycol(ncol,node) ncol= local column number
c     izcol: z position of column izcol(ncol,node) ncol= local column number

      module fft_data

      integer,allocatable,dimension(:,:) :: izcol_y,ixcol_y,ipe_y
      integer,allocatable,dimension(:,:) :: izcol_y2,ixcol_y2,ipe_y2
      integer,allocatable,dimension(:,:) :: ibca_y,iycol_z,ixcol_z
      integer,allocatable,dimension(:,:) :: ibca_y2,iycol_z2,ixcol_z2
      integer,allocatable,dimension(:,:) :: ipe_z,ibca_z,ivecadd
      integer,allocatable,dimension(:,:) :: ipe_z2,ibca_z2,ivecadd2
      integer,allocatable,dimension(:)   :: izchb,ixch,ichw,icount
      integer,allocatable,dimension(:)   :: izchb2,ixch2,ichw2,icount2

      real*8,allocatable,dimension(:)   :: tabnr1fw,tabnr2fw,tabnr3fw
      real*8,allocatable,dimension(:)   :: tabnr1in,tabnr2in,tabnr3in
      real*8,allocatable,dimension(:)   :: tabnr1lfw,tabnr2lfw
      real*8,allocatable,dimension(:)   :: tabnr1lin,tabnr2lin,
     &   tabnr3lrc,tabnr3lcr

      real*8,allocatable,dimension(:,:) :: vec,vec2

      integer,allocatable,dimension(:) :: ivpacn1,ivunpn1,ivpacn2,
     &                                    ivunpn2
      integer,allocatable,dimension(:) :: ivpac1,ivunp1,ivpac2,ivunp2

! amc new
      integer,allocatable,dimension(:) :: ivpacn1l,ivunpn1l,ivpacn2l,
     &                                    ivunpn2l
      integer,allocatable,dimension(:) :: ivpac1l,ivunp1l,ivpac2l,
     &                                    ivunp2l
      integer,allocatable,dimension(:) :: k0npac,k0nunp
      integer,allocatable,dimension(:,:) :: k0iunp

! amc new

      
      integer ncolx,ncoly,ncolz,ichunk,mnrx,mgz,ntabnr1,ntabnr2,ntabnr3
      integer ncolx2,ncoly2,ncolz2,ichunk2,mnr2x,mgz2
     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccc the tables for the Large real fft
      integer,allocatable,dimension(:,:) :: izcol_y2L,ixcol_y2L,ipe_y2L
      integer,allocatable,dimension(:,:) :: ibca_y2L,iycol_z2L,ixcol_z2L
      integer,allocatable,dimension(:,:) :: ipe_z2L,ibca_z2L,ivecadd2L
      integer,allocatable,dimension(:)   :: izchb2L,ixch2L,ichw2L,
     &     icount2L

      real*8,allocatable,dimension(:)   :: tabnr1lfwL,tabnr2lfwL
      real*8,allocatable,dimension(:)   :: tabnr1linL,tabnr2linL,
     &   tabnr3lrcL,tabnr3lcrL

      real*8,allocatable,dimension(:,:) :: vec2L

! amc new
      integer,allocatable,dimension(:) :: ivpacn1lL,ivunpn1lL,ivpacn2lL,
     &                                    ivunpn2lL
      integer,allocatable,dimension(:) :: ivpac1lL,ivunp1lL,ivpac2lL,
     &                                    ivunp2lL
      integer,allocatable,dimension(:) :: k0npacL,k0nunpL
      integer,allocatable,dimension(:,:) :: k0iunpL

! amc new

      integer ntabnr1L,ntabnr2L,ntabnr3L
      integer ncolx2L,ncoly2L,ncolz2L,ichunk2L,mnr2xL,mgz2L
      integer ncolxL,ncolzL

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccc the tables for the Large real fft
      integer,allocatable,dimension(:,:) :: 
     &                 izcol_y2L2,ixcol_y2L2,ipe_y2L2
      integer,allocatable,dimension(:,:) :: 
     &                 ibca_y2L2,iycol_z2L2,ixcol_z2L2
      integer,allocatable,dimension(:,:) :: 
     &                 ipe_z2L2,ibca_z2L2,ivecadd2L2
      integer,allocatable,dimension(:)   :: izchb2L2,ixch2L2,ichw2L2,
     &     icount2L2

      real*8,allocatable,dimension(:)   :: tabnr1lfwL2,tabnr2lfwL2
      real*8,allocatable,dimension(:)   :: tabnr1linL2,tabnr2linL2,
     &   tabnr3lrcL2,tabnr3lcrL2

      real*8,allocatable,dimension(:,:) :: vec2L2

! amc new
      integer,allocatable,dimension(:) :: ivpacn1lL2,ivunpn1lL2,
     &                                 ivpacn2lL2,ivunpn2lL2
      integer,allocatable,dimension(:) :: ivpac1lL2,ivunp1lL2,
     &                                ivpac2lL2,ivunp2lL2
      integer,allocatable,dimension(:) :: k0npacL2,k0nunpL2
      integer,allocatable,dimension(:,:) :: k0iunpL2

! amc new

      integer ntabnr1L2,ntabnr2L2,ntabnr3L2
      integer ncolx2L2,ncoly2L2,ncolz2L2,ichunk2L2,mnr2xL2,mgz2L2
      integer ncolxL2,ncolzL2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     
      contains
      
      subroutine fft_allocate(n1,n2,n3,nnodes)
      implicit none
      
      integer n1,n2,n3,nr1x,nr2x,nr3x,nnodes

      nr1x=n1            ! here, nr1x,nr2x,nr3x are just local variables
      nr2x=n2
      nr3x=n3+2
      ntabnr1=20000+2.28*nr1x
      ntabnr2=20000+2.28*nr2x
      ntabnr3=20000+2.28*nr3x

      mnrx=nr1x*nr2x*nr3x/8
      mnr2x=nr1x*nr2x*nr3x

cccccccc I don't think the mnrx, and mnr2x are correct ! 
cccccccc How can they be divided by nnodes again !! Needs to be clear out !

      mnrx=8*(mnrx/nnodes+50)      ! added by Lin-Wang
      mnr2x=8*(mnr2x/nnodes+50)    ! the factor of 2 might be reduced to 1

      ncolx= nr2x*nr3x/nnodes
      ncolz=n1*n2/nnodes
      
c     Allocate the integer arrays
      
      allocate(izcol_y((nr1x*nr3x)/nnodes,nnodes))
      allocate(izcol_y2((nr1x*nr3x)/nnodes,nnodes))
      allocate(ixcol_y((nr1x*nr3x)/nnodes,nnodes))
      allocate(ixcol_y2((nr1x*nr3x)/nnodes,nnodes))
      allocate(ipe_y(nr1x,nr3x))
      allocate(ipe_y2(nr1x,nr3x))
      allocate(ibca_y(nr1x,nr3x))
      allocate(ibca_y2(nr1x,nr3x))
      allocate(iycol_z((nr1x*nr2x)/nnodes,nnodes))
      allocate(iycol_z2((nr1x*nr2x)/nnodes,nnodes))
      allocate(ixcol_z((nr1x*nr2x)/nnodes,nnodes))
      allocate(ixcol_z2((nr1x*nr2x)/nnodes,nnodes))
      allocate(ipe_z(nr2x,nr1x))
      allocate(ipe_z2(nr2x,nr1x))
      allocate(ibca_z(nr2x,nr1x))
      allocate(ibca_z2(nr2x,nr1x))
      allocate(ivecadd(mnrx/nnodes,nnodes))
      allocate(ivecadd2(mnr2x/nnodes,nnodes))
      allocate(icount(nnodes))
      allocate(icount2(nnodes))


c Arrays required for MPI version

      allocate(ivpacn1(nnodes))
      allocate(ivunpn1(nnodes))
      allocate(ivpacn2(nnodes))
      allocate(ivunpn2(nnodes))
      allocate(ivpac1(mnrx))
      allocate(ivunp1(mnrx))
      allocate(ivpac2(mnrx))
      allocate(ivunp2(mnrx))

! amc new
      allocate(ivpacn1l(nnodes))
      allocate(ivunpn1l(nnodes))
      allocate(ivpacn2l(nnodes))
      allocate(ivunpn2l(nnodes))
      allocate(ivpac1l(mnr2x))
      allocate(ivunp1l(mnr2x))
      allocate(ivpac2l(mnr2x))
      allocate(ivunp2l(mnr2x))

      allocate(k0npac(nnodes))
      allocate(k0nunp(nnodes))
cccc      allocate(k0iunp(nr1x*nr2x,nnodes))  ! allocate seperately, too large
!amc new stuff

      
c Arrays required for MPI version
c     Now the real*8 arrays
      
      allocate(tabnr1fw(ntabnr1))
      allocate(tabnr2fw(ntabnr2))
      allocate(tabnr3fw(ntabnr3))

      allocate(tabnr1in(ntabnr1))
      allocate(tabnr2in(ntabnr2))
      allocate(tabnr3in(ntabnr3))

      allocate(tabnr1lfw(ntabnr1))
      allocate(tabnr2lfw(ntabnr2))
      allocate(tabnr1lin(ntabnr1))
      allocate(tabnr2lin(ntabnr2))


      allocate(tabnr3lrc(ntabnr3))
      allocate(tabnr3lcr(ntabnr3))


      allocate(vec(mnrx/(nnodes),nnodes))
      allocate(vec2(mnr2x/(nnodes),nnodes))
            
      end subroutine fft_allocate

cccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine fft_allocate2(k0nunp_max,nnodes)
      implicit none
      integer k0nunp_max,nnodes
      allocate(k0iunp(k0nunp_max,nnodes)) 
      return
      end subroutine fft_allocate2



      subroutine fft_allocateL(n1L,n2L,n3L,nnodes,iflag_fft2L)
      implicit none
      
      integer nr1xL,nr2xL,nr3xL,n1L,n2L,n3L,nnodes
      integer iflag_fft2L

      nr1xL=n1L
      nr2xL=n2L
      nr3xL=n3L+2
      ntabnr1L=20000+2.28*nr1xL
      ntabnr2L=20000+2.28*nr2xL
      ntabnr3L=20000+2.28*nr3xL

      mnr2xL=nr1xL*nr2xL*nr3xL

cccccccc I don't think the mnrx, and mnr2x are correct ! 
cccccccc How can they be divided by nnodes again !! Needs to be clear out !

      mnr2xL=8*(mnr2xL/nnodes+50)    ! the factor of 2 might be reduced to 1

      ncolxL= nr2xL*nr3xL/nnodes
      ncolzL=n1L*n2L/nnodes
      
c     Allocate the integer arrays

      if(iflag_fft2L.eq.0) then
      allocate(iycol_z2L((nr1xL*nr2xL)/nnodes,nnodes))   ! needed for real space manipulate
      allocate(ixcol_z2L((nr1xL*nr2xL)/nnodes,nnodes))
      endif

      
      if(iflag_fft2L.eq.1) then
      allocate(izcol_y2L((nr1xL*nr3xL)/nnodes,nnodes))
      allocate(ixcol_y2L((nr1xL*nr3xL)/nnodes,nnodes))
      allocate(ipe_y2L(nr1xL,nr3xL))
      allocate(ibca_y2L(nr1xL,nr3xL))
      allocate(iycol_z2L((nr1xL*nr2xL)/nnodes,nnodes))
      allocate(ixcol_z2L((nr1xL*nr2xL)/nnodes,nnodes))
      allocate(ipe_z2L(nr2xL,nr1xL))
      allocate(ibca_z2L(nr2xL,nr1xL))
      allocate(ivecadd2L(mnr2xL/nnodes,nnodes))
      allocate(icount2L(nnodes))


c Arrays required for MPI version

! amc new
      allocate(ivpacn1lL(nnodes))
      allocate(ivunpn1lL(nnodes))
      allocate(ivpacn2lL(nnodes))
      allocate(ivunpn2lL(nnodes))
      allocate(ivpac1lL(mnr2xL))
      allocate(ivunp1lL(mnr2xL))
      allocate(ivpac2lL(mnr2xL))
      allocate(ivunp2lL(mnr2xL))

      allocate(k0npacL(nnodes))
      allocate(k0nunpL(nnodes))
cccc      allocate(k0iunpL(nr1xL*nr2xL,nnodes))  ! allocate seperately, too large
!amc new stuff

      
c Arrays required for MPI version
c     Now the real*8 arrays
      
      allocate(tabnr1lfwL(ntabnr1L))
      allocate(tabnr2lfwL(ntabnr2L))
      allocate(tabnr1linL(ntabnr1L))
      allocate(tabnr2linL(ntabnr2L))


      allocate(tabnr3lrcL(ntabnr3L))
      allocate(tabnr3lcrL(ntabnr3L))


      allocate(vec2L(mnr2xL/(nnodes),nnodes))

      endif
            
      end subroutine fft_allocateL


      subroutine fft_allocate2L(k0nunp_maxL,nnodes)
      implicit none
      integer k0nunp_maxL,nnodes
      allocate(k0iunpL(k0nunp_maxL,nnodes)) 
      return
      end subroutine fft_allocate2L
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine fft_allocateL2(n1L2,n2L2,n3L2,nnodes)
      implicit none
      
      integer nr1xL2,nr2xL2,nr3xL2,n1L2,n2L2,n3L2,nnodes

      nr1xL2=n1L2
      nr2xL2=n2L2
      nr3xL2=n3L2+2
      ntabnr1L2=20000+2.28*nr1xL2
      ntabnr2L2=20000+2.28*nr2xL2
      ntabnr3L2=20000+2.28*nr3xL2

      mnr2xL2=nr1xL2*nr2xL2*nr3xL2

cccccccc I don't think the mnrx, and mnr2x are correct ! 
cccccccc How can they be divided by nnodes again !! Needs to be clear out !

      mnr2xL2=8*(mnr2xL2/nnodes+50)    ! the factor of 2 might be reduced to 1

      ncolxL2= nr2xL2*nr3xL2/nnodes
      ncolzL2=n1L2*n2L2/nnodes
      
c     Allocate the integer arrays
      
      allocate(izcol_y2L2((nr1xL2*nr3xL2)/nnodes,nnodes))
      allocate(ixcol_y2L2((nr1xL2*nr3xL2)/nnodes,nnodes))
      allocate(ipe_y2L2(nr1xL2,nr3xL2))
      allocate(ibca_y2L2(nr1xL2,nr3xL2))
      allocate(iycol_z2L2((nr1xL2*nr2xL2)/nnodes,nnodes))
      allocate(ixcol_z2L2((nr1xL2*nr2xL2)/nnodes,nnodes))
      allocate(ipe_z2L2(nr2xL2,nr1xL2))
      allocate(ibca_z2L2(nr2xL2,nr1xL2))
      allocate(ivecadd2L2(mnr2xL2/nnodes,nnodes))
      allocate(icount2L2(nnodes))


c Arrays required for MPI version

! amc new
      allocate(ivpacn1lL2(nnodes))
      allocate(ivunpn1lL2(nnodes))
      allocate(ivpacn2lL2(nnodes))
      allocate(ivunpn2lL2(nnodes))
      allocate(ivpac1lL2(mnr2xL2))
      allocate(ivunp1lL2(mnr2xL2))
      allocate(ivpac2lL2(mnr2xL2))
      allocate(ivunp2lL2(mnr2xL2))

      allocate(k0npacL2(nnodes))
      allocate(k0nunpL2(nnodes))
cccc      allocate(k0iunpL2(nr1xL2*nr2xL2,nnodes))  ! allocate seperately, too large
!amc new stuff

      
c Arrays required for MPI version
c     Now the real*8 arrays
      
      allocate(tabnr1lfwL2(ntabnr1L2))
      allocate(tabnr2lfwL2(ntabnr2L2))
      allocate(tabnr1linL2(ntabnr1L2))
      allocate(tabnr2linL2(ntabnr2L2))


      allocate(tabnr3lrcL2(ntabnr3L2))
      allocate(tabnr3lcrL2(ntabnr3L2))


      allocate(vec2L2(mnr2xL2/(nnodes),nnodes))
            
      end subroutine fft_allocateL2


      subroutine fft_allocate2L2(k0nunp_maxL2,nnodes)
      implicit none
      integer k0nunp_maxL2,nnodes
      allocate(k0iunpL2(k0nunp_maxL2,nnodes)) 
      return
      end subroutine fft_allocate2L2



cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccc




      
      subroutine fft_deallocate()
          deallocate(izcol_y)
          deallocate(izcol_y2)
          deallocate(ixcol_y)
          deallocate(ixcol_y2)
          deallocate(ipe_y)
          deallocate(ipe_y2)
          deallocate(ibca_y)
          deallocate(ibca_y2)
          deallocate(iycol_z)
          deallocate(iycol_z2)
          deallocate(ixcol_z)
          deallocate(ixcol_z2)
          deallocate(ipe_z)
          deallocate(ipe_z2)
          deallocate(ibca_z)
          deallocate(ibca_z2)
          deallocate(izchb)
          deallocate(izchb2)
          deallocate(ixch)
          deallocate(ixch2)
          deallocate(ichw)
          deallocate(ichw2)
          deallocate(ivecadd)
          deallocate(ivecadd2)
          deallocate(icount)
          deallocate(icount2)
          
          deallocate(tabnr1fw)
          deallocate(tabnr2fw)
          deallocate(tabnr3fw)
          deallocate(tabnr1in)
          deallocate(tabnr2in)
          deallocate(tabnr3in)

          deallocate(tabnr1lfw)
          deallocate(tabnr2lfw)
          deallocate(tabnr1lin)
          deallocate(tabnr2lin)

          deallocate(tabnr3lrc)
          deallocate(tabnr3lcr)
          deallocate(vec)
          deallocate(vec2)

          deallocate(ivpacn1)
          deallocate(ivunpn1)
          deallocate(ivpacn2)
          deallocate(ivunpn2)
          deallocate(ivpac1)
          deallocate(ivunp1)
          deallocate(ivpac2)
          deallocate(ivunp2)
          deallocate(ivpacn1l)
          deallocate(ivunpn1l)
          deallocate(ivpacn2l)
          deallocate(ivunpn2l)
          deallocate(ivpac1l)
          deallocate(ivunp1l)
          deallocate(ivpac2l)
          deallocate(ivunp2l)

          deallocate(k0npac)
          deallocate(k0nunp)
          deallocate(k0iunp)

      end subroutine fft_deallocate

cccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccc

      subroutine fft_deallocateL()
      implicit none

      deallocate(izcol_y2L)
      deallocate(ixcol_y2L)
      deallocate(ipe_y2L)
      deallocate(ibca_y2L)
      deallocate(iycol_z2L)
      deallocate(ixcol_z2L)
      deallocate(ipe_z2L)
      deallocate(ibca_z2L)
      deallocate(ivecadd2L)
      deallocate(icount2L)

      deallocate(ivpacn1lL)
      deallocate(ivunpn1lL)
      deallocate(ivpacn2lL)
      deallocate(ivunpn2lL)
      deallocate(ivpac1lL)
      deallocate(ivunp1lL)
      deallocate(ivpac2lL)
      deallocate(ivunp2lL)

      deallocate(k0npacL)
      deallocate(k0nunpL)
      deallocate(tabnr1lfwL)
      deallocate(tabnr2lfwL)
      deallocate(tabnr1linL)
      deallocate(tabnr2linL)

      deallocate(tabnr3lrcL)
      deallocate(tabnr3lcrL)

      deallocate(vec2L)
      deallocate(k0iunpL)
            
      end subroutine fft_deallocateL

ccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccc

      subroutine fft_deallocateL2()
      implicit none

      deallocate(izcol_y2L2)
      deallocate(ixcol_y2L2)
      deallocate(ipe_y2L2)
      deallocate(ibca_y2L2)
      deallocate(iycol_z2L2)
      deallocate(ixcol_z2L2)
      deallocate(ipe_z2L2)
      deallocate(ibca_z2L2)
      deallocate(ivecadd2L2)
      deallocate(icount2L2)

      deallocate(ivpacn1lL2)
      deallocate(ivunpn1lL2)
      deallocate(ivpacn2lL2)
      deallocate(ivunpn2lL2)
      deallocate(ivpac1lL2)
      deallocate(ivunp1lL2)
      deallocate(ivpac2lL2)
      deallocate(ivunp2lL2)

      deallocate(k0npacL2)
      deallocate(k0nunpL2)
      deallocate(tabnr1lfwL2)
      deallocate(tabnr2lfwL2)
      deallocate(tabnr1linL2)
      deallocate(tabnr2linL2)

      deallocate(tabnr3lrcL2)
      deallocate(tabnr3lcrL2)

      deallocate(vec2L2)
      deallocate(k0iunpL2)
            
      end subroutine fft_deallocateL2


      end module fft_data
