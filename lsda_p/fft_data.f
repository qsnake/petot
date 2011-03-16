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
     
      contains
      
      subroutine fft_allocate(n1,n2,n3,nnodes)
      implicit none
      
      integer n1,n2,n3,nr1x,nr2x,nr3x,nnodes

      nr1x=n1
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


      subroutine fft_allocate2(k0nunp_max,nnodes)
      implicit none
      integer k0nunp_max,nnodes
      allocate(k0iunp(k0nunp_max,nnodes)) 
      return
      end subroutine fft_allocate2

      
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

      end module fft_data
