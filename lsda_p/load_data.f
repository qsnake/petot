      module load_data
      
ccccccccccccccccccccccccccccccccccccc
ccc Written by Andrew Canning, 2001
ccccccccccccccccccccccccccccccccccccccccccc

c     Global arrays for load balancing

      integer,allocatable,dimension(:)  :: ngcol,ngcol2,ngcol0,ngcol02
      integer,allocatable,dimension(:)  :: ngycol,ngycol2,ngzcol,ngzcol2
      integer,allocatable,dimension(:)  :: ncol,ncol2,ng,ngtotnod2
      integer,allocatable,dimension(:,:) :: ngtotnod
      integer,allocatable,dimension(:,:) :: jjcol,jjcol2,jjnode,jjnode2
      integer,allocatable,dimension(:,:) :: iycol,iycol2,izcol,izcol2
      integer,allocatable,dimension(:)   :: igstar_jjcol2,igfin_jjcol2
      integer,allocatable,dimension(:)   :: n1p_n,n2p_n,n3p_n
      integer,allocatable,dimension(:)   :: n1p2_n,n2p2_n,n3p2_n
      integer iorg(2),iorg2(2)           ! node and ig number of origin
      integer nr1x,nr2x,nr3x,ng2

cccccccccccccccccccccccccccccccccccccccccccccccccc

      integer,allocatable,dimension(:)  :: ngcol2L,ngcol02L
      integer,allocatable,dimension(:)  :: ngycol2L,ngzcol2L
      integer,allocatable,dimension(:)  :: ncol2L,ngtotnod2L
      integer,allocatable,dimension(:,:) :: jjcol2L,jjnode2L
      integer,allocatable,dimension(:,:) :: iycol2L,izcol2L
      integer,allocatable,dimension(:)   :: igstar_jjcol2L,igfin_jjcol2L
      integer,allocatable,dimension(:)   :: n1p2_nL,n2p2_nL,n3p2_nL
      integer iorg2L(2)           ! node and ig number of origin
      integer nr1xL,nr2xL,nr3xL,ng2L

ccccccccccccccccccccccccccccccccccccccccccccccccc
     
      integer,allocatable,dimension(:)  :: ngcol2L2,ngcol02L2
      integer,allocatable,dimension(:)  :: ngycol2L2,ngzcol2L2
      integer,allocatable,dimension(:)  :: ncol2L2,ngtotnod2L2
      integer,allocatable,dimension(:,:) :: jjcol2L2,jjnode2L2
      integer,allocatable,dimension(:,:) :: iycol2L2,izcol2L2
      integer,allocatable,dimension(:) :: igstar_jjcol2L2,igfin_jjcol2L2
      integer,allocatable,dimension(:)   :: n1p2_nL2,n2p2_nL2,n3p2_nL2
      integer iorg2L2(2)           ! node and ig number of origin
      integer nr1xL2,nr2xL2,nr3xL2,ng2L2

      contains

      subroutine ngtot_allocate(nnodes,nkpt)
      implicit none
      integer nnodes,nkpt
      allocate(ngtotnod(nnodes,nkpt))
      allocate(ng(nkpt))
      allocate(ngtotnod2(nnodes))
      allocate(ngtotnod2L(nnodes))
      allocate(ngtotnod2L2(nnodes))
      return
      end subroutine ngtot_allocate
      
      subroutine load_allocate(n1,n2,n3,ncolx,nnodes,
     &    mg_nx,mr_n)
      implicit none
      integer ncolx,nnodes,mg_nx,itemp,mr_n
      integer n1,n2,n3

      nr1x=n1           ! here, nr1x,nr2x,nr3x are global variables
      nr2x=n2
      nr3x=n3+2

      itemp=nr2x*(nr3x-2)
      
      allocate(ngcol(itemp))
      allocate(ngcol2(itemp))
      allocate(ngcol0(itemp))
      allocate(ngcol02(itemp))
      allocate(ngycol(itemp))
      allocate(ngycol2(itemp))
      allocate(ngzcol(itemp))
      allocate(ngzcol2(itemp))
      
      allocate(ncol(nnodes))
      allocate(ncol2(nnodes))
      allocate(jjcol(nr2x,nr3x))
      allocate(jjcol2(nr2x,nr3x))
      allocate(jjnode(nr2x,nr3x))
      allocate(jjnode2(nr2x,nr3x))
      allocate(iycol(ncolx,nnodes))
      allocate(iycol2(ncolx,nnodes))
      allocate(izcol(ncolx,nnodes))
      allocate(izcol2(ncolx,nnodes))
      allocate(igstar_jjcol2(ncolx))
      allocate(igfin_jjcol2(ncolx))

     
      allocate(n1p_n(mg_nx))
      allocate(n2p_n(mg_nx))
      allocate(n3p_n(mg_nx))
      allocate(n1p2_n(mr_n))
      allocate(n2p2_n(mr_n))
      allocate(n3p2_n(mr_n))
     
      end subroutine load_allocate

ccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine load_allocateL(n1L,n2L,n3L,ncolxL,
     &      nnodes,mr_nL,iflag_fft2L)

      implicit none
      integer ncolxL,nnodes,itemp,mr_nL
      integer n1L,n2L,n3L
      integer iflag_fft2L

      nr1xL=n1L          ! nr1xL,nr2xL,nr3xL are global variables. 
      nr2xL=n2L
      nr3xL=n3L+2

      itemp=nr2xL*(nr3xL-2)

      if(iflag_fft2L.eq.1) then

      allocate(ngcol2L(itemp))
      allocate(ngcol02L(itemp))
      allocate(ngycol2L(itemp))
      allocate(ngzcol2L(itemp))
      
      allocate(ncol2L(nnodes))
      allocate(jjcol2L(nr2xL,nr3xL))
      allocate(jjnode2L(nr2xL,nr3xL))
      allocate(iycol2L(ncolxL,nnodes))
      allocate(izcol2L(ncolxL,nnodes))
      allocate(igstar_jjcol2L(ncolxL))
      allocate(igfin_jjcol2L(ncolxL))

      allocate(n1p2_nL(mr_nL))
      allocate(n2p2_nL(mr_nL))
      allocate(n3p2_nL(mr_nL))

      endif

      if(iflag_fft2L.eq.0) then
      allocate(jjcol2L(nr2xL,nr3xL))
      allocate(jjnode2L(nr2xL,nr3xL))
      allocate(igstar_jjcol2L(ncolxL))
      allocate(igfin_jjcol2L(ncolxL))
      allocate(n1p2_nL(mr_nL))
      endif

     
      end subroutine load_allocateL

ccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine load_allocateL2(n1L2,n2L2,n3L2,
     &    ncolxL2,nnodes,mr_nL2)

      implicit none
      integer ncolxL2,nnodes,itemp,mr_nL2
      integer n1L2,n2L2,n3L2

      nr1xL2=n1L2          ! nr1xL,nr2xL,nr3xL are global variables. 
      nr2xL2=n2L2
      nr3xL2=n3L2+2

      itemp=nr2xL2*(nr3xL2-2)

      allocate(ngcol2L2(itemp))
      allocate(ngcol02L2(itemp))
      allocate(ngycol2L2(itemp))
      allocate(ngzcol2L2(itemp))
      
      allocate(ncol2L2(nnodes))
      allocate(jjcol2L2(nr2xL2,nr3xL2))
      allocate(jjnode2L2(nr2xL2,nr3xL2))
      allocate(iycol2L2(ncolxL2,nnodes))
      allocate(izcol2L2(ncolxL2,nnodes))
      allocate(igstar_jjcol2L2(ncolxL2))
      allocate(igfin_jjcol2L2(ncolxL2))

      allocate(n1p2_nL2(mr_nL2))
      allocate(n2p2_nL2(mr_nL2))
      allocate(n3p2_nL2(mr_nL2))
     
      end subroutine load_allocateL2

ccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine load_deallocateL()
      implicit none
      
      deallocate(ngcol2L)
      deallocate(ngcol02L)
      deallocate(ngycol2L)
      deallocate(ngzcol2L)
      
      deallocate(ncol2L)
      deallocate(jjcol2L)
      deallocate(jjnode2L)
      deallocate(iycol2L)
      deallocate(izcol2L)
      deallocate(igstar_jjcol2L)
      deallocate(igfin_jjcol2L)

      deallocate(n1p2_nL)
      deallocate(n2p2_nL)
      deallocate(n3p2_nL)

      end subroutine load_deallocateL
cccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine load_deallocateL2()
      implicit none
      
      deallocate(ngcol2L2)
      deallocate(ngcol02L2)
      deallocate(ngycol2L2)
      deallocate(ngzcol2L2)
      
      deallocate(ncol2L2)
      deallocate(jjcol2L2)
      deallocate(jjnode2L2)
      deallocate(iycol2L2)
      deallocate(izcol2L2)
      deallocate(igstar_jjcol2L2)
      deallocate(igfin_jjcol2L2)

      deallocate(n1p2_nL2)
      deallocate(n2p2_nL2)
      deallocate(n3p2_nL2)

      end subroutine load_deallocateL2

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
     

      
      subroutine load_deallocate()
      
      deallocate(ngcol)
      deallocate(ngcol2)
      deallocate(ngcol0)
      deallocate(ngcol02)
      deallocate(ngycol)
      deallocate(ngycol2)
      deallocate(ngzcol)
      deallocate(ngzcol2)
      deallocate(ncol)
      deallocate(ncol2)
      deallocate(ngtotnod)
      deallocate(ngtotnod2)
      deallocate(ng)
      deallocate(jjcol)
      deallocate(jjcol2)
      deallocate(jjnode)
      deallocate(jjnode2)
      deallocate(iycol)
      deallocate(iycol2)
      deallocate(izcol)
      deallocate(izcol2)
     
      deallocate(n1p_n)
      deallocate(n2p_n)
      deallocate(n3p_n)
      deallocate(n1p2_n)
      deallocate(n2p2_n)
      deallocate(n3p2_n)
     
      end subroutine load_deallocate

      end module load_data
