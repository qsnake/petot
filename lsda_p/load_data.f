      module load_data
      
c     Global arrays for load balancing

      integer,allocatable,dimension(:)  :: ngcol,ngcol2,ngcol0,ngcol02
      integer,allocatable,dimension(:)  :: ngycol,ngycol2,ngzcol,ngzcol2
      integer,allocatable,dimension(:)  :: ncol,ncol2,ng,ngtotnod2
      integer,allocatable,dimension(:,:) :: ngtotnod
      integer,allocatable,dimension(:,:) :: jjcol,jjcol2,jjnode,jjnode2
      integer,allocatable,dimension(:,:) :: iycol,iycol2,izcol,izcol2
      integer,allocatable,dimension(:)   :: igstar_jjcol2,igfin_jjcol2
     
c     Local arrays for load balancing

      integer,allocatable,dimension(:)   :: n1p_n,n2p_n,n3p_n
      integer,allocatable,dimension(:)   :: n1p2_n,n2p2_n,n3p2_n

      integer iorg(2),iorg2(2)           ! node and ig number of origin
      integer nr1x,nr2x,nr3x,ng2

      contains

      subroutine ngtot_allocate(nnodes,nkpt)
      implicit none
      integer nnodes,nkpt
      allocate(ngtotnod(nnodes,nkpt))
      allocate(ng(nkpt))
      allocate(ngtotnod2(nnodes))
      return
      end subroutine ngtot_allocate
      
      subroutine load_allocate(ncolx,nnodes,mg_nx,mr_n)
      implicit none
      integer ncolx,nnodes,mg_nx,itemp,mr_n
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
