c     This module cotains arrays that were present in the original
c     program in common blocks and now need to be in a module as the
c     arrays in them are dynamically allocated.
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
******************************************

      
      module data
           
      real*8,allocatable,dimension(:)        :: gkk2_n
      real*8,allocatable,dimension(:)        :: gkx2_n,gky2_n,gkz2_n
      integer,allocatable,dimension(:)  :: ig_star,ig_local,ig_lenstar
      real*8,allocatable,dimension(:,:)      :: gkk_n,wg_n
      real*8,allocatable,dimension(:,:)      :: gkx_n,gky_n,gkz_n
      real*8,allocatable,dimension(:)        :: akx,aky,akz,weighkpt
      real*8,allocatable,dimension(:)        :: vion_n,vionT_n,rhocr_n
      real*8,allocatable,dimension(:,:)        :: vr_in_n,vr_out_n
      real*8,allocatable,dimension(:,:)        :: rho_n
      real*8,allocatable,dimension(:,:,:)      :: dw,dR
      real*8,allocatable,dimension(:,:)        :: R0,w_in0
      complex*16,allocatable,dimension(:,:,:):: wqmask
      complex*16,allocatable,dimension(:,:)  :: ug_n

       
      real*8,allocatable,dimension(:)  :: wmask,xyzmap
      complex*16,allocatable,dimension(:) :: cphase
      integer,allocatable,dimension(:)  :: indm
      real*8,allocatable,dimension(:)   :: wmaskX

      
      integer,parameter :: matom=600
      integer,parameter :: nmax=20  

      integer mr_n,mg_nx,mst
      integer num_gstar, ng_ig_local
       

      contains

      subroutine data_allocate_akx(nkpt)
      implicit none
      integer nkpt
      allocate(akx(nkpt))     ! the k-points
      allocate(aky(nkpt))
      allocate(akz(nkpt))
      allocate(weighkpt(nkpt))
      return

      end subroutine data_allocate_akx
      
      subroutine data_allocate(mg_nx,mx,mr_n,ilocal,
     &  inode,nkpt,islda,natom)
      implicit none
      
      integer mg_nx,mx,mr_n,ilocal,inode,islda
      integer,parameter :: nmax=20  
      integer mr_n,mg_nx,nkpt,natom
      
      allocate(gkk_n(mg_nx,nkpt))
      allocate(wg_n(mg_nx,nkpt))
      allocate(gkx_n(mg_nx,nkpt))
      allocate(gky_n(mg_nx,nkpt))
      allocate(gkz_n(mg_nx,nkpt))

      allocate(gkk2_n(mr_n/2))
      allocate(gkx2_n(mr_n/2))
      allocate(gky2_n(mr_n/2))
      allocate(gkz2_n(mr_n/2))

      allocate(vion_n(mr_n))
      allocate(vionT_n(mr_n))
      allocate(rhocr_n(mr_n))
      allocate(vr_in_n(mr_n,islda))
      allocate(vr_out_n(mr_n,islda))
      allocate(rho_n(mr_n,islda))
      allocate(dw(mr_n,nmax,islda))    ! store them on disk later
      allocate(dR(mr_n,nmax,islda))    ! store them on disk later
      allocate(R0(mr_n,islda))
      allocate(w_in0(mr_n,islda))
      allocate(ug_n(mg_nx,mx))
cccc if real space NL is used, we don't need wqmask()
      if(ilocal.eq.3) then
      allocate(wqmask(9,mg_nx,natom))
      else
      allocate(wqmask(1,1,1))
      endif


c     allocate(pg(mg_nx))
c     allocate(xatom(3,matom))
c     allocate(dx(3,matom))
c     allocate(itype(matom))

      end subroutine data_allocate

      subroutine data_deallocate_akx()
      deallocate(akx)
      deallocate(aky)
      deallocate(akz)
      deallocate(weighkpt)
      end subroutine data_deallocate_akx
      
      subroutine data_deallocate()
      
      deallocate(gkk_n)
      deallocate(gkk2_n)
      deallocate(wg_n)
      deallocate(gkx2_n)
      deallocate(gky2_n)
      deallocate(gkz2_n)
      deallocate(gkx_n)
      deallocate(gky_n)
      deallocate(gkz_n)
      deallocate(vion_n)
      deallocate(vionT_n)
      deallocate(vr_in_n)
      deallocate(vr_out_n)
      deallocate(rho_n)
      deallocate(dw)
      deallocate(dR)
      deallocate(R0)
      deallocate(w_in0)
      deallocate(wqmask)
      deallocate(ug_n)
c     deallocate(pg)
c     deallocate(xatom)
c     deallocate(dx)
c     deallocate(itype)
      end subroutine data_deallocate

      end module data


