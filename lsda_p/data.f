c     This module cotains arrays that were present in the original
c     program in common blocks and now need to be in a module as the
c     arrays in them are dynamically allocated.
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************

      
      module data
           
      real*8,allocatable,dimension(:)   :: gkk2_n
      real*8,allocatable,dimension(:)   :: gkk2_nL,gkk2_nL2
      real*8,allocatable,dimension(:)   :: gkx2_n,gky2_n,gkz2_n
      real*8,allocatable,dimension(:)   :: gkx2_nL,gky2_nL,gkz2_nL
      real*8,allocatable,dimension(:)   :: gkx2_nL2,gky2_nL2,gkz2_nL2
      integer,allocatable,dimension(:)  :: map_StoL,map_LtoS
      integer,allocatable,dimension(:)  :: ig_star,ig_local,ig_lenstar
      integer,allocatable,dimension(:)  :: ig_star_stop
      real*8,allocatable,dimension(:,:)      :: gkk_n,wg_n
      real*8,allocatable,dimension(:,:)      :: gkx_n,gky_n,gkz_n
      real*8,allocatable,dimension(:)        :: akx,aky,akz,weighkpt
      real*8,allocatable,dimension(:)     :: vion_nL,vionT_nL,rhocr_nL
      real*8,allocatable,dimension(:)     :: vext_nL
      real*8,allocatable,dimension(:,:)   :: vr_in_nL,vr_out_nL,vr_n
      real*8,allocatable,dimension(:,:)   :: vcoul_nL2
      real*8,allocatable,dimension(:,:)        :: rho_n,rho_nL,vtot_nL
      real*8,allocatable,dimension(:,:,:)      :: dw,dR
      real*8,allocatable,dimension(:,:)        :: R0,w_in0
      complex*16,allocatable,dimension(:,:)    :: wqmask
      complex*16,allocatable,dimension(:,:)    :: beta_psi
      complex*16,allocatable,dimension(:,:)  :: ug_n,sug_n

       
      real*8,allocatable,dimension(:)  :: wmask,xyzmap
      complex*16,allocatable,dimension(:) :: cphase
      real*8,allocatable,dimension(:,:) :: YYMask
      integer,allocatable,dimension(:)  :: indm
      real*8,allocatable,dimension(:,:)   :: wmaskX

      real*8,allocatable,dimension(:,:)  :: wmask_q,wmask_q0
      real*8,allocatable,dimension(:,:,:)  :: wmask_dq,wmask_dq0
      integer,allocatable,dimension(:)   :: indm_q
      real*8,allocatable,dimension(:,:,:,:)  :: Dij

      
      integer,parameter :: matom=2000

      integer npulay_max
      integer mr_n,mg_nx,mst
      integer mr_nL,mr_nL2
      integer num_gstar, ng_ig_local
      integer nref_tot
      integer iref_start(matom)
       
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

      subroutine data_allocate_YYMAsk(mm,nkk)
      implicit none
      integer mm,nkk
      allocate(YYMAsk(mm,nkk))
      return
      end subroutine data_allocate_YYMask

      subroutine data_deallocate_YYMAsk()
      implicit none
      deallocate(YYMAsk)
      return
      end subroutine data_deallocate_YYMask
 


      subroutine data_allocate_nL2(mr_nL2,ntype)
      implicit none
      integer mr_nL2,ntype
      allocate(gkk2_nL2(mr_nL2/2))
      allocate(gkx2_nL2(mr_nL2/2))
      allocate(gky2_nL2(mr_nL2/2))
      allocate(gkz2_nL2(mr_nL2/2))
      allocate(vcoul_nL2(mr_nL2,0:ntype))
      return
      end subroutine data_allocate_nL2

      
      subroutine data_allocate(mg_nx,mx,mr_n,mr_nL,
     &  ilocal,inode,nkpt,islda,natom,nref_tot,
     & ipsp_all,nnodes,npulay_max)
      implicit none
      
      integer mg_nx,mx,ilocal,inode,islda,ipsp_all,nnodes,
     & mr_nL,npulay_max
      integer mr_n,nkpt,natom,nref_tot
      
      allocate(gkk_n(mg_nx,nkpt))
      allocate(wg_n(mg_nx,nkpt))
      allocate(gkx_n(mg_nx,nkpt))
      allocate(gky_n(mg_nx,nkpt))
      allocate(gkz_n(mg_nx,nkpt))

      allocate(gkk2_n(mr_n/2))
      allocate(gkx2_n(mr_n/2))
      allocate(gky2_n(mr_n/2))
      allocate(gkz2_n(mr_n/2))

      allocate(vion_nL(mr_nL))
      allocate(vext_nL(mr_nL))
      allocate(vionT_nL(mr_nL))
      allocate(rhocr_nL(mr_nL))
      allocate(vr_in_nL(mr_nL,islda))
      allocate(vr_out_nL(mr_nL,islda))
      allocate(vr_n(mr_n,islda))
      allocate(rho_n(mr_n,islda))
      allocate(rho_nL(mr_nL,islda))
      allocate(vtot_nL(mr_nL,islda))
      allocate(dw(mr_nL,npulay_max,islda))    ! store them on disk later
      allocate(dR(mr_nL,npulay_max,islda))    ! store them on disk later
      allocate(R0(mr_nL,islda))
      allocate(w_in0(mr_nL,islda))
      allocate(ug_n(mg_nx,mx))
      if(ipsp_all.eq.1) then    
      allocate(sug_n(1,mx))            ! dummy
      else
      allocate(sug_n(mg_nx,mx))        ! used for orthogonization, sug_n=S*ug_n
      endif
      allocate(gkk2_nL(mr_nL/2))
      allocate(gkx2_nL(mr_nL/2))
      allocate(gky2_nL(mr_nL/2))
      allocate(gkz2_nL(mr_nL/2))
      allocate(map_StoL(mr_nL/2))
      allocate(map_LtoS(mr_nL/2))
cccc if real space NL is used, we don't need wqmask()
      if(ilocal.eq.3) then
      allocate(wqmask(mg_nx,nref_tot))
      else
      allocate(wqmask(1,1))
      endif
      allocate(beta_psi(nref_tot,mx/nnodes+1))
      allocate(Dij(32,32,natom,islda))


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
      deallocate(vion_nL)
      deallocate(vext_nL)
      deallocate(vionT_nL)
      deallocate(rhocr_nL)
      deallocate(vr_in_nL)
      deallocate(vr_out_nL)
      deallocate(vr_n)
      deallocate(rho_n)
      deallocate(rho_nL)
      deallocate(vtot_nL)
      deallocate(dw)
      deallocate(dR)
      deallocate(R0)
      deallocate(w_in0)
      deallocate(wqmask)
      deallocate(beta_psi)
      deallocate(ug_n)
c     deallocate(pg)
c     deallocate(xatom)
c     deallocate(dx)
c     deallocate(itype)
      deallocate(Dij)
      end subroutine data_deallocate

      subroutine data_deallocate_nL2()
      implicit none
      deallocate(gkk2_nL2)
      deallocate(gkx2_nL2)
      deallocate(gky2_nL2)
      deallocate(gkz2_nL2)
      deallocate(vcoul_nL2)
      return
      end subroutine data_deallocate_nL2

      end module data


