      subroutine Etotcalc(xatom,fatom,workr_n,E_tot,
     & iforce_cal,ido_rho,ido_vr,tolug,tolE,niter,nline,
     &  iCGmth,iscfmth,FermidE,itypeFermi,mCGbad,E_st,err_st,AL,
     &  nkpt,ntype,convergE,islda,igga,iwg_out,fwg_out,
     &  ivr_out)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************

******************************************
cccccc this program reads the initial wavefunction from  ugIO(ug_n,kpt,2,0,iislda)
cccccc and write the converged wavefunction to           ugIO(ug_n,kpt,1,0,iislda)
******************************************

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'mpif.h'

      include 'param.escan_real'
******************************************
       real*8 AL(3,3),ALt(3,3)
********************************************
       real*8 E_st(mst,nkpt,islda),err_st(mst,nkpt,islda) 
       real*8 occ(mst,nkpt,islda)
       integer icoul
       real*8 xcoul(3)
       character*20 file_tmp,fwg_out(2)
      


       real*8,allocatable,dimension(:)  :: xyzmaptmp
       real*8,allocatable,dimension(:,:)  :: wmasktmp
       real*8,allocatable,dimension(:,:,:)  :: wmaskXtmp
       integer,allocatable,dimension(:) :: indmtmp
       integer,allocatable,dimension(:) :: iatsum_all,iatsum_all2
       integer,allocatable,dimension(:) :: nmap_tmp
       real*8,allocatable,dimension(:,:) :: fatom_tmp,fatom_tmp2 


       real*8 xatom(3,matom)
       real*8 fatom(3,matom)

       real*8 occ_t(mtype)
       integer iiatom(mtype),icore(mtype),numref(matom)
       integer is_ref(mtype),ip_ref(mtype),id_ref(mtype)

       integer nmap_q(matom)
       real*8 rcut_q1(mtype),rcut_q2(mtype),rcut_qm

       integer nmap(matom)

       real*8 zatom(mtype)

       integer iCGmth(100),iscfmth(100)
       real*8 FermidE(100)
       integer itypeFermi(100)
cccccccccccccccccccccccc

       integer iatom(matom),ityatom(matom)
       real*8  totNel
       integer smatr(3,3,48),nrot,ilocal
       real*8 Ealpha(mtype)
    
       real*8  Dij0(32,32,mtype),Qij(32,32,mtype)
       integer isNLa(9,matom),ipsp_type(mtype),ipsp_all
*************************************************
       character*20 vwr_atom(mtype),file_test
       character*20 fdens_in,fdens_out
*************************************************

       common /comNL2/occ_t,iiatom,icore,numref,ityatom
       common /comispd_ref/is_ref,ip_ref,id_ref
       common /comEk/Ek
       common /comzatom/zatom
       common /comnmap/nmap
       common /com_rcut_q/rcut_q1,rcut_q2,rcut_qm
       common /comisNLa/isNLa,Dij0,Qij,ipsp_all,ipsp_type
       common /comcoul/icoul,xcoul
       common /comVext/ivext_in

       common /comMainEtot/iatom,totNel,
     &     smatr,nrot,ilocal,Ealpha,Ealphat
c

c saved arrays for mch_pulay
c
       real*8,allocatable,dimension (:,:) :: AA

       complex*16 workr_n(mr_n)

**************************************************
c initialize mpi and number each node


       allocate(AA(npulay_max,npulay_max))
c


       pi=4.0d0*datan(1.0d0)
       nr_nL=n1L*n2L*n3L/nnodes

       E_NSC0 = 0.0d0  !amc
       TS0= 0.0d0
       n33 = 0 
       E_COR0 = 0.0d0 
       E_HXC0 = 0.0d0
       E_TOT0 = 0.0d0
       DVE0 = 0.0d0
       E_IVext=0.d0
       E_rhoVext=0.d0

       Ek=0.5d0
       occ=1.d0


cc------------------------------------
       if(ivr_out.ne.2) then
       if(icoul.eq.0) then
       call getewald(fatom,xatom,AL,ityatom,ewald)
       endif

       if(icoul.eq.1) then
       call getewald3D(fatom,xatom,AL,ityatom,ewald)
       endif

       if(icoul.eq.11.or.icoul.eq.12.or.icoul.eq.13) then
       call getewald2D(fatom,xatom,AL,ityatom,ewald)
       endif

       if(ivext_in.eq.1) then
       call getEextV(fatom,xatom,AL,ityatom,E_IVext)
       endif
       endif  ! ivr_out.ne.2
cccccccccccccccccccccccccccccccccccccccccccccccccc
       if(iforce_cal.eq.0) fatom=0.d0

       call getVrhoL(AL,vion_nL,vionT_nL,xatom,ntype,iatom,
     &  rhocr_nL,totNel,ido_rho,islda)

        if(ivext_in.eq.1) then
        vion_nL=vion_nL+vext_nL
        endif


        if(ivr_out.eq.2) goto 911     ! jump out the wq,or mask_wr preparation
*********************************************************
***  ilocal.eq.3, q space nolocal Kleimen_Bylander
*********************************************************


       if(ilocal.eq.3) then
       do kpt=1,nkpt
       if(kpt.ge.kpt_dis(1).and.kpt.le.kpt_dis(2)) then
       call getwq(AL,ntype,iatom,xatom,kpt)
       call wqIO(nkpt,kpt,1)     ! 1: write, 2: read
       endif
       enddo
       endif

*********************************************************
***  ilocal.eq.2, r space nolocal Kleimen_Bylander
*********************************************************
       if(ilocal.eq.2) then

ccccccc this formula sometime is not right, if many atoms 
ccccccc are located inside one processor

       mrb2=2*(4*pi/3*rcut**3)/(vol/(n1*n2*n3))

       mrb2_matom_node=mrb2*(natom*1.d0/nnodes+2.d0)

       allocate(wmask(20*mrb2_matom_node))
       allocate(xyzmap(3*mrb2_matom_node))
       allocate(cphase(mrb2_matom_node))
       allocate(indm(mrb2_matom_node))
       allocate(xyzmaptmp(3*mrb2))
       allocate(indmtmp(mrb2))

       allocate(wmasktmp(20,mrb2))
       nmap=0

       iatsum = 0
       iatsum2=0 
       do 20 ia=natom_dis(1),natom_dis(2)

       iitype=ityatom(ia)

       inew=0
       iend=0
       if(ia.eq.natom_dis(1)) then
       inew=1
       else
       if(iitype.ne.ityatom(ia-1)) inew=1
       endif
       if(ia.eq.natom_dis(2)) then
       iend=1
       else
       if(iitype.ne.ityatom(ia+1)) iend=1
       endif

       call getwmask(xatom(1,ia),nmap(ia),indmtmp,
     &  ityatom(ia),wmasktmp,xyzmaptmp,AL,workr_n,mrb2,nref,
     &  inew,iend)

        if(numref(ia).ne.nref) then
        write(6,*) "numref(ia).ne.nref after gewmask,stop", 
     &   numref(ia), nref
        stop
        endif
       if(iatsum+nmap(ia).gt.mrb2_matom_node) then
       write(6,*) "iatsum.gt.mrb2_matom_node, stop"
       call mpi_abort(MPI_COMM_WORLD,ierr)
       endif

       i00=3*iatsum
         do i=1,nmap(ia)
         do j=1,nref
         iatsum2=iatsum2+1
         wmask(iatsum2)=wmasktmp(j,i)
         enddo
         enddo
         do i=1,nmap(ia)
         indm(i+iatsum)=indmtmp(i)
         enddo
         do i=1,3*nmap(ia)
         xyzmap(i+i00)=xyzmaptmp(i)
         enddo
       iatsum = iatsum + nmap(ia)
       
20     continue
       deallocate(wmasktmp)
ccccccccccccccccccccccccccccccccccccccc
       call mpi_barrier(MPI_COMM_WORLD,ierr)

       allocate(nmap_tmp(natom))
       call mpi_allreduce(nmap,nmap_tmp,natom,        ! this is used to pass nmap to other icolor
     & MPI_INTEGER,MPI_SUM,MPI_COMM_N,ierr)
       nmap(1:natom)=nmap_tmp(1:natom)
       deallocate(nmap_tmp)
       

       allocate(iatsum_all(num_group)) 
       allocate(iatsum_all2(num_group)) 
        iatsum_all(icolor+1)=iatsum
        iatsum_all2(icolor+1)=iatsum2
       do ii=1,num_group
       call mpi_bcast(iatsum_all(ii),1,MPI_INTEGER,
     &   ii-1,MPI_COMM_N,ierr)
       call mpi_bcast(iatsum_all2(ii),1,MPI_INTEGER,
     &   ii-1,MPI_COMM_N,ierr)
       enddo

       iatsum_posit=0
       iatsum_posit2=0
       do ii=1,num_group
        if(icolor.eq.ii-1) then
        do i=iatsum_all2(ii),1,-1
        wmask(i+iatsum_posit2)=wmask(i)       ! shift wmask into the right position
        enddo
        do i=iatsum_all(ii),1,-1
        indm(i+iatsum_posit)=indm(i)
        enddo
        do i=3*iatsum_all(ii),1,-1
        xyzmap(i+iatsum_posit*3)=xyzmap(i)
        enddo
        endif
        iatsum_posit=iatsum_posit+iatsum_all(ii)
        iatsum_posit2=iatsum_posit2+iatsum_all2(ii)
       enddo

       iatsum_posit=1
       iatsum_posit2=1
       do ii=1,num_group
       call mpi_bcast(wmask(iatsum_posit2),iatsum_all2(ii),
     &  MPI_REAL8,ii-1,MPI_COMM_N,ierr)       ! bcast this  wmask to different color groups
       call mpi_bcast(indm(iatsum_posit),iatsum_all(ii),
     &  MPI_INTEGER,ii-1,MPI_COMM_N,ierr)
       call mpi_bcast(xyzmap(3*iatsum_posit-2),3*iatsum_all(ii),
     &  MPI_REAL8,ii-1,MPI_COMM_N,ierr)
        iatsum_posit=iatsum_posit+iatsum_all(ii)
        iatsum_posit2=iatsum_posit2+iatsum_all2(ii)
       enddo

       iatsum=iatsum_posit-1

       s = dfloat(iatsum)
       call global_sumr(s)
       avenmap = s/dfloat(natom)

       if(inode_tot.eq.1) then
       write(6,*) "ave nmap=",avenmap
       endif

        deallocate(iatsum_all)
        deallocate(iatsum_all2)


       endif !ilocal.eq.2

**************************************************

      if(ipsp_all.eq.2) then
       mrb2_q=2*(4*pi/3*rcut_qm**3)/(vol/(n1L*n2L*n3L))
       mrb2_matom_node_q=mrb2_q*(natom*1.d0/nnodes+2.d0)
       
       allocate(wmask_q(49,mrb2_matom_node_q))
       allocate(wmask_q0(12,mrb2_matom_node_q))
       allocate(indm_q(mrb2_matom_node_q))

       call getwmask_q(xatom,nmap_q,iatom,rcut_q1,rcut_q2,
     &  AL,mrb2_matom_node_q)

       endif


911    continue      ! jump out point for ivr_out.eq.2
**************************************************
**** initial potential vr_in
**************************************************
**** Note, input rho_nL is the full charge density on grid nr_nL
**************************************************


       if(ido_vr.eq.1) then
       if(islda.eq.1.and.igga.eq.0) then
       call getpot2L(rho_nL,vion_nL,rhocr_nL,vr_in_nL,
     &        v0,E_Hxc,E_coul,E_ion)
       endif
       if(islda.eq.2.and.igga.eq.0) then
       call getpot3L(rho_nL,vion_nL,rhocr_nL,vr_in_nL,v0,
     &           E_Hxc,E_coul,E_ion)
       endif

       if(islda.eq.1.and.igga.eq.1) then
       call getpot4L(rho_nL,vion_nL,rhocr_nL,vr_in_nL,v0,
     &           E_Hxc,E_coul,E_ion)
       endif

       if(islda.eq.2.and.igga.eq.1) then
       call getpot5L(rho_nL,vion_nL,rhocr_nL,vr_in_nL,v0,
     &           E_Hxc,E_coul,E_ion)
       endif

       endif

       

**************************************************
       if(ivr_out.eq.2) return    ! calculate vr_in_nL only, and return

**************************************************
**** end initialize 
**************************************************

       if(ipsp_all.eq.1) then
       do iislda=1,islda
       do ia=1,natom
       iitype=ityatom(ia)
       do iref2=1,numref(ia)
       do iref1=1,numref(ia)
       Dij(iref1,iref2,ia,iislda)=Dij0(iref1,iref2,iitype)
       enddo
       enddo
       enddo
       enddo

       else

       do iislda=1,islda
       call get_Dij(Dij(1,1,1,iislda),vr_in_nL(1,iislda),
     &   nmap_q,1.d0)
       enddo

       endif
******************************************************
       do iislda=1,islda
       call convert_SLvr(vr_n(1,iislda),
     &      vr_in_nL(1,iislda),-1)      ! from nr_nL to nr_n 
       enddo
******************************************************


       nscf=0
       nint_last=1
       do 2000 nint=1,niter

**** temperarily use vr_out as a work array
         if(nr_nL.ge.nr_n) then
         do iislda=1,islda
         do i=1,nr_nL
         vr_out_nL(i,iislda)=rho_nL(i,iislda)    ! assuming nr_nL.gt.nr_n
         enddo
         enddo
         endif
******************************************************
**** wavefunction update calculations
**** vr_in & mx, are used in CG_real, CG_new, diag_comp.
****  CG_real, CG_new output the rho(r) of mx states.
******************************************************
       do iislda=1,islda
       call mpi_bcast(vr_n(1,iislda),nr_n,MPI_REAL8,0,
     &  MPI_COMM_N,ierr)              ! making sure the vr_n are the same for different color group
       enddo
******************************************************
       ave_line=0.d0
       do 201 iislda=1,islda
       do 200 kpt=1,nkpt

        if((iislda-1)*nkpt+kpt.ge.kpt_slda_dis(1).and.
     &     (iislda-1)*nkpt+kpt.le.kpt_slda_dis(2)) then      ! the big division of work

       if(ilocal.eq.2) then
       call getcphase()
       endif


       call gen_G_comp(kpt,0) 
       call fftprep_comp(n1,n2,n3)

       call ugIO(ug_n,kpt,2,0,iislda)

       if(ilocal.eq.3) then
       call wqIO(nkpt,kpt,2)
       endif


       if(iCGmth(nint).eq.1) then
       call CG_comp(ilocal,nline,tolug,
     &    E_st(1,kpt,iislda),err_st(1,kpt,iislda),ave_linek,
     &    vr_n(1,iislda),workr_n,0,1.d0,kpt,iislda)
       else
       call CG_new(ilocal,nline,tolug,
     &    E_st(1,kpt,iislda),err_st(1,kpt,iislda),ave_linek,
     &    vr_n(1,iislda),workr_n,mCGbad,kpt,iislda)
       endif



       call ugIO(ug_n,kpt,1,0,iislda)

       if(ipsp_all.eq.2.and.(nkpt.gt.1.or.islda.gt.1)) then
       call beta_psiIO(beta_psi,kpt,1,0,iislda)        ! write the beta_psi
       endif

       ave_line=ave_line+ave_linek

        endif          ! the kpt_slda_dis(1),(2)
200    continue
201    continue

ccccccccccccccccccc   bcast the E_st,err_st

       do  iislda=1,islda
       do  kpt=1,nkpt
        ibcast=0
        if((iislda-1)*nkpt+kpt.ge.kpt_slda_dis(1).and.
     & (iislda-1)*nkpt+kpt.le.kpt_slda_dis(2)) ibcast=icolor
       
        call mpi_allreduce(ibcast,ibcast_tmp,1,MPI_INTEGER,
     &  MPI_SUM,MPI_COMM_N,ierr)
        ibcast=ibcast_tmp

        call mpi_bcast(E_st(1,kpt,iislda),mx,MPI_REAL8,ibcast,
     &  MPI_COMM_N,ierr)

        call mpi_bcast(err_st(1,kpt,iislda),mx,MPI_REAL8,ibcast,
     &  MPI_COMM_N,ierr)
        enddo
        enddo
cccccccccccccccccccccccccccccccccccccccccccccccccc


       call mpi_allreduce(ave_line,ave_line_tmp,1,MPI_REAL8,
     &   MPI_SUM,MPI_COMM_N,ierr) 

           ave_line=ave_line_tmp/nkpt/islda

           errmax=-100.d0
           errave=0.d0
           do iislda=1,islda
           do kpt=1,nkpt
           do m=1,mx
           if(err_st(m,kpt,iislda).gt.errmax) 
     &      errmax=err_st(m,kpt,iislda)
           errave=errave+err_st(m,kpt,iislda)
           enddo
           enddo
           enddo
           errave=errave/(mx*nkpt*islda)

********************************************************
         if(iscfmth(nint).le.0)  goto 99    ! non-self-consistent
           nscf=nscf+1

           call occup(FermidE(nint),itypeFermi(nint),
     &       E_st,totNel,nkpt,
     &       occ,E_NSC,TS,workr_n,Ef,islda)
          nint_last=nint
          if(inode_tot.eq.1) then
          write(6,*) "Ef=", Ef
          endif
    


         do iislda=1,islda
         call convert_SLvr(rho_n(1,iislda),
     &      rho_nL(1,iislda),1)    ! from nr_n to nr_nL
         enddo

ccccccccccccccccccccccccccccccccccccccccc

         E_dDrho=0.d0

         if(ipsp_all.eq.2) then
         call add_rho_beta(rho_nL,occ,islda,nkpt,nmap_q,
     &    vr_in_nL,E_dDrho)
         endif
ccccc E_dDrho= sum_i,j1,j2 <psi_i|beta_j1> delta_D_j1,j2 <beta_j2|psi_i>
ccccc Thus, it can remove the part of delta_D in \sum_i E_i, and keep only D_0 contrib. 


         if(nrot.gt.1) then
         do iislda=1,islda
         call symmop(rho_nL(1,iislda))
         enddo
         endif
***************************************************
*** E_dDrho=\int dDij*drho d^3r=<psi|V_NL(Dij)|psi>-<psi|V_NL(Dij0)|psi>
******************************************************
**************************************************
***** scf calc., updating vr_in 
**************************************************
***  calculate the total energy E (or say E_trial)
***  and dvE=\int (v_out-v_in)^2 d^3r
**************************************************
       E_psiV=0.d0

       do iislda=1,islda
       do i=1,nr_n
       E_psiV=E_psiV+rho_n(i,iislda)*vr_n(i,iislda)   ! rho_n is not symmetrize, but okay
       enddo
       enddo

       drho=0.d0
       E_rhoVext=0.d0

       do iislda=1,islda
       do i=1,nr_nL
       drho=drho+dabs(rho_nL(i,iislda)-vr_out_nL(i,iislda))
       E_rhoVext=E_rhoVext+rho_nL(i,iislda)*vext_nL(i)
       enddo
       enddo

       call global_sumr(E_psiV)
       call global_sumr(drho)
       call global_sumr(E_rhoVext)

cccccccccccccccccccccccccccccccc
       E_psiV=E_psiV*vol/nr      ! E_psiV is summed over nr, not nrL 
       drho=drho*vol/nrL
       E_rhoVext=E_rhoVext*vol/nrL


       if(islda.eq.1.and.igga.eq.0) then
       call getpot2L(rho_nL,vion_nL,rhocr_nL,vr_out_nL,
     &        v0,E_Hxc,E_coul,E_ion)
       endif
       if(islda.eq.2.and.igga.eq.0) then
       call getpot3L(rho_nL,vion_nL,rhocr_nL,vr_out_nL,v0,
     &           E_Hxc,E_coul,E_ion)
       endif

       if(islda.eq.1.and.igga.eq.1) then
       call getpot4L(rho_nL,vion_nL,rhocr_nL,vr_out_nL,v0,
     &           E_Hxc,E_coul,E_ion)
       endif

       if(islda.eq.2.and.igga.eq.1) then
       call getpot5L(rho_nL,vion_nL,rhocr_nL,vr_out_nL,v0,
     &           E_Hxc,E_coul,E_ion)
       endif



cccc Note: E_NSC-E_psiV-E_dDrho= E_kinetc+E_NL(D0)


       E_cor=E_ion-E_rhoVext-E_psiV-E_dDrho      ! E_ion=Vion*rho, includes E_rhoVext 

       if(icoul.eq.1.or.icoul.eq.11.or.icoul.eq.12.
     &                      or.icoul.eq.13) Ealphat=0.d0

******  E_extV is the energy due to all the external potential corrections

       E_extV=E_IVext+E_rhoVext    ! sum over ionic part and electronic part

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc   E_NSC+E_cor= E_kinetic+E_nonlocal+E_ion (E_ion=Vion(no Vext)*rho)
cccccccccccccccccccccccccccccccc
       E_tot=E_NSC+E_cor+E_Hxc+ewald+Ealphat-TS
       E_tot=E_tot+E_extV

       dvE=0.d0
       do iislda=1,islda
       do i=1,nr_nL
       dvE=dvE+(vr_in_nL(i,iislda)-vr_out_nL(i,iislda))**2
       enddo
       enddo

       call global_sumr(dvE)

       dvE=dvE/nrL
       dv_ave=dsqrt(dvE)

ccccccccccccccccccccccccccccccccc

      if(islda.eq.2) then
      scharge1=0.d0
      scharge2=0.d0
      do ikpt=1,nkpt
      do m=1,mx
      scharge1=scharge1+occ(m,ikpt,1)
      scharge2=scharge2+occ(m,ikpt,2)
      enddo
      enddo
      dcharge=0.d0
      do i=1,nr_nL
      dcharge=dcharge+dabs(rho_nL(i,1)-rho_nL(i,2))
      enddo
      call global_sumr(dcharge)
      dcharge=dcharge*vol/nrL
      endif


       if(inode_tot.eq.1) then

       write(6,*) "***************************************************"
       write(6,*) "iter=", nint
       write(6,*) "---------------------------------------------------"
       if(islda.eq.2) then
       write(6,411) scharge1,scharge2,dcharge
       write(6,*) "---------------------------------------------------"
       endif
       write(6,401) dvE, dvE-dvE0
       write(6,402) dv_ave,drho
       write(6,403) erraved,errmaxd
       write(6,404) errave,errmax
       write(6,*) "---------------------------------------------------"
       write(6,398) ewald
       write(6,399) Ealphat
       write(6,412) E_extVm E_extV-E_extV0
       write(6,405) E_NSC, E_NSC-E_NSC0
       write(6,406) E_cor, E_cor-E_cor0
       write(6,407) E_Hxc, E_Hxc-E_Hxc0
       write(6,408) -TS, -TS+TS0
       write(6,409) E_tot, E_tot-E_tot0
       write(6,*) "---------------------------------------------------"
       write(6,*) "E_ion,E_coul,E_xc",E_ion-E_rhoVext,
     &                 E_coul,E_Hxc-E_coul
       write(6,*) "ave(vtot)=ave_(V_xc)=: v0=",v0
       write(6,*) "Eigen energies are values after setting ave_vtot=0" 
       write(22,*) "-------------------------------------------"
       write(22,410) nint,ave_line
       write(22,415) Ef*27.211396d0
       write(22,416) errave,erraved
       write(22,402) dv_ave,drho
       if(islda.eq.2) then
       write(22,411) scharge1,scharge2,dcharge
       endif
       write(22,417) E_tot, E_tot-E_tot0
       call system_flush(22)
411   format(" spin_up;dn;loc_diff  = ", 3(f18.10,2x))
410   format(" iter=",i4,"     ave_lin=", f5.1)
416   format(" ugerr CG; ugerr Diag = ", 2(E10.4,1x))
417   format(" E_tot                = ",  E20.14,4x,E10.4)  
 
      endif


       convergE=dabs(E_tot-E_tot0)
       if(convergE.lt.tolE) goto 2001     ! finished

       if(nint.ne.niter) then
       E_NSC0=E_NSC
       E_cor0=E_cor
       E_extV0=E_extV
       E_Hxc0=E_Hxc
       TS0=TS
       E_tot0=E_tot
       dvE0=dvE
       endif
**************************************************
*** mch_pulay: input the current vr_in->vr_out,
*** using previous such relation, and make a linear combination
*** output: a new vr_in->vr_out, with smaller difference
**************************************************
*** mch_kerk, Thomas3: input the vr_in -> vr_out, 
*** Using non linear-combination prediction (preconditioning)
*** output: a new vr_in->vr_out, with smaller difference.
***  the vr_in is used in the next iteration
**************************************************

******************************************************

       call mch_pulay(vr_in_nL,vr_out_nL,nscf,AA,nreset,islda)

       do iislda=1,islda
       if(iscfmth(nint).eq.1) then
       call mch_kerk(vr_in_nL(1,iislda),vr_out_nL(1,iislda))  ! generate new vr_in_nL
       else
       call Thomas3(vr_in_nL(1,iislda),vr_out_nL(1,iislda),   ! to be installed, vionT_n
     &      Ef,totNel,islda,ipsp_all)
       endif
       enddo


***************************************************
*** end selfconsistent updating vr_in
***************************************************


******************************************************
       if(ipsp_all.eq.2) then
       do iislda=1,islda
       call get_Dij(Dij(1,1,1,iislda),vr_in_nL(1,iislda),
     &   nmap_q,1.d0)
       enddo
       endif
******************************************************


       do iislda=1,islda
       call convert_SLvr(vr_n(1,iislda),
     &      vr_in_nL(1,iislda),-1)      ! from nr_nL to nr_n 
       enddo

cccccccccccccccccccccccccccccccccccccccccccccc

99     continue      ! jumping point for non-self-consistency



       do 198 iislda=1,islda
       do 199 kpt=1,nkpt

       if((iislda-1)*nkpt+kpt.ge.kpt_slda_dis(1).and.
     &     (iislda-1)*nkpt+kpt.le.kpt_slda_dis(2)) then

       if(ilocal.eq.2) then
       call getcphase()
       endif

       call gen_G_comp(kpt,0) 
       call fftprep_comp(n1,n2,n3)

       call ugIO(ug_n,kpt,2,0,iislda)
 
       if(ilocal.eq.3) then
       call wqIO(nkpt,kpt,2)
       endif

       call diag_comp(ilocal,E_st(1,kpt,iislda),
     &  err_st(1,kpt,iislda),vr_n(1,iislda),workr_n,kpt,iislda)

       call ugIO(ug_n,kpt,1,0,iislda)

       endif       !  kpt_slda_dis(1),(2)

199    continue
198    continue

cccccccccccccccccccccccccccccccccccccccccccccccc
       do  iislda=1,islda
       do  kpt=1,nkpt
        ibcast=0
        if((iislda-1)*nkpt+kpt.ge.kpt_slda_dis(1).and.
     & (iislda-1)*nkpt+kpt.le.kpt_slda_dis(2)) ibcast=icolor
       
        call mpi_allreduce(ibcast,ibcast_tmp,1,MPI_INTEGER,
     &  MPI_SUM,MPI_COMM_N,ierr)
        ibcast=ibcast_tmp

        call mpi_bcast(E_st(1,kpt,iislda),mx,MPI_REAL8,ibcast,
     &  MPI_COMM_N,ierr)

        call mpi_bcast(err_st(1,kpt,iislda),mx,MPI_REAL8,ibcast,
     &  MPI_COMM_N,ierr)
        enddo
        enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

       if(iwg_out.eq.2.and.icolor.eq.0) then      ! write out the wavefunction in fwg_out for checkpoint
       call write_wg(fwg_out,AL,islda,nkpt)
       endif

          errmaxd=-100.d0
          erraved=0.d0
          do iislda=1,islda
          do kpt=1,nkpt
          do m=1,mx
          if(err_st(m,kpt,iislda).gt.errmaxd)
     &                 errmaxd=err_st(m,kpt,iislda)
          erraved=erraved+err_st(m,kpt,iislda)
          enddo
          enddo
          enddo
          erraved=erraved/(mx*nkpt*islda)


2000   continue  
******************************************************
2001   continue    ! jump out point for convergE.lt.tolE


      if(inode_tot.eq.1.and.nscf.gt.0) then
     
       write(22,*) "---------------------------------------------------"
       write(22,*) "E_Fermi(eV)=", Ef*27.211396d0
       write(22,*) "---------------------------------------------------"
       if(islda.eq.2) then
       write(22,411) scharge1, scharge2,dcharge
       write(22,*) "---------------------------------------------------"
       endif
       write(22,415) Ef*27.211396d0 
415   format(" Ef(eV)               = ",  E13.7)
       write(22,401) dvE, dvE-dvE0
401   format(" dvE, dvE(n)-dvE(n-1) = ", 2(E10.4,1x))
       write(22,402) dv_ave,drho
402   format(" dv_ave, drho_tot     = ", 2(E10.4,1x))
       write(22,403) erraved,errmaxd
403   format(" ug err,diag,[ave,max]= ", 2(E10.4,1x))
       write(22,404) errave,errmax
404   format(" ug err, CG ,[ave,max]= ", 2(E10.4,1x))
       write(22,*) "---------------------------------------------------"
       write(22,398) ewald
398   format(" Ewald        = ", E20.14)
       write(22,399) Ealphat
399   format(" Alpha        = ", E20.14)
       write(22,412) E_extV,E_extV-E_extV0
412   format(" E_extV       = ", E20.14,4x,E10.4)
       write(22,405) E_NSC, E_NSC-E_NSC0
405   format(" E_NSC        = ", E20.14,4x,E10.4)  
       write(22,406) E_cor, E_cor-E_cor0
406   format(" E[-rho*V_Hxc]= ", E20.14,4x,E10.4)  
       write(22,407) E_Hxc, E_Hxc-E_Hxc0
407   format(" E_Hxc        = ", E20.14,4x,E10.4)  
       write(22,408) -TS, -TS+TS0
408   format(" -TS          = ", E20.14,4x,E10.4)  
       write(22,409) E_tot, E_tot-E_tot0
409   format(" E_tot        = ", E20.14,4x,E10.4)  
       write(22,*) "---------------------------------------------------"
       if(itypeFermi(nint_last).eq.1) then
       N_tmp=0.d0
       else
       N_tmp=itypeFermi(nint_last)-20
       endif
       F_dE=E_tot
       E_dE=E_tot+TS
       E_d0=((N_tmp+1)*F_dE+E_dE)/(N_tmp+2)
       write(22,414) E_d0,N_tmp
414    format("Zero temp. E_tot =", E20.14, 
     &  " Using formula: E_tot(T)+TS/(N+2), N=", i5)
       write(22,*) "---------------------------------------------------"
       write(22,396) E_coul,E_Hxc-E_coul,E_ion-E_rhoVext
       write(22,413) E_rhoVext,E_IVext
       write(22,397) E_psiV,E_dDrho
       write(22,395) v0
       write(22,*)
     &  "ave(V_ion_s(or p,d))=ave(V_Hatree)=0; ave(Vtot)=ave(V_xc)=v0"
       write(22,*) "---------------------------------------------------"
396   format(" E_Hart,E_xc,E_ion =", 3(E16.10,2x))
397   format("    E_psiV,E_dDrho =", 2(E16.10,2x))
395   format("      ave(vtot):v0 =", E16.10)
413   format(" E_rhoVext,E_IVext     =", 2(E16.10,2x))


c       do kpt=1,nkpt
c       write(22,*) "kpt= ", kpt
c       write(22,*) "err of each states, A.U"
c       write(22,102) (err_st(i,kpt), i=1,mx)
c       write(22,*) "eigen energies, in eV"
c       write(22,103) (E_st(i,kpt)*27.211396d0, i=1,mx)
c       write(22,*) "*********************************"
c       enddo

        endif
101   format(5(i6,7x))
102   format(5(E10.4,3x))
103   format(5(f12.8,1x))



******************************************************
*** calculate forces on each atom
****************************************************
**** using vr_out_nL(i,1) and vionT_nL(i) as working array, 
**** their contents will be destroyed. Bad practice
******************************************************
      if(iforce_cal.eq.1) then

      if(islda.eq.1.and.igga.eq.0) then
      do i=1,nr_nL
      vr_out_nL(i,1)=rho_nL(i,1)
      vionT_nL(i)=UxcCA(rho_nL(i,1)+rhocr_nL(i),uxc2)  
      enddo
      endif

      if(islda.eq.2.and.igga.eq.0) then
      do i=1,nr_nL
      vr_out_nL(i,1)=rho_nL(i,1)+rho_nL(i,2)
      call UxcCA2(rho_nL(i,1)+rhocr_nL(i)*0.5d0,
     &      rho_nL(i,2)+rhocr_nL(i)*0.5d0,
     &    vxc1,vxc2,uxc1,uxc2)
      vionT_nL(i)=(vxc1+vxc2)/2  
      enddo
      endif

      if(islda.eq.1.and.igga.eq.1) then
      s=0.d0
      do i=1,nr_nL
      vr_out_nL(i,1)=rho_nL(i,1)
      s=s+dabs(rhocr_nL(i))
      enddo
      call global_sumr(s)
      s=s*vol/nrL
        if(s.gt.1.D-5) then
        call getpot4_force(rho_nL,vionT_nL,rhocr_nL)
        else
        vionT_nL=0.d0          
        endif
      endif
      

      if(islda.eq.2.and.igga.eq.1) then
      s=0.d0
      do i=1,nr_nL
      vr_out_nL(i,1)=rho_nL(i,1)+rho_nL(i,2)
      s=s+dabs(rhocr_nL(i))
      enddo
      call global_sumr(s)
      s=s*vol/nrL
        if(s.gt.1.D-5) then
        call getpot5_force(rho_nL,vionT_nL,rhocr_nL)
        else
        vionT_nL=0.d0
        endif
      endif


cccc vionT_nL is used for force due to possible nonlinear core correction term

        if(icoul.eq.0) then
        call forcLC(AL,vr_out_nL(1,1),vionT_nL,    
     &   xatom,ntype,iatom,fatom)
        endif

        if(icoul.eq.1.or.icoul.eq.11.or.icoul.eq.12.
     &      or.icoul.eq.13) then
        call forcLC2(AL,vr_out_nL(1,1),vionT_nL,    
     &   xatom,ntype,iatom,fatom)
        endif

cccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccc
      if(ipsp_all.eq.2) then

       allocate(wmask_dq(49,0:3,mrb2_matom_node_q))
       allocate(wmask_dq0(12,0:3,mrb2_matom_node_q))

       call getwmask_dq(xatom,nmap_q,iatom,rcut_q1,rcut_q2,
     &  AL,mrb2_matom_node_q)

       call get_VdqdR(vr_in_nL,occ,islda,nkpt,nmap_q,
     &   fatom)

       deallocate(wmask_dq)
       deallocate(wmask_dq0)
       endif
cccccccccccccccccccccccccccccccccccccc
       
      if(ilocal.eq.3) then

      allocate(fatom_tmp(3,natom))
      allocate(fatom_tmp2(3,natom))
      fatom_tmp=0.d0
      fatom_tmp2=0.d0
      do iislda=1,islda
      do kpt=1,nkpt

        if((iislda-1)*nkpt+kpt.ge.kpt_slda_dis(1).and.
     &     (iislda-1)*nkpt+kpt.le.kpt_slda_dis(2)) then
       
      call ugIO(ug_n,kpt,2,0,iislda)
      call wqIO(nkpt,kpt,2)

      call forcNLq(fatom_tmp,occ,kpt,nkpt,iislda,islda,E_st)
ccccccc the resulting force is only correct after symmforc, since
ccccccc only the unsymmetrized k-points are used
         endif
      enddo
      enddo

      call mpi_allreduce(fatom_tmp,fatom_tmp2,3*natom,
     & MPI_REAL8,MPI_SUM,MPI_COMM_N,ierr)
       do ia=1,natom
       fatom(1,ia)=fatom(1,ia)+fatom_tmp2(1,ia)
       fatom(2,ia)=fatom(2,ia)+fatom_tmp2(2,ia)
       fatom(3,ia)=fatom(3,ia)+fatom_tmp2(3,ia)
       enddo
      deallocate(fatom_tmp)
      deallocate(fatom_tmp2)

      endif     ! ilocal.eq.1

cccccccccccccccccccccccccccccccccccc

      if(ilocal.eq.2) then

      allocate(wmaskX(20*mrb2_matom_node,3))    
      allocate(wmaskXtmp(20,mrb2,3))

       nmap=0    ! for the convenience of passing nmap to other icolor
       iatsum2=0
       do ia=natom_dis(1),natom_dis(2)
       iitype=ityatom(ia)
       inew=0
       iend=0
       if(ia.eq.natom_dis(1)) then
       inew=1
       else
       if(iitype.ne.ityatom(ia-1)) inew=1
       endif
       if(ia.eq.natom_dis(2)) then
       iend=1
       else
       if(iitype.ne.ityatom(ia+1)) iend=1
       endif

         call getwmaskX(xatom(1,ia),nmap(ia),indmtmp,
     &   ityatom(ia),wmaskXtmp,AL,workr_n,mrb2,
     &    nref,inew,iend)

         if(nref.ne.numref(ia)) then
         write(6,*) "nref.ne.numref,getwmaskX,stop",nref,numref(ia)
         call mpi_abort(MPI_COMM_WORLD,ierr)
         endif

         do i=1,nmap(ia)
         do j=1,nref
         iatsum2=iatsum2+1
         wmaskX(iatsum2,1)=wmaskXtmp(j,i,1)
         wmaskX(iatsum2,2)=wmaskXtmp(j,i,2)
         wmaskX(iatsum2,3)=wmaskXtmp(j,i,3)
         enddo
         enddo
       enddo  ! natom 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       call mpi_barrier(MPI_COMM_WORLD,ierr)

       allocate(nmap_tmp(natom))
       call mpi_allreduce(nmap,nmap_tmp,natom,        ! this is used to pass nmap to other icolor
     & MPI_INTEGER,MPI_SUM,MPI_COMM_N,ierr)
       nmap(1:natom)=nmap_tmp(1:natom)
       deallocate(nmap_tmp)


       allocate(iatsum_all2(num_group))
        iatsum_all2(icolor+1)=iatsum2
       do ii=1,num_group
       call mpi_bcast(iatsum_all2(ii),1,MPI_INTEGER,
     &   ii-1,MPI_COMM_N,ierr)
       enddo

         iatsum_posit2=0
       do ii=1,num_group
        if(icolor.eq.ii-1) then
        do i=iatsum_all2(ii),1,-1
        wmaskX(i+iatsum_posit2,1)=wmaskX(i,1)       ! shift wmask into the right position
        wmaskX(i+iatsum_posit2,2)=wmaskX(i,2)       ! shift wmask into the right position
        wmaskX(i+iatsum_posit2,3)=wmaskX(i,3)       ! shift wmask into the right position
        enddo
        endif
        iatsum_posit2=iatsum_posit2+iatsum_all2(ii)
       enddo

       iatsum_posit2=1
       do ii=1,num_group
       call mpi_bcast(wmaskX(iatsum_posit2,1),iatsum_all2(ii),
     &  MPI_REAL8,ii-1,MPI_COMM_N,ierr)       ! bcast this  wmask to different color groups
       call mpi_bcast(wmaskX(iatsum_posit2,2),iatsum_all2(ii),
     &  MPI_REAL8,ii-1,MPI_COMM_N,ierr)       ! bcast this  wmask to different color groups
       call mpi_bcast(wmaskX(iatsum_posit2,3),iatsum_all2(ii),
     &  MPI_REAL8,ii-1,MPI_COMM_N,ierr)       ! bcast this  wmask to different color groups
        iatsum_posit2=iatsum_posit2+iatsum_all2(ii)
       enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      allocate(fatom_tmp(3,natom))
      allocate(fatom_tmp2(3,natom))
      fatom_tmp=0.d0
      fatom_tmp2=0.d0

       do 301 iislda=1,islda
       do 302 kpt=1,nkpt

        if((iislda-1)*nkpt+kpt.ge.kpt_slda_dis(1).and.
     &     (iislda-1)*nkpt+kpt.le.kpt_slda_dis(2)) then

       call getcphase()
       call gen_G_comp(kpt,0) 
       call fftprep_comp(n1,n2,n3)

       call ugIO(ug_n,kpt,2,0,iislda)

       call forcNLr(workr_n,fatom_tmp,
     &  occ,kpt,nkpt,iislda,islda,E_st)       ! rewrite forcNLr inside, blas2 zgemv

          endif

302    continue
301    continue

      call mpi_allreduce(fatom_tmp,fatom_tmp2,3*natom,
     & MPI_REAL8,MPI_SUM,MPI_COMM_N,ierr)

       do ia=1,natom
       fatom(1,ia)=fatom(1,ia)+fatom_tmp2(1,ia)
       fatom(2,ia)=fatom(2,ia)+fatom_tmp2(2,ia)
       fatom(3,ia)=fatom(3,ia)+fatom_tmp2(3,ia)
       enddo
      deallocate(fatom_tmp)
      deallocate(fatom_tmp2)

       deallocate(wmaskX)
       deallocate(wmaskXtmp)

      endif     ! ilocal.eq.2
cccccccccccccccccccccccccccccccccccc
       if(nrot.gt.1) then
       call symmopf(smatr,nrot,AL,fatom,xatom,iatom)
       endif


      endif   ! iforce_cal.eq.1
**************************************************

  

      if(ilocal.eq.2) then
      deallocate(wmask)
      deallocate(xyzmap)
      deallocate(cphase)
      deallocate(indm)
      deallocate(xyzmaptmp)
      deallocate(indmtmp)
      endif

      if(ipsp_all.eq.2) then
       deallocate(wmask_q)
       deallocate(wmask_q0)
       deallocate(indm_q)
       endif

      
      deallocate(AA)
      return

      contains
      
*****************************************************
*****************************************************
       subroutine getcphase()
       implicit double precision (a-h,o-z)

       complex*16 cai
    
       cai=dcmplx(0.d0,1.d0)

       ico1=0
       do ia=1,natom
        do i=1,nmap(ia)
        ico1=ico1+1
        cphase(ico1)=cdexp(cai*(xyzmap(ico1*3-2)*akx(kpt)+
     &   xyzmap(ico1*3-1)*aky(kpt)+xyzmap(ico1*3)*akz(kpt)))
        enddo
       enddo

       return
       end subroutine getcphase

*****************************************************

        end
      
