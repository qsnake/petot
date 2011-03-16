      subroutine Etotcalc(xatom,fatom,workr_n,E_tot,
     & iforce_cal,ido_rho,ido_vr,tolug,tolE,niter,nline,
     &  iCGmth,iscfmth,FermidE,mCGbad,E_st,err_st,AL,
     &  nkpt,ntype,convergE,islda,igga)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
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


       real*8,allocatable,dimension(:)  :: xyzmaptmp
       real*8,allocatable,dimension(:,:)  :: wmasktmp
       integer,allocatable,dimension(:) :: indmtmp

       real*8 xatom(3,matom)
       real*8 fatom(3,matom)

       real*8 occ_t(mtype)
       integer iiatom(mtype),icore(mtype),numref(matom)
       integer is_ref(mtype),ip_ref(mtype),id_ref(mtype)

       integer nmap(matom)

       real*8 zatom(mtype)

       integer iCGmth(100),iscfmth(100)
       real*8 FermidE(100)
cccccccccccccccccccccccc

       integer iatom(matom),ityatom(matom)
       real*8  totNel
       integer smatr(3,3,48),nrot,ilocal
       real*8 Ealpha(mtype)
*************************************************
       character*20 vwr_atom(mtype)
       character*20 fwg_in,fwg_out,fdens_in,fdens_out
*************************************************

       common /comNL2/occ_t,iiatom,icore,numref
       common /comispd_ref/is_ref,ip_ref,id_ref
       common /comEk/Ek
       common /comzatom/zatom
       common /comnmap/nmap

       common /comMainEtot/iatom,ityatom,totNel,
     &     smatr,nrot,ilocal,Ealpha,Ealphat
c

c saved arrays for mch_pulay
c
       real*8 AA(nmax,nmax)


       complex*16 workr_n(mr_n)

**************************************************
c initialize mpi and number each node
c
       pi=4.0d0*datan(1.0d0)

       E_NSC0 = 0.0d0  !amc
       TS0= 0.0d0
       n33 = 0 
       E_COR0 = 0.0d0 
       E_HXC0 = 0.0d0
       E_TOT0 = 0.0d0
       DVE0 = 0.0d0

       Ek=0.5d0
       occ=1.d0

       call getewald(fatom,xatom,AL,ityatom,ewald)
       if(iforce_cal.eq.0) fatom=0.d0

       call getVrho(AL,vion_n,vionT_n,xatom,ntype,iatom,
     &  rhocr_n,totNel,ido_rho,workr_n,islda)

       ntot=totNel*1.00001d0
*********************************************************
***  ilocal.eq.3, q space nolocal Kleimen_Bylander
*********************************************************
       if(inode.eq.1) write(6,*) "into getwmask, or getwq"

       if(ilocal.eq.3) then
       do kpt=1,nkpt
       call getwq(AL,ntype,iatom,ityatom,xatom,kpt)
       call wqIO(nkpt,kpt,1)     ! 1: write, 2: read
       enddo
       endif

*********************************************************
***  ilocal.eq.2, r space nolocal Kleimen_Bylander
*********************************************************
       if(ilocal.eq.2) then

ccccccc this formula sometime is not right, if many atoms 
ccccccc are located inside one processor

       mrb2=2*(4*pi/3*rcut**3)/(vol/(n1*n2*n3))

       mrb2_matom_node=mrb2*natom/nnodes

       allocate(wmask(9*mrb2_matom_node))
       allocate(xyzmap(3*mrb2_matom_node))
       allocate(cphase(mrb2_matom_node))
       allocate(indm(mrb2_matom_node))
       allocate(wmasktmp(9,mrb2))
       allocate(xyzmaptmp(3*mrb2))
       allocate(indmtmp(mrb2))

       iatsum = 0
       iatsum2=0 
       do 20 ia=1,natom
       iitype=ityatom(ia)
       call getwmask(xatom(1,ia),nmap(ia),indmtmp,
     &  ityatom(ia),wmasktmp,xyzmaptmp,AL,workr_n,mrb2,
     &  is_ref(iitype),ip_ref(iitype),id_ref(iitype),nref)
        numref(ia)=nref
       if(iatsum+nmap(ia).gt.mrb2_matom_node) then
       write(6,*) "iatsum.gt.mbr2_matom_node, stop"
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

       s = dfloat(iatsum)
       call global_sumr(s)
       avenmap = s/dfloat(natom)

       if(inode.eq.1) then
       write(6,*) "ave nmap=",avenmap
       endif

       endif !ilocal.eq.2


       call mpi_barrier(MPI_COMM_WORLD,ierr)

**************************************************
**** initial potential vr_in
**************************************************
       if(ido_vr.eq.1) then
       if(islda.eq.1.and.igga.eq.0) then
       call getpot2(rho_n(1,1),vion_n,rhocr_n,vr_in_n(1,1),
     &        v0,E_Hxc,workr_n,E_coul,E_ion)
       endif
       if(islda.eq.2.and.igga.eq.0) then
       call getpot3(rho_n,vion_n,rhocr_n,vr_in_n,v0,
     &           E_Hxc,workr_n,E_coul,E_ion)   
       endif

       if(islda.eq.1.and.igga.eq.1) then
       call getpot4(rho_n,vion_n,rhocr_n,vr_in_n,v0,
     &           E_Hxc,workr_n,E_coul,E_ion)   
       endif

       if(islda.eq.2.and.igga.eq.1) then
       call getpot5(rho_n,vion_n,rhocr_n,vr_in_n,v0,
     &           E_Hxc,workr_n,E_coul,E_ion)   
       endif
       
       endif
**************************************************
**** end initialize 
**************************************************

       if(inode.eq.1) write(6,*) "start iterations"

       nscf=0
       do 2000 nint=1,niter

******************************************************
**** tmperarily use vr_out as a work array
         do iislda=1,islda
         do i=1,nr_n
         vr_out_n(i,iislda)=rho_n(i,iislda)
         enddo
         enddo
******************************************************
**** wavefunction update calculations
**** vr_in & mx, are used in CG_real, CG_new, diag_comp.
****  CG_real, CG_new output the rho(r) of mx states.
******************************************************
       ave_line=0.d0
       do 201 iislda=1,islda
       do 200 kpt=1,nkpt

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
     &    vr_in_n(1,iislda),workr_n,0,1.d0,kpt)
       else
       call CG_new(ilocal,nline,tolug,
     &    E_st(1,kpt,iislda),err_st(1,kpt,iislda),ave_linek,
     &    vr_in_n(1,iislda),workr_n,mCGbad,kpt)
       endif

       call ugIO(ug_n,kpt,1,0,iislda)

       ave_line=ave_line+ave_linek
200    continue
201    continue

           ave_line=ave_line/nkpt/islda

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

           call occup(FermidE(nint),E_st,totNel,nkpt,
     &       occ,E_NSC,TS,workr_n,Ef,islda)
    
         if(nrot.gt.1) then
         do iislda=1,islda
         call symmop(rho_n(1,iislda),workr_n)
         enddo
         endif
******************************************************
**************************************************
***** scf calc., updating vr_in 
**************************************************
***  calculate the total energy E (or say E_trial)
***  and dvE=\int (v_out-v_in)^2 d^3r
**************************************************
       E_cor=0.d0
       drho=0.d0

       do iislda=1,islda
       do i=1,nr_n
       E_cor=E_cor+rho_n(i,iislda)*
     &              (vion_n(i)-vr_in_n(i,iislda))
       drho=drho+dabs(rho_n(i,iislda)-vr_out_n(i,iislda))
       enddo
       enddo

       call global_sumr(E_cor)
       call global_sumr(drho)

       E_cor=E_cor*vol/nr
       drho=drho*vol/nr

       if(islda.eq.1.and.igga.eq.0) then
       call getpot2(rho_n(1,1),vion_n,rhocr_n,vr_out_n(1,1),
     &       v0,E_Hxc,workr_n,E_coul,E_ion)
       endif
       if(islda.eq.2.and.igga.eq.0) then
       call getpot3(rho_n,vion_n,rhocr_n,vr_out_n,v0,
     &           E_Hxc,workr_n,E_coul,E_ion)   
       endif
       if(islda.eq.1.and.igga.eq.1) then
       call getpot4(rho_n,vion_n,rhocr_n,vr_out_n,v0,
     &           E_Hxc,workr_n,E_coul,E_ion)   
       endif
       if(islda.eq.2.and.igga.eq.1) then
       call getpot5(rho_n,vion_n,rhocr_n,vr_out_n,v0,
     &           E_Hxc,workr_n,E_coul,E_ion)   
       endif

       E_tot=E_NSC+E_cor+E_Hxc+ewald+Ealphat+TS

       dvE=0.d0
       do iislda=1,islda
       do i=1,nr_n
       dvE=dvE+(vr_in_n(i,iislda)-vr_out_n(i,iislda))**2
       enddo
       enddo

       call global_sumr(dvE)

       dvE=dvE/nr
       dv_ave=dsqrt(dvE)

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
      do i=1,nr_n
      dcharge=dcharge+dabs(rho_n(i,1)-rho_n(i,2))
      enddo
      call global_sumr(dcharge)
      dcharge=dcharge*vol/nr
      endif

       if(inode.eq.1) then

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
       write(6,405) E_NSC, E_NSC-E_NSC0
       write(6,406) E_cor, E_cor-E_cor0
       write(6,407) E_Hxc, E_Hxc-E_Hxc0
       write(6,408) TS, TS-TS0
       write(6,409) E_tot, E_tot-E_tot0
       write(22,*) "-------------------------------------------"
       write(22,410) nint,ave_line,errave,erraved
       if(islda.eq.2) then
       write(22,411) scharge1,scharge2,dcharge
       endif
       write(22,409) E_tot, E_tot-E_tot0
       call system_flush(22)
411   format("charges of spin up,down and loc_diff ", 3(f18.10,2x))
410   format("iter=",i4,"  ave_lin=", f5.1, "  ugerr CG=", E8.1,
     &         "  ugerr Diag=", E8.1 )
 
      endif

       convergE=dabs(E_tot-E_tot0)
       if(convergE.lt.tolE) goto 2001     ! finished

       if(nint.ne.niter) then
       E_NSC0=E_NSC
       E_cor0=E_cor
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
       call mch_pulay(vr_in_n,vr_out_n,nscf,AA,nreset,islda)

       do iislda=1,islda
       if(iscfmth(nint).eq.1) then
       call mch_kerk(vr_in_n(1,iislda),vr_out_n(1,iislda),
     &     workr_n)
       else
       call Thomas3(vr_in_n(1,iislda),vr_out_n(1,iislda),
     &      Ef,ntot,workr_n,islda)
       endif
       enddo

***************************************************
*** end selfconsistent updating vr_in
***************************************************

99     continue      ! jumping point for non-self-consistency



       do 198 iislda=1,islda
       do 199 kpt=1,nkpt

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
     &  err_st(1,kpt,iislda),vr_in_n(1,iislda),workr_n,kpt)

       call ugIO(ug_n,kpt,1,0,iislda)

199    continue
198    continue

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

      if(inode.eq.1.and.nscf.gt.0) then
     
       write(22,*) "---------------------------------------------------"
       write(22,*) "E_Fermi(eV)=", Ef*27.211396d0
       write(22,*) "---------------------------------------------------"
       if(islda.eq.2) then
       write(22,411) scharge1, scharge2,dcharge
       write(22,*) "---------------------------------------------------"
       endif
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
       write(22,405) E_NSC, E_NSC-E_NSC0
405   format(" E_NSC        = ", E20.14,4x,E10.4)  
       write(22,406) E_cor, E_cor-E_cor0
406   format(" E[-rho*V_Hxc]= ", E20.14,4x,E10.4)  
       write(22,407) E_Hxc, E_Hxc-E_Hxc0
407   format(" E_Hxc        = ", E20.14,4x,E10.4)  
       write(22,408) TS, TS-TS0
408   format(" TS           = ", E20.14,4x,E10.4)  
       write(22,409) E_tot, E_tot-E_tot0
409   format(" E_tot        = ", E20.14,4x,E10.4)  
       write(22,*) "---------------------------------------------------"
       write(22,396) E_coul,E_Hxc-E_coul,E_ion
       write(22,*) "---------------------------------------------------"
396   format(" E_Hart,E_xc,E_ion =", 3(E16.10,2x))

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
**** using vr_out_n(i,1) and vionT_n(i) as working array, 
**** their contents will be destroyed. Bad practice
******************************************************
      if(iforce_cal.eq.1) then

      if(islda.eq.1.and.igga.eq.0) then
      do i=1,nr_n
      vr_out_n(i,1)=rho_n(i,1)
      vionT_n(i)=UxcCA(rho_n(i,1)+rhocr_n(i),uxc2)
      enddo
      endif

      if(islda.eq.2.and.igga.eq.0) then
      do i=1,nr_n
      vr_out_n(i,1)=rho_n(i,1)+rho_n(i,2)
      call UxcCA2(rho_n(i,1)+rhocr_n(i)*0.5d0,
     &      rho_n(i,2)+rhocr_n(i)*0.5d0,
     &    vxc1,vxc2,uxc1,uxc2)
      vionT_n(i)=(vxc1+vxc2)/2
      enddo
      endif

      if(islda.eq.1.and.igga.eq.1) then
      s=0.d0
      do i=1,nr_n
      vr_out_n(i,1)=rho_n(i,1)
      s=s+dabs(rhocr_n(i))
      enddo
      call global_sumr(s)
      s=s*vol/nr
        if(s.gt.1.D-5) then
        call getpot4_force(rho_n,vionT_n,rhocr_n,workr_n)
        else
        vionT_n=0.d0
        endif
      endif
      

      if(islda.eq.2.and.igga.eq.1) then
      s=0.d0
      do i=1,nr_n
      vr_out_n(i,1)=rho_n(i,1)+rho_n(i,2)
      s=s+dabs(rhocr_n(i))
      enddo
      call global_sumr(s)
      s=s*vol/nr
        if(s.gt.1.D-5) then
        call getpot5_force(rho_n,vionT_n,rhocr_n,workr_n)
        else
        vionT_n=0.d0
        endif
      endif

        call forcLC(AL,vr_out_n(1,1),vionT_n,
     &   xatom,ntype,iatom,fatom,workr_n)


      if(ilocal.eq.3) then
      do iislda=1,islda
      do kpt=1,nkpt
       
      call ugIO(ug_n,kpt,2,0,iislda)
      call wqIO(nkpt,kpt,2)

      call forcNLq(fatom,occ(1,1,iislda),kpt,nkpt)
ccccccc the resulting force is only correct after symmforc, since
ccccccc only the unsymmetrized k-points are used

      enddo
      enddo
      endif

cccccccccccccccccccccccccccccccccccc
      if(ilocal.eq.2) then

      allocate(wmaskX(9*mrb2_matom_node))

      do 300 ixyz=1,3

       iatsum2=0
       do ia=1,natom
       iitype=ityatom(ia)
         call getwmaskX(xatom(1,ia),nmap(ia),indmtmp,
     &   ityatom(ia),wmasktmp,AL,workr_n,ixyz,mrb2,
     &    is_ref(iitype),ip_ref(iitype),id_ref(iitype),nref)

         if(nref.ne.numref(ia)) then
         write(6,*) "nref.ne.numref,getwmaskX,stop",nref,numref(ia)
         call mpi_abort(MPI_COMM_WORLD,ierr)
         endif

         do i=1,nmap(ia)
         do j=1,nref
         iatsum2=iatsum2+1
         wmaskX(iatsum2)=wmasktmp(j,i)
         enddo
         enddo
       enddo  ! natom 

       do 301 iislda=1,islda
       do 302 kpt=1,nkpt

       call getcphase()
       call gen_G_comp(kpt,0) 
       call fftprep_comp(n1,n2,n3)

       call ugIO(ug_n,kpt,2,0,iislda)

       call forcNLr(workr_n,fatom,
     &  occ(1,1,iislda),kpt,nkpt,ixyz)

302    continue
301    continue
300    continue
       deallocate(wmaskX)

      endif
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
      deallocate(wmasktmp)
      deallocate(xyzmaptmp)
      deallocate(indmtmp)
      endif

      
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
      
