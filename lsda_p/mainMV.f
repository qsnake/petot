******************************************
*******************************************
      program mainMV
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************


      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'mpif.h'

      include 'param.escan_real'

      integer status(MPI_STATUS_SIZE)

******************************************
       real*8 AL(3,3),ALt(3,3)
********************************************
       real*8,allocatable,dimension(:,:,:)  :: E_st,err_st
       real*8,allocatable,dimension(:)  :: wmasktmp,xyzmaptmp
       integer,allocatable,dimension(:) :: indmtmp
       real*8,allocatable,dimension(:) :: rhonew1_nL,rhonew2_nL

       complex*16,allocatable,dimension(:,:)  :: ug_n_tmp
       real*8,allocatable,dimension(:)    :: workr_n_tmp
       complex*16,allocatable,dimension(:)  :: workr_n

       real*8, allocatable, dimension(:)    :: vr_nL

       integer icoul
       real*8 xcoul(3)


       real*8 xatom(3,matom),xatom_old(3,matom)
       integer imov_at(3,matom)

       real*8 fatom(3,matom)
       real*8 occ_t(mtype)
       real*8 zatom(mtype),Ealpha(mtype)
       integer iiatom(mtype),icore(mtype)
       integer nmap(matom)
       integer iatom(matom),ityatom(matom),numref(matom)
       integer is_ref(mtype),ip_ref(mtype),id_ref(mtype),
     &   nref_type(mtype),ipsp_type(mtype)
       integer iCGmth0(100),iCGmth1(100),iscfmth0(100),iscfmth1(100)
       real*8  FermidE0(100),FermidE1(100),totNel
       integer  itypeFermi0(100),itypeFermi1(100)
       integer smatr(3,3,48),nrot
       integer niter0,nline0,niter1,nline1,mCGbad0,mCGbad1
       integer kpt_dens(2),ispin_dens(2),iw_dens(2)
       integer lll(8,mtype),nbeta(mtype)
       integer isNLa(9,matom)
       real*8 Dij0(32,32,mtype),Qij(32,32,mtype)
*************************************************
       character*20 vwr_atom(mtype),fforce_out,fdens_out
       character*20 fwg_in(2),fwg_out(2),frho_in(2),frho_out(2),
     & fvr_in(2),fvr_out(2),f_tmp,fxatom_out,f_xatom,fvext_in
       character*20 file_tmp

       character*60 message
*************************************************

       common /comNL2/occ_t,iiatom,icore,numref,ityatom
       common /comisNLa/isNLa,Dij0,Qij,ipsp_all,ipsp_type
       common /comEk/Ek
       common /comzatom/zatom
       common /comnmap/nmap
       common /comispd_ref/is_ref,ip_ref,id_ref

       common /comMainEtot/iatom,totNel,
     &     smatr,nrot,ilocal,Ealpha,Ealphat
       common /comlll/lll,nbeta
       common /comcoul/icoul,xcoul
       common /comVext/ivext_in
c
c saved arrays for mch_pulay
c

       complex*16,allocatable,dimension(:) :: ugtemp
       complex*16,allocatable,dimension(:) :: temp_array


**************************************************
c      atime00=mclock()/100.d0
c
c initialize mpi and number each node
c

       open(10,file="etot.input")
       rewind(10)
       read(10,*) i1, nnodes_k,num_group
       close(10)


       call mpi_init(ierr)
       call mpi_comm_rank(MPI_COMM_WORLD,inode_tot,ierr)
       call mpi_comm_size(MPI_COMM_WORLD,nnodes_tot,ierr)

       if(nnodes_k*num_group.ne.nnodes_tot) then
       write(6,*) "nnodes_k*num_group.ne.nnodes_tot,stop",
     &  nnodes_k,num_group,nnodes_tot
       stop
       endif

       icolor=inode_tot/nnodes_k
       ikey=inode_tot-icolor*nnodes_k
 
       call mpi_comm_split(MPI_COMM_WORLD,icolor,ikey,
     &       MPI_COMM_K,ierr)
cccccc  splits the processor into num_group groups, the MPI_COMM_K is the group communicator
cccccc  the ikey is the local inode

       inode_tot=inode_tot+1
       nnodes=nnodes_k
       inode=ikey+1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       call mpi_comm_split(MPI_COMM_WORLD,ikey,icolor,
     &       MPI_COMM_N,ierr)
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       pi=4.0d0*datan(1.0d0)

       E_NSC0 = 0.0d0  !amc
       TS0= 0.0d0
       n33 = 0 
       E_COR0 = 0.0d0 
       E_HXC0 = 0.0d0
       E_TOT0 = 0.0d0
       DVE0 = 0.0d0

       call input(AL,ilocal,tolug,tolE,
     &  niter0,niter1,nline0,nline1,iCGmth0,iCGmth1,
     &  iscfmth0,iscfmth1,FermidE0,FermidE1,
     &  itypeFermi0,itypeFermi1,mCGbad0,mCGbad1,
     &  islda,igga,totNel,ntype,xatom,iatom,iiatom,ityatom,
     &  fxatom_out,iwg_in,fwg_in,iwg_out,fwg_out,irho_in,frho_in,
     &  irho_out,frho_out,ivr_in,fvr_in,ivr_out,fvr_out,
     &  iforce,fforce_out,ipsp_type,ipsp_all,vwr_atom,nkpt,smatr,
     &  nrot,num_mov,tolforce,imov_at,dtstart,dd_limit,
     &  idens_out,kpt_dens,ispin_dens,iw_dens,fdens_out,
     &  f_xatom,numref,nref_type,ivext_in,fvext_in,imv_cont)

cccccccccccccccccccccccccccccccccccccccc
        n_tmp1=nkpt*islda/num_group
        n_tmp2=nkpt*islda-n_tmp1*num_group
        if(icolor+1.le.n_tmp2) then
        kpt_slda_dis(1)=icolor*(n_tmp1+1)+1
        kpt_slda_dis(2)=(icolor+1)*(n_tmp1+1)
        else
        kpt_slda_dis(1)=icolor*n_tmp1+n_tmp2+1
        kpt_slda_dis(2)=(icolor+1)*n_tmp1+n_tmp2
        endif
cccccccccccccccccccccccccccccccccccccccc
        n_tmp1=nkpt/num_group
        n_tmp2=nkpt-n_tmp1*num_group
        if(icolor+1.le.n_tmp2) then
        kpt_dis(1)=icolor*(n_tmp1+1)+1
        kpt_dis(2)=(icolor+1)*(n_tmp1+1)
        else
        kpt_dis(1)=icolor*n_tmp1+n_tmp2+1
        kpt_dis(2)=(icolor+1)*n_tmp1+n_tmp2
        endif
cccccccccccccccccccccccccccccccccccccccc
        n_tmp1=natom/num_group
        n_tmp2=natom-n_tmp1*num_group
        if(icolor+1.le.n_tmp2) then
        natom_dis(1)=icolor*(n_tmp1+1)+1
        natom_dis(2)=(icolor+1)*(n_tmp1+1)
        else
        natom_dis(1)=icolor*n_tmp1+n_tmp2+1
        natom_dis(2)=(icolor+1)*n_tmp1+n_tmp2
        endif
        
cccccccccccccccccccccccccccccccccccccccc
        nrL=n1L*n2L*n3L
        nr_nL=n1L*n2L*n3L/nnodes
cccccccccccccccccccccccccccccccccccccccc
cccc this is the correct estimation for n1,n2,n3 even for nonorthogonal AL(3,3)
ccccccccccc
       fackpt=2*dsqrt(2.d0*Ecut)/(4*datan(1.d0))
       dd1=fackpt*dsqrt(AL(1,1)**2+AL(2,1)**2+AL(3,1)**2)
       dd2=fackpt*dsqrt(AL(1,2)**2+AL(2,2)**2+AL(3,2)**2)
       dd3=fackpt*dsqrt(AL(1,3)**2+AL(2,3)**2+AL(3,3)**2)

       fackpt=dsqrt(2.d0*Ecut2)/(4*datan(1.d0))
       dd11=fackpt*dsqrt(AL(1,1)**2+AL(2,1)**2+AL(3,1)**2)
       dd22=fackpt*dsqrt(AL(1,2)**2+AL(2,2)**2+AL(3,2)**2)
       dd33=fackpt*dsqrt(AL(1,3)**2+AL(2,3)**2+AL(3,3)**2)

       fackpt=dsqrt(2.d0*Ecut2L)/(4*datan(1.d0))
       dd11L=fackpt*dsqrt(AL(1,1)**2+AL(2,1)**2+AL(3,1)**2)
       dd22L=fackpt*dsqrt(AL(1,2)**2+AL(2,2)**2+AL(3,2)**2)
       dd33L=fackpt*dsqrt(AL(1,3)**2+AL(2,3)**2+AL(3,3)**2)

       if(inode_tot.eq.1) then
       open(22,file='report')
       rewind(22)

       write(6,*) "************************"
       write(6,*) "  nnodes,num_group=", nnodes,num_group
       write(6,*) "************************"
       write(6,302)  dd1,dd2,dd3
       write(6,303)  dd11,dd22,dd33
       write(6,304)  n1,n2,n3
       write(6,3303) dd11L,dd22L,dd33L
       write(6,3304) n1L,n2L,n3L
       write(6,*) "************************"

       write(22,*) "***** ))) OUTPUT FILE FROM PETOT ((( *******"
       write(22,*) "for more inf. see file ****: etot.input"
       write(22,*) "atom config file       ****: ", f_xatom
       write(22,*) "************************"
       write(22,*) "  nnodes,num_group=", nnodes,num_group
       write(22,*) "************************"
       write(22,302)  dd1,dd2,dd3
       write(22,303)  dd11,dd22,dd33
       write(22,304)  n1,n2,n3
       write(22,3303) dd11L,dd22L,dd33L
       write(22,3304) n1L,n2L,n3L
       write(22,*) "************************"
       write(22,309) natom
       write(22,301) islda,igga
       write(22,305) Ecut*2,Ecut2*2,Ecut2L*2,Smth
       write(22,306) totNel,mx,tolug,tolE
       write(22,307) ilocal,rcut,ntype
       write(22,308) nkpt,nrot,num_mov,dtsart,dd_limit
       write(22,*) "*********************"
       write(22,*) "AL1,AL2,AL3 in (x,y,z)"
       write(22,133) AL(1,1),AL(2,1),AL(3,1)
       write(22,133) AL(1,2),AL(2,2),AL(3,2)
       write(22,133) AL(1,3),AL(2,3),AL(3,3)
       write(22,*) "************************"

       endif
302    format("recommended n1,n2,n3 from Ecut  ***: ",3(f6.1,1x))
303    format("recommended n1,n2,n3 from Ecut2 ***: ",3(f6.1,1x))
304    format("Actual      n1,n2,n3 used here  ***: ",3(i6,1x))
3303   format("recomm.   n1L,n2L,n3Lfrom Ecut2L***: ",3(f6.1,1x))
3304   format("Actual    n1L,n2L,n3L used here ***: ",3(i6,1x))
309    format(" natom=",i7)
301    format(" islda=",i7,    "     igga=", i7)
305    format("  Ecut=",f7.2,  "    Ecut2=",f7.2,"    Ecut2L=",f7.2,
     &   "   Smth=",f7.2)
306    format("totNel=",f7.2,  "       mx=",i7,  "     tolug=",E7.1,
     &   "   tolE=",E7.1)
307    format("ilocal=",i7,    "     rcut=",f7.2,"     ntype=",i7)
308    format("numkpt=",i7,    "  num-sym=",i7,  "  atom_mov=",i7,
     & " dtstart=",f9.4, " dd_limit=",f9.4)
133    format(3(f14.7,1x))

**********************************
ccccc  mx passed from input from param.escan_real
       mst=mx
        allocate(E_st(mst,nkpt,islda))
        allocate(err_st(mst,nkpt,islda))
**********************************
       Ek=0.5d0

       nr=n1*n2*n3
       nr_n = nr/nnodes

       mr=n1*n2*(n3+2)
       mr_n=mr/nnodes*(1+3.d0*nnodes/(n2*n3))  ! possible imbalance factor, some proc. has 2-3 more column
       mr=mr_n*nnodes
ccccc  (n3+2) is for the half plane in d3fft_real2
ccccc  fact (1+2.d0*nnodes/(n2*n3)) is for possible column load imbalance between diff. processors. 
ccccc  need to check the d3fft_real2 more carefully, and fftprep_real2.f

       nrL=n1L*n2L*n3L
       nr_nL=nrL/nnodes

       mrL=n1L*n2L*(n3L+2)
       mr_nL=mrL/nnodes*(1+3.d0*nnodes/(n2L*n3L))
       mrL=mr_nL*nnodes
ccccc  (n3L+2) is for the half plane in d3fft_real2L
ccccc  fact (1+2.d0*nnodes/(n2*n3)) is for possible column load imbalance between diff. processors. 

c     Calculate the approx. no. of g points, so we can dimension arrays
      volume=al(3,1)*(al(1,2)*al(2,3)-al(1,3)*al(2,2))
     &     +al(3,2)*(al(1,3)*al(2,1)-al(1,1)*al(2,3))
     &     +al(3,3)*(al(1,1)*al(2,2)-al(1,2)*al(2,1))
      volume = dabs(volume)
      delta_k=(2*pi)**3/volume
      totg=(4.0d0/3.0d0*pi*(dsqrt(2.0d0*Ecut))**3)/delta_k
      mg_nx=int(1.1d0*totg/nnodes)+100

c     Check that nr is a mulitple of the no. of nodes
      if (mod(nr,nnodes)/=0.or.mod(mr,nnodes)/=0) then
         write(6,*)'No. of grid points must be multiple of
     &        the no. of nodes. for both nr and mr'
         stop
      end if


      num_scf0=0
      do iter=1,niter0
      if(iscfmth0(iter).ne.0) num_scf0=num_scf0+1
      enddo
      num_scf1=0
      do iter=1,niter1
      if(iscfmth1(iter).ne.0) num_scf1=num_scf1+1
      enddo
      num_scf_max=num_scf0
      if(num_mov.gt.0.and.num_scf1.gt.num_scf0)
     &   num_scf_max=num_scf1

      if(num_scf_max.gt.0) then
      npulay_max=num_scf_max+2     ! the max number of pulay mixing
      else
      npulay_max=1
      endif

      if(ivr_out.eq.2.or.irho_out.eq.2) then
      ivr_rho_out=2 
      else
      ivr_rho_out=0
      endif

      if(ivr_rho_out.eq.2) then            ! just calculate vr_out, then stop
      npulay_max=1
      mx=1                             ! reduce the memory
      nref_tot=1
      endif
 

       call fft_allocate(n1,n2,n3,nnodes)
cccc  ncolx is defined inside fft_allocate
       call load_allocate(n1,n2,n3,ncolx,nnodes,mg_nx,mr_n)
       call data_allocate(mg_nx,mx,mr_n,mr_nL,ilocal,inode,nkpt,
     &   islda,natom,nref_tot,ipsp_all,nnodes,npulay_max)

       call fft_allocateL(n1L,n2L,n3L,nnodes,iflag_fft2L)
       call load_allocateL(n1L,n2L,n3L,ncolxL,nnodes,mr_nL,
     &       iflag_fft2L)


        allocate(workr_n(mr_n))

**********************************

       call get_ALI(AL,ALI)
       vol=AL(1,1)*(AL(2,2)*AL(3,3)-AL(3,2)*AL(2,3))+
     &        AL(2,1)*(AL(3,2)*AL(1,3)-AL(1,2)*AL(3,3))+
     &        AL(3,1)*(AL(1,2)*AL(2,3)-AL(2,2)*AL(1,3))
       vol=dabs(vol)

       call gen_G2_real()
       call fftprep_real2(n1,n2,n3,mr_n)

       if(iflag_fft2L.eq.1) then
       call gen_G2L_real()
       call fftprep_real2L(n1L,n2L,n3L,mr_nL)
       else
       gkk2_nL=gkk2_n
       gkx2_nL=gkx2_n
       gky2_nL=gky2_n
       gkz2_nL=gkz2_n
       iorg2L=iorg2
       ngtotnod2L=ngtotnod2
       ncolz2L=ncolz2
       ncolzL=ncolz2          ! ncolzL and ncolz2L are the same
       ixcol_z2L=ixcol_z2
       iycol_z2L=iycol_z2
       jjnode2L=jjnode2
       jjcol2L=jjcol2
       igstar_jjcol2L=igstar_jjcol2
       igfin_jjcol2L=igfin_jjcol2
       n1p2_nL=n1p2_n
       endif

       if(nrot.gt.1) then
       allocate(ig_star(mr_nL))       ! defined in data.f
       allocate(ig_star_stop(mr_nL))       ! defined in data.f
       allocate(ig_local(mr_nL))      ! not, mr_n/2, since -k might be in it.
       allocate(ig_lenstar(mr_nL))
       call symmcheck(smatr,AL,nrot,xatom,iatom)
       call gen_Gstar_ind(smatr,nrot,AL)
       endif

*********************************************
       if(icolor.eq.0) then    
       if(ivr_rho_out.ne.2) then    ! for ivr_out.eq.2, ug initi is not necessary
       do iislda=1,islda
       if(iwg_in.eq.1.and.inode==nnodes) then
       open(11,file=fwg_in(iislda),form='unformatted',
     &   status='old',action='read',iostat=ierr)
       if(ierr.ne.0) then
       write(6,*) "file ***",fwg_in(iislda),"*** does not exist, stop"
       call mpi_abort(MPI_COMM_WORLD,ierr)
       endif

       rewind(11)
       endif
 
       do kpt=1,nkpt
       call gen_G_comp(kpt,0)         ! need to be called for all icolor to get the kpt related quantities 
       call fftprep_comp(n1,n2,n3)
       iranm=-4513*(iislda-1)-1371*kpt-5616*inode

       call init_ug(AL,iwg_in,workr_n,kpt,iranm)
       call ugIO(ug_n,kpt,1,0,iislda)  ! iflag=1, write, iflag=2, read.

       enddo

       if(iwg_in.eq.1.and.inode==nnodes) close(11)
       enddo   ! iislda=1,islda
       endif   ! ivr_rho_out.ne.2

       else   ! icolor.eq.0,  only do above for icolor.eq.0 due to strange arrangement for iwg_in=1
       do kpt=1,nkpt
       call gen_G_comp(kpt,0)         ! need to for  all icolor to initialize the kpt related quantities 
       enddo
       endif   ! icolor.eq.0
*************************************************************
       
*********************************************************
       if(idens_out.eq.2.or.idens_out.eq.21.
     &           or.idens_out.eq.22)  then     ! output charge at the end of calculation

       if(kpt_dens(2).gt.nkpt.or.ispin_dens(2).gt.islda.or.
     &  iw_dens(2).gt.mx) then
       if(inode.eq.1) then
       write(6,*) "something wrong for output charge, skip"
       write(6,*) "kpt_dens(2).gt.nkpt", kpt_dens(2),nkpt
       write(6,*) "or:ispin_dens(2).gt.islda",ispin_dens(2),islda
       write(6,*) "or:iw_dens(2).gt.mx", iw_dens(2),mx
       endif
       goto 5000   
       endif

       if(kpt_dens(1).gt.kpt_dens(2).or.ispin_dens(1).gt.
     &  ispin_dens(2).or.iw_dens(1).gt.iw_dens(2)) then
       if(inode.eq.1) then
       write(6,*) "something wrong for output charge, skip"
       write(6,*) "kpt_dens(1).gt.kpt_dens(2)",kpt_dens(1),
     &  kpt_dens(2)
       write(6,*) "or:ispin_dens(1).gt.ispin_dens(2)",
     &    ispin_dens(1),ispin_dens(2)
       write(6,*) "or:iw_dens(1).gt.iw_dens(2)",
     &    iw_dens(1), iw_dens(2)
       endif
       goto 5000   
       endif

       if(idens_out.eq.2.and.icolor.eq.0) then
       call dens_out(AL,workr_n,kpt_dens,ispin_dens,iw_dens,
     &    fdens_out)
       endif
       if(idens_out.eq.21.or.idens_out.eq.22) then
       if(icolor.eq.0) then
       call densWr_out(AL,workr_n,kpt_dens,ispin_dens,iw_dens,
     &    fdens_out,idens_out)
       endif
       endif
       goto 5000     ! clean up files
       endif
*********************************************************

       call w_line(ilocal,ntype,vwr_atom,Ealpha)

       if(icoul.eq.0) goto 401
       if(inode_tot.eq.1) write(6,*) "enter getvcoul"

       if(icoul.ne.1.and.icoul.ne.11.and.icoul.ne.12.and.
     &  icoul.ne.13) then
       if(inode_tot.eq.1) write(6,*) 
     &       "input icoul not correct, stop", icoul
       call mpi_abort(MPI_COMM_K,ierr)
       endif

       if(icoul.eq.1)   then     ! cluster calculation, all three dimension increase x2
       n1L2=n1L*2
       n2L2=n2L*2
       n3L2=n3L*2
       AL2=AL*2
       endif
       if(icoul.eq.11) then      ! slab calculation along the first direction
       n1L2=n1L*2
       n2L2=n2L
       n3L2=n3L
       AL2(:,1)=2*AL(:,1)
       AL2(:,2)=AL(:,2)
       AL2(:,3)=AL(:,3)
       endif
       if(icoul.eq.12) then      ! slab calculation along the second direction
       n1L2=n1L
       n2L2=n2L*2
       n3L2=n3L
       AL2(:,1)=AL(:,1)
       AL2(:,2)=2*AL(:,2)
       AL2(:,3)=AL(:,3)
       endif
       if(icoul.eq.13) then      ! slab calculation along the third direction
       n1L2=n1L
       n2L2=n2L
       n3L2=n3L*2
       AL2(:,1)=AL(:,1)
       AL2(:,2)=AL(:,2)
       AL2(:,3)=2*AL(:,3)
       endif

       if(inode_tot.eq.1) then
       write(6,*) "n1L2,n2L2,n3L2=",n1L2,n2L2,n3L2
       endif

       nrL2=n1L2*n2L2*n3L2
       nr_nL2=nrL2/nnodes
       mrL2=n1L2*n2L2*(n3L2+2)
       mr_nL2=mrL2/nnodes*(1+3.d0*nnodes/(n2L2*n3L2))
       mrL2=mr_nL2*nnodes
ccccc  (n3L2+2) is for the half plane in d3fft_real2L2
ccccc  fact (1+3.d0*nnodes/(n2L2*n3L2)) is for possible column load imbalance between diff. processors. 

       vol2=AL2(1,1)*(AL2(2,2)*AL2(3,3)-AL2(3,2)*AL2(2,3))+
     &        AL2(2,1)*(AL2(3,2)*AL2(1,3)-AL2(1,2)*AL2(3,3))+
     &        AL2(3,1)*(AL2(1,2)*AL2(2,3)-AL2(2,2)*AL2(1,3))
       vol2=dabs(vol2)

       call fft_allocateL2(n1L2,n2L2,n3L2,nnodes)
       call load_allocateL2(n1L2,n2L2,n3L2,ncolxL2,
     &        nnodes,mr_nL2)
       call data_allocate_nL2(mr_nL2,ntype)
       call get_ALI(AL2,ALI2)
       call gen_G2L2_real()
       call fftprep_real2L2(n1L2,n2L2,n3L2,mr_nL2)
       call getvcoul(ntype,Ealpha,icoul) 

401    continue    ! jump out point for icoul
***************************************************************

	  Ealphat=0.d0   
	  ch=0.d0
	  do ia=1,natom
	  Ealphat=Ealphat+Ealpha(ityatom(ia))
	  ch=ch+zatom(ityatom(ia))
	  enddo
	  Ealphat=Ealphat*ch/vol

      
       if(irho_in.eq.1) then
       do iislda=1,islda
       call rhoIO(AL,rho_nL(1,iislda),mr_nL,2,
     &                 n1L,n2L,n3L,frho_in(iislda))     
       enddo
       ido_rho=0
       else
       ido_rho=1   
       endif 

       if(ivr_in.eq.1) then
       do iislda=1,islda
       call rhoIO(AL,vr_in_nL(1,iislda),mr_nL,2,
     &                   n1L,n2L,n3L,fvr_in(iislda))    
       enddo
       ido_vr=0
       else
       ido_vr=1
       endif

       if(ivext_in.eq.1) then
       call rhoIO(AL,vext_nL,mr_nL,2,
     &          n1L,n2L,n3L,fvext_in)
       else
       vext_nL=0.d0
       endif


       if(iforce.eq.0.and.num_mov.eq.0) then
       iforce_cal=0 
       else
       iforce_cal=1
       endif

cccccccc  ido_rho=0;1, initially, not calc. or calc. rho_n from atom
cccccccc  ido_vr=0;1,  initially, not calc. or calc. vr_in_n from rho_n
cccccccc  iforce_cal=0;1    not calc. or calc. force

cccccccc   one total energy calculation before moving the atoms

       if(inode.eq.1) then
c        tim=rtc()
        tim=system_time()
       endif 


       call mpi_barrier(MPI_COMM_WORLD,ierr)

       call Etotcalc(xatom,fatom,workr_n,Etot,
     &  iforce_cal,ido_rho,ido_vr,tolug,tolE,niter0,nline0,
     &  iCGmth0,iscfmth0,FermidE0,itypeFermi0,mCGbad0,E_st,err_st,AL,
     &  nkpt,ntype,convergE,islda,igga,iwg_out,fwg_out,ivr_rho_out)


cc-----------------------------------------------------------
       if(ivr_rho_out.eq.2.and.ivr_out.eq.2) then      ! output vr_in_nL, then stop
       do iislda=1,islda
       call rhoIO(AL,vr_in_nL(1,iislda),mr_nL,1,
     &                  n1L,n2L,n3L,fvr_out(iislda))
       enddo
       stop
       endif

       if(ivr_rho_out.eq.2.and.irho_out.eq.2) then    ! output rho_nL, then stop
       do iislda=1,islda
       call rhoIO(AL,rho_nL(1,iislda),mr_nL,1,
     &                 n1L,n2L,n3L,frho_out(iislda))    
       enddo
       stop
       endif
cc-----------------------------------------------------------


       mov=0
       call output_all()

       if(num_mov.eq.0) goto 2001    ! last report and end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc Now, begin the atomic movements
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       E_pred=0.d0
       dt=1.0d0*dtstart        ! first trial dt    
       dt_trial=dt
       dt_fact_nextlinemin=1.d0
       dd_max=0.d0
       ave_dx=0.d0
************************************************************
       if(inode_tot.eq.1) then
       write(6,*) "**** finished input atom config calc.  ***" 
       write(6,*) "**** following are atomic  relaxation  ***"
       write(22,*) "**"
       write(22,*) "**"
       write(22,996) 0, Etot
       write(22,*) "**** finished input atom config calc.  ***" 
       write(22,*) "**** following are atomic  relaxation  ***"
       write(22,*) "**"
       write(22,*) "**"
       endif

       do 2000 mov=1,num_mov       ! the big do loop for atomic movement

       if(inode_tot.eq.1) then
       write(6,*) "** atomic relaxation, atom_mov_step: ", mov
       write(22,*) "============================================="
       write(22,*) "** atomic relaxation, atom_mov_step: ", mov
       write(22,*) "============================================="
       endif
******************** update some quantities
       Etot_old=Etot      ! the _old quantities correspond to xatom_oldx
       convergE_old=convergE
******************** check whether to stop for small remaining force
         force_max=-100.d0
         do i=1,natom
         if(dabs(fatom(1,i)*imov_at(1,i)).gt.force_max)
     &    force_max=dabs(fatom(1,i)*imov_at(1,i))
         if(dabs(fatom(2,i)*imov_at(2,i)).gt.force_max)
     &    force_max=dabs(fatom(2,i)*imov_at(2,i))
         if(dabs(fatom(3,i)*imov_at(3,i)).gt.force_max)
     &    force_max=dabs(fatom(3,i)*imov_at(3,i))
         enddo

         if(force_max.lt.tolforce) then
             if(inode_tot.eq.1) then
             write(22,*) "force_max.lt.tolforce, finished",
     &                 force_max,tolforce
             endif
         goto 2001  
         endif


**************************** store the ug and rho in file100xx
       do iislda=1,islda
       do kpt=1,nkpt
        if((iislda-1)*nkpt+kpt.ge.kpt_slda_dis(1).and.
     &     (iislda-1)*nkpt+kpt.le.kpt_slda_dis(2)) then
         call ugIO(ug_n,kpt,2,0,iislda)   ! read out the ug_n from Etotcalc
         call ugIO(ug_n,kpt,1,1,iislda)   ! store it in istep=1 file
         endif
       enddo
       enddo

       do iislda=1,islda
       if(iislda.eq.1) f_tmp="wrhofile1000"
       if(iislda.eq.2) f_tmp="wrhofile2000"
       call rhoIO(AL,rho_nL(1,iislda),mr_nL,1,n1L,n2L,n3L,f_tmp)    ! write down the charge density
************************************************************
***** There is some jittering due to the fact that, we are using charge 
***** in the interpolation of atom movement, but using potential inside Etotcalc
***** for scf potential mixing. As a result, we might lost 1-2 iterations, because
***** the charge here, is not the properly mixed charge, but the bare charge from
***** the output wavefunction. 
******
****** But TEST shows that vr_in_n interpolating is much worse than rho_n interp.
***** So we need keep using rho_n interpolating.  The problem is that we 
***** don't have Pulay mixing on charge (due to negative charge problem). 
***** Thus, we have to sacrify the 1-2 iteration lose. 
***** For most time, when dt_trial.ne.dt_final, this is not really a problem. 
***** Lin-Wang Wang, Dec.20,2000
**************************************************************
       enddo
*********************************************************
*************** first, along the force direction, move a trial step dt

        line_step=0
        call atomMV(mov,line_step,xatom,fatom,Etot,E_pred,
     &  dt,istop_line,dd_max,dd_limit,AL,imov_at,convergE,
     &  message,dt_fact_nextlinemin,xatom_old,imv_cont)
ccccc at this atomMVm line_step, xatom_old is determined
********* input: mov, line_step; 
*********  output: line_step=line_step+1, istop_line(decision, 1:stop, 0, not stop)

        call mpi_barrier(MPI_COMM_WORLD,ierr)
cc        call mpi_barrier(MPI_COMM_WORD,ierr)      ! a old bug in previous release
****************************************************************
         if(inode_tot.eq.1) then
         write(22,*) "****************************************"
         write(22,*) "**** within atom_mov_step:", mov
         write(22,*) "****  trial line_min_step:",line_step 
         write(22,994) dt,dd_max,force_max
         write(22,*) "****************************************"
         call system_flush(22)
         write(6,*) "**** within atom_mov_step:", mov
         write(6,*) "****  trial line_min_step:",line_step 
         write(6,994) dt,dd_max,force_max
         endif

994      format("  dt= ",E10.4,",  dd_max=",E10.4,",  force_max=",E10.4)
****************************************************************
****************************************************************
****** do move_rho before 1000, for returning cases, the rho has already been updated.
******  move_rho correct the current rho, (correspond to xatom_old), to rho (correspond to xatom)
******  The resulting rho might be negative in some region, but that is okay for Etotcalc, and getpotX.f
******  The Coulomb energy is calculated correctly, and the dabs(rho+rho_c) is used for Uxc.
******  We should not make rho_n positive here, that could lost local charge balance.
        call move_rho(1)      
1000    continue


        call mpi_barrier(MPI_COMM_WORLD,ierr)
       ido_rho=0
       ido_vr=1        ! generate the potential from the charge density
       iforce_cal=0
       if(dabs(dd_max/dd_limit-1).lt.0.0001) iforce_cal=1    ! possiblly for one Etotcalc calc.
       call Etotcalc(xatom,fatom,workr_n,Etot,
     &  iforce_cal,ido_rho,ido_vr,tolug,tolE,niter1,nline1,
     &  iCGmth1,iscfmth1,FermidE1,itypeFermi1,mCGbad1,E_st,err_st,AL,
     &  nkpt,ntype,convergE,islda,igga,iwg_out,fwg_out,ivr_rho_out)

        call mpi_barrier(MPI_COMM_WORLD,ierr)

        dt_trial=dt  
        call atomMV(mov,line_step,xatom,fatom,Etot,E_pred,
     &  dt,istop_line,dd_max,dd_limit,AL,imov_at,convergE,
     &  message,dt_fact_nextlinemin,xatom_old,imv_cont)
ccccc xatom_old is not changed by this atomMV step

****************************************************************
         if(inode_tot.eq.1) then
         write(22,*) "****************************************"
         write(22,*) "**** within atom_mov_step:", mov
         write(22,*) "****  trial line_min_step:",line_step 
         write(22,994) dt,dd_max,force_max
         write(22,*) "****************************************"
         call system_flush(22)

         write(6,*) "**** within atom_mov_step:", mov
         write(6,*) "****  trial line_min_step:",line_step 
         write(6,994) dt,dd_max,force_max

         endif
****************************************************************
*****   the situations not to have a interpolation
       if(dabs(dt/dt_trial-1).lt.0.001d0.and.
     &  dabs(dd_max/dd_limit-1.d0).lt.0.0001) goto 5001

        interp=1
        if(dt/dt_trial.lt.0.2d0) interp=0      ! use the dt=0 values for ug, move atom for rho
        if(dt/dt_trial.gt.2.5d0) interp=0      ! don't use interp for large extrapolation
        if(convergE.gt.100*tolE.and.convergE.gt.
     &      100*convergE_old) interp=0
        if(dd_max.gt.0.15) interp=0
        if(dabs(dt/dt_trial-1.d0).lt.0.02) interp=1    ! for many special case where dt=dt_trial
                
                
        if(inode_tot.eq.1) then
        write(6,*) "**** interpolation 1:yes; 0:no =", interp
        write(22,*) "**** interpolation 1:yes; 0:no =", interp
        endif
                
       if(interp.eq.1) then
       call interpolation()
       else
       call move_rho(2)     ! ug_n not interpolated, but using the xatom_old values
       endif

****************************************************************
        if(istop_line.eq.0) then
         if(inode_tot.eq.1) then
       write(22,*) "** WARNING: line_step>2, line_min not work well**"
       write(6,*) "** WARNING: line_step>2, line_min not work well**"
         call system_flush(22)
         endif
        goto 1000
        endif
****************************************************************

       ido_rho=0
       ido_vr=1
       iforce_cal=1

        call mpi_barrier(MPI_COMM_WORLD,ierr)

       call Etotcalc(xatom,fatom,workr_n,Etot,
     &  iforce_cal,ido_rho,ido_vr,tolug,tolE,niter1,nline1,
     &  iCGmth1,iscfmth1,FermidE1,itypeFermi1,mCGbad1,E_st,err_st,AL,
     &  nkpt,ntype,convergE,islda,igga,iwg_out,fwg_out,ivr_rho_out)

5001   continue       ! skip the above Etotcalc
**********************************************************
************* calculate force_max just for the report
         force_max=-100.d0
         do i=1,natom
         if(dabs(fatom(1,i)*imov_at(1,i)).gt.force_max)
     &    force_max=dabs(fatom(1,i)*imov_at(1,i))
         if(dabs(fatom(2,i)*imov_at(2,i)).gt.force_max)
     &    force_max=dabs(fatom(2,i)*imov_at(2,i))
         if(dabs(fatom(3,i)*imov_at(3,i)).gt.force_max)
     &    force_max=dabs(fatom(3,i)*imov_at(3,i))
         enddo
****************************************************************

         if(inode_tot.eq.1) then
         write(22,*) "================================="
         write(22,*) "**** Final result for atom_mov_step: ",mov
         write(22,*) "**** lin_min message: ", message
         write(22,997) dt,line_step,dt_fact_nextlinemin
         write(22,998) dd_max,force_max
         write(22,995) Etot_old,E_pred,Etot
         write(22,996) mov, Etot
         write(22,*) "================================="
         call system_flush(22)

995      format(" **** E_prev,E_pred,E_final", 3x, 3(E14.8,1x))
997      format(" ****      dt= ", E9.3, ",   line_step= ",i3,
     &    ",  dt_incr_next_step= ", E9.3)
998      format(" ****  dd_max= ", E9.3, ",   max_force= ",E9.3)
996      format(" RESULT: atom_move_step, E_tot: ", i4,1x,E25.15)

         write(6,*) "================================="
         write(6,*) "**** Final energy for atom_mov_step: ",mov
         write(6,*) "**** lin_min message: ", message
         write(6,997) dt,line_step,dt_fact_nextlinemin
         write(6,995) Etot_old,E_pred,Etot
         write(6,996) mov, Etot
         write(6,*) "================================="

         endif
**********************************************************
********  output everything after each atomic movements

        if(Etot.gt.Etot_old) then
              if(inode_tot.eq.1) then
              write(22,*) "========================="
              write(22,*) "Etot.gt.Etot_old, stop atom move"
              write(22,*) "stopped before the output of this step"
              write(22,*) "========================="
              endif
         stop
         endif

        call output_all()
**********************************************************

         if(dabs(Etot-Etot_old).lt.max(convergE,tolE)) then
          if(inode_tot.eq.1) then
        write(22,*) 
     & "Stop atom_mov, because E-E_old.lt.(tolE or E_converg_err)" 
        write(22,999)  dabs(Etot-Etot_old),tolE,convergE
999     format("|Etot-Etot_old|, tolE, E_converg_err: ",3(E12.4,1x))
          endif
         goto 2001       ! finished
         endif

          
**********************************************************
2000     continue  
****************************************************************
*****   end of the atomic movement
****************************************************************


2001   continue     ! exit point

*******************************************************
       if(idens_out.eq.1.or.idens_out.eq.11.
     &           or.idens_out.eq.12)  then     ! output charge at the end of calculation

       if(kpt_dens(2).gt.nkpt.or.ispin_dens(2).gt.islda.or.
     &  iw_dens(2).gt.mx) then
       if(inode_tot.eq.1) then
       write(6,*) "something wrong for output charge, skip" 
       write(6,*) "kpt_dens(2).gt.nkpt", kpt_dens(2),nkpt
       write(6,*) "or:ispin_dens(2).gt.islda",ispin_dens(2),islda
       write(6,*) "or:iw_dens(2).gt.mx", iw_dens(2),mx
       endif
       goto 912
       endif

       if(kpt_dens(1).gt.kpt_dens(2).or.ispin_dens(1).gt.
     &  ispin_dens(2).or.iw_dens(1).gt.iw_dens(2)) then
       if(inode_tot.eq.1) then
       write(6,*) "something wrong for output charge, skip"
       write(6,*) "kpt_dens(1).gt.kpt_dens(2)",kpt_dens(1),
     &  kpt_dens(2)
       write(6,*) "or:ispin_dens(1).gt.ispin_dens(2)",
     &    ispin_dens(1),ispin_dens(2)
       write(6,*) "or:iw_dens(1).gt.iw_dens(2)",
     &    iw_dens(1), iw_dens(2)
       endif
       goto 912
       endif

       if(idens_out.eq.1) then
       call dens_out(AL,workr_n,kpt_dens,ispin_dens,iw_dens,
     &    fdens_out)
       endif

       if(idens_out.eq.11.or.idens_out.eq.12) then
       call densWr_out(AL,workr_n,kpt_dens,ispin_dens,iw_dens,
     &    fdens_out,idens_out)
       endif

912    continue

       endif
*******************************************************

       if(inode_tot.eq.1) then
        tim=(system_time()-tim)*2.222D-9
       write(6,*) "computational time=", tim 
       endif

       if(inode_tot.eq.1) then
       write(22,*) "computational time=", tim
       write(22,*) "last report from diag_real"
       write(22,*) "*********************************"
       write(22,*) 
     & "Eigen energies are values after setting ave(Vtot)=0"
       write(22,*) 
     & "For Vtot=V_ion+V_Hartree+V_xc, and"
       write(22,*) 
     & "ave(V_ion+V_Hatree)=0, ave(V_xc).ne.0:  E=E+v0"
       do iislda=1,islda
       do kpt=1,nkpt
       write(22,*) "*********************************"
       write(22,104) iislda, kpt
       write(22,*) "err of each states, A.U"
       write(22,102) (err_st(i,kpt,iislda), i=1,mx)
       write(22,*) "eigen energies, in eV "
       write(22,103) (E_st(i,kpt,iislda)*27.211396d0, i=1,mx)
       write(22,*) "*********************************"
       enddo
       enddo
       endif

102   format(5(E10.4,3x))
103   format(5(f12.8,1x))
104   format("iislda,kpt= ",i3,", ", i4)

5000  continue
      
      call deletefile()
      call data_deallocate()
      call load_deallocate()
      call fft_deallocate()
      deallocate(workr_n)

      if(icoul.eq.1) then
      call data_deallocate_nL2()
      endif

      call mpi_finalize(ierr)

      stop

      contains
      
*****************************************************
*****************************************************

      subroutine output_all()

      implicit double precision (a-h,o-z)


      if(iwg_out.eq.1.or.iwg_out.eq.2) then
      if(icolor.eq.0) then
      call write_wg(fwg_out,AL,islda,nkpt)
      endif
      endif

***********************************************

       if(iforce_cal.eq.1) then
       if(inode_tot.eq.1) then
       open(13,file=fforce_out)
       rewind(13)
       fx=0.d0
       fy=0.d0
       fz=0.d0
       do i=1,natom
       write(13,337) iatom(i),fatom(1,i),fatom(2,i),fatom(3,i)
       fx=fx+fatom(1,i)
       fy=fy+fatom(2,i)
       fz=fz+fatom(3,i)
       enddo
       write(13,*) "****** total force ******************"
       ia=0
       write(13,337) ia, fx,fy,fz
       write(13,*) "**** atom_move_step: ", mov
       close(13)
c337    format(i4,3x,3(E11.5,1x))
337    format(i4,3x,3(E16.10,1x))
       endif
       endif


******************************************************
       if(ivr_out.eq.1) then
       do iislda=1,islda
       call rhoIO(AL,vr_in_nL(1,iislda),mr_nL,1,
     &                     n1L,n2L,n3L,fvr_out(iislda))      ! rhoIO should be called by all icolor groups
       enddo
       endif

       if(irho_out.eq.1) then
       do iislda=1,islda
       call rhoIO(AL,rho_nL(1,iislda),mr_nL,1,
     &                 n1L,n2L,n3L,frho_out(iislda))    
       enddo
       endif

       if(num_mov.gt.0) then
       if(inode_tot.eq.1) then
       open(11,file=fxatom_out)
       rewind(11)
       write(11,322) natom,mov
       write(11,320) AL(1,1),AL(2,1),AL(3,1)
       write(11,320) AL(1,2),AL(2,2),AL(3,2)
       write(11,320) AL(1,3),AL(2,3),AL(3,3)
       do ia=1,natom
       write(11,321) iatom(ia),xatom(1,ia),xatom(2,ia),xatom(3,ia),
     &  imov_at(1,ia),imov_at(2,ia),imov_at(3,ia)
       enddo
       close(11)
320    format(3(E15.7,1x))
321    format(i4,2x,3(f14.9,1x),2x,3(i2,1x))
322    format(i6,"     | atom_move_step: ", i4)
       endif
       endif

      return
      end subroutine output_all

*****************************************************
      subroutine interpolation()
      implicit double precision (a-h,o-z)

      complex*16 cc

      fac1=1.d0-dt/dt_trial
      fac2=dt/dt_trial

       if(inode_tot.eq.1) then
       write(6,*) "extrapolating,fac1,fac2",fac1,fac2 
       write(22,*) "extrapolating,fac1,fac2",fac1,fac2
       endif

      if(dabs(fac1).gt.20.d0)  return      ! don't do extreme extrapolation 

      allocate(ug_n_tmp(mg_nx,mx))
      do iislda=1,islda
      do kpt=1,nkpt

      if((iislda-1)*nkpt+kpt.ge.kpt_slda_dis(1).and.
     &     (iislda-1)*nkpt+kpt.le.kpt_slda_dis(2)) then

      call ugIO(ug_n_tmp,kpt,2,1,iislda)
      call ugIO(ug_n,kpt,2,0,iislda)

      do m=1,mx
      cc=dcmplx(0.d0,0.d0)
      do i=1,ngtotnod(inode,kpt)
      cc=cc+dconjg(ug_n_tmp(i,m))*ug_n(i,m)
      enddo
      call global_sumc(cc)
      cc=cc*vol
      cc=cc/cdabs(cc)

      do i=1,ngtotnod(inode,kpt)
      ug_n(i,m)=cc*ug_n_tmp(i,m)*fac1+ug_n(i,m)*fac2     !  phase factors corrected
      enddo
      enddo

      call ugIO(ug_n,kpt,1,0,iislda)

      endif      ! kpt_slda_dis(1)(2)

      enddo
      enddo


      deallocate(ug_n_tmp)

      allocate(workr_n_tmp(mr_nL))

      ss=0.d0
      do 300 iislda=1,islda

      if(iislda.eq.1)  f_tmp="wrhofile1000"
      if(iislda.eq.2)  f_tmp="wrhofile2000"
      call rhoIO(AL,workr_n_tmp,mr_nL,2,n1L,n2L,n3L,f_tmp)   ! rhoIO should be called by all icolor group   

       do i=1,nr_nL
       rho_nL(i,iislda)=dabs(workr_n_tmp(i)*fac1+
     &                         rho_nL(i,iislda)*fac2)
       ss=ss+rho_nL(i,iislda)
       enddo
       call mpi_barrier(MPI_COMM_K,ierr)
300    continue

       ss=ss*vol/(n1L*n2L*n3L)
       call global_sumr(ss)

       ss=totNel/ss

       do iislda=1,islda
       do i=1,nr_nL
       rho_nL(i,iislda)=ss*rho_nL(i,iislda)
       enddo
       enddo

       deallocate(workr_n_tmp)

       return
       end subroutine interpolation


       subroutine move_rho(iflag)
       implicit double precision(a-h,o-z)
       integer iflag
ccccccccc   This subroutine move the atomic electron charge with the atom
***** iflag.eq.1, current rho_n correspond to xatom_old
***** iflag.eq.2, rho_n in wrhofile1000 corrspond to xatom_old, also using ug_n of xatom_old as current ug_n
*****  The resulting rho_nL could be negative at some points, but that is okay to calc. potential.
                                                                                                         
        allocate(rhonew1_nL(mr_nL))
                                                                                                         
          if(iflag.eq.2) then
        f_tmp="wrhofile1000"
        call rhoIO(AL,rho_nL(1,1),mr_nL,2,n1L,n2L,n3L,f_tmp)
        if(islda.eq.2) then
        f_tmp="wrhofile2000"
        call rhoIO(AL,rho_nL(1,2),mr_nL,2,n1L,n2L,n3L,f_tmp)
        endif
          endif
               
        if(islda.eq.1) then
        call getrho_only(AL,xatom_old,ntype,iatom,totNel,rhonew1_nL)   ! calc. rho from xatom
        do i=1,nr_nL
        rho_nL(i,1)=rho_nL(i,1)-rhonew1_nL(i)
        enddo
        call getrho_only(AL,xatom,ntype,iatom,totNel,rhonew1_nL)
        do i=1,nr_nL
        rho_nL(i,1)=rho_nL(i,1)+rhonew1_nL(i)
        enddo
        endif

        if(islda.eq.2) then
        allocate(rhonew2_nL(mr_nL))
        call getrho_only(AL,xatom_old,ntype,iatom,totNel,rhonew1_nL)
        call getrho_only(AL,xatom,ntype,iatom,totNel,rhonew2_nL)
        do i=1,nr_nL
        f1=rho_nL(i,1)/(rho_nL(i,1)+rho_nL(i,2))
        f2=1.d0-f1
        rho_nL(i,1)=rho_nL(i,1)+(rhonew2_nL(i)-rhonew1_nL(i))*f1
        rho_nL(i,2)=rho_nL(i,2)+(rhonew2_nL(i)-rhonew1_nL(i))*f2
        enddo
        deallocate(rhonew2_nL)
        endif
                             
             
        deallocate(rhonew1_nL)
                             
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if(iflag.eq.2) then
         do iislda=1,islda
         do kpt=1,nkpt
          if((iislda-1)*nkpt+kpt.ge.kpt_slda_dis(1).and.
     &     (iislda-1)*nkpt+kpt.le.kpt_slda_dis(2)) then
         call ugIO(ug_n,kpt,2,1,iislda)    ! read out the xatom_old values
         call ugIO(ug_n,kpt,1,0,iislda)    ! write down as the current ug_n values
           endif
         enddo
         enddo
        endif
              
        return
        end subroutine move_rho



ccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine deletefile()
       implicit double precision(a-h,o-z)
       character*8 fname
       character*9 fname9


       if(inode_tot.eq.1) then

       do istep=0,1
       do iislda=1,islda
       do kpt=1,nkpt

       fname="ugiofile"
       kpt1=mod(kpt,10)
       kpt2=mod((kpt-kpt1)/10,10)
       kpt3=mod((kpt-kpt1-kpt2*10)/100,10)
       kpt4=mod((kpt-kpt1-kpt2*10-kpt3*100)/1000,10)

       open(10,file=
     & fname//char(iislda+48)//char(istep+48)//char(kpt4+48)//
     &  char(kpt3+48)//char(kpt2+48)//char(kpt1+48),
     &  form="unformatted")
       close(10,status='delete')

       if(ipsp_all.eq.2) then
       fname9="bpsiiofil"
       kpt1=mod(kpt,10)
       kpt2=mod((kpt-kpt1)/10,10)
       kpt3=mod((kpt-kpt1-kpt2*10)/100,10)
       kpt4=mod((kpt-kpt1-kpt2*10-kpt3*100)/1000,10)

       open(10,file=
     &fname9//char(iislda+48)//char(istep+48)//char(kpt4+48)//
     &  char(kpt3+48)//char(kpt2+48)//char(kpt1+48),
     &       form="unformatted")
       close(10,status='delete')
       endif

       enddo
       enddo
       enddo


       if(ilocal.eq.3) then
       do kpt=1,nkpt

       fname="wqiofile"
       kpt1=mod(kpt,10)
       kpt2=mod((kpt-kpt1)/10,10)
       kpt3=mod((kpt-kpt1-kpt2*10)/100,10)
       kpt4=mod((kpt-kpt1-kpt2*10-kpt3*100)/1000,10)

       open(10,file=
     & fname//char(kpt4+48)//char(kpt3+48)//
     &  char(kpt2+48)//char(kpt1+48),
     &       form="unformatted")
       close(10,status='delete')
       enddo
       endif

       open(10,file='wrhofile1000')
       close(10,status='delete')
       open(10,file='wrhofile2000')
       close(10,status='delete')

       endif

       return
       end subroutine deletefile
       
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

        end
      


