      subroutine input(AL,ilocal,tolug,tolE,
     &  niter0,niter1,nline0,nline1,iCGmth0,iCGmth1,
     &  iscfmth0,iscfmth1,FermidE0,FermidE1,mCGbad0,mCGbad1,
     &  islda,igga,totNel,ntype,xatom,iatom,fxatom_out,
     &  iwg_in,fwg_in,iwg_out,fwg_out,irho_in,frho_in,
     &  irho_out,frho_out,ivr_in,fvr_in,ivr_out,fvr_out,
     &  iforce,fforce_out,vwr_atom,nkpt,smatr,nrot,
     &  num_mov,tolforce,imov_at,dtstart,
     &  idens_out,kpt_dens,ispin_dens,iw_i,iw_f,fdens_out,
     &  f_xatom)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
******************************************

****************************************
****  It stores the wavefunction in G space, only in half
****  of the E_cut sphere. 
******************************************

       use fft_data
       use load_data
       use data

      implicit double precision (a-h,o-z)

      include 'mpif.h'

      include 'param.escan_real'
******************************************
       real*8 AL(3,3),AL_t(3,3)
***********************************************
****  single precision for vr(i)
***********************************************
*************************************************
       real*8 xatom(3,matom)
       integer iatom(matom),imov_at(3,matom)
       integer smatr(3,3,48),nrot
       integer iCGmth0(100),iCGmth1(100),iscfmth0(100),iscfmth1(100)
       real*8 FermidE0(100),FermidE1(100),totNel
       integer niter0,nline0,niter1,nline1,mCGbad0,mCGbad1

**************************************************
       character*20 vwr_atom(mtype),fforce_out,fdens_out
       character*20 fwg_in(2),fwg_out(2),frho_in(2),frho_out(2),
     & fvr_in(2),fvr_out(2),f_tmp,fxatom_out

       character*20 f_xatom,sym_file,kpt_file

**************************************************

       open(9,file='etot.input',status='old',action='read',iostat=ierr) 
       if(ierr.ne.0) then
       if(inode.eq.1)  
     &  write(6,*) "file ***etot.input*** does not exist, stop"
       stop
       endif

       read(9,*,iostat=ierr) i1, f_xatom
       if(ierr.ne.0) call error_stop(1)
       call readf_xatom()
       read(9,*,iostat=ierr) i1, n1,n2,n3
       if(ierr.ne.0) call error_stop(2)
       read(9,*,iostat=ierr) i1, islda,igga
       if(ierr.ne.0) call error_stop(3)
       read(9,*,iostat=ierr) i1, Ecut,Ecut2,Smth
       if(ierr.ne.0) call error_stop(4)
       if(islda.ne.1.and.islda.ne.2.and.inode.eq.1) then
       write(6,*) "islda must be 1 (lda) or 2 (slda), stop", islda
       call mpi_abort(MPI_COMM_WORLD,ierr)
       endif
       if(igga.ne.0.and.igga.ne.1.and.inode.eq.1) then
       write(6,*) "igga must be 0 (no gga) or 1 (gga), stop", igga
       call mpi_abort(MPI_COMM_WORLD,ierr)
       endif
       read(9,*,iostat=ierr) i1, iwg_in,(fwg_in(i),i=1,islda)
       if(ierr.ne.0) call error_stop(5)
       read(9,*,iostat=ierr) i1, iwg_out,(fwg_out(i),i=1,islda)
       if(ierr.ne.0) call error_stop(6)
       read(9,*,iostat=ierr) i1, irho_in,(frho_in(i),i=1,islda)
       if(ierr.ne.0) call error_stop(7)
       read(9,*,iostat=ierr) i1, irho_out,(frho_out(i),i=1,islda)
       if(ierr.ne.0) call error_stop(8)
       read(9,*,iostat=ierr) i1, ivr_in,(fvr_in(i),i=1,islda)
       if(ierr.ne.0) call error_stop(9)
       read(9,*,iostat=ierr) i1, ivr_out,(fvr_out(i),i=1,islda)
       if(ierr.ne.0) call error_stop(10)
       read(9,*,iostat=ierr) i1,idens_out,kpt_dens,ispin_dens,
     &                    iw_i,iw_f,fdens_out
       if(ierr.ne.0) call error_stop(11)
       read(9,*,iostat=ierr) i1, iforce,fforce_out
       if(ierr.ne.0) call error_stop(12)
       read(9,*,iostat=ierr) i1, isym,sym_file
       if(ierr.ne.0) call error_stop(13)
       read(9,*,iostat=ierr) i1, ikpt_yno,kpt_file
       if(ierr.ne.0) call error_stop(14)
       call readkpt()
         
       read(9,*,iostat=ierr) i1, totNel,mx,tolug,tolE
       if(ierr.ne.0) call error_stop(15)
       read(9,*,iostat=ierr) i1, niter0,nline0,mCGbad0
       if(ierr.ne.0) call error_stop(16)
         do i=1,niter0
         read(9,*,iostat=ierr) iCGmth0(i),iscfmth0(i),FermidE0(i)
       if(ierr.ne.0) call error_stop(16)
         FermidE0(i)=FermidE0(i)/27.211396d0
	 enddo
       read(9,*,iostat=ierr) i1, num_mov,tolforce,dtstart,fxatom_out
       if(ierr.ne.0) call error_stop(17)
       read(9,*,iostat=ierr) i1, niter1,nline1,mCGbad1
       if(ierr.ne.0) call error_stop(18)
         do i=1,niter1
         read(9,*,iostat=ierr) iCGmth1(i),iscfmth1(i),FermidE1(i)
       if(ierr.ne.0) call error_stop(18)
         FermidE1(i)=FermidE1(i)/27.211396d0
	 enddo
       read(9,*,iostat=ierr) i1, ilocal
       if(ierr.ne.0) call error_stop(19)
       read(9,*,iostat=ierr) i1, rcut
       if(ierr.ne.0) call error_stop(20)
       read(9,*,iostat=ierr) i1, ntype
       if(ierr.ne.0) call error_stop(21)
       do ia=1,ntype
       read(9,*,iostat=ierr) i1, vwr_atom(ia)
       if(ierr.ne.0) call error_stop(21+ia)
       enddo
       close(9)

*************************************************
**** change Ecut,Eref to A.U
*************************************************
       Ecut=Ecut/2
       Ecut2=Ecut2/2
*************************************************
           if(isym.eq.0) then
           nrot=1
           else
       open(10,file=sym_file,status='old',action='read',iostat=ierr)
       if(ierr.ne.0) then
       if(inode.eq.1)
     &   write(6,*) "file ***", sym_file, "*** does not exist, stop"
       stop
       endif
           rewind(10)
           read(10,*) nrot
           do irot=1,nrot
           read(10,*)
           do j=1,3
           read(10,*) smatr(1,j,irot),smatr(2,j,irot),smatr(3,j,irot)
           enddo
           enddo
           close(10)
           endif
*************************************************

       nr=n1*n2*n3

       if(inode.eq.1) then
        write(6,*) "number of nonlocal atom=", natom
        write(6,*) "n1,n2,n3=", n1,n2,n3
        write(6,*) "islda,igga=", islda, igga
        write(6,*) "AL1,AL2,AL3 in (x,y,z) components and A.U "
        write(6,*) "each line is one supercell edge vector"
        write(6,3) AL(1,1),AL(2,1),AL(3,1)
        write(6,3) AL(1,2),AL(2,2),AL(3,2)
        write(6,3) AL(1,3),AL(2,3),AL(3,3)
       endif

3      format(3(f13.7,1x))
*************************************************
       nh1=n1/2+1

*************************************************
      return
      contains

*******************************************
      subroutine readf_xatom()
       open(10,file=f_xatom,status='old',action='read',iostat=ierr)
       if(ierr.ne.0) then
       if(inode.eq.1) 
     &  write(6,*) "file ***",f_xatom, " ***does not exist, stop"
       stop
       endif
       
       rewind(10)
       read(10,*) natom
      
       if(natom.gt.matom) then
       if(inode.eq.1) then
       write(6,*) "natom.gt.matom, increase matom in data.f, stop"
       endif
       stop
       endif

       read(10,*) (AL(i,1),i=1,3)
       read(10,*) (AL(i,2),i=1,3)
       read(10,*) (AL(i,3),i=1,3)
       do i=1,natom
       read(10,*) iatom(i),xatom(1,i),xatom(2,i),xatom(3,i),
     & imov_at(1,i),imov_at(2,i),imov_at(3,i)
       enddo
       close(10)
      return
      end subroutine readf_xatom
**************************************************


      subroutine readkpt()

      implicit double precision (a-h,o-z)
         real*8 AL_t(3,3),tmp(3)

       if(ikpt_yno.eq.0)  then 
       nkpt=1
       else
       open(12,file=kpt_file,status='old',action='read',iostat=ierr)
       if(ierr.ne.0) then
       if(inode.eq.1) 
     &  write(6,*) "file ***",kpt_file, "*** does not exist, stop"
       stop
       endif
       rewind(12)
       read(12,*) nkpt
       endif

       call data_allocate_akx(nkpt)
       call ngtot_allocate(nnodes,nkpt)

       if(ikpt_yno.eq.0) then
       akx(1)=0.d0
       aky(1)=0.d0
       akz(1)=0.d0
       weighkpt(1)=1.d0
       return
       else
       
         pi=4*datan(1.d0)

         read(12,*) iflag, ALxyz

         if(iflag.eq.1)  then
         if(inode.eq.1) write(6,*) "input kpts in Cartesian Coord"
         sumw=0.d0
         do kpt=1,nkpt
         read(12,*) akx(kpt),aky(kpt),akz(kpt),weighkpt(kpt)
         sumw=sumw+weighkpt(kpt)
         akx(kpt)=akx(kpt)*2*pi/ALxyz
         aky(kpt)=aky(kpt)*2*pi/ALxyz
         akz(kpt)=akz(kpt)*2*pi/ALxyz
         enddo
         endif

         if(iflag.eq.2) then
         if(inode.eq.1) write(6,*) "input kpts in primary cell unit"

         do i=1,3
         do j=1,3
         AL_t(j,i)=AL(i,j)
         enddo
         tmp(i)=1
         enddo
         call gaussj(AL_t,3,3,tmp,1,1)

         sumw=0.d0
         do kpt=1,nkpt
         read(12,*) ak1_t,ak2_t,ak3_t,weighkpt(kpt)
         sumw=sumw+weighkpt(kpt)
         akx(kpt)=2*pi*(AL_t(1,1)*ak1_t+AL_t(1,2)*ak2_t+
     &           AL_t(1,3)*ak3_t)
         aky(kpt)=2*pi*(AL_t(2,1)*ak1_t+AL_t(2,2)*ak2_t+
     &           AL_t(2,3)*ak3_t)
         akz(kpt)=2*pi*(AL_t(3,1)*ak1_t+AL_t(3,2)*ak2_t+
     &           AL_t(3,3)*ak3_t)
         enddo
         endif

         close(12)

       endif


       if(dabs(sumw-1.d0).gt.0.0000001d0) then
       if(inode.eq.1) then
       write(6,*) "**** Warning, sum of kpt weights not eq. 1", sumw
       write(6,*) "**** renormalize the kpt weight to 1"
       endif
         do kpt=1,nkpt
         weighkpt(kpt)=weighkpt(kpt)/sumw
         enddo
       endif

      return
      end subroutine readkpt

      subroutine error_stop(i1)
      implicit double precision(a-h,o-z)
      if(inode.eq.1) then
      write(6,*) "error in etot.input, line=",i1
      endif
      stop

      return
      end subroutine error_stop

      end
      
      

