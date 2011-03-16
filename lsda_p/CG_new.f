      subroutine CG_new(ilocal,nline,tol,
     &  E_st,err_st,ave_line,vr,workr_n,mCGbad,kpt,iislda)
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

      include 'param.escan_real'

***********************************************
      parameter (mline=20)
      parameter (lwork=6000)

      complex*16 workx(lwork)
      real*8 workrx(3*mline)

      complex*16 sumdum(32,matom)
      complex*16, allocatable,dimension(:,:,:) :: sumdum_all


      complex*16 ugh(mg_nx),pg(mg_nx)

      complex*16 Rg(mg_nx,mline+1),Rgh(mg_nx,mline+1)
      complex*16 SRg(mg_nx,mline+1)
      
**** first, use the stupid but easy to program way

       real*8 vr(mr_n),workr_n(mr_n)

       real*8 prec(mg_nx)

       real*8 Dij0(32,32,mtype),Qij(32,32,mtype)
       integer isNLa(9,matom),ipsp_type(mtype)


       complex*16 hh(mline,mline),Z(mline,mline)
       complex*16 hs(mline,mline),ss(mline,mline)
       real*8 EE(mline),EE1(mline),w(mline)
       real*8 aux(3*mline)

       real*8 occ_t(mtype)
       integer iiatom(mtype),icore(mtype),numref(matom),ityatom(mtype)
**********************************************
        integer lin_st(mst)
	real*8 E_st(mst),err_st(mst)
	complex*16 c,cc,c1,c2,c3


       common /comisNLa/isNLa,Dij0,Qij,ipsp_all,ipsp_type
       common /com123b/m1,m2,m3,ngb,nghb,ntype,rrcut,msb
       common /comEk/Ek
       common /comNL2/occ_t,iiatom,icore,numref,ityatom


       ng_n=ngtotnod(inode,kpt)

       if(nline.gt.mline) then
       write(6,*) "nline > mline, stop", nline
       stop
       endif

       do i=1,ng_n
       x=gkk_n(i,kpt)/Ek
       y=27.d0+x*(18.d0+x*(12.d0+x*8.d0))
       prec(i)=y/(y+16.d0*x**4)
       enddo

       allocate(sumdum_all(32,matom,nline))

       do 1000 m=1,mx
***************************************************
       do i=1,ng_n
       Rg(i,1)=ug_n(i,m)
       enddo
***************************************************
       do 200 nt=1,nline

ccccccc  the input Rg(1,nt) might not be normalized, since SRg is not known !

       call Hpsi_comp(Rg(1,nt),Rgh(1,nt),
     &    ilocal,vr,workr_n,kpt,1,SRg(1,nt),sumdum,iislda) 

ccccccc  now, normalize Rg(1,nt)

       call orth_comp(Rg(1,nt),Rg,0,1,kpt,
     &  SRg,SRg(1,nt),ipsp_all,1,1,fnorm)     ! normalize Rg and SRg
        do i=1,ng_n
        Rgh(i,nt)=fnorm*Rgh(i,nt)
        enddo
        sumdum=fnorm*sumdum
        sumdum_all(:,:,nt)=sumdum(:,:)
cccccccccccccccccccccccccccccccccccccc



       if(nt.eq.1) then
       s=0.d0
       do i=1,ng_n
       s=s+dreal(Rgh(i,nt)*dconjg(Rg(i,nt)))
       enddo
       
       call global_sumr(s)
       s=s*vol
       Eref=s
       endif



       do 50 nt1=1,nt
       c1=dcmplx(0.d0,0.d0)
       do i=1,ng_n
       c1=c1+Rgh(i,nt)*dconjg(Rgh(i,nt1))
       enddo

       c2=dcmplx(0.d0,0.d0)

       if(ipsp_all.eq.1) then
       do i=1,ng_n
       c2=c2+Rgh(i,nt)*dconjg(Rg(i,nt1))
       enddo
       else
       do i=1,ng_n
       c2=c2+Rgh(i,nt)*dconjg(SRg(i,nt1))+
     &       SRg(i,nt)*dconjg(Rgh(i,nt1)) 
       enddo
       c2=c2/2
       endif

       call global_sumc(c1)
       call global_sumc(c2)
       c1=c1*vol
       c2=c2*vol
       hh(nt,nt1)=dconjg(c1)
       hh(nt1,nt)=c1
       hs(nt,nt1)=dconjg(c2)
       hs(nt1,nt)=c2

       if(ipsp_all.eq.2) then
       c3=dcmplx(0.d0,0.d0)
       do i=1,ng_n
       c3=c3+SRg(i,nt)*dconjg(SRg(i,nt1))
       enddo
       call global_sumc(c3)
       c3=c3*vol
       ss(nt,nt1)=dconjg(c3)
       ss(nt1,nt)=c3
       endif

 50    continue
       
cccccccccccccccccccccccccccccccccccccccccccc


       do i1=1,nt
       do i2=1,nt
       z(i1,i2)=hh(i1,i2)-2*Eref*hs(i1,i2)
       if(ipsp_all.eq.2) then
       z(i1,i2)=z(i1,i2)+Eref**2*ss(i1,i2)
       endif
       enddo
       if(ipsp_all.eq.1) then
       z(i1,i1)=z(i1,i1)+Eref**2
       endif
       enddo


       call system_czheev('v','U',nt,z,mline,EE,workx,
     &  lwork,workrx,info)


       do i=1,ng_n
       c=dcmplx(0.d0,0.d0)
       c1=dcmplx(0.d0,0.d0)
       do nt1=1,nt
       c=c+Z(nt1,1)*Rg(i,nt1)
       c1=c1+Z(nt1,1)*Rgh(i,nt1)
       enddo
       Rg(i,nt+1)=c
       Rgh(i,nt+1)=c1
       enddo


       if(ipsp_all.eq.2) then
       do i=1,ng_n
       c=dcmplx(0.d0,0.d0)
       do nt1=1,nt
       c=c+Z(nt1,1)*SRg(i,nt1)
       enddo
       SRg(i,nt+1)=c
       enddo
       endif

************************************************
*** update E_ref according to Rayleigh quotient
************************************************
       s=0.d0
       do i=1,ng_n
       s=s+dreal(Rgh(i,nt+1)*dconjg(Rg(i,nt+1)))
       enddo

       call global_sumr(s)
       s=s*vol
       Eref=s
**************************************************

       err=dsqrt(dabs(EE(1)))

cccccccc test
c       if(inode.eq.1) then
c       write(6,*) "m,nt,err,Eref",m,nt,err,Eref
c       endif
ccccccccccccccccccccccc

       if(nt.eq.nline.or.err.lt.tol) then
       lin_st(m)=nt
       do i=1,ng_n
       ug_n(i,m)=Rg(i,nt+1)
       enddo

       if(ipsp_all.eq.2) then
       do i=1,ng_n 
       sug_n(i,m)=SRg(i,nt+1)
       enddo
       endif

       sumdum=dcmplx(0.d0,0.d0)
       do nt1=1,nt
       sumdum(:,:)=sumdum(:,:)+Z(nt1,1)*sumdum_all(:,:,nt1)
       enddo

       goto 201

       else

       call orth_comp(Rgh(1,nt+1),Rg,nt+1,3,kpt,
     &   SRg,SRg(1,nt+1),ipsp_all,0,-1,fnorm)    ! SRg(1,nt+1) not used

       do i=1,ng_n
       Rg(i,nt+1)=Rgh(i,nt+1)*prec(i)
       enddo

       call orth_comp(Rg(1,nt+1),Rg,nt,2,kpt,
     &   SRg,SRg(1,nt+1),ipsp_all,0,1,fnorm)     ! SRg(1,nt+1) not used, not normalized

       endif

200    continue

201    continue

       err_st(m)=err
       E_st(m)=Eref

cccccc  store sumdum in beta_psi
       mx_n=mx/nnodes+1
       if((m-1)/mx_n+1.eq.inode) then
       im=m-(inode-1)*mx_n
       iref_t=0
       do ia=1,natom
       do iref=1,numref(ia)
       iref_t=iref_t+1
       beta_psi(iref_t,im)=sumdum(iref,ia)
       enddo
       enddo
       endif

1000   continue

       deallocate(sumdum_all)


       ave_line=0
       do i=1,mx
       ave_line=ave_line+lin_st(i)
       enddo
       ave_line=ave_line/mx
*****************************************************
*****************************************************
***** orthogonalization and looking for the bad states
       mbad=0
       m=0
2000   continue      ! do loop
       m=m+1

       call orth_comp(ug_n(1,m),ug_n,0,1,kpt,
     &  sug_n,sug_n(1,m),ipsp_all,1,1,fnorm)     ! normalize ug_n(1,m) and sug_n(1,m)

       call orth_comp(ug_n(1,m),ug_n,m-1,1,kpt,
     &  sug_n,sug_n(1,m),ipsp_all,1,1,fnorm)     ! orth ug_n(1,m) to ug_n(1,m'), m'=1,m-1

ccccccc looking for bad states (converged into a same state
ccccccc as one below it)

       if(fnorm.gt.5.d0) then
       call Hpsi_comp(ug_n(1,m),Rgh(1,1),
     &    ilocal,vr,workr_n,kpt,1,sug_n(1,m),sumdum,iislda) 

       s=0.d0
       do i=1,ng_n
       s=s+dreal(Rgh(i,1)*dconjg(ug_n(i,m)))
       enddo
       call global_sumr(s)
       s=s*vol
       s_err=0.d0
       if(ipsp_all.eq.1) then
       do i=1,ng_n
       s_err=s_err+cdabs(Rgh(i,1)-s*ug_n(i,m))**2
       enddo
       else
       do i=1,ng_n
       s_err=s_err+cdabs(Rgh(i,1)-s*sug_n(i,m))**2
       enddo
       endif

       call global_sumr(s_err)
       s_err=dsqrt(s_err*vol)

       if(s_err.gt.10*err_st(m)) then   ! yes, it is a bad state
       mbad=mbad+1
       do i=1,ng_n
       Rgh(i,1)=ug_n(i,m)
       enddo
       if(ipsp.eq.2) then
       do i=1,ng_n
       SRg(i,1)=sug_n(i,m)
       enddo
       endif

       do mm=m,mx-1
       err_st(mm)=err_st(mm+1)
       E_st(mm)=E_st(mm+1)
       do i=1,ng_n
       ug_n(i,mm)=ug_n(i,mm+1)
       enddo
       if(ipsp.eq.2) then
       do i=1,ng_n
       sug_n(i,mm)=sug_n(i,mm+1)
       enddo
       endif
       enddo

       err_st(mx)=s_err
       E_st(mx)=s
       do i=1,ng_n
       ug_n(i,mx)=Rgh(i,1)
       enddo
       if(ipsp_all.eq.2) then
       do i=1,ng_n
       sug_n(i,mx)=SRg(i,1)
       enddo
       endif
       

       m=m-1

       endif      ! s_err.gt.10.err_st(m)
       endif      ! fnorm.gt.5

       if(m.lt.mx-mbad) goto 2000
*****************************************************
**** note, the last mbad states are not orthogonalize to previous
**** states. But, they will be orthogonalized in CG_comp
*****************************************************
       if(mbad.gt.0) then
       if(inode.eq.1) then
       write(6,*) "mbad.gt.0, recalc. the last mbad states,kpt", 
     &   mbad,kpt
       endif

       ave_err=0.d0
       do m=1,mx-mbad
       ave_err=ave_err+err_st(m)
       enddo
       ave_err=ave_err/(mx-mbad)

       call CG_comp(ilocal,10,tol,E_st,err_st,
     &    ave_line_tmp,vr,workr_n,mbad,ave_err,kpt,iislda)
       endif

       if(mCGbad.gt.0) then
       if(inode.eq.1) then
       write(6,*) "mCGbad.gt.0, recalc. the last mCGbad states,kpt", 
     &    mCGbad,kpt
       endif
       ave_err=0.d0
       do m=mx-mCGbad+1,mx
       ave_err=ave_err+err_st(m)
       enddo
       ave_err=ave_err/mCGbad/3

       if(inode.eq.1) then
       write(6,*) "kpt,ave_err=",kpt,ave_err
       endif

       call CG_comp(ilocal,10,tol,E_st,err_st,
     &    ave_line_tmp,vr,workr_n,mCGbad,ave_err,kpt,iislda)

        endif

*****************************************************
       Esum=0.d0
       do i=1,mx
       Esum=Esum+E_st(i)
       enddo
*****************************************************
       if(inode_tot.eq.1) then
       write(6,*) "*********************************"
       write(6,*) "****** kpt=",kpt
       write(6,*) "err of each states, A.U"
       write(6,102) (err_st(i), i=1,mx)
       write(6,*) "eigen energies, in eV"
       write(6,103) (E_st(i)*27.211396d0, i=1,mx)
       write(6,*) "*********************************"
       write(6,*) "sum eigen(dbl occ, Ryd)=", Esum*2*2
       write(6,*) "*********************************"
       endif

101   format(5(i6,7x))
102   format(5(E10.4,3x))
103   format(5(f13.8,1x))
***********************************************
*** for metal only
***********************************************


      return
      end

