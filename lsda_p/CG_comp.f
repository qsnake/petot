      subroutine CG_comp(ilocal,nline,tol,
     &  E_st,err_st,ave_line,vr,workr_n,mbad,ave_err,kpt,iislda)
****************************************
cc     Written by Lin-Wang Wang, March 30, 2001. 
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************
****************************************
****   The main conjugate gradient routine
****   mbad=0, normal iCG=1 run, ave_err not used
****   mbad>0, called from CG_new, fix the states: mx-mbad,mx to ave_err
******************************************
******************************************

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'param.escan_real'
***********************************************
       complex*16 pg(mg_nx),ugh(mg_nx),pgh(mg_nx)
       complex*16 pg_old(mg_nx),ughh_old(mg_nx)
       complex*16,allocatable,dimension(:)  ::  spg

       real*8 Dij0(32,32,mtype),Qij(32,32,mtype)
       integer isNLa(9,matom),ipsp_type(mtype)

       real*8 vr(mr_n)
       real*8 prec(mg_nx)
       complex*16 workr_n(mr_n)
       complex*16 cc
**********************************************
**** if iopt=0 is used, pghh_old can be deleted
**********************************************
        integer lin_st(mst)
	real*8 E_st(mst),err_st(mst),err2_st(mst)
	real*8 Ef
        complex*16 sumdum(32,matom),sumdum2(32,matom)
        real*8 occ_t(mtype)
       integer iiatom(mtype),icore(mtype),numref(matom),ityatom(mtype)

       common /com123b/m1,m2,m3,ngb,nghb,ntype,rrcut,msb
       common /comEk/Ek
       common /comisNLa/isNLa,Dij0,Qij,ipsp_all,ipsp_type
       common /comNL2/occ_t,iiatom,icore,numref,ityatom

       ng_n=ngtotnod(inode,kpt)

       if(ipsp_all.eq.1) then
       allocate(spg(1))
       else
       allocate(spg(mg_nx))
       endif


       err2=1.d0
       Eref=0.d0
       ughh_old = (0.0d0,0.0d0)
       pg_old = (0.0d0,0.0d0)
       iopt=0

      if(mbad.eq.0) mstart=1
      if(mbad.gt.0) mstart=mx-mbad+1

      do 2000 m=mstart,mx

      nltot=0
555   continue   ! turning back point for mbad.gt.0 run


       rr0=1.d+40


       do i=1,ng_n
       x=gkk_n(i,kpt)/Ek
       y=27.d0+x*(18.d0+x*(12.d0+x*8.d0))
       prec(i)=y/(y+16.d0*x**4)
       enddo

      call orth_comp(ug_n(1,m),ug_n,m-1,2,kpt,
     &  sug_n,sug_n(1,m),ipsp_all,0,1,fnorm)      ! no normalization, since sug_n(1,m) is not known yet

        
        

      do 3000 nint2=1,nline
      nltot=nltot+1

************************************************
**** ugh = (H-Eref) * ug
************************************************
      if(nint2.eq.1) then
        call Hpsi_comp(ug_n(1,m),ugh,ilocal,vr,workr_n,kpt,
     &   1,sug_n(1,m),sumdum,iislda)
        call orth_comp(ug_n(1,m),ug_n,0,1,kpt,
     &  sug_n,sug_n(1,m),ipsp_all,1,1,fnorm)     ! normalize ug_n and sug_n
        ugh=fnorm*ugh
        sumdum=fnorm*sumdum

      else
        do i=1,ng_n
        ugh(i)=ugh(i)*cos_th+pgh(i)*sin_th
        enddo
        sumdum=sumdum*cos_th+sumdum2*sin_th
      endif

************************************************
**** E_ug = < ug(i) | (H-Eref)^2 | ug(i) >
****      = < ugh(i) | ugh(i) >
************************************************
**** orthogonalization to ug_n(i,m) is done here outside, instead of 
**** inside orth_com, since we want E_ug, and err 

      s=0.d0


      do i=1,ng_n
      s=s+dreal(ugh(i)*dconjg(ug_n(i,m)))
      enddo

      call global_sumr(s)

      s=s*vol
      E_ug=s
      eigen=s

      err=0.d0

      if(ipsp_all.eq.1) then
      do i=1,ng_n
      err=err+cdabs(ugh(i)-s*ug_n(i,m))**2
      pg(i)=ugh(i)-s*ug_n(i,m)
      enddo
      else
      do i=1,ng_n
      err=err+cdabs(ugh(i)-s*sug_n(i,m))**2
      pg(i)=ugh(i)-s*sug_n(i,m)
      enddo
      endif

      call global_sumr(err)

      err=dsqrt(dabs(err*vol))

************************************************
****  orthogonalize pg to ug, 
****  pg=pg-s*sug
****  Note, the resulting pg*ug=0 (not pg*sug=0). 
************************************************
**** note, iss=-1, not 1 in orth_comp
************************************************
      call orth_comp(pg,ug_n,m-1,2,kpt,
     & sug_n,spg,ipsp_all,0,-1,fnorm)
************************************************
      err2=0.d0
      do i=1,ng_n
      err2=err2+cdabs(pg(i))**2
      enddo
      call global_sumr(err2)
      err2=dsqrt(dabs(err2*vol))

************************************************
****** check convergency
************************************************

c      if(err.lt.tol.or.nint2.eq.nline) then
      if(err2.lt.tol.or.nint2.eq.nline) then
      lin_st(m)=nint2 
      E_st(m)=eigen
      err_st(m)=err
      err2_st(m)=err2
      goto 3001
      endif
**********************************************

************************************************
**** begin conjugate gradient
**** should I use pg**2 formula, or pg*s*pg formula here !!??
**** This is still not clear. 
************************************************
      rr1=0.d0
      rr00=0.d0
      do i=1,ng_n
      rr00=rr00+cdabs(pg(i))**2*dble(prec(i))
      rr1=rr1+dble(prec(i))*(pg(i)-iopt*ughh_old(i))*
     &      dconjg(pg(i))
      ughh_old(i)=pg(i)
      enddo

      call global_sumr(rr00)
      call global_sumr(rr1)
**********************************************
**** here, pg == ughh
**********************************************

      beta=rr1/rr0
      rr0=rr00

**********************************************
***** pg(i) is the line minimization direction
**********************************************
***** preconditioning should be done before this
**********************************************

       do i=1,ng_n
       pg(i)=-pg(i)*dble(prec(i))+beta*pg_old(i)
       enddo
**********************************************
************************************************
      call orth_comp(pg,ug_n,m,2,kpt,
     & sug_n,spg,ipsp_all,0,1,fnorm)

       do i=1,ng_n
       pg_old(i)=pg(i)
       enddo
**********************************************
***** pgh = (H-Eref) * pg
**********************************************

      call Hpsi_comp(pg,pgh,ilocal,vr,workr_n,kpt,1,spg,
     & sumdum2,iislda)

      call orth_comp(pg,ug_n,0,1,kpt,       ! now we can normalize pg
     & sug_n,spg,ipsp_all,1,1,fnorm)   

      do i=1,ng_n
      pgh(i)=fnorm*pgh(i)         ! also normalize pgh
      enddo
      sumdum2=fnorm*sumdum2

**********************************************

      E_upg=0.d0
      E_pg=0.d0

        do i=1,ng_n
        E_upg=E_upg+dreal(pg(i)*dconjg(ugh(i)))
        E_pg=E_pg+dreal(pg(i)*dconjg(pgh(i)))
        enddo

      call global_sumr(E_upg)
      call global_sumr(E_pg)

      E_upg=E_upg*vol
      E_pg=E_pg*vol
**********************************************
**** calculate the theta from E_ug,E_pg,E_upg:
**** ug_new = ug*cos(th)+pg*sin(th)
**** E=E_ug*cos(th)^2+2*E_upg*cos(th)*sin(th)+E_pg*sin(th)^2
**********************************************

      theta=0.5d0*dabs(datan(2.d0*E_upg/(E_ug-E_pg)))
      cos_th=dcos(theta)
      sin_th=dsin(theta)

      pred_E=E_ug*cos_th**2+2*E_upg*cos_th*sin_th+
     &       E_pg*sin_th**2

**********************************************
**** update ug using theta
**********************************************
      do i=1,ng_n
      ug_n(i,m)=ug_n(i,m)*cos_th+pg(i)*sin_th
      enddo

      if(ipsp_all.eq.2) then     ! also update sug_n in the same way
      do i=1,ng_n
      sug_n(i,m)=sug_n(i,m)*cos_th+spg(i)*sin_th
      enddo
      endif

**********************************************
***** debuging:
**********************************************
c      if(inode.eq.1) then
c      write(6,*) nint2,E_ug,pred_E,rr1,rr0
c      write(6,700) nint2,E_ug,pred_E,err,err2,
c     &  E_upg,E_pg,theta
c700   format(i2,1x,2(f17.12,1x),2(E8.2,1x),3(E10.4,1x))
c      endif
**********************************************
**** do 3000, is for the nline line minimization
**********************************************
3000  continue
3001  continue

ccccccc using err2 to stop, instead of err. 
ccccccc The closeness to minimum for this conjugate gradient 
ccccccc depends on err2, not err. Err depends on the accuracy of previous states. 

c      if(mbad.gt.0.and.err2.gt.ave_err.and.
c     &   err_st(m).ge.tol.and.nltot.lt.40) goto 555
      if(mbad.gt.0.and.err2.gt.ave_err.and.
     &   err2.ge.tol.and.nltot.lt.40) goto 555

      if(mbad.gt.0.and.inode_tot.eq.1) then
      write(6,*) "nltot,err2,err_st", nltot,err2,err_st(m)
      endif

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
**********************************************
****  do 2000, is for all the wavefunctions m
**********************************************
2000  continue

      deallocate(spg)

      if(mbad.gt.0) return

      Esum=0.d0
      ave_line=0.d0
      ave_err2=0.d0
      do i=1,mx
      Esum=Esum+E_st(i)
      ave_line=ave_line+lin_st(i)
      ave_err2=ave_err2+err2_st(i)
      enddo

      ave_line=ave_line/mx
      ave_err2=ave_err2/mx

      if(nline.ne.0.and.inode_tot.eq.1) then
      write(6,*) "*********************************"
      write(6,*) "**** kpt=", kpt
      write(6,*) "the number of line minimization used"
      write(6,101) (lin_st(i),i=1,mx)
      write(6,104) ave_err2
104   format("the err of each states, A.U.    Aver err2=", E10.4)
      write(6,102) (err_st(i), i=1,mx)
      write(6,*) "the eigen energies, in eV"
      write(6,103) (E_st(i)*27.211396d0, i=1,mx)
      write(6,*) "*********************************"
      write(6,*) "sum eigen(dbl occ, Ryd)=", Esum*2*2
      write(6,*) "*********************************"
      endif
101   format(5(i6,7x))
102   format(5(E10.4,3x))
103   format(5(f13.8,1x))

***********************************************

      return
      end

