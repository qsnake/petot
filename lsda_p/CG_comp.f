      subroutine CG_comp(ilocal,nline,tol,
     &  E_st,err_st,ave_line,vr,workr_n,mbad,ave_err,kpt)
****************************************
cc     Written by Lin-Wang Wang, March 30, 2001. 
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work 
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

       real*8 vr(mr_n)
       real*8 prec(mg_nx)
       complex*16 workr_n(mg_nx)
**********************************************
**** if iopt=0 is used, pghh_old can be deleted
**********************************************
        integer lin_st(mst)
	real*8 E_st(mst),err_st(mst)
	real*8 Ef,occ(mst)

       common /com123b/m1,m2,m3,ngb,nghb,ntype,rrcut,msb
       common /comEk/Ek

       ng_n=ngtotnod(inode,kpt)


       err2=1.d0
       Eref=0.d0
       ughh_old = (0.0d0,0.0d0)
       pg_old = (0.0d0,0.0d0)
       iopt=0

      if(mbad.eq.0) mstart=1
      if(mbad.gt.0) mstart=mx-mbad+1

      do m=mstart,mx
      call orth_comp(ug_n(1,m),ug_n,m-1,1,kpt)
      enddo

      do 2000 m=mstart,mx

      nltot=0
555   continue   ! turning back point for mbad.gt.0 run


       rr0=1.d+40


       do i=1,ng_n
       x=gkk_n(i,kpt)/Ek
       y=27.d0+x*(18.d0+x*(12.d0+x*8.d0))
       prec(i)=y/(y+16.d0*x**4)
       enddo

      call orth_comp(ug_n(1,m),ug_n,m-1,1,kpt)

      do 3000 nint2=1,nline
      nltot=nltot+1

************************************************
**** ugh = (H-Eref) * ug
************************************************
      if(nint2.eq.1) then
        call Hpsi_comp(ug_n(1,m),ugh,ilocal,vr,workr_n,kpt)
      else
        do i=1,ng_n
        ugh(i)=ugh(i)*cos_th+pgh(i)*sin_th
        enddo
      endif

************************************************
**** E_ug = < ug(i) | (H-Eref)^2 | ug(i) >
****      = < ugh(i) | ugh(i) >
************************************************
      s=0.d0

      do i=1,ng_n
      s=s+dreal(ugh(i)*dconjg(ug_n(i,m)))
      enddo

      call global_sumr(s)

      s=s*vol
      E_ug=s
      eigen=s

      err=0.d0

      do i=1,ng_n
      err=err+cdabs(ugh(i)-s*ug_n(i,m))**2
      enddo

      call global_sumr(err)

      err=dsqrt(dabs(err*vol))

************************************************
****** check convergency
************************************************

      if(err.lt.tol.or.nint2.eq.nline) then
      lin_st(m)=nint2 
      E_st(m)=eigen
      err_st(m)=err
      goto 3001
      endif
**********************************************
**** here, pg == ughh
************************************************
        do i=1,ng_n
        pg(i)=ugh(i)
        enddo
**********************************************
****  orthogonalize ughh to ug
************************************************
      call orth_comp(pg,ug_n,m,2,kpt)
************************************************
      err2=0.d0
      do i=1,ng_n
      err2=err2+cdabs(pg(i))**2
      enddo
      call global_sumr(err2)
      err2=dsqrt(dabs(err2*vol))

************************************************
**** begin conjugate gradient
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
       do i=1,ng_n
       pg(i)=-pg(i)*dble(prec(i))+beta*pg_old(i)
       enddo
**********************************************
****  another option: orth_real(pg,ug,mx,2) 
************************************************
      call orth_comp(pg,ug_n,m,2,kpt)

      s=0.d0
      do i=1,ng_n
      s=s+cdabs(pg(i))**2
      enddo
      
      call global_sumr(s)

      s=1/dsqrt(s*vol)

      do i=1,ng_n
      pg_old(i)=pg(i)
      pg(i)=s*pg(i)
      enddo
**********************************************
***** pgh = (H-Eref) * pg
**********************************************

      call Hpsi_comp(pg,pgh,ilocal,vr,workr_n,kpt)

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

**********************************************
***** debuging:
**********************************************
*      write(6,*) nint2,E_ug,pred_E
**********************************************
**** do 3000, is for the nline line minimization
**********************************************
3000  continue
3001  continue
      if(mbad.gt.0.and.err2.gt.ave_err.and.
     &   err_st(m).ge.tol.and.nltot.lt.40) goto 555

      if(mbad.gt.0.and.inode.eq.1) then
      write(6,*) "nltot,err2,err_st", nltot,err2,err_st(m)
      endif



**********************************************
****  do 2000, is for all the wavefunctions m
**********************************************
2000  continue

      if(mbad.gt.0) return

      Esum=0.d0
      ave_line=0.d0
      do i=1,mx
      Esum=Esum+E_st(i)
      ave_line=ave_line+lin_st(i)
      enddo

      ave_line=ave_line/mx

      if(nline.ne.0.and.inode.eq.1) then
      write(6,*) "*********************************"
      write(6,*) "the number of line minimization used"
      write(6,101) (lin_st(i),i=1,mx)
      write(6,*) "the err of each states, A.U"
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

