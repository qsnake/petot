      subroutine CG_new(ilocal,nline,tol,
     &  E_st,err_st,ave_line,vr,workr_n,mCGbad,kpt)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
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

      complex*16 ugh(mg_nx),pg(mg_nx)

      complex*16 Rg(mg_nx,mline+1),Rgh(mg_nx,mline+1)
**** first, use the stupid but easy to program way

       real*8 vr(mr_n),workr_n(mr_n)

       real*8 prec(mg_nx)

       complex*16 hh(mline,mline),Z(mline,mline)
       complex*16 h(mline,mline)
       real*8 EE(mline),EE1(mline),w(mline)
       real*8 aux(3*mline)
**********************************************
        integer lin_st(mst)
	real*8 E_st(mst),err_st(mst)
	complex*16 c,cc,c1,c2


       common /com123b/m1,m2,m3,ngb,nghb,ntype,rrcut,msb
       common /comEk/Ek


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

       do 1000 m=1,mx
***************************************************
       do i=1,ng_n
       Rg(i,1)=ug_n(i,m)
       enddo
***************************************************
       do 200 nt=1,nline

       call Hpsi_comp(Rg(1,nt),Rgh(1,nt),
     &        ilocal,vr,workr_n,kpt)

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
       c2=dcmplx(0.d0,0.d0)
       do i=1,ng_n
       c1=c1+Rgh(i,nt)*dconjg(Rgh(i,nt1))
       c2=c2+Rgh(i,nt)*dconjg(Rg(i,nt1))
       enddo

       call global_sumc(c1)
       call global_sumc(c2)
       c1=c1*vol
       c2=c2*vol
       hh(nt,nt1)=dconjg(c1)
       hh(nt1,nt)=c1
       h(nt,nt1)=dconjg(c2)
       h(nt1,nt)=c2
 50    continue


       do i1=1,nt
       do i2=1,nt
       z(i1,i2)=hh(i1,i2)-2*Eref*h(i1,i2)
       enddo
       z(i1,i1)=z(i1,i1)+Eref**2
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

       if(nt.eq.nline.or.err.lt.tol) then
       lin_st(m)=nt
       do i=1,ng_n
       ug_n(i,m)=Rg(i,nt+1)
       enddo
       goto 201

       else

       call orth_comp(Rgh(1,nt+1),Rg,nt+1,3,kpt)

       do i=1,ng_n
       Rg(i,nt+1)=Rgh(i,nt+1)*prec(i)
       enddo

       call orth_comp(Rg(1,nt+1),Rg,nt,1,kpt)

       endif

200    continue

201    continue

       err_st(m)=err
       E_st(m)=Eref

1000   continue


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

       s=0.d0
       do i=1,ng_n
       s=s+cdabs(ug_n(i,m))**2
       enddo
       call global_sumr(s)
       s=1./dsqrt(s*vol)

       do i=1,ng_n
       ug_n(i,m)=s*ug_n(i,m)
       enddo

       call orth_comp(ug_n(1,m),ug_n,m-1,2,kpt)

       s=0.d0
       do i=1,ng_n
       s=s+cdabs(ug_n(i,m))**2
       enddo
       call global_sumr(s)
       s=1./dsqrt(s*vol)
       do i=1,ng_n
       ug_n(i,m)=s*ug_n(i,m)
       enddo

ccccccc looking for bad states (converged into a same state
ccccccc as one below it)

       if(s.gt.5.d0) then
       call Hpsi_comp(ug_n(1,m),Rgh(1,1),
     &        ilocal,vr,workr_n,kpt)
       s=0.d0
       do i=1,ng_n
       s=s+dreal(Rgh(i,1)*dconjg(ug_n(i,m)))
       enddo
       call global_sumr(s)
       s=s*vol
       s_err=0.d0
       do i=1,ng_n
       s_err=s_err+cdabs(Rgh(i,1)-s*ug_n(i,m))**2
       enddo
       call global_sumr(s_err)
       s_err=dsqrt(s_err*vol)

       if(s_err.gt.10*err_st(m)) then   ! yes, it is a bad state
       mbad=mbad+1
       do i=1,ng_n
       Rgh(i,1)=ug_n(i,m)
       enddo

       do mm=m,mx-1
       err_st(mm)=err_st(mm+1)
       E_st(mm)=E_st(mm+1)
       do i=1,ng_n
       ug_n(i,mm)=ug_n(i,mm+1)
       enddo
       enddo

       err_st(mx)=s_err
       E_st(mx)=s
       do i=1,ng_n
       ug_n(i,mx)=Rgh(i,1)
       enddo

       m=m-1

       endif
       endif

       if(m.lt.mx-mbad) goto 2000
*****************************************************
**** note, the last mbad states are not orthogonalize to previous
**** states. But, they will be orthogonalized in CG_comp
*****************************************************
       if(mbad.gt.0) then
       if(inode.eq.1) then
       write(6,*) "mbad.gt.0, recalc. the last mbad states", mbad
       endif
       ave_err=0.d0
       do m=1,mx-mbad
       ave_err=ave_err+err_st(m)
       enddo
       ave_err=ave_err/(mx-mbad)

       call CG_comp(ilocal,10,tol,E_st,err_st,
     &    ave_line_tmp,vr,workr_n,mbad,ave_err,kpt)
       endif

       if(mCGbad.gt.0) then
       if(inode.eq.1) then
       write(6,*) "mCGbad.gt.0, recalc. the last mCGbad states", mCGbad
       endif
       ave_err=0.d0
       do m=mx-mCGbad+1,mx
       ave_err=ave_err+err_st(m)
       enddo
       ave_err=ave_err/mCGbad/3

       if(inode.eq.1) then
       write(6,*) "ave_err=",ave_err
       endif

       call CG_comp(ilocal,10,tol,E_st,err_st,
     &     ave_line_tmp,vr,workr_n,mCGbad,ave_err,kpt)
        endif

*****************************************************
       Esum=0.d0
       do i=1,mx
       Esum=Esum+E_st(i)
       enddo
*****************************************************
       if(inode.eq.1) then
       write(6,*) "*********************************"
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

