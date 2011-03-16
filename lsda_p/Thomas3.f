******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Currently, we don't have a version which deal with islda=2 for both spin together. 
ccc The current code deal sperately v_in,v_out up and v_in,v_out down both with the total charge density. 
ccc As a result, for islda=2, there could be a oscillation between the up and down spins,
cccc thus Thomas3 might not be ideal for islda=2, but it can still be helpful.
cccccccccccccccccccccccccccccccccccccc

*******************************************
*** This program using (v_in, rho) -> (v_out) to fixed the
*** Hamiltonian: V_ion: (so V_hxc(rho)+V_ion=v_out)
***                       Note, this V_ion is not the true V_ion, because v_out 
****                      is the pulay mixed result. 
***              dvp: (so v_in & rho, Ef0 satisfy the solution eq.)
***  Then, it resolve rho', so that v_in'=v_out'
***  Finally, the v_in' is used for small k 
***  and kerker mixing is used for large k
****************************************************************
*** The final equation for psi=dsqrt(rho)
*** -0.5 \nabla^2 psi +( rho^(2/3)+ v_in + vionT_nL)*psi + dvp = Ef*psi
***  and v_in  = V_ion + v_Hxc(rho)
*** vionT_nL is the Thomas-Fermi corr. for the nonlocal part of the potential
*** Note, below, in this program vion = V_ion + vionT_nL
****************************************************************
*** For ipsp_all.eq.2 (ultrasoft), vionT_nL=0 (since there is no vs,vp,vd)
*** we will generate vionT_nL assuming dvp=0, then for rho < 10^-3 rho_ave region, 
*** set vionT_nL=0, then generate dvp
****************************************************************

      subroutine Thomas3(v_in,v_out,Ef0,totNel,islda,ipsp_all)

*****************************************************************
*** in return, v_in will be replaced by new v_in
*** rho, v_out, dv are not changed
*****************************************************************

      use fft_data
      use load_data
      use data

      implicit double precision(a-h,o-z)
      include 'param.escan_real'


      real*8 v_in(mr_nL),v_out(mr_nL)
********* inputs
********* vion, dvp are derived from input v_in,v_out,rho
********* vion = V_ion + vionT_nL 
********* V[rho]=V_HXC[rho] + V_ion = V_HXC[rho] + vion - vionT_nL

      real*8,allocatable,dimension(:) :: workr_n,rho_0,ugr,
     &   ughr,pgr,pgr_old,vion,dvp

      allocate(workr_n(mr_nL))
      allocate(rho_0(mr_nL))
      allocate(ugr(mr_nL))
      allocate(ughr(mr_nL))
      allocate(pgr(mr_nL))
      allocate(pgr_old(mr_nL))
      allocate(vion(mr_nL))
      allocate(dvp(mr_nL))

      if(inode.eq.1) then
      write(6,*) "from Thomas3, ipsp_all=",ipsp_all
      endif


      ng2_nL=ngtotnod2L(inode)

      if(inode.eq.1) then
        write(6,*) "Ef0=", Ef0
      endif

      pi=4*datan(1.d0)
      alpha=1.d0
      pgr_old = 0.0d0

      rho_0=0.d0

      do iislda=1,islda
      do i=1,nr_nL
      rho_0(i) = rho_0(i)+ rho_nL(i,iislda)
      enddo
      enddo
********************************************
*** make rho(i) = workr(i)**2 and workr(i) is in the Ecut2
********************************************
      s1=0.d0
 
      do i=1,nr_nL
       s1=s1+rho_0(i)
       workr_n(i)=dsqrt(dabs(rho_0(i)))
      enddo

      call global_sumr(s1)

      call d3fft_real2L(pgr,workr_n,1,1)
      call d3fft_real2L(pgr,workr_n,-1,1)

      s2=0.d0

      do i=1,nr_nL
        rho_0(i)=workr_n(i)**2
        s2=s2+rho_0(i)
      enddo

      call global_sumr(s2)

      s=s1/s2
      do i=1,nr_nL
        rho_0(i)=rho_0(i)*s
      enddo

*********************************************
*** here, use pgr,pgi as work functions
*********************************************

      fa=0.5d0*(3*pi**2)**(2.d0/3.d0)
ccccccccccccccccccccccccccccccccccccc


      if(ipsp_all.eq.1) then
      do i=1,nr_nL
       pgr(i)=fa*rho_0(i)**(2.d0/3.d0)+(v_in(i)+vionT_nL(i))
      enddo
      else
       vionT_nL=0.d0
      do i=1,nr_nL
       pgr(i)=fa*rho_0(i)**(2.d0/3.d0)+v_in(i)
      enddo
      endif
     

ccccccccccccccccccccccccccccccccccccccccccccccc

      do i=1,nr_nL
      workr_n(i)=dsqrt(rho_0(i))
      enddo


      call d3fft_real2L(ugr,workr_n,1,1)

      do i=1,ng2_nL
        ugr(i*2)=alpha*gkk2_nL(i)*ugr(i*2)
        ugr(i*2-1)=alpha*gkk2_nL(i)*ugr(i*2-1)
      enddo


      call d3fft_real2L(ugr,workr_n,-1,1)

      ugr = workr_n

      do i=1,nr_nL
       dvp(i)=(Ef0-pgr(i))*dsqrt(rho_0(i))-ugr(i)
       dvp(i)=dvp(i)/dsqrt(totNel)
      enddo

      rho_ave=totNel/vol

      if(ipsp_all.eq.2) then
      do i=1,nr_nL
       vionT_nL(i)=Ef0-pgr(i)-ugr(i)/dsqrt(rho_0(i))
       if(rho_0(i).lt.rho_ave*1.D-3) vionT_nL(i)=0.d0
       dvp(i)=dvp(i)-vionT_nL(i)*dsqrt(rho_0(i))/dsqrt(totNel)
      enddo
      endif
cccccccccccccccc finished generate dvp(i)


      do i=1,nr_nL
        workr_n(i) = rho_0(i)
      enddo

      call d3fft_real2L(pgr,workr_n,1,1)

      do i=1,ng2_nL

       if(inode.eq.iorg2L(1).and.i.eq.iorg2L(2)) then

        pgr(i*2)=0.d0
        pgr(i*2-1)=0.d0
        else
        pgr(i*2)=pgr(i*2)*2*pi/gkk2_nL(i)
        pgr(i*2-1)=pgr(i*2-1)*2*pi/gkk2_nL(i)
       endif

      enddo


      call d3fft_real2L(pgr,workr_n,-1,1)

      pgr = workr_n


      do i=1,nr_nL
       pgr(i)=pgr(i)+UxcCA(rho_0(i)+rhocr_nL(i),uxc2)
       vion(i)=(v_out(i)-pgr(i))+vionT_nL(i)
      enddo

cccccccccccc  finished generate vion


ccccccccccccccccccccccccccccccccccccccccccccccc
ccc finished obtaining vion(r) [=V_ion + vionT_nL] and dvp(r)
ccc Now, the Hamiltonian has be determined, so just need
ccc to solve the Hamiltonian selfconsistently
ccccccccccccccccccccccccccccccccccccccccccccccc
      s1=0.d0
      s2=0.d0
      s3=0.d0
      s4=0.d0

      do i=1,nr_nL
      s1=s1+dvp(i)**2
      s2=s2+rho_0(i)
      s3=s3+dabs(v_in(i)+vionT_nL(i))
      s4=s4+vion(i)**2
      enddo

      call global_sumr(s1)
      call global_sumr(s2)
      call global_sumr(s3)
      call global_sumr(s4)

      s1=dsqrt(s1/s2)
      s3=s3/nrL
      s4=dsqrt(s4/nrL)

      if(inode.eq.1) then
        write(6,*) "average dvp (Hartree)", s1
        write(6,*) "average v_in+vionT", s3
        write(6,*) "average vion", s4
      endif
******************************************************
*****, Now, the potential vion(i) and dvp(i) has been calculated
******************************************************

       do i=1,nr_nL
       workr_n(i)=dsqrt(rho_0(i)/totNel)
       enddo

       call d3fft_real2L(ugr,workr_n,1,1)

*****************************************************
*****************************************************
ccccc ugr,ugi in G space

       nmax2=30

       do 1000 nint=1,nmax2

       if(nint.eq.1) rr0=1.D+40

       s=0.d0

       do i=1,ng2_nL
         s=s+ugr(i*2-1)**2+ugr(i*2)**2
       enddo

       call global_sumr(s)

       s=s*vol*2
       s=dsqrt(1.d0/s)

       do i=1,ng2_nL
        ugr(i*2)=s*ugr(i*2)
        ugr(i*2-1)=s*ugr(i*2-1)
       enddo

******************
ccccc calculate H*u and Etot
       call Hpsi3(ugr,ughr,vion,dvp,Etot,alpha,totNel,workr_n)
******************

      do i=1,ng2_nL
       pgr(i*2)=ughr(i*2)
       pgr(i*2-1)=ughr(i*2-1)
      enddo

      Ef=0.d0

      do i=1,ng2_nL
       Ef=Ef+pgr(i*2-1)*ugr(i*2-1)+pgr(i*2)*ugr(i*2)
      enddo

      call global_sumr(Ef)

      Ef=Ef*vol*2

      err=0.d0

      do i=1,ng2_nL
        pgr(i*2)=pgr(i*2)-Ef*ugr(i*2)
        pgr(i*2-1)=pgr(i*2-1)-Ef*ugr(i*2-1)
        err=err+pgr(i*2)**2+pgr(i*2-1)**2
      enddo

      call global_sumr(err)

      err=dsqrt(err*vol*2)

      if((nint.eq.1.or.nint.eq.nmax2).and.inode.eq.1) then
      write(6,*) "nint,Ef,err",nint,Ef,err
      endif

      Ek=0.5d0
      rr1=0.d0

      do i=1,ng2_nL
        x=gkk2_nL(i)/Ek
        y=1/(1+x**4)**0.25d0
        rr1=rr1+y*(pgr(i*2-1)**2+pgr(i*2)**2)
      enddo

      call global_sumr(rr1)

      rr1=rr1*vol*2

      beta=rr1/rr0

      rr0=rr1

      do i=1,ng2_nL
       x=gkk2_nL(i)/Ek
       y=1/(1+x**4)**0.25d0
       pgr(i*2)=-pgr(i*2)*y+beta*pgr_old(i*2)
       pgr(i*2-1)=-pgr(i*2-1)*y+beta*pgr_old(i*2-1)
      enddo

      s=0.d0

      do i=1,ng2_nL
      s=s+pgr(i*2)*ugr(i*2)+pgr(i*2-1)*ugr(i*2-1)
      enddo

      call global_sumr(s)

      s=s*vol*2

      do i=1,ng2_nL
        pgr(i*2)=pgr(i*2)-s*ugr(i*2)
        pgr(i*2-1)=pgr(i*2-1)-s*ugr(i*2-1)
      enddo

      s=0.d0

      do i=1,ng2_nL
       s=s+pgr(i*2)**2+pgr(i*2-1)**2
      enddo

      call global_sumr(s)

      s=1.d0/dsqrt(s*vol*2)

      do i=1,ng2_nL

        pgr_old(i*2)=pgr(i*2)
        pgr_old(i*2-1)=pgr(i*2-1)

        pgr(i*2)=pgr(i*2)*s
        pgr(i*2-1)=pgr(i*2-1)*s

      enddo

**************************************************

      E_upg=0.d0

      do i=1,ng2_nL
      E_upg=E_upg+pgr(i*2)*ughr(i*2)+pgr(i*2-1)*ughr(i*2-1)
      enddo

      call global_sumr(E_upg)

      E_upg=E_upg*vol*2
      E_upg=2*E_upg*totNel

**************************************************
ccccc  theta0=0.2D-3, might cause slightly different results for
ccccc  parallel and serie code due to round off error.
cccccc for debug, use theta0=0.2D-2,but for actual run, use 0.2D-03.
**************************************************
      theta0=0.2D-03
c      theta0=0.2D-02

      cos_th=dcos(theta0)
      sin_th=dsin(theta0)

      do i=1,ng2_nL
        ugr(i*2)=ugr(i*2)*cos_th+pgr(i*2)*sin_th
        ugr(i*2-1)=ugr(i*2-1)*cos_th+pgr(i*2-1)*sin_th
      enddo

       call Hpsi3(ugr,ughr,vion,dvp,Etot1,alpha,totNel,workr_n)

      do i=1,ng2_nL
       ugr(i*2)=(ugr(i*2)-pgr(i*2)*sin_th)/cos_th
       ugr(i*2-1)=(ugr(i*2-1)-pgr(i*2-1)*sin_th)/cos_th
      enddo


      E_thth=2*(Etot1-Etot-E_upg*theta0)/theta0**2

      theta=0.5d0*dabs(datan(E_upg/E_thth))
      cos_th=dcos(theta)
      sin_th=dsin(theta)

      predE=0.5d0*E_upg*dsin(2*theta)-
     &  0.25d0*E_thth*dcos(2*theta)+
     &  (Etot+0.25d0*E_thth)


c       if(inode.eq.1) then
c      write(6,*) "E_upg,dE",E_upg,(Etot1-Etot)/theta0
c      write(6,*) "Etot,predE,err,theta",Etot,predE,err,theta 
c      write(6,500) Etot,predE,err,theta 
c        endif
500   format("Etot,predE,err,theta", 
     &   2(f20.15,1x), 2(E10.5,1x))
**************************************************

      do i=1,ng2_nL
      ugr(i*2)=ugr(i*2)*cos_th+pgr(i*2)*sin_th
      ugr(i*2-1)=ugr(i*2-1)*cos_th+pgr(i*2-1)*sin_th
      enddo

1000  continue


      call d3fft_real2L(ugr,workr_n,-1,1)

      do i=1,nr_nL
        workr_n(i)=totNel*workr_n(i)**2
      enddo

*******************************************

      do i=1,nr_nL
      pgr(i)=vion(i)-vionT_nL(i)+UxcCA(workr_n(i)+rhocr_nL(i),uxc2)
      enddo
*** vion-vionT_nL=V_ion
*******************************************

      call d3fft_real2L(ugr,workr_n,1,1)

      do i=1,ng2_nL

       if(inode.eq.iorg2L(1).and.i.eq.iorg2L(2)) then

        ugr(i*2)=0.d0
        ugr(i*2-1)=0.d0
        else
        ugr(i*2)=ugr(i*2)*2*pi/gkk2_nL(i)
        ugr(i*2-1)=ugr(i*2-1)*2*pi/gkk2_nL(i)
       endif

      enddo


      call d3fft_real2L(ugr,workr_n,-1,1)

*******************************************
      s=0.d0

      do i=1,nr_nL
      pgr(i)=pgr(i)+workr_n(i)
      s=s+pgr(i)
      enddo

      call global_sumr(s)

      s=s/nrL

      do i=1,nr_nL
      pgr(i)=pgr(i)-s
      enddo

*******************************************
*** Now, the pgr is the new v_in potential
*******************************************

      do i=1,nr_nL
      ugr(i)=v_out(i)
      ughr(i)=v_in(i)
      enddo

      workr_n = pgr
      call d3fft_real2L(pgr,workr_n,1,1)
      workr_n = ughr
      call d3fft_real2L(ughr,workr_n,1,1)
      workr_n = ugr
      call d3fft_real2L(ugr,workr_n,1,1)


      Ecut0=0.5d0

      do i=1,ng2_nL
        x=0.8d0*gkk2_nL(i)/(gkk2_nL(i)+0.5d0)
        y=dexp(-gkk2_nL(i)/Ecut0)
        pgr(i*2)=pgr(i*2)*y+(1-y)*(ughr(i*2)*(1-x)+ugr(i*2)*x)
     &    -0.8d0*ugr(i*2)-0.2d0*ughr(i*2)   
        pgr(i*2-1)=pgr(i*2-1)*y+(1-y)*(ughr(i*2-1)*(1-x)+ugr(i*2-1)*x)
     &    -0.8d0*ugr(i*2-1)-0.2d0*ughr(i*2-1)   
      enddo

      call d3fft_real2L(pgr,workr_n,-1,1)

      s1=0.d0
      s2=0.d0

      do i=1,nr_nL
        workr_n(i)=workr_n(i)+0.2d0*v_in(i)+0.8d0*v_out(i)
        s1=s1+dabs(workr_n(i)-v_in(i))
        s2=s2+dabs(workr_n(i)-v_out(i))
        v_in(i)=workr_n(i)
**** the 0.2*v_in+0.8*v_out is used for charge mixing for gkk2_nL > Ecut2L
**** for the parts withing Ecut2L, they are cancelled by -0.8*ugr-0.2*ughr above
      enddo

      call global_sumr(s1)
      call global_sumr(s2)

      s1=s1/nrL
      s2=s2/nrL

      if(inode.eq.1) then 
        write(6,*) "ave, v_new-v_in and v_new-v_out from Thomas"
        write(6,*) s1,s2
      endif

      deallocate(workr_n)
      deallocate(rho_0)
      deallocate(ugr)
      deallocate(ughr)
      deallocate(pgr)
      deallocate(pgr_old)
      deallocate(vion)
      deallocate(dvp)

      return
      end

**************************************************
**************************************************

      subroutine Hpsi3(ugr,ur,vion,dvp,Etot,alpha,totNel,workr_n)

      use fft_data
      use load_data
      use data

      implicit double precision(a-h,o-z)
      include 'param.escan_real'

      real*8 rho_d(mr_nL)

      real*8 workr_n(mr_nL)

      real*8 ugr(mr_nL)

      real*8 ur(mr_nL)

      real*8 vion(mr_nL),dvp(mr_nL)


      ng2_nL=ngtotnod2L(inode)

      pi=4*datan(1.d0)
      fa=0.5d0*(3*pi**2)**(2.d0/3.d0)

      call d3fft_real2L(ugr,workr_n,-1,1)

      ur = workr_n

***************************************

        do i=1,nr_nL
         workr_n(i)=totNel*ur(i)**2
        enddo


        call d3fft_real2L(rho_d,workr_n,1,1)

      do i=1,ng2_nL

       if(inode.eq.iorg2L(1).and.i.eq.iorg2L(2)) then
        rho_d(i*2)=0.d0
        rho_d(i*2-1)=0.d0
       else
        rho_d(i*2)=rho_d(i*2)*2*pi/gkk2_nL(i)
        rho_d(i*2-1)=rho_d(i*2-1)*2*pi/gkk2_nL(i)
       endif

      enddo

      call d3fft_real2L(rho_d,workr_n,-1,1)

      rho_d = workr_n

***************************************

      Etot=0.d0

      do i=1,nr_nL

      d=totNel*ur(i)**2
      vtot=vion(i)+UxcCA(d+rhocr_nL(i),uxc2)+fa*d**(2.d0/3.d0)+
     &  rho_d(i)

      Etot=Etot+vion(i)*d+uxc2*(d+rhocr_nL(i))+fa*(3.d0/5.d0)*
     &  d**(5.d0/3.d0)+0.5d0*rho_d(i)*d
     & +2*dvp(i)*totNel*ur(i)

      workr_n(i)=vtot*ur(i)+dvp(i)
      enddo

      call global_sumr(Etot)

      Etot=Etot*vol/nrL


      call d3fft_real2L(ur,workr_n,1,1)

      Ekk=0.d0

      do i=1,ng2_nL

        ur(i*2)=ur(i*2)+alpha*gkk2_nL(i)*ugr(i*2)
        ur(i*2-1)=ur(i*2-1)+alpha*gkk2_nL(i)*ugr(i*2-1)

        Ekk=Ekk+alpha*gkk2_nL(i)*(ugr(i*2-1)**2+ugr(i*2)**2)

      enddo

      call global_sumr(Ekk)

      Etot=Etot+Ekk*vol*totNel*2



      return
      end


