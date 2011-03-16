******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
******************************************

*******************************************
*** This program using (v_in+dv, rho) -> (v_out+dv) to fixed the
*** Hamiltonian: vion (so V_hxc(rho)+vion=v_out+dv)
***              dvp: (so v_in+dv & rho, Ef0 satisfy the solution eq.)
***  Then, it resolve rho', so that v_in'=v_out'
***  Finally, the v_in' is used for small k 
***  and kerker mixing is used for large k
****************************************************************
*** The final equation for psi=dsqrt(rho)
*** -0.5 \nabla^2 psi +( rho^(2/3)+ v_in)*psi + dvp = Ef*psi
***  and v_in= vion + v_Hxc(rho)
****************************************************************

      subroutine Thomas3(v_in,v_out,Ef0,ntot,workr_n,islda)

*****************************************************************
*** in return, v_in will be replaced by new v_in
*** rho, v_out, dv are not changed
*****************************************************************

      use fft_data
      use load_data
      use data

      implicit double precision(a-h,o-z)
      include 'param.escan_real'

      real*8 workr_n(mr_n)

      real*8 rho_0(mr_n)
      real*8 v_in(mr_n),v_out(mr_n)
********* inputs
      real*8 vion(mr_n),dvp(mr_n)
********* vion, dvp are derived from input v_in,v_out,rho

      real*8 ugr(mr_n)
      real*8 ughr(mr_n)
      real*8 pgr(mr_n)
      real*8 pgr_old(mr_n)

      ng2_n=ngtotnod2(inode)

      if(inode.eq.1) then
        write(6,*) "Ef0=", Ef0
      endif

      pi=4*datan(1.d0)
      alpha=1.d0
      pgr_old = 0.0d0

      rho_0=0.d0

      do iislda=1,islda
      do i=1,nr_n
      rho_0(i) = rho_0(i)+ rho_n(i,iislda)
      enddo
      enddo
********************************************
*** make rho(i) = workr(i)**2 and workr(i) is in the Ecut2
********************************************
      s1=0.d0
 
      do i=1,nr_n
       s1=s1+rho_0(i)
       workr_n(i)=dsqrt(rho_0(i))
      enddo

      call global_sumr(s1)

      call d3fft_real2(pgr,workr_n,1,1)
      call d3fft_real2(pgr,workr_n,-1,1)

      s2=0.d0

      do i=1,nr_n
        rho_0(i)=workr_n(i)**2
        s2=s2+rho_0(i)
      enddo

      call global_sumr(s2)

      s=s1/s2
      do i=1,nr_n
        rho_0(i)=rho_0(i)*s
      enddo

*********************************************
*** here, use pgr,pgi as work functions
*********************************************

      do i=1,nr_n
        workr_n(i) = rho_0(i)
      enddo

      call d3fft_real2(pgr,workr_n,1,1)

      do i=1,ng2_n

       if(inode.eq.iorg2(1).and.i.eq.iorg2(2)) then

        pgr(i*2)=0.d0
        pgr(i*2-1)=0.d0
        else
        pgr(i*2)=pgr(i*2)*2*pi/gkk2_n(i)
        pgr(i*2-1)=pgr(i*2-1)*2*pi/gkk2_n(i)
       endif

      enddo


      call d3fft_real2(pgr,workr_n,-1,1)

      pgr = workr_n

      fa=0.5d0*(3*pi**2)**(2.d0/3.d0)

      do i=1,nr_n
c      pgr(i)=pgr(i)+UxcMik(rho_0(i),uxc2)
       pgr(i)=pgr(i)+UxcCA(rho_0(i)+rhocr_n(i),uxc2)
       vion(i)=(v_out(i)+vionT_n(i))-pgr(i)
      enddo

      do i=1,nr_n
       pgr(i)=fa*rho_0(i)**(2.d0/3.d0)+(v_in(i)+vionT_n(i))
      enddo

ccccccccccccccccccccccccccccccccccccccccccccccc

      do i=1,nr_n
      workr_n(i)=dsqrt(rho_0(i))
      enddo


      call d3fft_real2(ugr,workr_n,1,1)

      do i=1,ng2_n
        ugr(i*2)=alpha*gkk2_n(i)*ugr(i*2)
        ugr(i*2-1)=alpha*gkk2_n(i)*ugr(i*2-1)
      enddo


      call d3fft_real2(ugr,workr_n,-1,1)

      ugr = workr_n

      do i=1,nr_n
       dvp(i)=(Ef0-pgr(i))*dsqrt(rho_0(i))-ugr(i)
       dvp(i)=dvp(i)/dsqrt(1.d0*ntot)
      enddo

ccccccccccccccccccccccccccccccccccccccccccccccc
ccc finished obtaining vion(r) and dvp(r)
ccc Now, the Hamiltonian has be determined, so just need
ccc to solve the Hamiltonian selfconsistently
ccccccccccccccccccccccccccccccccccccccccccccccc
      s1=0.d0
      s2=0.d0
      s3=0.d0

      do i=1,nr_n
      s1=s1+dvp(i)**2
      s2=s2+rho_0(i)
      s3=s3+dabs(v_in(i)+vionT_n(i))
      enddo

      call global_sumr(s1)
      call global_sumr(s2)
      call global_sumr(s3)

      s1=dsqrt(s1/s2)
      s3=s3/nr

      if(inode.eq.1) then
        write(6,*) "average dvp (Hartree)", s1
        write(6,*) "average v_in+dv", s3
      endif
******************************************************
*****, Now, the potential vion(i) and dvp(i) has been calculated
******************************************************

       do i=1,nr_n
       workr_n(i)=dsqrt(rho_0(i)/ntot)
       enddo

       call d3fft_real2(ugr,workr_n,1,1)

*****************************************************
*****************************************************
ccccc ugr,ugi in G space

       do 1000 nint=1,nmax

       if(nint.eq.1) rr0=1.D+40

       s=0.d0

       do i=1,ng2_n
         s=s+ugr(i*2-1)**2+ugr(i*2)**2
       enddo

       call global_sumr(s)

       s=s*vol*2
       s=dsqrt(1.d0/s)

       do i=1,ng2_n
        ugr(i*2)=s*ugr(i*2)
        ugr(i*2-1)=s*ugr(i*2-1)
       enddo

******************
ccccc calculate H*u and Etot
       call Hpsi3(ugr,ughr,vion,dvp,Etot,alpha,ntot,workr_n)
******************

      do i=1,ng2_n
       pgr(i*2)=ughr(i*2)
       pgr(i*2-1)=ughr(i*2-1)
      enddo

      Ef=0.d0

      do i=1,ng2_n
       Ef=Ef+pgr(i*2-1)*ugr(i*2-1)+pgr(i*2)*ugr(i*2)
      enddo

      call global_sumr(Ef)

      Ef=Ef*vol*2

      err=0.d0

      do i=1,ng2_n
        pgr(i*2)=pgr(i*2)-Ef*ugr(i*2)
        pgr(i*2-1)=pgr(i*2-1)-Ef*ugr(i*2-1)
        err=err+pgr(i*2)**2+pgr(i*2-1)**2
      enddo

      call global_sumr(err)

      err=dsqrt(err*vol*2)

      if((nint.eq.1.or.nint.eq.nmax).and.inode.eq.1) then
      write(6,*) "Ef,err",Ef,err
      endif

      Ek=0.5d0
      rr1=0.d0

      do i=1,ng2_n
        x=gkk2_n(i)/Ek
        y=1/(1+x**4)**0.25d0
        rr1=rr1+y*(pgr(i*2-1)**2+pgr(i*2)**2)
      enddo

      call global_sumr(rr1)

      rr1=rr1*vol*2

      beta=rr1/rr0

      rr0=rr1

      do i=1,ng2_n
       x=gkk2_n(i)/Ek
       y=1/(1+x**4)**0.25d0
       pgr(i*2)=-pgr(i*2)*y+beta*pgr_old(i*2)
       pgr(i*2-1)=-pgr(i*2-1)*y+beta*pgr_old(i*2-1)
      enddo

      s=0.d0

      do i=1,ng2_n
      s=s+pgr(i*2)*ugr(i*2)+pgr(i*2-1)*ugr(i*2-1)
      enddo

      call global_sumr(s)

      s=s*vol*2

      do i=1,ng2_n
        pgr(i*2)=pgr(i*2)-s*ugr(i*2)
        pgr(i*2-1)=pgr(i*2-1)-s*ugr(i*2-1)
      enddo

      s=0.d0

      do i=1,ng2_n
       s=s+pgr(i*2)**2+pgr(i*2-1)**2
      enddo

      call global_sumr(s)

      s=1.d0/dsqrt(s*vol*2)

      do i=1,ng2_n

        pgr_old(i*2)=pgr(i*2)
        pgr_old(i*2-1)=pgr(i*2-1)

        pgr(i*2)=pgr(i*2)*s
        pgr(i*2-1)=pgr(i*2-1)*s

      enddo

**************************************************

      E_upg=0.d0

      do i=1,ng2_n
      E_upg=E_upg+pgr(i*2)*ughr(i*2)+pgr(i*2-1)*ughr(i*2-1)
      enddo

      call global_sumr(E_upg)

      E_upg=E_upg*vol*2
      E_upg=2*E_upg*ntot

**************************************************
ccccc  theta0=0.2D-3, might cause slightly different results for
ccccc  parallel and serie code due to round off error.
cccccc for debug, use theta0=0.2D-2,but for actual run, use 0.2D-03.
**************************************************
      theta0=0.2D-03
c      theta0=0.2D-02

      cos_th=dcos(theta0)
      sin_th=dsin(theta0)

      do i=1,ng2_n
        ugr(i*2)=ugr(i*2)*cos_th+pgr(i*2)*sin_th
        ugr(i*2-1)=ugr(i*2-1)*cos_th+pgr(i*2-1)*sin_th
      enddo

       call Hpsi3(ugr,ughr,vion,dvp,Etot1,alpha,ntot,workr_n)

      do i=1,ng2_n
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
500   format(2(f20.15,1x), 2(E10.5,1x))
**************************************************

      do i=1,ng2_n
      ugr(i*2)=ugr(i*2)*cos_th+pgr(i*2)*sin_th
      ugr(i*2-1)=ugr(i*2-1)*cos_th+pgr(i*2-1)*sin_th
      enddo

1000  continue


      call d3fft_real2(ugr,workr_n,-1,1)

      do i=1,nr_n
        workr_n(i)=ntot*workr_n(i)**2
      enddo

*******************************************

      do i=1,nr_n
c      pgr(i)=vion(i)+UxcMik(workr_n(i),uxc2)
      pgr(i)=vion(i)+UxcCA(workr_n(i)+rhocr_n(i),uxc2)
      enddo

*******************************************

      call d3fft_real2(ugr,workr_n,1,1)

      do i=1,ng2_n

       if(inode.eq.iorg2(1).and.i.eq.iorg2(2)) then

        ugr(i*2)=0.d0
        ugr(i*2-1)=0.d0
        else
        ugr(i*2)=ugr(i*2)*2*pi/gkk2_n(i)
        ugr(i*2-1)=ugr(i*2-1)*2*pi/gkk2_n(i)
       endif

      enddo


      call d3fft_real2(ugr,workr_n,-1,1)

*******************************************
      s=0.d0

      do i=1,nr_n
      pgr(i)=pgr(i)+workr_n(i)-vionT_n(i)
      s=s+pgr(i)
      enddo

      call global_sumr(s)

      s=s/nr

      do i=1,nr_n
      pgr(i)=pgr(i)-s
      enddo

*******************************************
*** Now, the pgr is the new v_in potential
*******************************************

      do i=1,nr_n
      ugr(i)=v_out(i)
      ughr(i)=v_in(i)
      enddo

      workr_n = pgr
      call d3fft_real2(pgr,workr_n,1,1)
      workr_n = ughr
      call d3fft_real2(ughr,workr_n,1,1)
      workr_n = ugr
      call d3fft_real2(ugr,workr_n,1,1)


      Ecut0=0.5d0

      do i=1,ng2_n
        x=0.8d0*gkk2_n(i)/(gkk2_n(i)+0.5d0)
        y=dexp(-gkk2_n(i)/Ecut0)
        pgr(i*2)=pgr(i*2)*y+(1-y)*(ughr(i*2)*(1-x)+ugr(i*2)*x)
        pgr(i*2-1)=pgr(i*2-1)*y+(1-y)*(ughr(i*2-1)*(1-x)+ugr(i*2-1)*x)
      enddo

      call d3fft_real2(pgr,workr_n,-1,1)

      s1=0.d0
      s2=0.d0

      do i=1,nr_n
        s1=s1+dabs(workr_n(i)-v_in(i))
        s2=s2+dabs(workr_n(i)-v_out(i))
        v_in(i)=workr_n(i)
      enddo

      call global_sumr(s1)
      call global_sumr(s2)

      s1=s1/nr
      s2=s2/nr

      if(inode.eq.1) then 
        write(6,*) "ave, v_new-v_in and v_new-v_out from Thomas"
        write(6,*) s1,s2
      endif

      return
      end

**************************************************
**************************************************

      subroutine Hpsi3(ugr,ur,vion,dvp,Etot,alpha,ntot,workr_n)

      use fft_data
      use load_data
      use data

      implicit double precision(a-h,o-z)
      include 'param.escan_real'

      real*8 rho_d(mr_n)

      real*8 workr_n(mr_n)

      real*8 ugr(mr_n)

      real*8 ur(mr_n)

      real*8 vion(mr_n),dvp(mr_n)


      ng2_n=ngtotnod2(inode)

      pi=4*datan(1.d0)
      fa=0.5d0*(3*pi**2)**(2.d0/3.d0)

      call d3fft_real2(ugr,workr_n,-1,1)

      ur = workr_n

***************************************

        do i=1,nr_n
         workr_n(i)=ntot*ur(i)**2
        enddo


        call d3fft_real2(rho_d,workr_n,1,1)

      do i=1,ng2_n

       if(inode.eq.iorg2(1).and.i.eq.iorg2(2)) then
        rho_d(i*2)=0.d0
        rho_d(i*2-1)=0.d0
       else
        rho_d(i*2)=rho_d(i*2)*2*pi/gkk2_n(i)
        rho_d(i*2-1)=rho_d(i*2-1)*2*pi/gkk2_n(i)
       endif

      enddo

      call d3fft_real2(rho_d,workr_n,-1,1)

      rho_d = workr_n

***************************************

      Etot=0.d0

      do i=1,nr_n

      d=ntot*ur(i)**2
c      vtot=vion(i)+UxcMik(d,uxc2)+fa*d**(2.d0/3.d0)+
c      &  rho_d(i)
      vtot=vion(i)+UxcCA(d+rhocr_n(i),uxc2)+fa*d**(2.d0/3.d0)+
     &  rho_d(i)

      Etot=Etot+vion(i)*d+uxc2*(d+rhocr_n(i))+fa*(3.d0/5.d0)*
     &  d**(5.d0/3.d0)+0.5d0*rho_d(i)*d
     & +2*dvp(i)*ntot*ur(i)

      workr_n(i)=vtot*ur(i)+dvp(i)
      enddo

      call global_sumr(Etot)

      Etot=Etot*vol/nr


      call d3fft_real2(ur,workr_n,1,1)

      Ekk=0.d0

      do i=1,ng2_n

        ur(i*2)=ur(i*2)+alpha*gkk2_n(i)*ugr(i*2)
        ur(i*2-1)=ur(i*2-1)+alpha*gkk2_n(i)*ugr(i*2-1)

        Ekk=Ekk+alpha*gkk2_n(i)*(ugr(i*2-1)**2+ugr(i*2)**2)

      enddo

      call global_sumr(Ekk)

      Etot=Etot+Ekk*vol*ntot*2


      return
      end


