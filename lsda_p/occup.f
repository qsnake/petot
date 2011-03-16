      subroutine occup(dE,itypeFermi,E_st,totNel,nkpt,
     & occ,E_NSC,TS,workr_n,Ef,islda)
****************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************

******************************************
      use fft_data
      use load_data
      use data
      implicit double precision (a-h,o-z)

      include 'mpif.h'
      include 'param.escan_real'
      
      real*8 Ef,occ(mst,nkpt,islda),E_st(mst,nkpt,islda)
      complex*16 workr_n(mr_n)
      real*8, allocatable, dimension (:) :: rho_n_tmp
******************************************

     
      n=totNel/2+0.1d0
      if(n.gt.mx) then
      write(6,*) "n.gt.mx, stop", n,mx,mst
      stop
      endif

      Ef=E_st(n,1,1)
      if(mx.gt.n) Ef=(E_st(n,1,1)+E_st(n+1,1,1))/2

****** first, using three point Emax, Emin, then try (Emax+Emin)/2
****** only at the final steps, using the following derivatives.

      Emax=-10000.d0
      Emin=10000.d0
      do iislda=1,islda
      do kpt=1,nkpt
      do m=1,mx
      if(E_st(m,kpt,iislda).gt.Emax) Emax=E_st(m,kpt,iislda)
      if(E_st(m,kpt,iislda).lt.Emin) Emin=E_st(m,kpt,iislda)
      enddo
      enddo
      enddo
      Emax=Emax+20*dE
      Emin=Emin-20*dE

      s1=0.d0
      s2=0.d0
      do iislda=1,islda
      do kpt=1,nkpt
      do m=1,mx
      y1=(E_st(m,kpt,iislda)-Emax)/dE
      call Fermi(f_occ,S_occ,y1,itypeFermi)
      s1=s1+f_occ*weighkpt(kpt)

      y2=(E_st(m,kpt,iislda)-Emin)/dE
      call Fermi(f_occ,S_occ,y2,itypeFermi)
      s2=s2+f_occ*weighkpt(kpt)

      enddo
      enddo
      enddo
      s1=s1*2/islda
      s2=s2*2/islda

      if(s1.le.totNel) then
      Ef=Emax
      s=s1
      goto 90
      endif

      if(s2.ge.totNel) then
      Ef=Emin
      s=s2
      goto 90
      endif

      do it=1,100
      Ef=(Emax+Emin)/2
      s=0.d0
      do iislda=1,islda
      do kpt=1,nkpt
      do m=1,mx
      y=(E_st(m,kpt,iislda)-Ef)/dE
      call Fermi(f_occ,S_occ,y,itypeFermi)
      s=s+f_occ*weighkpt(kpt)
      enddo
      enddo
      enddo
      s=s*2/islda
      if(s.ge.totNel) then
      Emax=Ef
      else
      Emin=Ef
      endif
      enddo

      Ef=(Emin+Emax)/2

90    continue
cc      write(6,*) "Ef,s=", Ef*27.211396,s
      occ_charge=s

ccccccc TS is the kT*entropy term for free energy, F=E-TS
ccccccc This term is needed in order for the total
ccccccc energy (free energy) be the local minimont at
ccccccc the Schrodinger's equation.
cccccccccccc The current formula is only correct for Fermi-Dirac distribution
      TS=0.d0
      do iislda=1,islda
      do kpt=1,nkpt
      do m=1,mx
      y=(E_st(m,kpt,iislda)-Ef)/dE
      call Fermi(yy,S_occ,y,itypeFermi)
      occ(m,kpt,iislda)=2.d0/islda*yy*weighkpt(kpt)
      TS=TS+2.d0/islda*weighkpt(kpt)*dE*S_occ
      enddo
      enddo
      enddo
***********************************************
cccccccccccccccccccccccccccccccccccccccccccccccccccc
***********************************************

*** get the charge density
      do 100 iislda=1,islda
      do i=1,nr_n
      rho_n(i,iislda)=0.d0
      enddo

      do kpt=1,nkpt

       if((iislda-1)*nkpt+kpt.ge.kpt_slda_dis(1).and.
     &     (iislda-1)*nkpt+kpt.le.kpt_slda_dis(2)) then      ! the big division of work

      call gen_G_comp(kpt,0)
      call fftprep_comp(n1,n2,n3)

      call ugIO(ug_n,kpt,2,0,iislda)

      do m=1,mx

      call d3fft_comp(ug_n(1,m),workr_n,-1,kpt)

      do i=1,nr_n
      rho_n(i,iislda)=rho_n(i,iislda)+
     &     occ(m,kpt,iislda)*cdabs(workr_n(i))**2      ! rho_n could be negative for itypeFermi=21,22,etc
      enddo

      enddo      ! m=1,mx

       endif     ! kpt_slda_dis(1),(2)
      enddo
100   continue

cccccccccccccccccccccccccccccccccccccccccccccccc
      allocate(rho_n_tmp(nr_n))
      do iislda=1,islda
      call mpi_allreduce(rho_n(1,iislda),rho_n_tmp,nr_n,
     & MPI_REAL8,MPI_SUM,MPI_COMM_N,ierr)
      do i=1,nr_n
      rho_n(i,iislda)=rho_n_tmp(i)
      enddo
      enddo
      deallocate(rho_n_tmp)
ccccccccccccccccccccccccccccccccccccccccccccccccc


      E_NSC=0.d0
      do iislda=1,islda
      do kpt=1,nkpt
      do m=1,mx
      E_NSC=E_NSC+occ(m,kpt,iislda)*E_st(m,kpt,iislda)
      enddo
      enddo
      enddo
*********************************************

      return
      end


      
      subroutine Fermi(f_occ,S_occ,x,itypeFermi)
      implicit double precision(a-h,o-z)

cccccccccc   According to: 
cccc  (1) M. Methfessel, A.T. Paxton, Phys. Rev. B 40, 3616 (1989)
cccc  (2) G. Kresse, J. Furthmuller, Comp. Mat. Sci. 6, 15 (1996). 

      real*8 h(0:30)


      y=x

      if(itypeFermi.eq.1) then       !  Fermi-Diract smearing
      if(y.gt.100.d0) y=100.d0
      if(y.lt.-100.d0) y=-100.d0
      f_occ=1.d0/(dexp(y)+1.d0)
      S_occ=-((dabs(1.d0-f_occ)+1.D-30)*dlog(dabs(1.d0-f_occ)+
     & 1.D-30)+(dabs(f_occ)+1.D-30)*dlog(dabs(f_occ)+1.D-30))
      return
      endif

ccc Using M. Methfessel, A.T. Paxton, P.R.B, 40, 3616 (1989). 
      if(itypeFermi.ge.20.and.itypeFermi.lt.30) then 
      n=itypeFermi-20    ! the order of S_n function
      if(y.gt.9.d0) y=9.d0 
      if(y.lt.-9.d0) y=-9.d0
      f_occ=0.5d0*(1-erf(y))
      fg=exp(-y**2)

      h(0)=1.d0
      h(1)=2*y
      do i=1,2*n
      h(i+1)=2*y*h(i)-2*i*h(i-1)
      enddo

      a=1.d0/dsqrt(4*datan(1.d0))

      do i=1,n
      a=-a/(i*4)
      f_occ=f_occ+a*h(2*i-1)*fg
      enddo
      S_occ=0.5d0*a*h(2*n)*fg

      if(y.lt.-8.99) then
      f_occ=1.d0
      S_occ=0.d0
      endif

      return
      endif
ccccccccccccccccccccccccccccccccccccc

      write(6,*) "itypeFermi not defined, stop",itypeFermi
      stop
      end
      
      
      
