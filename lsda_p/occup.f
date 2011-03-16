      subroutine occup(dE,E_st,totNel,nkpt,
     & occ,E_NSC,TS,workr_n,Ef,islda)
****************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
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
      y2=(E_st(m,kpt,iislda)-Emin)/dE
      if(y1.gt.100.d0) y1=100.d0
      if(y1.lt.-100.d0) y1=-100.d0
      if(y2.gt.100.d0) y2=100.d0
      if(y2.lt.-100.d0) y2=-100.d0
      s1=s1+1.d0/(dexp(y1)+1.d0)*weighkpt(kpt)
      s2=s2+1.d0/(dexp(y2)+1.d0)*weighkpt(kpt)
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
      if(y.gt.100.d0) y=100.d0
      if(y.lt.-100.d0) y=-100.d0
      s=s+1.d0/(dexp(y)+1.d0)*weighkpt(kpt)
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

ccccccc TS is the kT*entropy term for free energy
ccccccc This term is needed in order for the total
ccccccc energy (free energy) be the local minimont at
ccccccc the Schrodinger's equation.
      TS=0.d0
      do iislda=1,islda
      do kpt=1,nkpt
      do m=1,mx
      y=(E_st(m,kpt,iislda)-Ef)/dE
      if(y.gt.100.d0) y=100.d0
      if(y.lt.-100.d0) y=-100.d0
      yy=1.d0/(dexp(y)+1.d0)
      occ(m,kpt,iislda)=2.d0/islda*yy*weighkpt(kpt)
      TS=TS+2.d0/islda*weighkpt(kpt)*dE*
     & ((dabs(1.d0-yy)+1.D-30)*dlog(dabs(1.d0-yy)+1.D-30)+
     & (dabs(yy)+1.D-30)*dlog(dabs(yy)+1.D-30))
      enddo
      enddo
      enddo
***********************************************
***********************************************

*** get the charge density
      do 100 iislda=1,islda
      do i=1,nr_n
      rho_n(i,iislda)=0.d0
      enddo

      do kpt=1,nkpt

      call gen_G_comp(kpt,0)
      call fftprep_comp(n1,n2,n3)

      call ugIO(ug_n,kpt,2,0,iislda)

      do m=1,mx

      call d3fft_comp(ug_n(1,m),workr_n,-1,kpt)

      do i=1,nr_n
      rho_n(i,iislda)=rho_n(i,iislda)+
     &     occ(m,kpt,iislda)*cdabs(workr_n(i))**2
      enddo

      enddo
      enddo
100   continue

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



      
