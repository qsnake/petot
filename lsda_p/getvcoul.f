      subroutine getvcoul(ntype,Ealpha,icoul)
***************************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************
****** This program calculate the modified v(r) for a cluster of slab.
****** For icoul=1, it does not have periodicity at all three directions. 
****** For icoul=11, it does not have periodicity at the first n1 direction
****** For icoul=12, it does not have periodicity at the second n2 direction
****** For icoul=13, it does not have periodicity at the third n3 direction
************************************************
******  For icoul=1, it takes the periodic v'(r), then subtract \sum'_R 1/(r-R), 
******  sum' means R not equal 0 (the original one). So what left v(r)=v'(r)-\sum'_R 1/(r-R) 
******  should be similar to 1/r, but the values are smoothed near r=0 from the Fourier 
******  space representation. 
******  To do \sum'_R 1/(r-R), we have used the Ewald's method. Basically, 
******  \sum'_R 1/(r-R)= \sum'_R (1/(r-R)-G(r-R)/(r-R))+ \sum_R G(r-R)/(r-R) -G(r)/r
******  here G(r-R)/(r-R) denotes the potential due to a spherical Gaussian charge at R. 
*****   The first term is sum over real space \sum'_R, the second term \sum_R G(r-R)/(r-R) 
*****   is calculated in Fourier space, and the last term G(r)/r is calculated in real space. 
*******************************************
*****   For icoul=11,12,13, first, take v(r) from the icoul=1 calculation, then
*****   v_sl(r)=v(r)+\sum'_Rp 1/(r-Rp).   Here Rp are the periodic points in the 
*****   two dimensional plane. \sum' means Rp.ne.0.    Again, the above sum is done using 
*****   Ewald like technique. The details follow the formula in: 
*****   D.M. Heyes, M. Barber, J.H.R. Clarke, J. Chem. Soc. Faraday Trans. 2: 73, 1485(1977). 
***********************************************

       use fft_data
       use load_data
       use data

      implicit double precision (a-h,o-z)
      include 'param.escan_real'
      include "mpif.h"

      real*8 zatom(mtype),dv_corr(0:mtype)
      integer, allocatable, dimension(:) :: ncount_tmp

      real*8, allocatable, dimension(:)  :: workr_nL2

      real*8 qi2(mnq),vq(mnq,mtype),rhoq(mnq,mtype)
      real*8 vqT(mnq,mtype),rhocq(mnq,mtype)
      real*8 Ealpha(mtype)
      real*8 AL2D(3,2),ALI2D(3,2)

      character*20 file_tmp


      common /comVrho/qi2,vq,rhoq,vqT,rhocq
      common /comzatom/zatom


cccccccccccccccccccccccccccccccccccccccccccccccccc
       allocate(workr_nL2(mr_nL2))
       vins2=1.d0/vol2
       ng2_nL2=ngtotnod2L2(inode)
       pi=4*datan(1.d0)


 
       vcoul_nL2=0.d0
       do 100  itype=0,ntype      ! itype=0, plain coulomb, 4pi/q**2

       do 10 i=1,ng2_nL2

      q=dsqrt(gkx2_nL2(i)**2+gky2_nL2(i)**2+gkz2_nL2(i)**2)

ccccccc vq(iq),qi2(iq) not defined for 0.5*q**2 > Ecut2L

      if(0.5*q**2.ge.Ecut2L) goto 10  

      iq=1+q*(mnq-1.d0)/qi2(mnq)

      x=(q-qi2(iq))/(qi2(iq+1)-qi2(iq))    ! assuming equal distance grid
 
      f1=1-x-0.5d0*x*(1-x)
      f2=x+x*(1-x)
      f3=-0.5d0*x*(1-x)       ! using quadratic interpolation based on vq*q**2

      amask=1.d0
      if(0.5d0*q**2.gt.Ecut2L*0.9d0) then          
      x=(0.5d0*q**2-Ecut2L*0.9d0)/(0.1d0*Ecut2L)
      amask=(dcos(x*pi)+1.d0)/2.d0              ! provide a smooth cutoff around Ecut2L, so the real space v(r) should be smooth, approach to 1/r at large distance
      endif

      if(itype.ge.1) then
      y=vq(iq,itype)*qi2(iq)**2*f1+vq(iq+1,itype)*qi2(iq+1)**2*f2+
     &   vq(iq+2,itype)*qi2(iq+2)**2*f3
      else
      y=4*pi       ! itype=0, v(q)=4pi/q**2
      endif

 

      if(q.lt.1.D-6) then
      y=0.d0       ! v(q=0)=0, the constant to be settled later. 
      else
      y=y/q**2
      endif

      y=y*amask

      i2=i*2


      vcoul_nL2(i2-1,itype)=vcoul_nL2(i2-1,itype)+y*vins2
      vcoul_nL2(i2,itype)=vcoul_nL2(i2,itype)+0.d0          ! it is real in G space, sph. symmetric


 10   continue

      call d3fft_real2L2(vcoul_nL2(1,itype),workr_nL2,-1,0)
      vcoul_nL2(:,itype) = workr_nL2(:)      ! this is v(r) with the atom at the origin

100   continue

***************************************************
***  The vcoul_nL2 obtained above is the results of periodic images by the AL2 box. 
***  Now, for icoul=1, we want to get ride of this periodicity, by subtract out all the
***  1/r potentials from all the image charges. The sum of all these image charges are 
***  calculated from Ewald sum method. 

ccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccc

      d1=1.d0/dsqrt(ALI2(1,1)**2+ALI2(2,1)**2+ALI2(3,1)**2)
      d2=1.d0/dsqrt(ALI2(1,2)**2+ALI2(2,2)**2+ALI2(3,2)**2)
      d3=1.d0/dsqrt(ALI2(1,3)**2+ALI2(2,3)**2+ALI2(3,3)**2)

      d1q=1.d0/dsqrt(AL2(1,1)**2+AL2(2,1)**2+AL2(3,1)**2)
      d2q=1.d0/dsqrt(AL2(1,2)**2+AL2(2,2)**2+AL2(3,2)**2)
      d3q=1.d0/dsqrt(AL2(1,3)**2+AL2(2,3)**2+AL2(3,3)**2)

ccccc This beta balances the real and G space calculations.
ccccc The computational time in both real and G space are O(N^2).
ccccc In both real and G space, the calc. propto ~ 1000*natom**2

      beta=1.2*pi**0.5*
     &     (d1q*d2q*d3q/(d1*d2*d3))**0.16666666

ccccc this beta is to make nc1*nc2*nc3=nq1*nq2*nq3

c      cut=5.5d0     !  exp(-cut**2)
      cut=5.0d0     !  exp(-cut**2)

      cut_gkk=4*(cut*beta)**2/(2*pi)**2
      cut_rr=(cut/beta)**2
      fac_gkk=(2*pi)**2/(4*beta**2)

      nc1=dsqrt(cut_rr)/d1+0.25 
      nc2=dsqrt(cut_rr)/d2+0.25 
      nc3=dsqrt(cut_rr)/d3+0.25 

      nq1=dsqrt(cut_gkk)/d1q  
      nq2=dsqrt(cut_gkk)/d2q 
      nq3=dsqrt(cut_gkk)/d3q  

      dv=pi/beta**2/vol2


      if(inode_tot.eq.1) then
      write(6,*) "nc1,2,3, nq1,2,3,beta:", 
     &  nc1,nc2,nc3,nq1,nq2,nq3,beta
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccc

      do 2000 ii=1,nr_nL2 
      jj=ii+(inode-1)*nr_nL2
      i1=(jj-1)/(n2L2*n3L2)+1
      j1=(jj-1-(i1-1)*n2L2*n3L2)/n3L2+1
      k1=jj-(i1-1)*n2L2*n3L2-(j1-1)*n3L2

      if(i1.le.n1L2/2+1) then
      idx1=i1-1
      else
      idx1=i1-n1L2-1
      endif

      if(j1.le.n2L2/2+1) then
      idx2=j1-1
      else
      idx2=j1-n2L2-1
      endif

      if(k1.le.n3L2/2+1) then
      idx3=k1-1
      else
      idx3=k1-n3L2-1
      endif

      x0=AL2(1,1)*idx1/n1L2+AL2(1,2)*idx2/n2L2+
     &     AL2(1,3)*idx3/n3L2
      y0=AL2(2,1)*idx1/n1L2+AL2(2,2)*idx2/n2L2+
     &     AL2(2,3)*idx3/n3L2
      z0=AL2(3,1)*idx1/n1L2+AL2(3,2)*idx2/n2L2+
     &     AL2(3,3)*idx3/n3L2

      sr=0.d0
      do 200 i=-nc1,nc1     ! this is the real space summation, rho(r)=gauss(r-R)-delta(r-R)
      do 200 j=-nc2,nc2
      do 200 k=-nc3,nc3
      x=x0+AL2(1,1)*i+AL2(1,2)*j+AL2(1,3)*k
      y=y0+AL2(2,1)*i+AL2(2,2)*j+AL2(2,3)*k
      z=z0+AL2(3,1)*i+AL2(3,2)*j+AL2(3,3)*k
      rr=x**2+y**2+z**2
      if((i.ne.0.or.j.ne.0.or.k.ne.0).and.rr.gt.cut_rr) goto 200

      r=dsqrt(rr)

      if(i.eq.0.and.j.eq.0.and.k.eq.0) then  
***  special treatment for i,j,k=0, the point itself,not the images
***  for this point, only add the potential from gauss(r) [(erfc(beta*r)-1)/r], not the delta(r-R)
       if(r.lt.1.D-6) r=1.D-6
       if(rr.gt.cut_rr) then
       sr=sr-1.d0/r
       else
       sr=sr+(erfc(beta*r)-1.d0)/r
       endif
      else
**** for all the other image points, add the potential from gauss(r-R)-delta(r-R) [erfc(beta*r)/r]
      sr=sr+erfc(beta*r)/r
      endif

200   continue

      sr=sr-dv    ! the important constant, to correct the v(q=0)=0 overall shifting

      sq=0.d0
      do 300 i=0,nq1     ! this is q space summation, subtract out the lattice gauss(r-R) potential
      do 300 j=-nq2,nq2
      do 300 k=-nq3,nq3
cccc the summation over nq1,nq2,nq3 is done for each grid points, so no FFT etc. (It actually 
cccc can be done by FFT !) 

      if(i.eq.0.and.j.lt.0) goto 300
      if(i.eq.0.and.j.eq.0.and.k.lt.0) goto 300
      if(i.eq.0.and.j.eq.0.and.k.eq.0) goto 300

      gkx=ALI2(1,1)*i+ALI2(1,2)*j+ALI2(1,3)*k
      gky=ALI2(2,1)*i+ALI2(2,2)*j+ALI2(2,3)*k
      gkz=ALI2(3,1)*i+ALI2(3,2)*j+ALI2(3,3)*k
     
      gkk=gkx**2+gky**2+gkz**2        ! gkk,gkx,gky,gkz missing 2*pi

      if(gkk.gt.cut_gkk) goto 300

      ph=(x0*gkx+y0*gky+z0*gkz)*2*pi

      sq=sq+dexp(-gkk*fac_gkk)*dcos(ph)/gkk

300   continue

      sq=sq*4*pi/vol2*2.d0/(2*pi)**2

      vcoul_nL2(ii,0)=vcoul_nL2(ii,0)-sr-sq  

      do itype=1,ntype
      vcoul_nL2(ii,itype)=vcoul_nL2(ii,itype)+
     &                    zatom(itype)*(sr+sq)  
      enddo

2000  continue

ccccccccccc  after this dv_corr, the vcoul_nL2(ii,0) is like 1/r (no shift)
ccccccccccc                      the vcoul_nL2(ii,itype) is like v_loc(r,itype) (no shift)
ccccccccccc                      There is no need for the original Ealpha

       dv_corr(0)=0.d0
       do itype=1,ntype
       dv_corr(itype)=Ealpha(itype)/vol2    
       enddo

       do itype=0,ntype
       vcoul_nL2(:,itype)=vcoul_nL2(:,itype)+dv_corr(itype)
       enddo 

       if(icoul.eq.1) then
cccccccccc  provide a buffer cutoff surface. Actually not really needed. 
cccccccccc  But maybe it is a good idea to keep the symmetric of the potential
       do ii=1,nr_nL2
       jj=ii+(inode-1)*nr_nL2
       i1=(jj-1)/(n2L2*n3L2)+1
       j1=(jj-1-(i1-1)*n2L2*n3L2)/n3L2+1
       k1=jj-(i1-1)*n2L2*n3L2-(j1-1)*n3L2
                                                                                             
      amask1=1.d0
      if(mod(n1L2,2).eq.0.and.i1.eq.n1L2/2+1) amask1=0.d0
      amask2=1.d0
      if(mod(n2L2,2).eq.0.and.j1.eq.n2L2/2+1) amask2=0.d0
      amask3=1.d0
      if(mod(n3L2,2).eq.0.and.k1.eq.n3L2/2+1) amask3=0.d0
                                                                                             
      amask=amask1*amask2*amask3
                                                                                             
ccccccc not really necessary for the symmetry

cc      do itype=0,ntype
cc      vcoul_nL2(ii,itype)=vcoul_nL2(ii,itype)*amask
cc      enddo
      enddo
       
      goto 4000  ! don't do the following, which is for slab calculation
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
****  Now, we need to add back the sum of 2D image potentials !
****  This again, it done by a Ewald like summation. 
****  First, get the two 2D lattice vectors: 

      if(icoul.eq.11) then
      AL2D(:,1)=AL2(:,2)
      AL2D(:,2)=AL2(:,3)
      endif
      if(icoul.eq.12) then
      AL2D(:,1)=AL2(:,1)
      AL2D(:,2)=AL2(:,3)
      endif
      if(icoul.eq.13) then
      AL2D(:,1)=AL2(:,1)
      AL2D(:,2)=AL2(:,2)
      endif

cccccccccccccccccccccccccccccccccccccc
      d1q_2D=1.d0/dsqrt(AL2D(1,1)**2+AL2D(2,1)**2+AL2D(3,1)**2)
      d2q_2D=1.d0/dsqrt(AL2D(1,2)**2+AL2D(2,2)**2+AL2D(3,2)**2)

      dcosth=AL2D(1,1)*AL2D(1,2)+AL2D(2,1)*AL2D(2,2)+
     &         AL2D(3,1)*AL2D(3,2)
      dcosth=dcosth*d1q_2D*d2q_2D
      ALI2D(:,1)=AL2D(:,1)-AL2D(:,2)*d2q_2D/d1q_2D*dcosth
      ALI2D(:,2)=AL2D(:,2)-AL2D(:,1)*d1q_2D/d2q_2D*dcosth

      S_plane=dsqrt(ALI2D(1,1)**2+ALI2D(2,1)**2+ALI2D(3,1)**2)/d2q_2D

      sum1=ALI2D(1,1)*AL2D(1,1)+ALI2D(2,1)*AL2D(2,1)+
     &         ALI2D(3,1)*AL2D(3,1)
      sum2=ALI2D(1,2)*AL2D(1,2)+ALI2D(2,2)*AL2D(2,2)+
     &         ALI2D(3,2)*AL2D(3,2)
      ALI2D(:,1)=ALI2D(:,1)/sum1    ! ALI2D is the reciprocal of AL2D, on the same plane
      ALI2D(:,2)=ALI2D(:,2)/sum2
      
      d1_2D=1.d0/dsqrt(ALI2D(1,1)**2+ALI2D(2,1)**2+ALI2D(3,1)**2)
      d2_2D=1.d0/dsqrt(ALI2D(1,2)**2+ALI2D(2,2)**2+ALI2D(3,2)**2)

      nc1_2D=dsqrt(cut_rr)/d1_2D+0.25   ! need to look at more carefully
      nc2_2D=dsqrt(cut_rr)/d2_2D+0.25   ! need to look at more carefully

      nq1_2D=dsqrt(cut_gkk)/d1q_2D  
      nq2_2D=dsqrt(cut_gkk)/d2q_2D 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(inode_tot.eq.1) then
      write(6,*) "d1,2_2D",d1_2D,d2_2D
      write(6,*) "nc1,2_2D",nc1_2D,nc2_2D
      write(6,*) "nq1,2_2D",nq1_2D,nq2_2D
      write(6,*) "AL2D",AL2D(1,1),AL2D(2,1),AL2D(3,1)
      write(6,*) "AL2D",AL2D(1,2),AL2D(2,2),AL2D(3,2)
      write(6,*) "ALI2D",ALI2D(1,1),ALI2D(2,1),ALI2D(3,1)
      write(6,*) "ALI2D",ALI2D(1,2),ALI2D(2,2),ALI2D(3,2)
      endif
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 2009 ii=1,nr_nL2 
      jj=ii+(inode-1)*nr_nL2
      i1=(jj-1)/(n2L2*n3L2)+1
      j1=(jj-1-(i1-1)*n2L2*n3L2)/n3L2+1
      k1=jj-(i1-1)*n2L2*n3L2-(j1-1)*n3L2

      if(i1.le.n1L2/2+1) then
      idx1=i1-1
      else
      idx1=i1-n1L2-1
      endif

      if(j1.le.n2L2/2+1) then
      idx2=j1-1
      else
      idx2=j1-n2L2-1
      endif

      if(k1.le.n3L2/2+1) then
      idx3=k1-1
      else
      idx3=k1-n3L2-1
      endif

      x0=AL2(1,1)*idx1/n1L2+AL2(1,2)*idx2/n2L2+
     &     AL2(1,3)*idx3/n3L2
      y0=AL2(2,1)*idx1/n1L2+AL2(2,2)*idx2/n2L2+
     &     AL2(2,3)*idx3/n3L2
      z0=AL2(3,1)*idx1/n1L2+AL2(3,2)*idx2/n2L2+
     &     AL2(3,3)*idx3/n3L2

      x11=ALI2D(1,1)*x0+ALI2D(2,1)*y0+ALI2D(3,1)*z0
      x22=ALI2D(1,2)*x0+ALI2D(2,2)*y0+ALI2D(3,2)*z0

      dz_x0=x0-AL2D(1,1)*x11-AL2D(1,2)*x22
      dz_y0=y0-AL2D(2,1)*x11-AL2D(2,2)*x22
      dz_z0=z0-AL2D(3,1)*x11-AL2D(3,2)*x22

      dz=dsqrt(dz_x0**2+dz_y0**2+dz_z0**2)

      sr=0.d0
      do 209 i=-nc1_2D,nc1_2D     ! this is the real space summation, rho(r)=gauss(r-R)-delta(r-R)
      do 209 j=-nc2_2D,nc2_2D
      x=x0+AL2D(1,1)*i+AL2D(1,2)*j
      y=y0+AL2D(2,1)*i+AL2D(2,2)*j
      z=z0+AL2D(3,1)*i+AL2D(3,2)*j
      rr=x**2+y**2+z**2
      if((i.ne.0.or.j.ne.0).and.rr.gt.cut_rr) goto 209

      r=dsqrt(rr)

      if(i.eq.0.and.j.eq.0) then  
***  special treatment for i,j=0, the point itself,not the images
***  for this point, only add the potential from gauss(r) [(erfc(beta*r)-1)/r], not the delta(r-R)
       if(r.lt.1.D-6) r=1.D-6
       if(rr.gt.cut_rr) then
       sr=sr-1.d0/r
       else
       sr=sr+(erfc(beta*r)-1.d0)/r
       endif
      else
**** for all the other image points, add the potential from gauss(r-R)-delta(r-R) [erfc(beta*r)/r]
      sr=sr+erfc(beta*r)/r
      endif

209   continue

cccccc not so sure about the zero point
ccc    sr=sr-dv    ! the important constant, to correct the v(q=0)=0 overall shifting

      sq=0.d0
      do 309 i=0,nq1_2D     ! this is q space summation, subtract out the lattice gauss(r-R) potential
      do 309 j=-nq2_2D,nq2_2D
cccc the summation over nq1,nq2,nq3 is done for each grid points, so no FFT etc. (It actually 
cccc can be done by FFT !) 

      if(i.eq.0.and.j.lt.0) goto 309       ! only half space need to be summed
      if(i.eq.0.and.j.eq.0) goto 309

      gkx=ALI2D(1,1)*i+ALI2D(1,2)*j
      gky=ALI2D(2,1)*i+ALI2D(2,2)*j
      gkz=ALI2D(3,1)*i+ALI2D(3,2)*j
     
      gkk=gkx**2+gky**2+gkz**2        ! gkk,gkx,gky,gkz missing 2*pi

      if(gkk.gt.cut_gkk) goto 309       ! this same criterion is okay for 2D q sum

      gk=dsqrt(gkk)*2*pi

      ph=(x0*gkx+y0*gky+z0*gkz)*2*pi

cccccc the Fortran erfc(x) is good for -infity < x < 10.d0
cccccc! make sure dexp(dz*gk) not overflow, erfc(y) not be zero
      y=dz*beta+gk/beta/2
      if(y.lt.9.5) then    
      sq=sq+dcos(ph)/gk*(
     &  dexp(dz*gk)*erfc(dz*beta+gk/beta/2)+dexp(-dz*gk)*
     &  erfc(-dz*beta+gk/beta/2))
      else
      f22=1.d0     ! f22=exp(y**2)*erfc(y)
      term_iter=1.d0
      do m_iter=1,10    ! serial expansion for exp(y**2)*erfc(y). 
      term_iter=-term_iter*(2*m_iter-1)/(2*y**2)
      f22=f22+term_iter
      enddo
      f22=f22/y/dsqrt(pi)
      sq=sq+dcos(ph)/gk*(f22*dexp(-(dz*beta)**2-(gk/beta/2)**2)
     &  +dexp(-dz*gk)*erfc(-dz*beta+gk/beta/2))
      endif

309   continue

      sq=sq*pi/S_plane*2    ! fact of 2 because  i=0,nq1_2D, only sum half of the sphere
ccccc THIS GAUSSIAN PART LOOKS TOO SMALL !!!   IS THIS REALLY CORRECT ! fact of 10000 !

      sq=sq-2*pi/S_plane/beta*(dexp(-(dz*beta)**2)/dsqrt(pi)+
     &   dz*beta*erf(dz*beta))       ! this formula is okay numerically 

***** erf(-x)=-erf(x), checked 

****** The constant is not added yet. 

cccccccccccc opposite sign as before
      vcoul_nL2(ii,0)=vcoul_nL2(ii,0)+sr+sq  

      do itype=1,ntype
      vcoul_nL2(ii,itype)=vcoul_nL2(ii,itype)-
     &                    zatom(itype)*(sr+sq)  
      enddo

2009  continue


cccccccccccccccccccccccccccccccccccccccccccccc
cccccccccc  provide a buffer cutoff surface. Actually not really needed. 
cccccccccc  But maybe it is a good idea to keep the symmetric of the potential
      do ii=1,nr_nL2
       jj=ii+(inode-1)*nr_nL2
       i1=(jj-1)/(n2L2*n3L2)+1
       j1=(jj-1-(i1-1)*n2L2*n3L2)/n3L2+1
       k1=jj-(i1-1)*n2L2*n3L2-(j1-1)*n3L2
                                                                                             
      amask=1.d0
      if(icoul.eq.11.and.mod(n1L2,2).eq.0.and.i1.eq.n1L2/2+1) 
     &    amask=0.d0
      if(icoul.eq.12.and.mod(n2L2,2).eq.0.and.j1.eq.n2L2/2+1) 
     &    amask=0.d0
      if(icoul.eq.13.and.mod(n3L2,2).eq.0.and.k1.eq.n3L2/2+1) 
     &    amask=0.d0
                                                                                             
ccccccc  not really necessary for the symmetry
c      do itype=0,ntype
c      vcoul_nL2(ii,itype)=vcoul_nL2(ii,itype)*amask
c      enddo
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  vcoul_nL2 is like 1/r
ccccccccccccccccccccccccccccccccccccccccccccccccccc

4000  continue     ! jumping point for icoul=1


      do itype=0,ntype
      workr_nL2(:)= vcoul_nL2(:,itype)
      call d3fft_real2L2(vcoul_nL2(1,itype),workr_nL2,1,0)
      enddo

      vcoul_nL2(:,0)=vcoul_nL2(:,0)*vol2

ccccccccc  Now, vcoul_nL2 is in q space. 

      deallocate(workr_nL2)

      return
      end



  
