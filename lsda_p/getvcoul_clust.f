      subroutine getvcoul_clust(ntype,Ealpha)
***************************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************
****** This program calculate the modified v(r), either without periodicity

ccccc This program needs some fine tunning.

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

      character*20 file_tmp


      common /comVrho/qi2,vq,rhoq,vqT,rhocq
      common /comzatom/zatom


cccccccccccccccccccccccccccccccccccccccccccccccccc
       allocate(workr_nL2(mr_nL2))
       vins2=1.d0/vol2
       ng2_nL2=ngtotnod2L2(inode)
       pi=4*datan(1.d0)


       vcoul_nL2=0.d0
       do 100  itype=0,ntype

       do 10 i=1,ng2_nL2

      q=dsqrt(gkx2_nL2(i)**2+gky2_nL2(i)**2+gkz2_nL2(i)**2)


ccccccc vq(iq),qi2(iq) not defined for 0.5*q**2 > Ecut2L

      if(0.5*q**2.ge.Ecut2L) goto 10  

      iq=1+q*(mnq-1.d0)/qi2(mnq)

      x=(q-qi2(iq))/(qi2(iq+1)-qi2(iq))    ! assuming equal distance grid
 
      f1=1-x-0.5d0*x*(1-x)
      f2=x+x*(1-x)
      f3=-0.5d0*x*(1-x)       ! using quadratic interpolation

      amask=1.d0
      if(0.5d0*q**2.gt.Ecut2L*0.9d0) then
      x=(0.5d0*q**2-Ecut2L*0.9d0)/(0.1d0*Ecut2L)
      amask=(dcos(x*pi)+1.d0)/2.d0
      endif

      if(itype.ge.1) then
      y=vq(iq,itype)*qi2(iq)**2*f1+vq(iq+1,itype)*qi2(iq+1)**2*f2+
     &   vq(iq+2,itype)*qi2(iq+2)**2*f3
      else
      y=4*pi
      endif

 

      if(q.lt.1.D-6) then
      y=0.d0
      else
      y=y/q**2
      endif

      y=y*amask

      i2=i*2


      vcoul_nL2(i2-1,itype)=vcoul_nL2(i2-1,itype)+y*vins2
      vcoul_nL2(i2,itype)=vcoul_nL2(i2,itype)+0.d0


 10   continue

      call d3fft_real2L2(vcoul_nL2(1,itype),workr_nL2,-1,0)
      vcoul_nL2(:,itype) = workr_nL2(:)

100   continue

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


      if(inode.eq.1) then
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
      if(mod(n1L2,2).eq.0.and.i1.eq.n1L2/2+1) goto 2000

      if(j1.le.n2L2/2+1) then
      idx2=j1-1
      else
      idx2=j1-n2L2-1
      endif
      if(mod(n2L2,2).eq.0.and.j1.eq.n2L2/2+1) goto 2000

      if(k1.le.n3L2/2+1) then
      idx3=k1-1
      else
      idx3=k1-n3L2-1
      endif
      if(mod(n3L2,2).eq.0.and.k1.eq.n3L2/2+1) goto 2000

      x0=AL2(1,1)*idx1/n1L2+AL2(1,2)*idx2/n2L2+
     &     AL2(1,3)*idx3/n3L2
      y0=AL2(2,1)*idx1/n1L2+AL2(2,2)*idx2/n2L2+
     &     AL2(2,3)*idx3/n3L2
      z0=AL2(3,1)*idx1/n1L2+AL2(3,2)*idx2/n2L2+
     &     AL2(3,3)*idx3/n3L2

      sr=0.d0
      do 200 i=-nc1,nc1
      do 200 j=-nc2,nc2
      do 200 k=-nc3,nc3
      x=x0+AL2(1,1)*i+AL2(1,2)*j+AL2(1,3)*k
      y=y0+AL2(2,1)*i+AL2(2,2)*j+AL2(2,3)*k
      z=z0+AL2(3,1)*i+AL2(3,2)*j+AL2(3,3)*k
      rr=x**2+y**2+z**2
      if((i.ne.0.or.j.ne.0.or.k.ne.0).and.rr.gt.cut_rr) goto 200

      r=dsqrt(rr)

      if(i.eq.0.and.j.eq.0.and.k.eq.0) then
       if(r.lt.1.D-6) r=1.D-6
       if(rr.gt.cut_rr) then
       sr=sr-1.d0/r
       else
       sr=sr+(erfc(beta*r)-1.d0)/r
       endif
      else
      sr=sr+erfc(beta*r)/r
      endif

200   continue

      sr=sr-dv

      sq=0.d0
      do 300 i=0,nq1
      do 300 j=-nq2,nq2
      do 300 k=-nq3,nq3

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

ccccccccccccccccccccccccccccccccccccccccccccccccccc

       do 3000 ii=1,nr_nL2 
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

      do itype=0,ntype
      vcoul_nL2(ii,itype)=(vcoul_nL2(ii,itype)+
     &        dv_corr(itype))*amask
      enddo

3000   continue




      do itype=0,ntype
      workr_nL2(:)= vcoul_nL2(:,itype)
      call d3fft_real2L2(vcoul_nL2(1,itype),workr_nL2,1,0)
      enddo


      vcoul_nL2(:,0)=vcoul_nL2(:,0)*vol2

      deallocate(workr_nL2)

      return
      end



  
