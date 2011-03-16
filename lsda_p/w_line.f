      subroutine w_line(ilocal,ntype,vwr_atom,Ealpha)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
******************************************


      use data

      implicit double precision(a-h,o-z)


      include "mpif.h"
      include "param.escan_real"
      real*8 r(2000),ws(2000),wp(2000),wd(2000)
      real*8 vloc(2000),vs(2000),vp(2000),vd(2000)
      real*8 vlocT(2000),rhoc(2000)
      real*8 vw(2000),vw2(2000),vw4(2000)

ccccccc  assume ws,wp and vs,vp,vd has the same 
ccccccc  r(i)
ccccccc  wq2() is for the old fashion G space implementation
ccccccc  wq4() is for TBinit_ug.f 
      real*8 qi2(mnq),vq(mnq,mtype),rhoq(mnq,mtype)
      real*8 rhocq(mnq,mtype)
      real*8 vqT(mnq,mtype)

cccccccc wq is for real space implementation with 1/rmask
cccccccc wq2 is for q space implementation without 1/rmask

      real*8 qi(mnq),wq(mnq,3,mtype),wq2(mnq,3,mtype)
      real*8 wq4(mnq,3,mtype)
      real*8 ri(201),amr(201)
      real*8 qi3(mnq),qi4(mnq)

      integer  iiatom(mtype),numref(matom)
      real*8  zatom(mtype),Ealpha(mtype)
      real*8  occ_s(mtype),occ_p(mtype),occ_d(mtype),occ_t(mtype)
      integer isNL(3,mtype),icore(mtype)
      integer is_ref(mtype),ip_ref(mtype),id_ref(mtype)
      integer is_TB(mtype),ip_TB(mtype),id_TB(mtype)


      character*20  vwr_atom(mtype)

      common /comVrho/qi2,vq,rhoq,vqT,rhocq
      common /comline/qi,wq,ri,amr
      common /comline2/qi3,wq2
      common /comline4/qi4,wq4
      common /comNL2/occ_t,iiatom,icore,numref
      common /comzatom/zatom
      common /comisNL/isNL
      common /comispd_ref/is_ref,ip_ref,id_ref

      if(ilocal.eq.2) then
      open(10,file='maskr',status='old',action='read',iostat=ierr)
      if(ierr.ne.0) then
      if(inode.eq.1) 
     & write(6,*) "file ***maskr*** is needed for ilocal=2,stop"
      stop
      endif
      
      rewind(10)
      do i=1,201
      read(10,*) ri(i),amr(i)
      enddo
      close(10)
      else
      do i=1,201
      ri(i)=10.d0*(i-1)/200.d0
      amr(i)=1.d0
      enddo
      endif

      do 2000 ia=1,ntype

      open(10,file=vwr_atom(ia),status='old',action='read',iostat=ierr)
      if(ierr.ne.0) then
      if(inode.eq.1)
     & write(6,*) "file ***",vwr_atom(ia),"*** does not exist, stop"
      stop
      endif
      read(10,*) nrr,ic,iiatom(ia),zatom(ia),iloc,occ_s(ia),
     &  occ_p(ia),occ_d(ia)
      read(10,*) is_ref(ia),ip_ref(ia),id_ref(ia),
     &  is_TB(ia),ip_TB(ia),id_TB(ia)

      if(iloc.eq.1) is_ref(ia)=0
      if(iloc.eq.2) ip_ref(ia)=0
      if(iloc.eq.3) id_ref(ia)=0

      icore(ia)=ic


         if(nrr.gt.2000) then
         write(6,*) "nrr > 2000, in vwr.atom, cut nrr, stop"
         stop
         endif

         do i=1,nrr
         rhoc(i)=0.d0
         enddo


      if(iloc.ne.0.and.ic.eq.0) then
      do i=1,nrr
      read(10,*) r(i),vs(i),vp(i),vd(i),ws(i),wp(i),wd(i)
******* vs,vp,vd are in Hartee
      enddo
      endif

      if(iloc.eq.0.and.ic.eq.0) then
      do i=1,nrr
      read(10,*) r(i),vs(i),vp(i),vd(i),ws(i),wp(i),wd(i),vloc(i)
      enddo
      endif

      if(iloc.ne.0.and.ic.eq.1) then
      do i=1,nrr
      read(10,*) r(i),vs(i),vp(i),vd(i),ws(i),wp(i),wd(i),rhoc(i)
      enddo
      endif

      if(iloc.eq.0.and.ic.eq.1) then
      do i=1,nrr
      read(10,*) r(i),vs(i),vp(i),vd(i),ws(i),wp(i),wd(i),
     &      vloc(i),rhoc(i)
      enddo
      endif

      close(10)

************************************************
**** take the local potential
      do i=1,nrr
      if(iloc.eq.0) vloc(i)=vloc(i)
      if(iloc.eq.1) vloc(i)=vs(i)
      if(iloc.eq.2) vloc(i)=vp(i)
      if(iloc.eq.3) vloc(i)=vd(i)
      if(iloc.eq.12) vloc(i)=(vs(i)+vp(i))/2
      if(iloc.eq.13) vloc(i)=(vs(i)+vd(i))/2
      if(iloc.eq.23) vloc(i)=(vp(i)+vd(i))/2
      enddo
************************************************
      s=0.d0
      ch=zatom(ia)
      pi=4*datan(1.d0)
      do i=2,nrr-1
      if(r(i).lt.15.d0) then
      s=s+(ch*r(i)+vloc(i)*r(i)**2)*(r(i+1)-r(i-1))/2
      endif
      enddo
      Ealpha(ia)=s*4*pi

************************************************
*** vlocT is for the use of Thomas procedure
************************************************
      do i=1,nrr
      if(r(i).lt.4.d0) then
      a=vs(i)*occ_s(ia)*ws(i)**2+vp(i)*occ_p(ia)*wp(i)**2
     &  +vd(i)*occ_d(ia)*wd(i)**2
      b=occ_s(ia)*ws(i)**2+occ_p(ia)*wp(i)**2+
     &   occ_d(ia)*wd(i)**2
      vlocT(i)=a/b
      else
      vlocT(i)=vp(i)
      endif
      enddo

************************************************
   
      qmx=dsqrt(2*Ecut2)*1.2d0

      pi=4*datan(1.d0)
***********************************************************
**** for the ionic potential and atomic charge density
**** using vs as the local potential for III-V
***********************************************************
      if(occ_d(ia).gt.0.d0.and.inode.eq.1) then
      write(6,*) "warning, occ_d > 0, wd might not be normalized"
      endif

      occ_t(ia)=occ_s(ia)+occ_p(ia)+occ_d(ia)

      z=zatom(ia)
      rst=10.d0

      do 101 iq=1,mnq
      g1=(iq-1)*qmx/(mnq-1.d0)
      if(g1.lt.1.D-3) g1=1.D-3

      s=0.d0
      s1=0.d0
      s2=0.d0
      s3=0.d0

      do i=1,nrr-1
      if(r(i).gt.rst) then
      rst1=r(i)
      goto 97
      endif
      r1=r(i)
      r2=r(i+1)
      x1=g1*r1
      x2=g1*r2

      b1=vloc(i)*r1
      b2=(vloc(i+1)*r2-b1)/(x2-x1)

      c1=(dcos(x1)-dcos(x2))*b1+(x1-x2)*b2*dcos(x2)+
     &  b2*(dsin(x2)-dsin(x1))
      s=s+c1

      rhoi=(occ_s(ia)*ws(i)**2+occ_p(ia)*wp(i)**2+
     & occ_d(ia)*wd(i)**2)/(4*pi)
      rhoi1=(occ_s(ia)*ws(i+1)**2+occ_p(ia)*wp(i+1)**2+
     & occ_d(ia)*wd(i+1)**2)/(4*pi)

      b11=rhoi*r1
      b22=(rhoi1*r2-b11)/(x2-x1)
      c11=(dcos(x1)-dcos(x2))*b11+(x1-x2)*b22*dcos(x2)+
     &  b22*(dsin(x2)-dsin(x1))
      s1=s1+c11

      b1=vlocT(i)*r1
      b2=(vlocT(i+1)*r2-b1)/(x2-x1)

      c1=(dcos(x1)-dcos(x2))*b1+(x1-x2)*b2*dcos(x2)+
     &  b2*(dsin(x2)-dsin(x1))
      s2=s2+c1

      b1=rhoc(i)*r1
      b2=(rhoc(i+1)*r2-b1)/(x2-x1)
      c1=(dcos(x1)-dcos(x2))*b1+(x1-x2)*b2*dcos(x2)+
     &  b2*(dsin(x2)-dsin(x1))
      s3=s3+c1

      enddo
97    continue

      s=s/g1**2
      s=s*4*pi
      s=s-z*4*pi*dcos(g1*rst1)/g1**2

      s1=s1/g1**2
      s1=s1*4*pi

      s2=s2/g1**2
      s2=s2*4*pi
      s2=s2-z*4*pi*dcos(g1*rst1)/g1**2

      s3=s3/g1**2
      s3=s3*4*pi

      qi2(iq)=g1
      vq(iq,ia)=s
      rhoq(iq,ia)=s1
      vqT(iq,ia)=s2
      rhocq(iq,ia)=s3
101   continue

********************************************************
*** for the nonlocal potentail Kleiman-Bylander ref. dv*psi
*** using vs as the local potential
********************************************************
*** isp=1, s state; 2, p state; 3, d state.
********************************************************
      rst=rcut/1.4
      do 1000 isp=1,3

      s=0.d0
      do i=1,nrr

      if(r(i).gt.rst) then
      amrI=0.d0
      else
      r2=r(i)/rcut
      ir=1+r2*200.d0
      f1=(ri(ir+1)-r2)/(ri(ir+1)-ri(ir))
      f2=(r2-ri(ir))/(ri(ir+1)-ri(ir))
      amrI=amr(ir)*f1+amr(ir+1)*f2
      amrI=1.d0/amrI
      endif

      if(isp.eq.1) then
      vw(i)=(vs(i)-vloc(i))*ws(i)
      vw4(i)=ws(i)
        if(r(i).lt.rst.and.i.gt.1) then
        s=s+ws(i)**2*(vs(i)-vloc(i))*(r(i+1)-r(i-1))
     &      /2*r(i)**2 
        endif
      endif

      if(isp.eq.2) then
      vw(i)=(vp(i)-vloc(i))*wp(i)
      vw4(i)=wp(i)
        if(r(i).lt.rst.and.i.gt.1) then
        s=s+wp(i)**2*(vp(i)-vloc(i))*(r(i+1)-r(i-1))
     &      /2*r(i)**2 
        endif
      endif

      if(isp.eq.3) then
       vw(i)=(vd(i)-vloc(i))*wd(i)
       vw4(i)=wd(i)
        if(r(i).lt.rst.and.i.gt.1) then
        s=s+wd(i)**2*(vd(i)-vloc(i))*(r(i+1)-r(i-1))
     &      /2*r(i)**2 
        endif
      endif

      vw2(i)=vw(i)
      vw(i)=vw(i)*amrI

      enddo

      s=4*pi*s

      if(isp.eq.iloc.or.s.eq.0.d0) then
      scale=0.d0
      else
      scale=1/dsqrt(dabs(s))
      endif


      if(s.ge.0.d0) then
      isNL(isp,ia)=1
      else
      isNL(isp,ia)=-1
      endif

      if(inode.eq.1) then
      write(6,*) "ia,iloc,isp,scale", ia,iloc,isp,scale*isNL(isp,ia)
      endif
******************************
      if(isp.eq.1) then

      do 202 iq=1,mnq
      g1=(iq-1)*qmx/(mnq-1.d0)

      if(g1.lt.1.D-3) g1=1.D-3

      qi(iq)=g1

      s=0.d0
      s1=0.d0
      s2=0.d0
      do i=2,nrr-1
      if(r(i).gt.rst) goto 96
      x=r(i)*g1
      if(x.gt.0.01d0) then
      aj0=dsin(x)/x
      else
      aj0=1.d0-x**2/6.d0+x**4/120.d0-x**6/5040.d0
      endif
      s=s+aj0*r(i)**2*vw(i)*(r(i+1)-r(i-1))/2
      s1=s1+aj0*r(i)**2*vw2(i)*(r(i+1)-r(i-1))/2
      s2=s2+aj0*r(i)**2*vw4(i)*(r(i+1)-r(i-1))/2
      enddo
96    continue
      s=s*4*pi
      s=s*scale

      s1=s1*4*pi
      s1=s1*scale

      s2=s2*4*pi

      qi(iq)=g1
      qi3(iq)=g1
      qi4(iq)=g1
      wq(iq,isp,ia)=s
      wq2(iq,isp,ia)=s1
      wq4(iq,isp,ia)=s2
202   enddo
      endif


      if(isp.eq.2) then

      do 200 iq=1,mnq
      g1=(iq-1)*qmx/(mnq-1.d0)

      if(g1.lt.1.D-3) g1=1.D-3

      qi(iq)=g1

      s=0.d0
      s1=0.d0
      s2=0.d0
      do i=2,nrr-1
      if(r(i).gt.rst) goto 99
      x=r(i)*g1
      if(x.gt.0.1) then
      aj1=dsin(x)/x**2-dcos(x)/x
      else
      aj1=x/3.d0-x**3/30.d0+x**5/840.d0-x**7/45360.d0
      endif
      s=s+aj1*r(i)**2*vw(i)*(r(i+1)-r(i-1))/2
      s1=s1+aj1*r(i)**2*vw2(i)*(r(i+1)-r(i-1))/2
      s2=s2+aj1*r(i)**2*vw4(i)*(r(i+1)-r(i-1))/2
      enddo
99    continue
      s=s*4*pi
      s=s*scale

      s1=s1*4*pi
      s1=s1*scale

      s2=s2*4*pi

      qi(iq)=g1
      qi3(iq)=g1
      wq(iq,isp,ia)=s
      wq2(iq,isp,ia)=s1
      wq4(iq,isp,ia)=s2
200   enddo
      endif


      if(isp.eq.3) then

      do 201 iq=1,mnq
      g1=(iq-1)*qmx/(mnq-1.d0)

      if(g1.lt.1.D-3) g1=1.D-3

      qi(iq)=g1

      s=0.d0
      s1=0.d0
      s2=0.d0
      do i=2,nrr-1
      if(r(i).gt.rst) goto 98
      x=r(i)*g1

      if(x.gt.0.2d0) then
      aj2=(3/x**3-1/x)*dsin(x)-3*dcos(x)/x**2
      else
      aj2=x**2/15.d0-x**4/210.d0+x**6/7560.d0-x**8/498960.d0 
      endif

      s=s+aj2*r(i)**2*vw(i)*(r(i+1)-r(i-1))/2
      s1=s1+aj2*r(i)**2*vw2(i)*(r(i+1)-r(i-1))/2
      s2=s2+aj2*r(i)**2*vw4(i)*(r(i+1)-r(i-1))/2
      enddo
98    continue
      s=s*4*pi
      s=s*scale

      s1=s1*4*pi
      s1=s1*scale

      s2=s2*4*pi

      qi(iq)=g1
      qi3(iq)=g1
      wq(iq,isp,ia)=s
      wq2(iq,isp,ia)=s1
      wq4(iq,isp,ia)=s2
201   enddo
      endif

1000  continue


c     if(ia.eq.2) then
c     open(11,file="graph2")
c     rewind(11)
c     do i=1,mnq
c     write(11,300) qi(i),wq2(i,1,ia),wq2(i,2,ia),wq2(i,3,ia),
c    &     vqT(i,ia)-vq(i,ia),rhoq(i,ia)
c     enddo
c     close(11)
c     endif


300   format(6(E12.4,1x))

2000  continue

c      open(11,file="graph.vq")
c      rewind(11)
c      do i=1,mnq
c      write(11,300) qi(i),vq(i,1),vq(i,2)
c      enddo
c      close(11)

      return
      end

