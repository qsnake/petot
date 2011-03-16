      subroutine w_line_vwr(vwr_atom,zatom,iiatom,
     &  Ealpha,occ_t,isNL,Dij0,Qij,icore,is_ref,
     &  ip_ref,id_ref,qi,wq,qi2,vq,
     &  rhoq,rhocq,vqT,ri,amr,lll,nbeta)

******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************
      implicit double precision(a-h,o-z)

      include "param.escan_real"

      real*8 r(2000),ws(2000),wp(2000),wd(2000)
      real*8 vloc(2000),vs(2000),vp(2000),vd(2000)
      real*8 vlocT(2000),rhoc(2000)
      real*8 vw(2000),vw2(2000)

ccccccc  assume ws,wp and vs,vp,vd has the same r(i)
      real*8 qi2(mnq),vq(mnq),rhoq(mnq)
      real*8 rhocq(mnq)
      real*8 vqT(mnq)
      real*8 Dij0(32,32),Qij(32,32)

cccccccc wq has 1/rmask for real space implementation
cccccccc    and has no 1/rmask for q space implementation

      real*8 qi(mnq),wq(mnq,8)
      real*8 ri(201),amr(201)

      integer  iiatom

      real*8  zatom,Ealpha
      real*8  occ_s,occ_p,occ_d,occ_t
      integer isNL(3),icore
      integer is_ref,ip_ref,id_ref
      integer is_TB,ip_TB,id_TB
      integer lll(8),nbeta


      character*20  vwr_atom


      open(10,file=vwr_atom,status='old',action='read',iostat=ierr)
      if(ierr.ne.0) then
      if(inode.eq.1)
     & write(6,*) "file ***",vwr_atom,"*** does not exist, stop"
      stop
      endif
      read(10,*) nrr,ic,iiatom,zatom,iloc,occ_s,
     &  occ_p,occ_d
      read(10,*) is_ref,ip_ref,id_ref,
     &  is_TB,ip_TB,id_TB

ccccccccccc turn off this dangerous option
      is_ref=1
      ip_ref=1
      id_ref=1

      if(iloc.eq.1) is_ref=0
      if(iloc.eq.2) ip_ref=0
      if(iloc.eq.3) id_ref=0

      nbeta=0
      if(is_ref.ne.0) then
      nbeta=nbeta+1
      lll(nbeta)=0
      endif
      if(ip_ref.ne.0) then
      nbeta=nbeta+1
      lll(nbeta)=1
      endif
      if(id_ref.ne.0) then
      nbeta=nbeta+1
      lll(nbeta)=2
      endif

      icore=ic


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
      ch=zatom
      pi=4*datan(1.d0)
      do i=2,nrr-1
      if(r(i).lt.15.d0) then
      s=s+(ch*r(i)+vloc(i)*r(i)**2)*(r(i+1)-r(i-1))/2
      endif
      enddo
      Ealpha=s*4*pi

************************************************
*** vlocT is for the use of Thomas procedure
************************************************
      do i=1,nrr
      if(r(i).lt.4.d0) then
      a=vs(i)*occ_s*ws(i)**2+vp(i)*occ_p*wp(i)**2
     &  +vd(i)*occ_d*wd(i)**2
      b=occ_s*ws(i)**2+occ_p*wp(i)**2+
     &   occ_d*wd(i)**2
      vlocT(i)=a/b
      else
      vlocT(i)=vp(i)
      endif
      enddo

************************************************
   
      qmx=dsqrt(2*Ecut2L)*1.2d0

      pi=4*datan(1.d0)
***********************************************************
**** for the ionic potential and atomic charge density
**** using vs as the local potential for III-V
***********************************************************
      if(occ_d.gt.0.d0.and.inode_tot.eq.1) then
      write(6,*) "warning, occ_d > 0, wd might not be normalized"
      endif

      occ_t=occ_s+occ_p+occ_d

      z=zatom
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

      rhoi=(occ_s*ws(i)**2+occ_p*wp(i)**2+
     & occ_d*wd(i)**2)/(4*pi)
      rhoi1=(occ_s*ws(i+1)**2+occ_p*wp(i+1)**2+
     & occ_d*wd(i+1)**2)/(4*pi)

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

      s2=s2-s

      s=s/g1**2
      s=s*4*pi
      s=s-z*4*pi*dcos(g1*rst1)/g1**2

      s1=s1/g1**2
      s1=s1*4*pi

      s2=s2/g1**2
      s2=s2*4*pi

      s3=s3/g1**2
      s3=s3*4*pi

      qi2(iq)=g1
      vq(iq)=s
      rhoq(iq)=s1
      vqT(iq)=s2
      rhocq(iq)=s3
101   continue

********************************************************
*** for the nonlocal potentail Kleiman-Bylander ref. dv*psi
*** using vs as the local potential
********************************************************
*** isp=1, s state; 2, p state; 3, d state.
********************************************************
      rst=rcut/1.03
      iref=0
      Dij0=0.d0
      Qij=0.d0
      do 1000 ibeta=1,nbeta

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

      if(lll(ibeta).eq.0) then
      vw(i)=(vs(i)-vloc(i))*ws(i)
        if(r(i).lt.rst.and.i.gt.1) then
        s=s+ws(i)**2*(vs(i)-vloc(i))*(r(i+1)-r(i-1))
     &      /2*r(i)**2 
        endif
      endif

      if(lll(ibeta).eq.1) then
      vw(i)=(vp(i)-vloc(i))*wp(i)
        if(r(i).lt.rst.and.i.gt.1) then
        s=s+wp(i)**2*(vp(i)-vloc(i))*(r(i+1)-r(i-1))
     &      /2*r(i)**2 
        endif
      endif

      if(lll(ibeta).eq.2) then
       vw(i)=(vd(i)-vloc(i))*wd(i)
        if(r(i).lt.rst.and.i.gt.1) then
        s=s+wd(i)**2*(vd(i)-vloc(i))*(r(i+1)-r(i-1))
     &      /2*r(i)**2 
        endif
      endif

      vw(i)=vw(i)*amrI

      enddo

      s=4*pi*s

      if(s.eq.0.d0) then
      scale=0.d0
      else
      scale=1/dsqrt(dabs(s))
      endif

      do ll=1,2*lll(ibeta)+1
      iref=iref+1
      if(s.ge.0.d0) then
      Dij0(iref,iref)=1.d0
      else
      Dij0(iref,iref)=-1.d0
      endif
      enddo

      if(s.ge.0.d0) then
      isNL(ibeta)=1
      else
      isNL(ibeta)=-1
      endif

      if(inode_tot.eq.1) then
      write(6,*) "iloc,isp,scale", iloc,isp,scale*isNL(ibeta)
      endif
******************************
      if(lll(ibeta).eq.0) then

      do 202 iq=1,mnq
      g1=(iq-1)*qmx/(mnq-1.d0)

      if(g1.lt.1.D-3) g1=1.D-3

      qi(iq)=g1

      s=0.d0
      do i=2,nrr-1
      if(r(i).gt.rst) goto 96
      x=r(i)*g1
      if(x.gt.0.01d0) then
      aj0=dsin(x)/x
      else
      aj0=1.d0-x**2/6.d0+x**4/120.d0-x**6/5040.d0
      endif
      s=s+aj0*r(i)**2*vw(i)*(r(i+1)-r(i-1))/2
      enddo
96    continue
      s=s*4*pi
      s=s*scale

      qi(iq)=g1
      wq(iq,ibeta)=s
202   enddo
      endif


      if(lll(ibeta).eq.1) then

      do 200 iq=1,mnq
      g1=(iq-1)*qmx/(mnq-1.d0)

      if(g1.lt.1.D-3) g1=1.D-3

      qi(iq)=g1

      s=0.d0
      do i=2,nrr-1
      if(r(i).gt.rst) goto 99
      x=r(i)*g1
      if(x.gt.0.1) then
      aj1=dsin(x)/x**2-dcos(x)/x
      else
      aj1=x/3.d0-x**3/30.d0+x**5/840.d0-x**7/45360.d0
      endif
      s=s+aj1*r(i)**2*vw(i)*(r(i+1)-r(i-1))/2
      enddo
99    continue
      s=s*4*pi
      s=s*scale

      qi(iq)=g1
      wq(iq,ibeta)=s
200   enddo
      endif


      if(lll(ibeta).eq.2) then

      do 201 iq=1,mnq
      g1=(iq-1)*qmx/(mnq-1.d0)

      if(g1.lt.1.D-3) g1=1.D-3

      qi(iq)=g1

      s=0.d0
      do i=2,nrr-1
      if(r(i).gt.rst) goto 98
      x=r(i)*g1

      if(x.gt.0.2d0) then
      aj2=(3/x**3-1/x)*dsin(x)-3*dcos(x)/x**2
      else
      aj2=x**2/15.d0-x**4/210.d0+x**6/7560.d0-x**8/498960.d0 
      endif

      s=s+aj2*r(i)**2*vw(i)*(r(i+1)-r(i-1))/2
      enddo
98    continue
      s=s*4*pi
      s=s*scale


      qi(iq)=g1
      wq(iq,ibeta)=s
201   enddo
      endif

1000  continue


300   format(6(E12.4,1x))


      return
      end

