      subroutine getewald2D(fatom,xatom,AL,
     &    ityatom,ewald)
***************************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************

ccccccccccccccccccccccccccccccccccccccc
ccccc  final result is fatom=+dE/dR

ccccc This program needs some fine tunning.

       use fft_data
       use load_data
       use data

      implicit double precision (a-h,o-z)
      include 'param.escan_real'
      include "mpif.h"

      real*8 AL(3,3)   ! ALI is passed in from param.escan_real
      real*8 AL2D(3,2),ALI2D(3,2)
      real*8 xatom(3,matom),fatom(3,matom),fatomtmp(3,matom)
      real*8 zatom(mtype)
      integer ityatom(matom)
      integer, allocatable, dimension(:) :: ncount_tmp
      integer icoul
      real*8 xcoul(3)
                                                                                                    
      common /comzatom/zatom
      common /comcoul/icoul,xcoul

cccccccccccccccccccccccccccccccccccccccccccccc
ccccc assigne atoms to each node
      allocate(ncount_tmp(nnodes))
      ncount_tmp=0
      jnode=0
      do ia=1,natom
      jnode=jnode+1
      if(jnode.gt.nnodes) jnode=jnode-nnodes
      ncount_tmp(jnode)=ncount_tmp(jnode)+1
      enddo
      natom_st=0
      do jnode=1,inode-1
      natom_st=natom_st+ncount_tmp(jnode)
      enddo
      natom_st=natom_st+1
      natom_fn=natom_st+ncount_tmp(inode)-1
      deallocate(ncount_tmp)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      pi=4*datan(1.d0)
*******************************************
      if(icoul.eq.11) then
      AL2D(:,1)=AL(:,2)
      AL2D(:,2)=AL(:,3)
      endif
      if(icoul.eq.12) then
      AL2D(:,1)=AL(:,1)
      AL2D(:,2)=AL(:,3)
      endif
      if(icoul.eq.13) then
      AL2D(:,1)=AL(:,1)
      AL2D(:,2)=AL(:,2)
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

ccccccccccccccccccccccccccccccccccccccccc
      beta=1.2*pi**0.5*
     &     (d1q_2D*d2q_2D/(d1_2D*d2_2D))**0.25d0

      cut=5.0d0

      cut_gkk=4*(cut*beta)**2/(2*pi)**2
      cut_rr=(cut/beta)**2
      fac_gkk=(2*pi)**2/(4*beta**2)
       
      nc1_2D=dsqrt(cut_rr)/d1_2D+0.25+1   ! need to look at more carefully
      nc2_2D=dsqrt(cut_rr)/d2_2D+0.25+1   ! need to look at more carefully
      nq1_2D=dsqrt(cut_gkk)/d1q_2D+1
      nq2_2D=dsqrt(cut_gkk)/d2q_2D+1

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ewald=0.d0
      do 1000 ia1=1,natom
      fatom(1,ia1)=0.d0
      fatom(2,ia1)=0.d0   
      fatom(3,ia1)=0.d0  

      x1=xatom(1,ia1)
      if(x1.ge.1.d0) x1=x1-1.d0
      if(x1.lt.0.d0) x1=x1+1.d0
      if(x1.gt.xcoul(1)) x1=x1-1.d0
      y1=xatom(2,ia1)
      if(y1.ge.1.d0) y1=y1-1.d0
      if(y1.lt.0.d0) y1=y1+1.d0
      if(y1.gt.xcoul(2)) y1=y1-1.d0
      z1=xatom(3,ia1)
      if(z1.ge.1.d0) z1=z1-1.d0
      if(z1.lt.0.d0) z1=z1+1.d0
      if(z1.gt.xcoul(3)) z1=z1-1.d0


      do 1000 ia2=natom_st,natom_fn

      x2=xatom(1,ia2)
      if(x2.ge.1.d0) x2=x2-1.d0
      if(x2.lt.0.d0) x2=x2+1.d0
      if(x2.gt.xcoul(1)) x2=x2-1.d0
      y2=xatom(2,ia2)
      if(y2.ge.1.d0) y2=y2-1.d0
      if(y2.lt.0.d0) y2=y2+1.d0
      if(y2.gt.xcoul(2)) y2=y2-1.d0
      z2=xatom(3,ia2)
      if(z2.ge.1.d0) z2=z2-1.d0
      if(z2.lt.0.d0) z2=z2+1.d0
      if(z2.gt.xcoul(3)) z2=z2-1.d0

      dx1=x1-x2    ! the definition of dx1 is different from getewald, minus sign
      dx2=y1-y2
      dx3=z1-z2

      x0=AL(1,1)*dx1+AL(1,2)*dx2+AL(1,3)*dx3
      y0=AL(2,1)*dx1+AL(2,2)*dx2+AL(2,3)*dx3
      z0=AL(3,1)*dx1+AL(3,2)*dx2+AL(3,3)*dx3

      x11=ALI2D(1,1)*x0+ALI2D(2,1)*y0+ALI2D(3,1)*z0
      x22=ALI2D(1,2)*x0+ALI2D(2,2)*y0+ALI2D(3,2)*z0
       
      dz_x0=x0-AL2D(1,1)*x11-AL2D(1,2)*x22
      dz_y0=y0-AL2D(2,1)*x11-AL2D(2,2)*x22
      dz_z0=z0-AL2D(3,1)*x11-AL2D(3,2)*x22
ccccccc dz_x0,dz_y0,dz_z0 are the perpendicular z direction of (x0,y0,z0)

      dz=dsqrt(dz_x0**2+dz_y0**2+dz_z0**2)
      dztmp=dz
      if(dztmp.lt.1.D-8) dztmp=1.D-8

      ch1=zatom(ityatom(ia1))
      ch2=zatom(ityatom(ia2))

      sr=0.d0
      srx=0.d0
      sry=0.d0
      srz=0.d0
      do 200 i=-nc1_2D,nc1_2D
      do 200 j=-nc2_2D,nc2_2D
      x=x0+AL2D(1,1)*i+AL2D(1,2)*j
      y=y0+AL2D(2,1)*i+AL2D(2,2)*j
      z=z0+AL2D(3,1)*i+AL2D(3,2)*j
      

      rr=x**2+y**2+z**2

      if(rr.gt.cut_rr) goto 200
      r=dsqrt(rr)

ccccc do not calc. force for ia1.eq.ia2, they cancels out anyway for ia1.eq.ia2 image atoms
      if(ia2.ne.ia1) then    ! r cannot be zero
      derfrdr=(erfc(beta*(r+2.d-5))/(r+2.d-5)-
     &  erfc(beta*r)/r)/2.d-5     !  The analytical expression for this is available

      srx=srx-derfrdr*x/r
      sry=sry-derfrdr*y/r
      srz=srz-derfrdr*z/r
      endif

      if(ia1.eq.ia2.and.i.eq.0.and.
     &      j.eq.0) goto 200
      sr=sr+erfc(beta*r)/r

200   continue

cccc  this term is for the r=0 potential due to the onsite Gauss. density
      if(ia1.eq.ia2) sr=sr-2*beta/dsqrt(pi)  

      ewald=ewald+sr*ch1*ch2*0.5d0

      fatom(1,ia1)=fatom(1,ia1)+srx*ch1*ch2
      fatom(2,ia1)=fatom(2,ia1)+sry*ch1*ch2
      fatom(3,ia1)=fatom(3,ia1)+srz*ch1*ch2


      sq=0.d0
      sqx=0.d0
      sqy=0.d0
      sqz=0.d0
      do 309 i=0,nq1_2D    
      do 309 j=-nq2_2D,nq2_2D
       
      if(i.eq.0.and.j.lt.0) goto 309       ! only half space need to be summed
      if(i.eq.0.and.j.eq.0) goto 309       ! don't do the G=0 point
       
      gkx=ALI2D(1,1)*i+ALI2D(1,2)*j
      gky=ALI2D(2,1)*i+ALI2D(2,2)*j
      gkz=ALI2D(3,1)*i+ALI2D(3,2)*j
       
      gkk=gkx**2+gky**2+gkz**2        ! gkk,gkx,gky,gkz missing 2*pi
       
      if(gkk.gt.cut_gkk) goto 309       ! this same criterion is okay for 2D q sum

      gk=dsqrt(gkk)*2*pi
      gkx=gkx*2*pi
      gky=gky*2*pi
      gkz=gkz*2*pi
       
       
      ph=(x0*gkx+y0*gky+z0*gkz)

cccccc the Fortran erfc(x) is good for -infity < x < 10.d0
cccccc! make sure dexp(dz*gk) not overflow, erfc(y) not be zero
      y=dz*beta+gk/beta/2
      if(y.lt.9.5) then
      f33=dexp(dz*gk)*erfc(dz*beta+gk/beta/2)
      else
      f22=1.d0     ! f22=exp(y**2)*erfc(y)
      term_iter=1.d0
      do m_iter=1,10    ! serial expansion for exp(y**2)*erfc(y).
      term_iter=-term_iter*(2*m_iter-1)/(2*y**2)
      f22=f22+term_iter
      enddo
      f22=f22/y/dsqrt(pi)
      f33=f22*dexp(-(dz*beta)**2-(gk/beta/2)**2)
      endif
      f44=dexp(-dz*gk)*erfc(-dz*beta+gk/beta/2)
      ff=f33+f44
      sq=sq+dcos(ph)/gk*ff
      sqx=sqx+dsin(ph)/gk*ff*gkx
      sqy=sqy+dsin(ph)/gk*ff*gky
      sqz=sqz+dsin(ph)/gk*ff*gkz

      ff2=dcos(ph)*(f33-f44)
      sqx=sqx-ff2*dz_x0/dztmp      ! these are the dz direction force
      sqy=sqy-ff2*dz_y0/dztmp      ! the paper has a minus sign error. 
      sqz=sqz-ff2*dz_z0/dztmp

309   continue

ccccccccccc  the dz direction force is still not right, it depends on beta !!!

      ff2=2*0.5d0*erf(dz*beta)     ! 0.5, since later, it will be multiplied by one extra 2
      sqx=sqx+ff2*dz_x0/dztmp      ! these are the dz direction force
      sqy=sqy+ff2*dz_y0/dztmp      
      sqz=sqz+ff2*dz_z0/dztmp

      fact=ch1*ch2*pi*2/S_plane    ! fact of 2 because  i=0,nq1_2D, only sum half of the spheree

      sq=sq-1.d0/beta*(dexp(-(dz*beta)**2)/dsqrt(pi)+
     &   dz*beta*erf(dz*beta))       ! this is the z direction term for total energy

      ewald=ewald+sq*fact*0.5d0     ! 0.5d0 due to the double counting
      fatom(1,ia1)=fatom(1,ia1)+sqx*fact
      fatom(2,ia1)=fatom(2,ia1)+sqy*fact
      fatom(3,ia1)=fatom(3,ia1)+sqz*fact

1000  continue

      call mpi_allreduce(ewald,ewaldtmp,1,MPI_REAL8,
     & MPI_SUM,MPI_COMM_K,ierr)
      ewald = ewaldtmp
      call mpi_allreduce(fatom,fatomtmp,3*natom,MPI_REAL8,
     & MPI_SUM,MPI_COMM_K,ierr)
      fatom = fatomtmp

ccccc above formula is fatom=- dE/dR
ccccc we need to change it to fatom=dE/dR

      fatom=-fatom  

ccccccccc  THERE MIGHT BE A CONSTANT WHICH IS MISSING IN EWALD !!

      return
      end



  
