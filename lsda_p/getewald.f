      subroutine getewald(fatom,xatom,AL,
     &    ityatom,ewald)
***************************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
******************************************

ccccc This program needs some fine tunning.

       use fft_data
       use load_data
       use data

      implicit double precision (a-h,o-z)
      include 'param.escan_real'
      include "mpif.h"

      real*8 AL(3,3)
      real*8 xatom(3,matom),fatom(3,matom)
      real*8 zatom(mtype)
      integer ityatom(matom)
      integer, allocatable, dimension(:) :: ncount_tmp

      common /comzatom/zatom

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


      pi=4*datan(1.d0)
      d1=dsqrt(AL(1,1)**2+AL(2,1)**2+AL(3,1)**2)
      d2=dsqrt(AL(1,2)**2+AL(2,2)**2+AL(3,2)**2)
      d3=dsqrt(AL(1,3)**2+AL(2,3)**2+AL(3,3)**2)

ccccc This beta balances the real and G space calculations.
ccccc The computational time in both real and G space are O(N^2).
ccccc In both real and G space, the calc. propto ~ 1000*natom**2

      dave=(d1*d2*d3)**0.333333d0

      if(dave.le.16.d0) then 
      beta=3.0d0/dave
      else
      beta=5.d0/dave
      endif


      nc1=16.d0/beta/d1+1
      nc2=16.d0/beta/d2+1
      nc3=16.d0/beta/d3+1


      ewald=0.d0
      do 1000 ia1=1,natom
      fatom(1,ia1)=0.d0
      fatom(2,ia1)=0.d0   
      fatom(3,ia1)=0.d0  

      do 1000 ia2=natom_st,natom_fn
      x1=xatom(1,ia2)-xatom(1,ia1)
      y1=xatom(2,ia2)-xatom(2,ia1)
      z1=xatom(3,ia2)-xatom(3,ia1)
      x0=AL(1,1)*x1+AL(1,2)*y1+AL(1,3)*z1
      y0=AL(2,1)*x1+AL(2,2)*y1+AL(2,3)*z1
      z0=AL(3,1)*x1+AL(3,2)*y1+AL(3,3)*z1

      ch1=zatom(ityatom(ia1))
      ch2=zatom(ityatom(ia2))

      sr=0.d0
      srx=0.d0
      sry=0.d0
      srz=0.d0
      do 200 i=-nc1,nc1
      do 200 j=-nc2,nc2
      do 200 k=-nc3,nc3
      x=x0+AL(1,1)*i+AL(1,2)*j+AL(1,3)*k
      y=y0+AL(2,1)*i+AL(2,2)*j+AL(2,3)*k
      z=z0+AL(3,1)*i+AL(3,2)*j+AL(3,3)*k
      r=dsqrt(x**2+y**2+z**2)

      if(beta*r.gt.8.d0) goto 200

      if(ia2.ne.ia1) then
      derfrdr=(erfc(beta*(r+2.d-5))/(r+2.d-5)-
     &  erfc(beta*r)/r)/2.d-5

      srx=srx-derfrdr*x/r
      sry=sry-derfrdr*y/r
      srz=srz-derfrdr*z/r
      endif

      if(ia1.eq.ia2.and.i.eq.0.and.
     &      j.eq.0.and.k.eq.0) goto 200
      sr=sr+erfc(beta*r)/r

200   continue

      if(ia1.eq.ia2) sr=sr-2*beta/dsqrt(pi)

ccccc factor 2 for the sum 1,ng2, is only half the sphere

      ewald=ewald+sr*ch1*ch2*0.5d0

      fatom(1,ia1)=fatom(1,ia1)+srx*ch1*ch2
      fatom(2,ia1)=fatom(2,ia1)+sry*ch1*ch2
      fatom(3,ia1)=fatom(3,ia1)+srz*ch1*ch2
1000  continue

      ng2_n=ngtotnod2(inode)


      sq=0.d0
      do 300 i=1,ng2_n
      gkk=gkx2_n(i)**2+gky2_n(i)**2+gkz2_n(i)**2
      yy=gkk/(4*beta**2)

      if(gkk.le.1.D-10.or.yy.gt.50.d0) goto 300

      ff1=dexp(-yy)
      ff2=ff1*4*pi/vol*2.d0/gkk


      s_cos=0.d0
      do 2002 ia1=1,natom
      s_sin=0.d0
      do 2000 ia2=1,natom

      x1=xatom(1,ia2)-xatom(1,ia1)
      y1=xatom(2,ia2)-xatom(2,ia1)
      z1=xatom(3,ia2)-xatom(3,ia1)

      x0=AL(1,1)*x1+AL(1,2)*y1+AL(1,3)*z1
      y0=AL(2,1)*x1+AL(2,2)*y1+AL(2,3)*z1
      z0=AL(3,1)*x1+AL(3,2)*y1+AL(3,3)*z1

      ph=x0*gkx2_n(i)+y0*gky2_n(i)+z0*gkz2_n(i)

      ch=zatom(ityatom(ia1))*zatom(ityatom(ia2))

      s_cos=s_cos+ch*dcos(ph)
      s_sin=s_sin+ch*dsin(ph)
2000  continue

      s_sin=s_sin*ff2
      fatom(1,ia1)=fatom(1,ia1)+s_sin*gkx2_n(i)
      fatom(2,ia1)=fatom(2,ia1)+s_sin*gky2_n(i)
      fatom(3,ia1)=fatom(3,ia1)+s_sin*gkz2_n(i)
2002  continue

      sq=sq+ff1*s_cos/gkk

300   continue
ccccc factor 2 for the sum 1,ng2, is only half the sphere

      sq=sq*4*pi/vol*2.d0

      ch=0.d0
      do ia1=1,natom
      do ia2=1,natom
      ch=ch+zatom(ityatom(ia1))*zatom(ityatom(ia2))
      enddo
      enddo



      call mpi_allreduce(sq,sq,1,MPI_REAL8,
     & MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(ewald,ewald,1,MPI_REAL8,
     & MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(fatom,fatom,3*natom,MPI_REAL8,
     & MPI_SUM,MPI_COMM_WORLD,ierr)


      sq=sq-ch*pi/beta**2/vol
      ewald=ewald+sq*0.5d0

      return
      end



  
