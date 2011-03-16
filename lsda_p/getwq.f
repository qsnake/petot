      subroutine getwq(AL,ntype,iatom,xatom,kpt)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'mpif.h'
      include 'param.escan_real'

      real*8 xatom(3,matom)

      real*8 AL(3,3)

      real*8 qi(mnq),wq(mnq,8,mtype)
      real*8 ri(201),amr(201)

      real*8 occ_t(mtype)
      integer iiatom(mtype),iatom(matom),icore(mtype),numref(matom),
     & ityatom(matom)
      integer is_ref(mtype),ip_ref(mtype),id_ref(mtype)
      integer lll(8,mtype),nbeta(mtype)


      complex*16 cc,cai

***************************************************
****  xatom(1),xatom(2),xatom(3) are the coord in unit of AL(3,3)
****  supercell edges
***************************************************

      common /comline/qi,wq,ri,amr
      common /comNL2/occ_t,iiatom,icore,numref,ityatom
      common /comispd_ref/is_ref,ip_ref,id_ref
      common /comlll/lll,nbeta

*******************************************************
**** generate the Kleiman-Bylander reference wavefunction
*******************************************************
      ng_n=ngtotnod(inode,kpt)

      vins=1.d0/vol
      cai=dcmplx(0.d0,1.d0)

      do 20 ia=1,natom
      iref_start2=iref_start(ia)
      

      iitype=ityatom(ia)

      do 10 i=1,ng_n

        x1=xatom(1,ia)
        y1=xatom(2,ia)
        z1=xatom(3,ia)
      
        x11=AL(1,1)*x1+AL(1,2)*y1+AL(1,3)*z1
        y11=AL(2,1)*x1+AL(2,2)*y1+AL(2,3)*z1
        z11=AL(3,1)*x1+AL(3,2)*y1+AL(3,3)*z1

        ph=gkx_n(i,kpt)*x11+gky_n(i,kpt)*y11+gkz_n(i,kpt)*z11
        cc=cdexp(dcmplx(0.d0,ph))

      q=dsqrt(gkx_n(i,kpt)**2+gky_n(i,kpt)**2+gkz_n(i,kpt)**2)

      iq=1+q*(mnq-1.d0)/qi(mnq)

      x=(q-qi(iq))/(qi(iq+1)-qi(iq))

      f1=1-x-0.5d0*x*(1-x)
      f2=x+x*(1-x)
      f3=-0.5d0*x*(1-x)

c-----------------------------------------

       if(q.lt.1.D-6) then

       kk=0
       do 301 ibeta=1,nbeta(iitype)

      yspd=wq(iq,ibeta,iitype)*f1+
     &  wq(iq+1,ibeta,iitype)*f2+wq(iq+2,ibeta,iitype)*f3

       if(lll(ibeta,iitype).eq.0) then
       wqmask(i,iref_start2+kk+1)=yspd*cc*vins
       kk=kk+1
       endif
       if(lll(ibeta,iitype).eq.1) then
       wqmask(i,iref_start2+kk+1)=0.d0
       wqmask(i,iref_start2+kk+2)=0.d0
       wqmask(i,iref_start2+kk+3)=0.d0
       kk=kk+3
       endif
       if(lll(ibeta,iitype).eq.2) then
       wqmask(i,iref_start2+kk+1)=0.d0
       wqmask(i,iref_start2+kk+2)=0.d0
       wqmask(i,iref_start2+kk+3)=0.d0
       wqmask(i,iref_start2+kk+4)=0.d0
       wqmask(i,iref_start2+kk+5)=0.d0
       kk=kk+5
       endif
       if(lll(ibeta,iitype).gt.2) then
       write(6,*) "lll.gt.2, not programed, stop",
     & ibeta,iitype,lll(ibeta,iitype)
       endif

301    continue

c-------------------------------------------------------
       else   !  q.lt.1.D-6

       kk=0
       do 302 ibeta=1,nbeta(iitype)

       yspd=wq(iq,ibeta,iitype)*f1+
     &  wq(iq+1,ibeta,iitype)*f2+wq(iq+2,ibeta,iitype)*f3

       if(lll(ibeta,iitype).eq.0) then
       wqmask(i,iref_start2+kk+1)=yspd*cc*vins
       kk=kk+1
       endif

       if(lll(ibeta,iitype).eq.1) then
c       wqmask(i,iref_start2+kk+1)=
c     &            dsqrt(3.d0)*cai*gkx_n(i,kpt)/q*yspd*cc*vins
c       wqmask(i,iref_start2+kk+2)=
c     &            dsqrt(3.d0)*cai*gky_n(i,kpt)/q*yspd*cc*vins
c       wqmask(i,iref_start2+kk+3)=
c     &            dsqrt(3.d0)*cai*gkz_n(i,kpt)/q*yspd*cc*vins

cccccc  Yl0, Ylm+=[Ylm+(-)^m Yl-m]/sqrt(2), Ylm=[Ylm-(-)^m Yl-m]/i sqrt(2)  
cccccc cos(phi)=x/sqrt(x**2+y**2)  

       wqmask(i,iref_start2+kk+1)=
     &            dsqrt(3.d0)*cai*gkz_n(i,kpt)/q*yspd*cc*vins    ! Y10
       wqmask(i,iref_start2+kk+2)=
     &           -dsqrt(3.d0)*cai*gkx_n(i,kpt)/q*yspd*cc*vins    ! Y11+
       wqmask(i,iref_start2+kk+3)=
     &           -dsqrt(3.d0)*cai*gky_n(i,kpt)/q*yspd*cc*vins    ! Y11-
       kk=kk+3
       endif

       if(lll(ibeta,iitype).eq.2) then
c       wqmask(i,iref_start2+kk+1)=dsqrt(5.d0)*yspd*cc*vins*(
c     & dsqrt(3.d0)*gkx_n(i,kpt)*gky_n(i,kpt)/q**2)
c       wqmask(i,iref_start2+kk+2)=dsqrt(5.d0)*yspd*cc*vins*(
c     & dsqrt(3.d0)*gkx_n(i,kpt)*gkz_n(i,kpt)/q**2)
c       wqmask(i,iref_start2+kk+3)=dsqrt(5.d0)*yspd*cc*vins*(
c     & dsqrt(3.d0)*gky_n(i,kpt)*gkz_n(i,kpt)/q**2)
c       wqmask(i,iref_start2+kk+4)=dsqrt(5.d0)*yspd*cc*vins*(
c     & dsqrt(3.d0)/2*(gkx_n(i,kpt)**2-gky_n(i,kpt)**2)/q**2)
c       wqmask(i,iref_start2+kk+5)=dsqrt(5.d0)*yspd*cc*vins*(
c     & 1.d0/2*(gkx_n(i,kpt)**2+
c     &              gky_n(i,kpt)**2-2*gkz_n(i,kpt)**2)/q**2)
       wqmask(i,iref_start2+kk+1)=-dsqrt(5.d0)*yspd*cc*vins*(
     & 1.d0/2*(2*gkz_n(i,kpt)**2-gkx_n(i,kpt)**2-
     &              gky_n(i,kpt)**2)/q**2)                       ! Y20
       wqmask(i,iref_start2+kk+2)=dsqrt(5.d0)*yspd*cc*vins*(     ! Y21+
     & dsqrt(3.d0)*gkx_n(i,kpt)*gkz_n(i,kpt)/q**2)
       wqmask(i,iref_start2+kk+3)=dsqrt(5.d0)*yspd*cc*vins*(     ! Y21-
     & dsqrt(3.d0)*gky_n(i,kpt)*gkz_n(i,kpt)/q**2)
       wqmask(i,iref_start2+kk+4)=-dsqrt(5.d0)*yspd*cc*vins*(    ! Y22+
     & dsqrt(3.d0)/2*(gkx_n(i,kpt)**2-gky_n(i,kpt)**2)/q**2)
       wqmask(i,iref_start2+kk+5)=-dsqrt(5.d0)*yspd*cc*vins*(    ! Y22-
     & dsqrt(3.d0)*gkx_n(i,kpt)*gky_n(i,kpt)/q**2)

       kk=kk+5
       endif
       if(lll(ibeta,iitype).gt.2) then
       write(6,*) "lll.gt.2, not programed, stop",
     & ibeta,iitype,lll(ibeta,iitype)
       endif
302    continue

       endif  !  q.lt.1.D-6

 10   continue    ! i=1,ng
 20   continue    ! ia=1,natom

******************************************************************
      call mpi_barrier(MPI_COMM_K,ierr)
      return
      end

