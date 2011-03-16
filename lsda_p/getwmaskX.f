      subroutine getwmaskX(xatom,nmap,indmtmp,iitype,wmaskXtmp,
     &  AL,workr_n,mrb2,nref,inew,iend)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************

*****************************************************
*****************************************************

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include "mpif.h"
      include 'param.escan_real'

      real*8 xatom(3)
      real*8 AL(3,3),ALntmp(3,3)
      real*8 workr_n(mr_n)     ! use only half of the space

      real*8 dx(3)
      real*8 wmaskXtmp(20,mrb2,3)
      integer indmtmp(mrb2)
      real*8,allocatable,dimension(:)   :: fr,fr1,fr2,fr3
      real*8,allocatable,dimension(:)   :: ymask_tmp,ymask_tmp1,
     & ymask_tmp2,ymask_tmp3

      real*8 qi(mnq),wq(mnq,8,mtype)
      real*8 ri(201),amr(201)
      integer lll(8,mtype),nbeta(mtype)


      complex*16 cc,cYY,cc2,cai,cc_tmp

***************************************************
****  xatom(1),xatom(2),xatom(3) are the coord in unit of AL(3,3)
****  supercell edges
***************************************************

      common /comline/qi,wq,ri,amr
      common /comlll/lll,nbeta

      ng2_n=ngtotnod2(inode)
      cai=dcmplx(0.d0,1.d0)
      f_sqrt3=dsqrt(3.d0)
      f_sqrt5=dsqrt(5.d0)

      allocate(fr(mr_n))
      allocate(fr1(mr_n))
      allocate(fr2(mr_n))
      allocate(fr3(mr_n))
      allocate(ymask_tmp(mrb2))
      allocate(ymask_tmp1(mrb2))
      allocate(ymask_tmp2(mrb2))
      allocate(ymask_tmp3(mrb2))
*******************************************************
**** generate the Kleiman-Bylander reference wavefunction
*******************************************************

      nh1=(n1+1)/2+1

      vins=1.d0/vol
      nref=0
      nmap=0

cccccccccccccccccccc
      if(inew.eq.1) then      ! calculate YYMask, save it for later atoms with the same atom type
      kk=0
      do ibeta=1,nbeta(iitype)
      do lm=1,2*lll(ibeta,iitype)+1
      kk=kk+1
      enddo
      enddo
      nkk=kk


      call data_allocate_YYMask(ng2_n,nkk)

      kk=0
      do 11 ibeta=1,nbeta(iitype)
      ltmp=lll(ibeta,iitype)
      do 11 lm=1,2*ltmp+1
      kk=kk+1

      do 10 i=1,ng2_n

      q=dsqrt(gkx2_n(i)**2+gky2_n(i)**2+gkz2_n(i)**2)

      iq=1+q*(mnq-1.d0)/qi(mnq)

      x=(q-qi(iq))/(qi(iq+1)-qi(iq))

      f1=1-x-0.5d0*x*(1-x)
      f2=x+x*(1-x)
      f3=-0.5d0*x*(1-x)

      y=wq(iq,ibeta,iitype)*f1+wq(iq+1,ibeta,iitype)*f2+
     &  wq(iq+2,ibeta,iitype)*f3


      if(ltmp.eq.0) then
      YY=1.d0
      else

      if(q.lt.1.D-6) then
      YY=0.d0
      else
       if(ltmp.eq.1.and.lm.eq.1) then
       YY=f_sqrt3*gkz2_n(i)/q                 ! Y10
       endif
       if(ltmp.eq.1.and.lm.eq.2) then
       YY=-f_sqrt3*gkx2_n(i)/q                ! Y11+
       endif
       if(ltmp.eq.1.and.lm.eq.3) then
       YY=-f_sqrt3*gky2_n(i)/q                ! Y11-
       endif
       if(ltmp.eq.2.and.lm.eq.1) then
       YY=-f_sqrt5/2*
     &       (2*gkz2_n(i)**2-gkx2_n(i)**2-gky2_n(i)**2)/q**2    ! Y20
       endif
       if(ltmp.eq.2.and.lm.eq.2) then
       YY=f_sqrt5*f_sqrt3*gkx2_n(i)*gkz2_n(i)/q**2     ! Y21+
       endif
       if(ltmp.eq.2.and.lm.eq.3) then
       YY=f_sqrt5*f_sqrt3*gky2_n(i)*gkz2_n(i)/q**2     ! Y21-
       endif
       if(ltmp.eq.2.and.lm.eq.4) then
       YY=-f_sqrt5*f_sqrt3/2*
     &       (gkx2_n(i)**2-gky2_n(i)**2)/q**2                   ! Y22+
       endif
       if(ltmp.eq.2.and.lm.eq.5) then
       YY=-f_sqrt5*f_sqrt3*gkx2_n(i)*gky2_n(i)/q**2    ! Y22-
       endif

      endif
      endif

      YYMask(i,kk)=y*YY*vins

 10   continue
 11   continue
      endif     ! inew=1
cc  inew=1, just to get cYYMask, store everything, cYYMask will be saved, using Module
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      x1=xatom(1)*n1
      y1=xatom(2)*n2
      z1=xatom(3)*n3
      if(x1.lt.0.d0) x1=x1+n1
      if(y1.lt.0.d0) y1=y1+n2
      if(z1.lt.0.d0) z1=z1+n3
      if(x1.gt.n1) x1=x1-n1
      if(y1.gt.n2) y1=y1-n2
      if(z1.gt.n3) z1=z1-n3

      
      x11=AL(1,1)*x1/n1+AL(1,2)*y1/n2+AL(1,3)*z1/n3
      y11=AL(2,1)*x1/n1+AL(2,2)*y1/n2+AL(2,3)*z1/n3
      z11=AL(3,1)*x1/n1+AL(3,2)*y1/n2+AL(3,3)*z1/n3

      kk=0
      do 1000 ibeta=1,nbeta(iitype)
      ltmp=lll(ibeta,iitype)
      cc_tmp=dcmplx(1.d0,0.d0)
      if(ltmp.eq.1) cc_tmp=dcmplx(0.d0,1.d0)
               
      do 1000 lm=1,2*ltmp+1
      kk=kk+1

      fr=0.d0
      fr1=0.d0
      fr2=0.d0
      fr3=0.d0
      do i=1,ng2_n
      ph=gkx2_n(i)*x11+gky2_n(i)*y11+gkz2_n(i)*z11
      cc=cdexp(dcmplx(0.d0,ph))*YYMask(i,kk)*cc_tmp


      fr(i*2-1)=fr(i*2-1)+dreal(cc)
      fr(i*2)=fr(i*2)+dimag(cc)

      cc2=cc*cai*gkx2_n(i)

      fr1(i*2-1)=fr1(i*2-1)+dreal(cc2)
      fr1(i*2)=fr1(i*2)+dimag(cc2)

      cc2=cc*cai*gky2_n(i)

      fr2(i*2-1)=fr2(i*2-1)+dreal(cc2)
      fr2(i*2)=fr2(i*2)+dimag(cc2)

      cc2=cc*cai*gkz2_n(i)

      fr3(i*2-1)=fr3(i*2-1)+dreal(cc2)
      fr3(i*2)=fr3(i*2)+dimag(cc2)
      enddo

      call d3fft_real2(fr,workr_n,-1,0)
      fr = workr_n

      call d3fft_real2(fr1,workr_n,-1,0)
      fr1 = workr_n

      call d3fft_real2(fr2,workr_n,-1,0)
      fr2 = workr_n

      call d3fft_real2(fr3,workr_n,-1,0)
      fr3 = workr_n

******************************************************************
******************************************************
**** x1,y1,z1 are real number grid indexes of the atom
******************************************************
*******************************************************
**** The actual mapping
**** msb is a shifting param, so that all points inside rrcut
**** is in the do 50 loop
*******************************************************
      if(kk.eq.1) then   ! get indmtmp(imap) and ymask_tmp(imap)
      imap=0
      ALntmp(:,1)=AL(:,1)/n1
      ALntmp(:,2)=AL(:,2)/n2
      ALntmp(:,3)=AL(:,3)/n3
               
      do 50 ii=1,nr_n
      jj=ii+(inode-1)*nr_n
      i=(jj-1)/(n2*n3)
      j=(jj-1-i*n2*n3)/n3
      k=jj-i*n2*n3-j*n3-1

      x10=i-x1
      if(dabs(x10-n1).lt.dabs(x10)) x10=x10-n1
      if(dabs(x10+n1).lt.dabs(x10)) x10=x10+n1
      y10=j-y1
      if(dabs(y10-n2).lt.dabs(y10)) y10=y10-n2
      if(dabs(y10+n2).lt.dabs(y10)) y10=y10+n2
      z10=k-z1
      if(dabs(z10-n3).lt.dabs(z10)) z10=z10-n3
      if(dabs(z10+n3).lt.dabs(z10)) z10=z10+n3

************************************************
****  rr the distance from the atom to the grid point
************************************************
      xt=ALntmp(1,1)*x10+ALntmp(1,2)*y10+ALntmp(1,3)*z10
      yt=ALntmp(2,1)*x10+ALntmp(2,2)*y10+ALntmp(2,3)*z10
      zt=ALntmp(3,1)*x10+ALntmp(3,2)*y10+ALntmp(3,3)*z10

      rr=xt**2+yt**2+zt**2

      
      if(rr.ge.rcut**2-1.D-6) goto 50
************************************************
************************************************
      r=dsqrt(rr)
      imap=imap+1

      indmtmp(imap)=ii
      r2=r/rcut
      ir=1+r2*200.d0
      f1=(ri(ir+1)-r2)/(ri(ir+1)-ri(ir))
      f2=(r2-ri(ir))/(ri(ir+1)-ri(ir))
      rmask=amr(ir)*f1+amr(ir+1)*f2

      if(ir.eq.1) ir=ir+1
      if(ir.ge.200) ir=200-1
      damr0=(amr(ir+1)-amr(ir-1))/(ri(ir+1)-ri(ir-1))
      damr1=(amr(ir+2)-amr(ir))/(ri(ir+2)-ri(ir))
      drmask=(damr0*f1+damr1*f2)/rcut

      ymask_tmp(imap)=rmask
      if(r.lt.1.E-9) then
      ymask_tmp1(imap)=0.d0
      ymask_tmp2(imap)=0.d0
      ymask_tmp3(imap)=0.d0
      else
      ymask_tmp1(imap)=-drmask*xt/r
      ymask_tmp2(imap)=-drmask*yt/r
      ymask_tmp3(imap)=-drmask*zt/r
      endif


 50   continue
      nmap=imap
        if(nmap.gt.mrb2) then
        write(6,*) "nmap > mrb2, stop", nmap,mrb2
        call mpi_abort(MPI_COMM_WORLD,ierr)
        endif
      endif       ! kk=1, save the ymask_tmp for other kk
ccccccccccccccccccccccccccccccccccccccccccccccccc

      do imap=1,nmap
      wmaskXtmp(kk,imap,1)=fr(indmtmp(imap))*
     &  ymask_tmp1(imap)+
     &     fr1(indmtmp(imap))*ymask_tmp(imap)

      wmaskXtmp(kk,imap,2)=fr(indmtmp(imap))*
     &  ymask_tmp2(imap)+
     &     fr2(indmtmp(imap))*ymask_tmp(imap)

      wmaskXtmp(kk,imap,3)=fr(indmtmp(imap))*
     &  ymask_tmp3(imap)+
     &     fr3(indmtmp(imap))*ymask_tmp(imap)
      enddo

1000  continue

      nref=kk

      deallocate(fr)
      deallocate(fr1)
      deallocate(fr2)
      deallocate(fr3)
      deallocate(ymask_tmp)
      deallocate(ymask_tmp1)
      deallocate(ymask_tmp2)
      deallocate(ymask_tmp3)

      if(iend.eq.1) then
      call data_deallocate_YYMask()
      endif

      return
      end
	 

