      subroutine getwmask(xatom,nmap,indmtmp,iitype,wmasktmp,
     &  xyzmaptmp,AL,workr_n,mrb2,nref,inew,iend)
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
ccc inew=1, new type (iitype) of atoms compared to previous one
ccc iend=1, the last atom of this type (iitype)

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'mpif.h'
      include 'param.escan_real'

      real*8 xatom(3)
      real*8 AL(3,3),ALntmp(3,3)
      real*8 workr_n(mr_n)     ! use only half of the space

      real*8 dx(3)
      real*8 wmasktmp(20,mrb2)   ! should be (32,mrb2) if we have double (s,p,d,f) ref.
      real*8 xyzmaptmp(3,mrb2)
      integer indmtmp(mrb2)
      real*8,allocatable,dimension(:)    :: fr  
      real*8  ymask_tmp(mrb2)

      real*8 qi(mnq),wq(mnq,8,mtype)
      real*8 ri(201),amr(201)

      integer lll(8,mtype),nbeta(mtype)


      complex*16 cc,cYY,cc_tmp

***************************************************
****  xatom(1),xatom(2),xatom(3) are the coord in unit of AL(3,3)
****  supercell edges
***************************************************

      common /comline/qi,wq,ri,amr
      common /comlll/lll,nbeta

      ng2_n=ngtotnod2(inode)
      f_sqrt3=dsqrt(3.d0)
      f_sqrt5=dsqrt(5.d0)

      allocate(fr(mr_n))
*******************************************************
**** generate the Kleiman-Bylander reference wavefunction
*******************************************************
      
      nh1=(n1+1)/2+1

      vins=1.d0/vol
cccccccccccc  
      if(inew.eq.1) then
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

ccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc Yl0, Ylm+=[Ylm+(-)^m Yl-m]/sqrt(2), Ylm=[Ylm-(-)^m Yl-m]/i sqrt(2)
ccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(ltmp.eq.0) then
      YY=1.d0
      else
      if(q.lt.1.D-6) then
      YY=0.d0
      else
       if(ltmp.eq.1.and.lm.eq.1) then
       YY=f_sqrt3*gkz2_n(i)/q    ! Y10
       endif
       if(ltmp.eq.1.and.lm.eq.2) then
       YY=-f_sqrt3*gkx2_n(i)/q   ! Y11+
       endif
       if(ltmp.eq.1.and.lm.eq.3) then
       YY=-f_sqrt3*gky2_n(i)/q   ! Y11-
       endif
       if(ltmp.eq.2.and.lm.eq.1) then
       YY=-f_sqrt5/2*
     &       (2*gkz2_n(i)**2-gkx2_n(i)**2-gky2_n(i)**2)/q**2   ! Y20
       endif
       if(ltmp.eq.2.and.lm.eq.2) then
       YY=f_sqrt5*f_sqrt3*gkx2_n(i)*gkz2_n(i)/q**2    ! Y21+
       endif
       if(ltmp.eq.2.and.lm.eq.3) then
       YY=f_sqrt5*f_sqrt3*gky2_n(i)*gkz2_n(i)/q**2    ! Y21-
       endif
       if(ltmp.eq.2.and.lm.eq.4) then
       YY=-f_sqrt5*f_sqrt3/2*
     &       (gkx2_n(i)**2-gky2_n(i)**2)/q**2                  ! Y22+
       endif
       if(ltmp.eq.2.and.lm.eq.5) then                                      
       YY=-f_sqrt5*f_sqrt3*gkx2_n(i)*gky2_n(i)/q**2   ! Y22-
       endif

      endif
      endif

      YYMask(i,kk)=y*YY*vins

 10   continue
 11   continue
      endif     
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
      do i=1,ng2_n
      ph=gkx2_n(i)*x11+gky2_n(i)*y11+gkz2_n(i)*z11
      cc=cdexp(dcmplx(0.d0,ph))*YYMask(i,kk)*cc_tmp
      fr(i*2-1)=fr(i*2-1)+dreal(cc)
      fr(i*2)=fr(i*2)+dimag(cc)
      enddo

      call d3fft_real2(fr,workr_n,-1,0)

      fr = workr_n
       
******************************************************
**** x1,y1,z1 are real number grid indexes of the atom
******************************************************
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
      if(imap.gt.mrb2) then
      write(6,*) "imap > mrb2, stop", imap,mrb2
      call mpi_abort(MPI_COMM_WORLD,ierr)
      endif
      indmtmp(imap)=ii
      r2=r/rcut
      ir=1+r2*200.d0
      f1=(ri(ir+1)-r2)/(ri(ir+1)-ri(ir))
      f2=(r2-ri(ir))/(ri(ir+1)-ri(ir))
      y=amr(ir)*f1+amr(ir+1)*f2
      ymask_tmp(imap)=y
      xyzmaptmp(1,imap)=xt
      xyzmaptmp(2,imap)=yt
      xyzmaptmp(3,imap)=zt
50    continue
      nmap=imap
      endif       ! kk=1
ccccccccccccccccccccccccccccccccccccccccccc

      do imap=1,nmap
      wmasktmp(kk,imap)=fr(indmtmp(imap))*ymask_tmp(imap)
      enddo
      
1000  continue

      nref=kk

      deallocate(fr)

      if(iend.eq.1) then
      call data_deallocate_YYMask()
      endif

      return
      end
	 

