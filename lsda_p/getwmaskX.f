      subroutine getwmaskX(xatom,nmap,indmtmp,iitype,wmaskXtmp,
     &  AL,workr_n,ixyz,mrb2,is_ref,ip_ref,id_ref,nref)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
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
      real*8 AL(3,3)
      real*8 workr_n(mr_n)     ! use only half of the space

      real*8 dx(3)
      real*8 wmaskXtmp(9,mrb2)
      integer indmtmp(mrb2)
      real*8,allocatable,dimension(:)   :: fr,fr1

      real*8 qi(mnq),wq(mnq,3,mtype)
      real*8 ri(201),amr(201)


      complex*16 cc,cYY,cc2,cai

***************************************************
****  xatom(1),xatom(2),xatom(3) are the coord in unit of AL(3,3)
****  supercell edges
***************************************************

      common /comline/qi,wq,ri,amr

      ng2_n=ngtotnod2(inode)
      cai=dcmplx(0.d0,1.d0)

      allocate(fr(mr_n))
      allocate(fr1(mr_n))
*******************************************************
**** generate the Kleiman-Bylander reference wavefunction
*******************************************************
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

      nh1=(n1+1)/2+1

      vins=1.d0/vol
      nref=0

      do 1000 iref=1,9

      if(iref.eq.1) isp=1
      if(iref.ge.2.and.iref.le.4) isp=2
      if(iref.ge.5) isp=3

      if(iref.eq.1.and.is_ref.eq.0) goto 1000
      if(iref.ge.2.and.iref.le.4.and.ip_ref.eq.0) goto 1000
      if(iref.ge.5.and.id_ref.eq.0) goto 1000
      nref=nref+1


      fr=0.d0
      fr1=0.d0

      do 10 i=1,ng2_n
      ph=gkx2_n(i)*x11+gky2_n(i)*y11+gkz2_n(i)*z11

      cc=cdexp(dcmplx(0.d0,ph))

      q=dsqrt(gkx2_n(i)**2+gky2_n(i)**2+gkz2_n(i)**2)

      iq=1+q*(mnq-1.d0)/qi(mnq)

      x=(q-qi(iq))/(qi(iq+1)-qi(iq))

      f1=1-x-0.5d0*x*(1-x)
      f2=x+x*(1-x)
      f3=-0.5d0*x*(1-x)

      y=wq(iq,isp,iitype)*f1+wq(iq+1,isp,iitype)*f2+
     &  wq(iq+2,isp,iitype)*f3


      if(iref.eq.1) then
      cYY=dcmplx(1.d0,0.d0)
      else

      if(q.lt.1.D-6) then
      cYY=dcmplx(0.d0,0.d0) 
      else
       if(iref.eq.2) then
       cYY=dsqrt(3.d0)*dcmplx(0.d0,gkx2_n(i)/q)
       endif
       if(iref.eq.3) then
       cYY=dsqrt(3.d0)*dcmplx(0.d0,gky2_n(i)/q)
       endif
       if(iref.eq.4) then
       cYY=dsqrt(3.d0)*dcmplx(0.d0,gkz2_n(i)/q)
       endif
       if(iref.eq.5) then
       cYY=dsqrt(5.d0)*dsqrt(3.d0)*gkx2_n(i)*gky2_n(i)/q**2
       endif
       if(iref.eq.6) then
       cYY=dsqrt(5.d0)*dsqrt(3.d0)*gkx2_n(i)*gkz2_n(i)/q**2
       endif
       if(iref.eq.7) then
       cYY=dsqrt(5.d0)*dsqrt(3.d0)*gky2_n(i)*gkz2_n(i)/q**2
       endif
       if(iref.eq.8) then
       cYY=dsqrt(5.d0)*dsqrt(3.d0)/2*
     &       (gkx2_n(i)**2-gky2_n(i)**2)/q**2
       endif
       if(iref.eq.9) then
       cYY=dsqrt(5.d0)/2*
     &       (gkx2_n(i)**2+gky2_n(i)**2-2*gkz2_n(i)**2)/q**2
       endif

      endif
      endif

      fr(i*2-1)=fr(i*2-1)+y*dreal(cYY*cc)*vins
      fr(i*2)=fr(i*2)+y*dimag(cYY*cc)*vins

      if(ixyz.eq.1) cc2=cYY*cc*cai*gkx2_n(i)
      if(ixyz.eq.2) cc2=cYY*cc*cai*gky2_n(i)
      if(ixyz.eq.3) cc2=cYY*cc*cai*gkz2_n(i)

      fr1(i*2-1)=fr1(i*2-1)+y*dreal(cc2)*vins
      fr1(i*2)=fr1(i*2)+y*dimag(cc2)*vins


 10   continue

      call d3fft_real2(fr,workr_n,-1,0)
      fr = workr_n

      call d3fft_real2(fr1,workr_n,-1,0)
      fr1 = workr_n
      

******************************************************************

******************************************************
**** x1,y1,z1 are real number grid indexes of the atom
******************************************************
*******************************************************
**** The actual mapping
**** msb is a shifting param, so that all points inside rrcut
**** is in the do 50 loop
*******************************************************
cccccc This is a over kill, might be reduced if it is too time
cccccc consuming.
***************************************
c in parallel code each PE will go through grid points that it holds
c each has ncolz  z cols. 
c

      imap=0

      ii1 = 0

      do 50 ico = 1,ncolz

      i  = ixcol_z(ico,inode) - 1
      x10=i-x1
      if(dabs(x10-n1).lt.dabs(x10)) x10=x10-n1
      if(dabs(x10+n1).lt.dabs(x10)) x10=x10+n1


      j  = iycol_z(ico,inode) - 1
      y10=j-y1
      if(dabs(y10-n2).lt.dabs(y10)) y10=y10-n2
      if(dabs(y10+n2).lt.dabs(y10)) y10=y10+n2

      do 50 k=0,n3-1
      ii1 = ii1 + 1
      z10=k-z1
      if(dabs(z10-n3).lt.dabs(z10)) z10=z10-n3
      if(dabs(z10+n3).lt.dabs(z10)) z10=z10+n3

************************************************
****  rr the distance from the atom to the grid point
************************************************
      xt=AL(1,1)*x10/n1+AL(1,2)*y10/n2+AL(1,3)*z10/n3
      yt=AL(2,1)*x10/n1+AL(2,2)*y10/n2+AL(2,3)*z10/n3
      zt=AL(3,1)*x10/n1+AL(3,2)*y10/n2+AL(3,3)*z10/n3

      rr=xt**2+yt**2+zt**2

      r=dsqrt(rr)
      
      if(r.ge.rcut-1.D-6) goto 50
************************************************
************************************************
      imap=imap+1

      indmtmp(imap)=ii1

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
      if(r.lt.1.E-9) then
       drmask=0.d0
      else
      if(ixyz.eq.1) drmask=-drmask*xt/r
      if(ixyz.eq.2) drmask=-drmask*yt/r
      if(ixyz.eq.3) drmask=-drmask*zt/r
      endif

      wmaskXtmp(nref,imap)=fr(ii1)*drmask+fr1(ii1)*rmask

 50   continue

      nmap=imap

      if(nmap.gt.mrb2) then
      write(6,*) "nmap > mrb2, stop", nmap,mrb2
      call mpi_abort(MPI_COMM_WORLD,ierr)
      endif

1000  continue


      deallocate(fr)
      deallocate(fr1)

      return
      end
	 

