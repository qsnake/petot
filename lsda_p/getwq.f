      subroutine getwq(AL,ntype,iatom,ityatom,xatom,kpt)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
******************************************


      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'param.escan_real'

      real*8 xatom(3,matom)

      real*8 AL(3,3)

      real*8 qi3(mnq),wq2(mnq,3,mtype)

      real*8 occ_t(mtype)
      integer iiatom(mtype),iatom(matom),icore(mtype),numref(matom)
      integer is_ref(mtype),ip_ref(mtype),id_ref(mtype)
      integer ityatom(matom)


      complex*16 cc,cai

***************************************************
****  xatom(1),xatom(2),xatom(3) are the coord in unit of AL(3,3)
****  supercell edges
***************************************************

      common /comline2/qi3,wq2
      common /comNL2/occ_t,iiatom,icore,numref
      common /comispd_ref/is_ref,ip_ref,id_ref

*******************************************************
**** generate the Kleiman-Bylander reference wavefunction
*******************************************************
      ng_n=ngtotnod(inode,kpt)

      vins=1.d0/vol
      cai=dcmplx(0.d0,1.d0)

      do 20 ia=1,natom

      iitype=ityatom(ia)

      is_ref1=is_ref(iitype)
      ip_ref1=ip_ref(iitype)
      id_ref1=id_ref(iitype)
      kk=0
      if(is_ref1.eq.1) kk=kk+1
      if(ip_ref1.eq.1) kk=kk+3
      if(id_ref1.eq.1) kk=kk+5
      numref(ia)=kk


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

      iq=1+q*(mnq-1.d0)/qi3(mnq)

      x=(q-qi3(iq))/(qi3(iq+1)-qi3(iq))

      f1=1-x-0.5d0*x*(1-x)
      f2=x+x*(1-x)
      f3=-0.5d0*x*(1-x)

      ys=wq2(iq,1,iitype)*f1+wq2(iq+1,1,iitype)*f2+
     &   wq2(iq+2,1,iitype)*f3
      yp=wq2(iq,2,iitype)*f1+wq2(iq+1,2,iitype)*f2+
     &   wq2(iq+2,2,iitype)*f3
      yd=wq2(iq,3,iitype)*f1+wq2(iq+1,3,iitype)*f2+
     &   wq2(iq+2,3,iitype)*f3
       

       if(q.lt.1.D-6) then
       kk=0
       if(is_ref1.eq.1) then
       wqmask(kk+1,i,ia)=ys*cc*vins
       kk=kk+1
       endif
       if(ip_ref1.eq.1) then
       wqmask(kk+1,i,ia)=0.d0
       wqmask(kk+2,i,ia)=0.d0
       wqmask(kk+3,i,ia)=0.d0
       kk=kk+3
       endif
       if(id_ref1.eq.1) then
       wqmask(kk+1,i,ia)=0.d0
       wqmask(kk+2,i,ia)=0.d0
       wqmask(kk+3,i,ia)=0.d0
       wqmask(kk+4,i,ia)=0.d0
       wqmask(kk+5,i,ia)=0.d0
       kk=kk+5
       endif
       else

       kk=0
       if(is_ref1.eq.1) then
       wqmask(kk+1,i,ia)=ys*cc*vins
       kk=kk+1
       endif

       if(ip_ref1.eq.1) then
       wqmask(kk+1,i,ia)=dsqrt(3.d0)*cai*gkx_n(i,kpt)/q*yp*cc*vins
       wqmask(kk+2,i,ia)=dsqrt(3.d0)*cai*gky_n(i,kpt)/q*yp*cc*vins
       wqmask(kk+3,i,ia)=dsqrt(3.d0)*cai*gkz_n(i,kpt)/q*yp*cc*vins
       kk=kk+3
       endif

       if(id_ref1.eq.1) then
       wqmask(kk+1,i,ia)=dsqrt(5.d0)*yd*cc*vins*(
     & dsqrt(3.d0)*gkx_n(i,kpt)*gky_n(i,kpt)/q**2)
       wqmask(kk+2,i,ia)=dsqrt(5.d0)*yd*cc*vins*(
     & dsqrt(3.d0)*gkx_n(i,kpt)*gkz_n(i,kpt)/q**2)
       wqmask(kk+3,i,ia)=dsqrt(5.d0)*yd*cc*vins*(
     & dsqrt(3.d0)*gky_n(i,kpt)*gkz_n(i,kpt)/q**2)
       wqmask(kk+4,i,ia)=dsqrt(5.d0)*yd*cc*vins*(
     & dsqrt(3.d0)/2*(gkx_n(i,kpt)**2-gky_n(i,kpt)**2)/q**2)
       wqmask(kk+5,i,ia)=dsqrt(5.d0)*yd*cc*vins*(
     & 1.d0/2*(gkx_n(i,kpt)**2+
     &              gky_n(i,kpt)**2-2*gkz_n(i,kpt)**2)/q**2)
       endif
       endif

 10   continue
 20   continue

******************************************************************
      return
      end

