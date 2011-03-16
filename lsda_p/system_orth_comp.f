      subroutine orth_comp(wg,ug_nt,mt,isign,kpt,
     &    sug_nt,swg,ipsp_all,isug,iss,fnorm)
*************************************************************************
*** Written by Lin-Wang Wang, 2001
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

cccccccccccccccccccccccccccccccccccccccccccc
cccc system dependent part, blas: cgemv or zgemv
cccc Note, there are other subroutines contain zgemv!
cccc use "grep zgemv *.f" to find them out and change them into cgemv for Cray T3E machine
cccccccccccccccccccccccccccccccccccccccccccc
**** isign=1: orth wg to ug_n(1,m) for m=1,mt, then normalize
**** isign=2: orth wg to ug_n(1,m) for m=1,mt, without normalization
**** isign=3: orth wg just to ug(1,mt)
**************************************************
****  for ipsp_all.eq.2, it could, or could not provide swg
****  isug=0.and.ipsp_all=2, do not update swg
****  isug=1.and.ipsp_all=2, do update swg
****  iss=1.and.ipsp_all=2, wg=wg-c*ug  and resulting wg*sug=0
****  iss=-1.and.ipsp_all=2, wg=wg-c*sug and resulting wg*ug=0
****  ipsp_all=2 and isign=1 and iss=-1 is not defined
****  ipsp_all=2 and iss=-1 and isug=1  is not defined
**************************************************

      use data
      use load_data 

      implicit double precision (a-h,o-z)

      include 'param.escan_real'
      include "mpif.h"

      complex*16 wg(mg_nx),swg(1)
      complex*16 sumdumc(mt),sumdumctmp(mt)
      complex*16 ug_nt(mg_nx,mt),sug_nt(mg_nx,mt)     ! mt could be zero. It is okay ?!
cccccc  Note, when ipsp_all=1, the outside sug_nt is (1,mt). But as long as we
cccccc  don't touch sug_nt, it should be okay, and there is no new memory being allocated. 
      complex*16 cc,ccc,cc1,cc2


      call mpi_barrier(MPI_COMM_K,ierr)

       if(ipsp_all.eq.2.and.isign.eq.1.and.iss.eq.-1) then
       write(6,*) "orth: isign.eq.1.and.iss.eq.-1,not defined, stop"
       call mpi_abort(MPI_COMM_WORLD,ierr)
       endif

       if(ipsp_all.eq.2.and.isug.eq.1.and.iss.eq.-1) then
       write(6,*) "orth: isug.eq.1.and.iss.eq.-1,not defined, stop"
       call mpi_abort(MPI_COMM_WORLD,ierr)
       endif


      ng_n = ngtotnod(inode,kpt)

      sum_norm=0.d0
      if(isign.eq.1.and.ipsp_all.eq.2) then
        cc=dcmplx(0.d0,0.d0)
        do i=1,ng_n
        cc=cc+wg(i)*dconjg(swg(i))
        enddo
   
        call global_sumc(cc)
        cc=cc*vol
        sum_norm=cdabs(cc)
      endif
***************************************
      
      if(isign.eq.3.or.mt.le.1) then  

      mi=mt
      if(mt.le.0) goto 25

      do 24 m=mi,mt
      s=0.d0
      cc=dcmplx(0.d0,0.d0)

      if(ipsp_all.eq.1.or.iss.eq.-1) then
      do i=1,ng_n
      cc=cc+wg(i)*dconjg(ug_nt(i,m))
      enddo
      else
      do i=1,ng_n
      cc=cc+wg(i)*dconjg(sug_nt(i,m))
      enddo
      endif

      call global_sumc(cc)
      cc=cc*vol

      sum_norm=sum_norm-cdabs(cc)**2

      if(ipsp_all.eq.1.or.iss.eq.1) then
      do i=1,ng_n
      wg(i)=wg(i)-cc*ug_nt(i,m)
      enddo
      else
      do i=1,ng_n
      wg(i)=wg(i)-cc*sug_nt(i,m)
      enddo
      endif

       if(isug.eq.1.and.ipsp_all.eq.2) then
       do i=1,ng_n
       swg(i)=swg(i)-cc*sug_nt(i,m)
       enddo
       endif

24    continue
25    continue

      else

      cc1=dcmplx(1.d0,0.d0)
      cc2=dcmplx(0.d0,0.d0)

      if(ipsp_all.eq.1.or.iss.eq.-1) then
cccccc   T3E, blas, only 
c      call cgemv('c',ng_n,mt,cc1,ug_nt,mg_nx,wg,1,cc2,sumdumc,1)
cccccccc  IBM SP and all other machine, blas
      call zgemv('C',ng_n,mt,cc1,ug_nt,mg_nx,wg,1,cc2,sumdumc,1)
ccccccccccccccccccccccccccccc
      else
cccccc   T3E, blas, only
c      call cgemv('c',ng_n,mt,cc1,sug_nt,mg_nx,wg,1,cc2,sumdumc,1)
cccccccc  IBM SP, and all other machine blas
      call zgemv('C',ng_n,mt,cc1,sug_nt,mg_nx,wg,1,cc2,sumdumc,1)
ccccccccccccccccccccccccccccc
        
      endif

      call mpi_allreduce(sumdumc,sumdumctmp,mt,
     $     MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_K,ierr)
      sumdumc = sumdumctmp

      do i = 1,mt
       sumdumc(i) = sumdumc(i)*vol
       sum_norm=sum_norm-cdabs(sumdumc(i))**2
      enddo

      cc1=dcmplx(-1.d0,0.d0)
      cc2=dcmplx(1.d0,0.d0)

      if(ipsp_all.eq.1.or.iss.eq.1) then
cccccccc T3E, blas
c      call cgemv('n',ng_n,mt,cc1,ug_nt,mg_nx,sumdumc,1,cc2,wg,1)
cccccccc IBM SP, blas
      call zgemv('N',ng_n,mt,cc1,ug_nt,mg_nx,sumdumc,1,cc2,wg,1)
cccccccccccccccccccccc
      else
cccccccc T3E, blas
c      call cgemv('n',ng_n,mt,cc1,sug_nt,mg_nx,sumdumc,1,cc2,wg,1)
cccccccc IBM SP, blas
      call zgemv('N',ng_n,mt,cc1,sug_nt,mg_nx,sumdumc,1,cc2,wg,1)
      endif


          if(isug.eq.1.and.ipsp_all.eq.2) then
cccccccc T3E, blas only, Note there are other subroutines contains
cccccccc zgemv. Use  "grep zgemv *.f" to find them out and change them for cgemv
c      call cgemv('n',ng_n,mt,cc1,sug_nt,mg_nx,sumdumc,1,cc2,swg,1)
cccccccc IBM SP, and all other machines
      call zgemv('N',ng_n,mt,cc1,sug_nt,mg_nx,sumdumc,1,cc2,swg,1)
cccccccccccccccccccccc
          endif

      endif 

**************************************************
****  renormalize
**************************************************
      if(isign.eq.1.and.ipsp_all.eq.1) then
      s=0.d0

      do i=1,ng_n
      s=s+cdabs(wg(i))**2
      enddo

      call global_sumr(s)

      if(s.lt.1.D-200) then
      write(6,*) "test, warning, s=", s
      endif

      s=dsqrt(1.d0/(s*vol))

      do i=1,ng_n
      wg(i)=s*wg(i)
      enddo
      fnorm=s
      endif
*************************

      if(isign.eq.1.and.ipsp_all.eq.2) then
      if(sum_norm.lt.1.D-200) then
      write(6,*) "sum_norm.lt.1.D-200 in orth, stop",sum_norm
      call mpi_abort(MPI_COMM_WORLD,ierr)
      endif
      s=dsqrt(1.d0/sum_norm)
      do i=1,ng_n
      wg(i)=s*wg(i)
      enddo
        if(isug.eq.1) then
        do i=1,ng_n
        swg(i)=s*swg(i)
        enddo
        endif
      fnorm=s
      endif


      return
      end








