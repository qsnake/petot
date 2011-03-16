      subroutine orth_comp(wg,ug_nt,mt,isign,kpt)
**************************************************
**** isign=1: orth wg to ug_n(1,m) for m=1,mt, then normalize
**** isign=2: orth wg to ug_n(1,m) for m=1,mt, without normalization
**** isign=3: orth wg just to ug(1,mt)
cccccccccccccccccccccccccccccccccccccccccccc
cccc system dependent part, blas: cgemv or zgemv
cccccccccccccccccccccccccccccccccccccccccccc

      use data
      use load_data 

      implicit double precision (a-h,o-z)

      include 'param.escan_real'
      include "mpif.h"

      complex*16 wg(mg_nx)
      complex*16 sumdumc(mt)
      complex*16 ug_nt(mg_nx,mx)
      complex*16 cc,ccc,cc1,cc2

c     integer iorg(2) ! node and ig number of origin


      call mpi_barrier(MPI_COMM_WORLD,ierr)

      ng_n = ngtotnod(inode,kpt)
      
      if(isign.eq.3.or.mt.le.1) then  

      mi=mt
      if(mt.le.0) goto 25

      do 24 m=mi,mt
      s=0.d0
      cc=dcmplx(0.d0,0.d0)

      do i=1,ng_n
      cc=cc+wg(i)*dconjg(ug_nt(i,m))
      enddo

      call global_sumc(cc)
      cc=cc*vol

      do i=1,ng_n
      wg(i)=wg(i)-cc*ug_nt(i,m)
      enddo

24    continue
25    continue

      else

      cc1=dcmplx(1.d0,0.d0)
      cc2=dcmplx(0.d0,0.d0)

cccccc   T3E, blas, 
      call cgemv('c',ng_n,mt,cc1,ug_nt,mg_nx,wg,1,cc2,sumdumc,1)
cccccccc  IBM SP, blas
c      call zgemv('C',ng_n,mt,cc1,ug_nt,mg_nx,wg,1,cc2,sumdumc,1)
ccccccccccccccccccccccccccccc

      call mpi_allreduce(sumdumc,sumdumc,mt,
     $     MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_WORLD,ierr)

      do i = 1,mt
       sumdumc(i) = sumdumc(i)*vol
      enddo

      cc1=dcmplx(-1.d0,0.d0)
      cc2=dcmplx(1.d0,0.d0)

cccccccc T3E, blas
      call cgemv('n',ng_n,mt,cc1,ug_nt,mg_nx,sumdumc,1,cc2,wg,1)
cccccccc IBM SP, blas
cc      call zgemv('N',ng_n,mt,cc1,ug_nt,mg_nx,sumdumc,1,cc2,wg,1)
cccccccccccccccccccccc

      endif 

**************************************************
****  renormalize
**************************************************
      if(isign.eq.1) then
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
      endif

      return
      end








