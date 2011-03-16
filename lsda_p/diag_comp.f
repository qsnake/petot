      subroutine diag_comp(ilocal,E_st,err_st,vr,workr_n,kpt,iislda)
*************************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************
***  Note for subspace diagonalization, the purpose is to make:
****  ug_m1 * H * ug_m2 = delta_m1,m2. 
****  as a result, the S operation is not really needed here. 
****  We will diagonalize h(m1,m2), find U.  We don't need to call the
****  general eigen vector solver for h*psi=e s*psi


      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      parameter (lwork=6000)

      include 'param.escan_real'

      complex*16 ugh(mg_nx)
      complex*16 h(mst,mst)
      real*8 E_st(mst),err_st(mst)
      complex*16 workx(lwork)
      real*8 workrx(3*mst)

      real*8 vr(mr_n)
      complex*16 workr_n(mr_n)

      real*8 Dij0(32,32,mtype),Qij(32,32,mtype)
      integer isNLa(9,matom),ipsp_type(mtype)
      complex*16 sumdum(32,matom)



******* extra due to djacobi()
      integer ind(mst)


      complex*16 ctmp(mst),c

      common /comisNLa/isNLa,Dij0,Qij,ipsp_all,ipsp_type


**********************************************
**** prepare the matrix
**********************************************
      ng_n=ngtotnod(inode,kpt)

      do 100 m1=1,mx

      call Hpsi_comp(ug_n(1,m1),ugh,ilocal,vr,workr_n,kpt,1,
     &  sug_n(1,m1),sumdum,iislda)


      do 100 m2=1,m1

      c=dcmplx(0.d0,0.d0)
      do i=1,ng_n
      c=c+ugh(i)*dconjg(ug_n(i,m2))
      enddo

      call global_sumc(c)
      c=c*vol

      h(m1,m2)=dconjg(c)
      h(m2,m1)=c
100   continue

********************************************
********************************************
      call system_czheev('V','U',mx,h,mst,E_st,workx,
     &  lwork,workrx,info)
************************************
********* the unitary rotation

      do 200 i=1,ng_n

      do m=1,mx
       c=dcmplx(0.d0,0.d0)

       do m1=1,mx
       c=c+h(m1,m)*ug_n(i,m1)    ! terrible operation
       enddo

       ctmp(m)=c

      enddo

      do m=1,mx
      ug_n(i,m)=ctmp(m)       ! terrible operation
      enddo
200   continue
*********************************
      if(ipsp_all.eq.2) then
      do 201 i=1,ng_n

      do m=1,mx
       c=dcmplx(0.d0,0.d0)

       do m1=1,mx
       c=c+h(m1,m)*sug_n(i,m1)    ! terrible operation
       enddo

       ctmp(m)=c

      enddo

      do m=1,mx
      sug_n(i,m)=ctmp(m)         ! terrible operation
      enddo
201   continue
      endif
*********************************
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc   test   ! we should not do this one !!!
c      do m=1,mx
c
c      call Hpsi_comp(ug_n(1,m),ugh,ilocal,vr,workr_n,kpt,1,sug_n(1,m),iislda)
c
c      s=0.d0
c
c      if(ipsp_all.eq.1) then
c      do i=1,ng_n
c      s=s+cdabs(ugh(i)-E_st(m)*ug_n(i,m))**2
c      enddo
c       else
c      do i=1,ng_n
c      s=s+cdabs(ugh(i)-E_st(m)*sug_n(i,m))**2
c      enddo
c      endif
c
c
c      call global_sumr(s)
c
c      s=s*vol
c      err_st(m)=dsqrt(s)
c      enddo
***********************************
      if(inode_tot.eq.1) then
       write(6,*) "*********************************"
       write(6,*) "*****kpt= ", kpt
       write(6,*) "report from diag_real"
c       write(6,*) "err of each states, A.U"
c       write(6,102) (err_st(i), i=1,mx)
       write(6,*) "eigen energies, in eV"
       write(6,103) (E_st(i)*27.211396d0, i=1,mx)
       write(6,*) "*********************************"
      endif
101   format(5(i6,7x))
102   format(5(E10.4,3x))
103   format(5(f12.8,1x))

      return
      end



