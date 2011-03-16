      subroutine w_line(ilocal,ntype,vwr_atom,Ealpha)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************


      use data

      implicit double precision(a-h,o-z)


      include "mpif.h"
      include "param.escan_real"
      real*8 r(2000),ws(2000),wp(2000),wd(2000)
      real*8 vloc(2000),vs(2000),vp(2000),vd(2000)
      real*8 vlocT(2000),rhoc(2000)
      real*8 vw(2000),vw2(2000)
      real*8 qfuncLM0(12,2000,mtype),r_at(2000,mtype)
      real*8 a_r(mtype),b_r(mtype)

ccccccc  assume ws,wp and vs,vp,vd has the same r(i)
      real*8 qi2(mnq),vq(mnq,mtype),rhoq(mnq,mtype)
      real*8 rhocq(mnq,mtype)
      real*8 vqT(mnq,mtype)
      real*8 rcut_q1(mtype),rcut_q2(mtype),rcut_qm

cccccccc wq has 1/rmask for real space implementation
cccccccc    and has no 1/rmask for q space implementation

      real*8 qi(mnq),wq(mnq,8,mtype)
      real*8 ri(201),amr(201)

      integer  iiatom(mtype),numref(matom),ipsp_type(mtype),
     &   ityatom(matom)
      integer  lll(8,mtype),nbeta(mtype)

      real*8  zatom(mtype),Ealpha(mtype)
      real*8  occ_s(mtype),occ_p(mtype),occ_d(mtype),occ_t(mtype)
      integer isNL(3,mtype),icore(mtype)
      integer is_ref(mtype),ip_ref(mtype),id_ref(mtype)
      integer is_TB(mtype),ip_TB(mtype),id_TB(mtype)
      integer isNLa(9,matom)
      real*8  Dij0(32,32,mtype),Qij(32,32,mtype)
      real*8  qijrad(0:6,32,32,mtype)

      real*8  q_LMnm(32,32,49,mtype),q_LMnm0(32,32,12,mtype)


      character*20  vwr_atom(mtype)

      common /comVrho/qi2,vq,rhoq,vqT,rhocq
      common /comline/qi,wq,ri,amr
      common /comNL2/occ_t,iiatom,icore,numref,ityatom
      common /comzatom/zatom
      common /comisNL/isNL
      common /comisNLa/isNLa,Dij0,Qij,ipsp_all,ipsp_type
      common /com_qLMnm/q_LMnm,q_LMnm0
      common /comispd_ref/is_ref,ip_ref,id_ref
      common /comlll/lll,nbeta
      common /com_rcut_q/rcut_q1,rcut_q2,rcut_qm
      common /com_qfuncLM0/qfuncLM0,r_at,a_r,b_r


      if(ilocal.eq.2) then
      open(10,file='maskr',status='old',action='read',iostat=ierr)
      if(ierr.ne.0) then
      if(inode.eq.1) 
     & write(6,*) "file ***maskr*** is needed for ilocal=2,stop"
      stop
      endif

      rewind(10)
      do i=1,201
      read(10,*) ri(i),amr(i)
      enddo
      close(10)
      else
      do i=1,201
      ri(i)=10.d0*(i-1)/200.d0
      amr(i)=1.d0             ! for ilocal.eq.3, amr=1, no mask, cannot to have amr
      enddo
      endif


      rcut_qm=0.d0
      do 2000 ia=1,ntype
      if(ipsp_type(ia).eq.1) then

      call w_line_vwr(vwr_atom(ia),zatom(ia),iiatom(ia),
     &  Ealpha(ia),occ_t(ia),isNL(1,ia),Dij0(1,1,ia),
     &  Qij(1,1,ia),icore(ia),
     &  is_ref(ia),ip_ref(ia),id_ref(ia),qi,wq(1,1,ia),qi2,vq(1,ia),
     &  rhoq(1,ia),rhocq(1,ia),vqT(1,ia),ri,amr,lll(1,ia),nbeta(ia))

       rcut_q1(ia)=0.d0
       rcut_q2(ia)=0.d0

      else 

      call w_line_usp(vwr_atom(ia),zatom(ia),iiatom(ia),
     &  Ealpha(ia),occ_t(ia),isNL(1,ia),Dij0(1,1,ia),
     &  Qij(1,1,ia),qijrad(0,1,1,ia),icore(ia),
     &  is_ref(ia),ip_ref(ia),id_ref(ia),qi,wq(1,1,ia),qi2,vq(1,ia),
     &  rhoq(1,ia),rhocq(1,ia),vqT(1,ia),ri,amr,lll(1,ia),nbeta(ia),
     &  rcut_q1(ia),rcut_q2(ia),qfuncLM0(1,1,ia),
     &  r_at(1,ia),a_r(ia),b_r(ia))

       call getq_LMnm(q_LMnm(1,1,1,ia),q_LMnm0(1,1,1,ia),
     &              qijrad(0,1,1,ia),ia)
cccc WARNING, we might need to have different q_LMnm(L,M=0,0) from different
cccc n,m (l). They are very different. Correspondingly, we need different
cccc wmask_q for these different (L,M=0,0,and l). 

      endif
       if(rcut_q2(ia).gt.rcut_qm) rcut_qm=rcut_q2(ia)
2000  continue

      return

      contains

      subroutine getq_LMnm(q_LMnmtt,q_LMnm0,qijradt,ia)
      implicit double precision(a-h,o-z)
      real*8  qijradt(0:6,32,32)
      complex*16  q_LMnmt(32,32,49),q_LMnmtmp(32,32,49)
      integer indec(-6:6,0:6)
      complex*16 cai
      real*8 q_LMnmtt(32,32,49),q_LMnm0(32,32,12)
      integer beta_ind1(20),beta_ind2(20)

      cai=dcmplx(0.d0,1.d0)
      pi=4*datan(1.d0)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc first, get q_LMnm0
      LM0=0
      do ibeta1=1,nbeta(ia)
      do ibeta2=1,ibeta1
      if(lll(ibeta1,ia).eq.lll(ibeta2,ia)) then
      LM0=LM0+1
      beta_ind1(LM0)=ibeta1
      beta_ind2(LM0)=ibeta2
      endif
      enddo
      enddo
      if(LM0.gt.12) then
      write(6,*) "LM0.gt.12, stop"
      stop
      endif
cccccccccccccccccccccccccc

      af4pi=1.d0/dsqrt(4*pi)

      q_LMnm0=0.d0

      kk1=0
      do 301 ibeta1=1,nbeta(ia)
      do 301 m1=1,2*lll(ibeta1,ia)+1
      kk1=kk1+1

      kk2=0
      do 201 ibeta2=1,nbeta(ia) 
      do 201 m2=1,2*lll(ibeta2,ia)+1
      kk2=kk2+1

      if(lll(ibeta1,ia).eq.lll(ibeta2,ia).and.m1.eq.m2) then
       LM00=0
       do LM0=1,12
       if((ibeta1.eq.beta_ind1(LM0).and.
     &           ibeta2.eq.beta_ind2(LM0)).or.
     &     (ibeta1.eq.beta_ind2(LM0).and.
     &           ibeta2.eq.beta_ind1(LM0))) LM00=LM0
        enddo
       if(LM00.eq.0) then
       write(6,*) "something wrong,stop"
       stop
       endif
      q_LMnm0(kk1,kk2,LM00)=qijradt(0,kk1,kk2)*af4pi
      endif
201   continue
301   continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      q_LMnmtmp=dcmplx(0.d0,0.d0)

      kk1=0
      do 300 ibeta1=1,nbeta(ia)
      do 300 m1=-lll(ibeta1,ia),lll(ibeta1,ia)
      kk1=kk1+1
      
      kk2=0
      do 200 ibeta2=1,nbeta(ia) 
      do 200 m2=-lll(ibeta2,ia),lll(ibeta2,ia)
      kk2=kk2+1

      kk3=0
      do Lt=0,6
      do Mt=-Lt,Lt
      kk3=kk3+1
      c_3j_1=clebsch_gordan(2*Lt,-2*Mt,2*lll(ibeta1,ia),-2*m1,
     & 2*lll(ibeta2,ia),-2*m2)/dsqrt(2*lll(ibeta2,ia)+1.d0)*
     & (-1)**(Lt-lll(ibeta1,ia)-m2)
      c_3j_2=clebsch_gordan(2*Lt,0,2*lll(ibeta1,ia),0,
     & 2*lll(ibeta2,ia),0)/dsqrt(2*lll(ibeta2,ia)+1.d0)*
     & (-1)**(Lt-lll(ibeta1,ia))
      c_LMnm=dsqrt((2*lll(ibeta1,ia)+1.d0)*(2*lll(ibeta2,ia)+1.d0)*
     & (2*Lt+1.d0)/(4*pi))*(-1)**Mt*c_3j_1*c_3j_2

      q_LMnmtmp(kk1,kk2,kk3)=qijradt(Lt,kk1,kk2)*
     &  (-1)**m1*c_LMnm


cccc  Note, there is a dconj() on Y_l1,m1, which causes (-1)**m1,and -m1 
cccc  qijradt(Lt,kk1,kk2) only depends on l1,l2 of kk1,kk2, not on m1,m2
cccc  in the clebsch_gordan function, entries are 2*(l,m)
      enddo
      enddo

200   continue
300   continue


cccc  some unitary transformation within kk1,kk2,and kk3
cccc  The C_CG coeff is for Ylm, Yl-m. Change them into 
cccc  Ylm+=[Ylm+(-)^m Yl-m]/sqrt(2) 
cccc  Ylm-=[Ylm-(-)^m Yl-m]/i sqrt(2) 
cccc  Ylm =[Ylm+ + iYlm-]/sqrt(2)
cccc  Yl-m=[Ylm+ - iYlm-] (-)^m/sqrt(2)

cccc  Unitary transformation about M,-M
cccc  current: q_LMnmtmp(:,:,kk3)*(YLM,YL-M)
cccc  change that to: q_LMnmt(:,:,kk3)*(YLM+,YLM-)

      kk3=0
      do Lt=0,6
      do Mt=-Lt,Lt
      kk3=kk3+1
      indec(Mt,Lt)=kk3
      enddo
      enddo

      sqr2I=1.d0/dsqrt(2.d0)

      q_LMnmt=dcmplx(0.d0,0.d0)

      kk3=0
      do Lt=0,6
      do Mt=0,Lt
      kk3=kk3+1
      if(Mt.eq.0) then
      q_LMnmt(:,:,kk3)=q_LMnmtmp(:,:,indec(Mt,Lt))
      else
      q_LMnmt(:,:,kk3)=(q_LMnmtmp(:,:,indec(Mt,Lt))+(-1)**Mt*
     &                q_LMnmtmp(:,:,indec(-Mt,Lt)))*sqr2I
      kk3=kk3+1
      q_LMnmt(:,:,kk3)=cai*(q_LMnmtmp(:,:,indec(Mt,Lt))-(-1)**Mt*
     &                q_LMnmtmp(:,:,indec(-Mt,Lt)))*sqr2I
      endif
      enddo
      enddo

      

ccccccccccccccccccccccccccccccccccccccccccccc
cccc  unitary transform q_LMnmt(kk1,:,:)
cccc  current: q_LMnmt(kk1,:,:)*(Yl1m1,Yl1-m1)
cccc  change that to: q_LMnmtmp(kk1,:,:)*(Yl1m1+,Yl1m1-)
      
      indec=0

      kk1=0
      do ibeta1=1,nbeta(ia)
      do m1=-lll(ibeta1,ia),lll(ibeta1,ia)
      kk1=kk1+1
      indec(m1,ibeta1)=kk1
      enddo
      enddo

      q_LMnmtmp=dcmplx(0.d0,0.d0)

      kk1=0
      do ibeta1=1,nbeta(ia)
      do m1=0,lll(ibeta1,ia)
      kk1=kk1+1
      if(m1.eq.0) then
      q_LMnmtmp(kk1,:,:)=q_LMnmt(indec(m1,ibeta1),:,:)
      else
      q_LMnmtmp(kk1,:,:)=(q_LMnmt(indec(m1,ibeta1),:,:)+
     &    (-1)**m1*q_LMnmt(indec(-m1,ibeta1),:,:))*sqr2I
      kk1=kk1+1
      q_LMnmtmp(kk1,:,:)=cai*(q_LMnmt(indec(m1,ibeta1),:,:)-
     &    (-1)**m1*q_LMnmt(indec(-m1,ibeta1),:,:))*sqr2I
      endif
      enddo
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc  unitary transform q_LMnmt(:,kk2,:)
cccc  current: q_LMnmtmp(:,kk2,:)*(Yl2m2^*,Yl2-m2^*)
cccc  change that to: q_LMnmt(:,kk2,:)*(Yl2m2+^*+,Yl2m2-^*)

      indec=0

      kk2=0
      do ibeta2=1,nbeta(ia)
      do m2=-lll(ibeta2,ia),lll(ibeta2,ia)
      kk2=kk2+1
      indec(m2,ibeta2)=kk2
      enddo
      enddo

      q_LMnmt=dcmplx(0.d0,0.d0)

      kk2=0
      do ibeta2=1,nbeta(ia)
      do m2=0,lll(ibeta2,ia)
      kk2=kk2+1
      if(m2.eq.0) then
      q_LMnmt(:,kk2,:)=q_LMnmtmp(:,indec(m2,ibeta2),:)
      else
      q_LMnmt(:,kk2,:)=(q_LMnmtmp(:,indec(m2,ibeta2),:)+
     &    (-1)**m2*q_LMnmtmp(:,indec(-m2,ibeta2),:))*sqr2I
      kk2=kk2+1
      q_LMnmt(:,kk2,:)=-cai*(q_LMnmtmp(:,indec(m2,ibeta2),:)-    ! extra - sign due to Ylm^*
     &    (-1)**m2*q_LMnmtmp(:,indec(-m2,ibeta2),:))*sqr2I
      endif
      enddo
      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc q_LMnm is the final result, to be used in add_rho_beta
ccccccccccccccccccccccc
cccc Note, by now, the basis for q_LMnm(n,m,LM) are all real 
cccc in terms of LM+,LM-, ln+,ln-,lm+,lm-. Thus, the q_LMnm is real

      q_LMnmtt=dreal(q_LMnmt)

      return
      end subroutine getq_LMnm
      end
      
      
      


