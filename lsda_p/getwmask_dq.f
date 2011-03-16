      subroutine getwmask_dq(xatom,nmap_q,iatom,rcut_q1,
     &  rcut_q2,AL,mrb2_matom_node_q)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************


*****************************************************
cc This subroutine calculates  d[qmask_q(R)]/dR = qmask_dq
*****************************************************

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'mpif.h'
      include 'param.escan_real'

      real*8 xatom(3,matom)
      real*8 AL(3,3),ALIt(3,3),ALntmp(3,3)

      integer nmap_q(matom),iatom(matom)
      real*8 rcut_q1(mtype),rcut_q2(mtype)

      real*8 tmp_dq(49,0:3),tmp_dq0(12,0:3)

      real*8 Plm(-6:6,0:6)
      real*8 gl_tmp(0:6),gl(1001,0:6),rgl(1001)
      real*8 occ_t(mtype)
      real*8 sum_all(49),sum_alltmp(49),sum_all2(12,0:3)
      integer iiatom(mtype),icore(mtype),numref(matom),ityatom(matom)

      real*8 qfuncLM0(12,2000,mtype),r_at(2000,mtype)
      real*8 a_r(mtype),b_r(mtype)


      integer isNLa(9,matom),ipsp_type(mtype)
      real*8  Dij0(32,32,mtype),Qij(32,32,mtype)
      common /comisNLa/isNLa,Dij0,Qij,ipsp_all,ipsp_type
      common /comNL2/occ_t,iiatom,icore,numref,ityatom
      common /com_qfuncLM0/qfuncLM0,r_at,a_r,b_r
*****************************************************
      if(inode_tot.eq.1) then
      open(15,file='graph.j',status='old',action='read',iostat=ierr)
      if(ierr.ne.0) then
      write(6,*) "file **** graph.j **** does not exist, stop"
      call mpi_abort(MPI_COMM_WORLD,ierr)
      endif
      do i=1,1001
      read(15,*) rgl(i),(gl(i,l),l=0,6)
      enddo
      close(15)
      endif
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call mpi_bcast(rgl,1001,MPI_REAL8,
     &   0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gl,1001*7,MPI_REAL8,
     &   0,MPI_COMM_WORLD,ierr)


      call get_ALI(AL,ALIt)
      pi=4*atan(1.d0)
      sqr2=dsqrt(2.d0)
      nr_nL=n1L*n2L*n3L/nnodes

      a1I=sqrt(ALIt(1,1)**2 + ALIt(2,1)**2 + ALIt(3,1)**2)
      a2I=sqrt(ALIt(1,2)**2 + ALIt(2,2)**2 + ALIt(3,2)**2)
      a3I=sqrt(ALIt(1,3)**2 + ALIt(2,3)**2 + ALIt(3,3)**2)

      ALntmp(:,1)=AL(:,1)/n1L
      ALntmp(:,2)=AL(:,2)/n2L
      ALntmp(:,3)=AL(:,3)/n2L


c      ddx=0.02
      ddx=0.00001d0
      ddx_i=1/ddx

      wmask_dq=0.d0

      imap=0
      ico2=0
      do 200 ia=1,natom


      itype=ityatom(ia)

      if(ipsp_type(itype).eq.1) then    ! not ultrasoft potential
      nmap_q(ia)=0
      goto 200
      endif

      rcut_tmpL=rcut_q2(itype)
      rcut_tmp=rcut_q1(itype)

      
      scale0=1.d0/rcut_tmp**3
      scale1=scale0/rcut_tmp
      scale2=scale1/rcut_tmp
      scale3=scale2/rcut_tmp
      scale4=scale3/rcut_tmp
      scale5=scale4/rcut_tmp
      scale6=scale5/rcut_tmp

      di_r=rcut_tmp*n1L*a1I+1.5d0
      dj_r=rcut_tmp*n2L*a2I+1.5d0
      dk_r=rcut_tmp*n3L*a3I+1.5d0
      if(di_r.gt.n1L) di_r=n1L
      if(dj_r.gt.n2L) dj_r=n2L
      if(dk_r.gt.n3L) dk_r=n3L

      x1=xatom(1,ia)*n1L
      x2=xatom(2,ia)*n2L
      x3=xatom(3,ia)*n3L
      if(x1.lt.0.d0) x1=x1+n1L
      if(x2.lt.0.d0) x2=x2+n2L
      if(x3.lt.0.d0) x3=x3+n3L
      if(x1.gt.n1L) x1=x1-n1L
      if(x2.gt.n2L) x2=x2-n2L
      if(x3.gt.n3L) x3=x3-n3L

      
      imap_iat=0
ccccccccccccccccccccccccccccccccccccc
c      do 100 ii=1,nr_nL
c       jj=ii+(inode-1)*nr_nL
c       i=(jj-1)/(n2L*n3L)
c       j=(jj-1-i*n2L*n3L)/n3L
c       k=jj-i*n2L*n3L-j*n3L-1        
ccccc I believe these two ways to do ii is equivalent
cccccc do 100 ii=1,nr_nL => do 100 ico=1,ncolzL/do 101 k=0,n3L-1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do 100 ico = 1,ncolz2L       ! actually ncolzL and ncolz2L are the same
      i  = ixcol_z2L(ico,inode) - 1
         
      dx1=i-x1
      if(dabs(dx1-n1L).lt.dabs(dx1)) dx1=dx1-n1L
      if(dabs(dx1+n1L).lt.dabs(dx1)) dx1=dx1+n1L
      if(dabs(dx1).gt.di_r) goto 100

 
      j  = iycol_z2L(ico,inode) - 1

      dx2=j-x2
      if(dabs(dx2-n2L).lt.dabs(dx2)) dx2=dx2-n2L
      if(dabs(dx2+n2L).lt.dabs(dx2)) dx2=dx2+n2L
      if(dabs(dx2).gt.dj_r) goto 100


      do 101 k=0,n3L-1
      ii=(ico-1)*n3L+k+1
ccccccc the equivalent way to do ii loop

      dx3=k-x3
      if(dabs(dx3-n3L).lt.dabs(dx3)) dx3=dx3-n3L
      if(dabs(dx3+n3L).lt.dabs(dx3)) dx3=dx3+n3L
      if(dabs(dx3).gt.dk_r) goto 101

      xt=ALntmp(1,1)*dx1+ALntmp(1,2)*dx2+ALntmp(1,3)*dx3
      yt=ALntmp(2,1)*dx1+ALntmp(2,2)*dx2+ALntmp(2,3)*dx3
      zt=ALntmp(3,1)*dx1+ALntmp(3,2)*dx2+ALntmp(3,3)*dx3

      rr=xt**2+yt**2+zt**2


      if(rr.ge.rcut_tmpL**2-1.D-5) goto 101

      r=dsqrt(rr)

      imap=imap+1
      imap_iat=imap_iat+1
      if(imap.gt.mrb2_matom_node_q) then
      write(6,*) "imap > mrb2_matom_node_q, stop", inode,imap,
     &   mrb2_matom_node_q
      call mpi_abort(MPI_COMM_WORLD,ierr)
      endif
      indm_q(imap)=ii

ccccccccccccccccccccccccccccc
cccc   rcut_tmp=0.8*rcut_tmpL 
cccc   rcut_tmp is used for qwmask_q(iref,imap)
cccc   rcut_tmpL is used for qwmask_q0(LM0,imap)


      if(r.gt.rcut_tmp-1.D-5) goto 333     ! wmask_q=0 for these imap

      r2=r/rcut_tmp
      ir=1+r2*1000.d0
      if(ir.gt.1000) ir=1000

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 334 ixyz=0,3
      xtn=xt
      ytn=yt
      ztn=zt
      if(ixyz.eq.1) xtn=xt+ddx
      if(ixyz.eq.2) ytn=yt+ddx 
      if(ixyz.eq.3) ztn=zt+ddx
      rrn=xtn**2+ytn**2+ztn**2
      rn=dsqrt(rrn)
      r2n=rn/rcut_tmp

      f1=(rgl(ir+1)-r2n)/(rgl(ir+1)-rgl(ir))
      f2=(r2n-rgl(ir))/(rgl(ir+1)-rgl(ir))

      gl_tmp(0)=(gl(ir,0)*f1+gl(ir+1,0)*f2)*scale0   
      gl_tmp(1)=(gl(ir,1)*f1+gl(ir+1,1)*f2)*scale1
      gl_tmp(2)=(gl(ir,2)*f1+gl(ir+1,2)*f2)*scale2
      gl_tmp(3)=(gl(ir,3)*f1+gl(ir+1,3)*f2)*scale3
      gl_tmp(4)=(gl(ir,4)*f1+gl(ir+1,4)*f2)*scale4
      gl_tmp(5)=(gl(ir,5)*f1+gl(ir+1,5)*f2)*scale5
      gl_tmp(6)=(gl(ir,6)*f1+gl(ir+1,6)*f2)*scale6

      if(rn.lt.1.D-10) then
      dcosth=0.d0
      else
      dcosth=ztn/rn
      endif

      call LegendreSP(Plm,dcosth)

      rtn=dsqrt(xtn**2+ytn**2)
      if(rtn.lt.1.D-10) then
      dsinphi=0.d0
      else
      dsinphi=ytn/rtn
      endif
      phi=asin(dsinphi)
      if(xtn.lt.0.d0) then
      if(ytn.ge.0.d0) then
      phi=pi-phi
      else
      phi=-pi-phi
      endif
      endif

      iref=0
      do l=0,6
      do m=0,l
      iref=iref+1

      if(m.eq.0) then
      tmp_dq(iref,ixyz)=gl_tmp(l)*Plm(m,l)
      else
      tmp_dq(iref,ixyz)=gl_tmp(l)*Plm(m,l)*cos(m*phi)*sqr2   !Ylm+=[Ylm+(-)^m Yl-m]/sqrt(2)
      iref=iref+1
      tmp_dq(iref,ixyz)=gl_tmp(l)*Plm(m,l)*sin(m*phi)*sqr2   !Ylm-=[Ylm-(-)^m Yl-m]/i sqrt(2)
      endif
      enddo
      enddo
334   continue

      do 335 ixyz=0,3
      iref=0
      do l=0,6
      do m=0,l
      iref=iref+1
      if(m.eq.0) then
      wmask_dq(iref,ixyz,imap)=tmp_dq(iref,ixyz)
      else
      wmask_dq(iref,ixyz,imap)=tmp_dq(iref,ixyz)
      iref=iref+1
      wmask_dq(iref,ixyz,imap)=tmp_dq(iref,ixyz)
      endif
      enddo
      enddo
335   continue
      

333   continue
ccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccc  calculate  wmask_q0(LM0,imap)

      ir=log(r/a_r(itype)+1.d0)/b_r(itype)+1.00001d0      ! logrithmic grid

      do ixyz=0,3
      xtn=xt
      ytn=yt
      ztn=zt
      if(ixyz.eq.1) xtn=xt+ddx
      if(ixyz.eq.2) ytn=yt+ddx 
      if(ixyz.eq.3) ztn=zt+ddx
      rrn=xtn**2+ytn**2+ztn**2
      rn=dsqrt(rrn)

      f1=(r_at(ir+1,itype)-rn)/(r_at(ir+1,itype)-r_at(ir,itype))
      f2=1.d0-f1

      do LM0=1,12
       tmp_dq0(LM0,ixyz)=f1*qfuncLM0(LM0,ir,itype)+
     &         f2*qfuncLM0(LM0,ir+1,itype)
      enddo
      enddo

      do ixyz=0,3
      do LM0=1,12
      wmask_dq0(LM0,ixyz,imap)=tmp_dq0(LM0,ixyz)
      enddo
      enddo



101   continue
100   continue
      nmap_q(ia)=imap_iat

200   continue
cccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccc
 
      do ixyz=0,3
      ico=0
      ico2=0
      do ia=1,natom


      sum_all=0.d0
      do i=1,nmap_q(ia)
      ico=ico+1
      do LM0=1,12
      sum_all(LM0)=sum_all(LM0)+wmask_dq0(LM0,ixyz,ico)
      enddo
      enddo
      call mpi_allreduce(sum_all,sum_alltmp,12,MPI_REAL8,
     & MPI_SUM,MPI_COMM_K,ierr)
      sum_all = sum_alltmp
      do LM0=1,12
      sum_all(LM0)=sum_all(LM0)*vol/(n1L*n2L*n3L)
      if(dabs(sum_all(LM0)).gt.1.D-13) then
      sum_all(LM0)=dsqrt(4*pi)/sum_all(LM0)
      endif
      enddo

      do i=1,nmap_q(ia)
      ico2=ico2+1
      do LM0=1,12
      wmask_dq0(LM0,ixyz,ico2)=wmask_dq0(LM0,ixyz,ico2)*sum_all(LM0)      ! sum equals sqrt(4*pi)
      enddo
      enddo

      enddo   ! ia=1,natom
      enddo   ! ixyz
cccccccccccccccccccccccccccccc

      do ixyz=0,3
      ico=0
      ico2=0
      do ia=1,natom            ! maybe should divid this among icolor

      sum_all=0.d0
      do i=1,nmap_q(ia)
      ico=ico+1
      do LM=1,49
      sum_all(LM)=sum_all(LM)+wmask_dq(LM,ixyz,ico)
      enddo
      enddo
      call mpi_allreduce(sum_all,sum_alltmp,49,MPI_REAL8,
     & MPI_SUM,MPI_COMM_K,ierr)
      sum_all = sum_alltmp
      sum=sum_all(1)*vol/(n1L*n2L*n3L)
      sum=dsqrt(4*pi)/sum

      nmap_q_total=nmap_q(ia)
      call mpi_allreduce(nmap_q_total,nmap_q_tmp,1,MPI_INTEGER,
     & MPI_SUM,MPI_COMM_K,ierr)
      nmap_q_total = nmap_q_tmp
      do LM=2,49
      sum_all(LM)=sum_all(LM)/nmap_q_total
      enddo

      do i=1,nmap_q(ia)
      ico2=ico2+1

      wmask_dq(1,ixyz,ico2)=wmask_dq(1,ixyz,ico2)*sum      ! sum equals sqrt(4*pi)

      do LM=2,49
      wmask_dq(LM,ixyz,ico2)=wmask_dq(LM,ixyz,ico2)-sum_all(LM)    ! enforce sum=0.d0
      enddo
      enddo

      enddo   ! atom
      enddo   ! ixyz

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc   calculate the derivative

      ico=0
      do ia=1,natom
      do i=1,nmap_q(ia)
      ico=ico+1
      do ixyz=1,3

      do LM0=1,12
      wmask_dq0(LM0,ixyz,ico)=(wmask_dq0(LM0,ixyz,ico)-
     &         wmask_dq0(LM0,0,ico))*ddx_i
      enddo

      do LM=1,49
      wmask_dq(LM,ixyz,ico)=(wmask_dq(LM,ixyz,ico)-
     &         wmask_dq(LM,0,ico))*ddx_i
      enddo

ccccccccccc  TEST TEST
c      do LM0=1,12
c      wmask_dq0(LM0,ixyz,ico)=wmask_dq(1,ixyz,ico)
c      enddo
cccccccccccccccccccccc

      enddo
      enddo
      enddo
cccccccccccccccccccccccccccccccccccccccc





      return
      end
	 

