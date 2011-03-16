       subroutine forcNLq(fatom,occ,kpt,nkpt)
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
       include "mpif.h"

       real*8 occ_t(mtype)
       integer iiatom(mtype),icore(mtype),numref(matom)

       real*8 fatom(3,matom)
       real*8 occ(mst,nkpt)
       complex*16 s(9),sx(9),sy(9),sz(9)
       complex*16 cc,cx,cy,cz,cai
*************************************************
       integer isNLa(9,matom)
       common /comisNLa/isNLa
       common /comNL2/occ_t,iiatom,icore,numref


       cai=dcmplx(0.d0,1.d0)
       s=0.d0
       sx=0.d0
       sy=0.d0
       sz=0.d0

       ng_n=ngtotnod(inode,kpt)

       do 40 ia=1,natom
       nref=numref(ia)

       do 40 m=1,mx

	do j=1,nref
	s(j)=dcmplx(0.d0,0.d0)
	sx(j)=dcmplx(0.d0,0.d0)
	sy(j)=dcmplx(0.d0,0.d0)
	sz(j)=dcmplx(0.d0,0.d0)
	enddo

	do i=1,ng_n
	cc=dconjg(ug_n(i,m))
	cx=cai*gkx_n(i,kpt)*cc
	cy=cai*gky_n(i,kpt)*cc
	cz=cai*gkz_n(i,kpt)*cc
          do j=1,nref
          s(j)=s(j)+wqmask(j,i,ia)*cc
          sx(j)=sx(j)+wqmask(j,i,ia)*cx
          sy(j)=sy(j)+wqmask(j,i,ia)*cy
          sz(j)=sz(j)+wqmask(j,i,ia)*cz
	  enddo
	enddo

        call mpi_allreduce(s,s,9,MPI_DOUBLE_COMPLEX,
     & MPI_SUM,MPI_COMM_WORLD,ierr)
        call mpi_allreduce(sx,sx,9,MPI_DOUBLE_COMPLEX,
     & MPI_SUM,MPI_COMM_WORLD,ierr)
        call mpi_allreduce(sy,sy,9,MPI_DOUBLE_COMPLEX,
     & MPI_SUM,MPI_COMM_WORLD,ierr)
        call mpi_allreduce(sz,sz,9,MPI_DOUBLE_COMPLEX,
     & MPI_SUM,MPI_COMM_WORLD,ierr)


	 fx=0.d0
	 fy=0.d0
	 fz=0.d0
	 do j=1,nref
	 fx=fx+dreal(dconjg(s(j))*sx(j))*isNLa(j,ia)
	 fy=fy+dreal(dconjg(s(j))*sy(j))*isNLa(j,ia)
	 fz=fz+dreal(dconjg(s(j))*sz(j))*isNLa(j,ia)
	 enddo

	 fatom(1,ia)=fatom(1,ia)+2*fx*vol**2*occ(m,kpt)
	 fatom(2,ia)=fatom(2,ia)+2*fy*vol**2*occ(m,kpt)
	 fatom(3,ia)=fatom(3,ia)+2*fz*vol**2*occ(m,kpt)
 40     continue
*****************************************
*** The factors are right. 
*** The 2 is from the derivative of (\psi*\psi)^2
*****************************************


       return
       end
