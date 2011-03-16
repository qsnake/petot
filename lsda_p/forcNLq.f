       subroutine forcNLq(fatom,occ,kpt,nkpt,iislda,islda,E_st)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************

cccc  Add Qij term also in this program !!


       use fft_data
       use load_data
       use data

       implicit double precision (a-h,o-z)

       include 'param.escan_real'
       include "mpif.h"

       real*8 occ_t(mtype)
       integer iiatom(mtype),icore(mtype),numref(matom),
     &   ityatom(matom)

       real*8 fatom(3,matom)
       real*8 occ(mst,nkpt,islda),E_st(mst,nkpt,islda)
       complex*16 s(32),sx(32),sy(32),sz(32),stmp(32)
       complex*16 ssx(32,32),ssy(32,32),ssz(32,32)
       complex*16 ssx2(32,32),ssy2(32,32),ssz2(32,32)
       complex*16 cc,cx,cy,cz,cai
*************************************************
       integer isNLa(9,matom),ipsp_type(mtype)
       real*8 Dij0(32,32,mtype),Qij(32,32,mtype)
       common /comisNLa/isNLa,Dij0,Qij,ipsp_all,ipsp_type
       common /comNL2/occ_t,iiatom,icore,numref,ityatom


       cai=dcmplx(0.d0,1.d0)
       s=0.d0
       sx=0.d0
       sy=0.d0
       sz=0.d0

       ng_n=ngtotnod(inode,kpt)

       do 50 ia=1,natom
       nref=numref(ia)
       iref_start2=iref_start(ia)
       iitype=ityatom(ia)


         if(ipsp_type(iitype).eq.2) then
       ssx=dcmplx(0.d0,0.d0)
       ssy=dcmplx(0.d0,0.d0)
       ssz=dcmplx(0.d0,0.d0)
       ssx2=dcmplx(0.d0,0.d0)
       ssy2=dcmplx(0.d0,0.d0)
       ssz2=dcmplx(0.d0,0.d0)
         else
       fx=0.d0
       fy=0.d0
       fz=0.d0
         endif

       do 40 m=1,mx

	s=dcmplx(0.d0,0.d0)
	sx=dcmplx(0.d0,0.d0)
	sy=dcmplx(0.d0,0.d0)
	sz=dcmplx(0.d0,0.d0)

        do j=1,nref
        iref=iref_start2+j
          do i=1,ng_n
	  cc=dconjg(ug_n(i,m))
          s(j)=s(j)+wqmask(i,iref)*cc
          sx(j)=sx(j)+wqmask(i,iref)*gkx_n(i,kpt)*cc
          sy(j)=sy(j)+wqmask(i,iref)*gky_n(i,kpt)*cc
          sz(j)=sz(j)+wqmask(i,iref)*gkz_n(i,kpt)*cc
	  enddo
        sx(j)=sx(j)*cai
        sy(j)=sy(j)*cai
        sz(j)=sz(j)*cai
       enddo

*****************************
**** it used: mpi_allreduce(s,s,9,MPI_DOUBLE_COMPLEX,...)
**** That causes overflow, since the s(nref+1) might be used
**** in the previous atoms, and is not set to zero at this atom (in the previous version). 
**** Each time mpi_allreduce is called, it has been multiplied by node. 
**** And it will be done for mx times for each atom !

        call mpi_allreduce(s,stmp,nref,MPI_DOUBLE_COMPLEX,
     & MPI_SUM,MPI_COMM_K,ierr)
        s = stmp
        call mpi_allreduce(sx,stmp,nref,MPI_DOUBLE_COMPLEX,
     & MPI_SUM,MPI_COMM_K,ierr)
        sx = stmp
        call mpi_allreduce(sy,stmp,nref,MPI_DOUBLE_COMPLEX,
     & MPI_SUM,MPI_COMM_K,ierr)
        sy = stmp
        call mpi_allreduce(sz,stmp,nref,MPI_DOUBLE_COMPLEX,
     & MPI_SUM,MPI_COMM_K,ierr)
        sz = stmp


         if(ipsp_type(iitype).eq.2) then
         do j2=1,nref
         do j1=1,nref
         ssx(j1,j2)=ssx(j1,j2)+dconjg(s(j2))*sx(j1)*occ(m,kpt,iislda)
         ssy(j1,j2)=ssy(j1,j2)+dconjg(s(j2))*sy(j1)*occ(m,kpt,iislda)
         ssz(j1,j2)=ssz(j1,j2)+dconjg(s(j2))*sz(j1)*occ(m,kpt,iislda)
         ssx2(j1,j2)=ssx2(j1,j2)+dconjg(s(j2))*sx(j1)*
     &       occ(m,kpt,iislda)*E_st(m,kpt,iislda)
         ssy2(j1,j2)=ssy2(j1,j2)+dconjg(s(j2))*sy(j1)*
     &       occ(m,kpt,iislda)*E_st(m,kpt,iislda)
         ssz2(j1,j2)=ssz2(j1,j2)+dconjg(s(j2))*sz(j1)*
     &       occ(m,kpt,iislda)*E_st(m,kpt,iislda)
         enddo
         enddo

         else

         do j1=1,nref
         fx=fx+Dij(j1,j1,ia,iislda)*dconjg(s(j1))*
     &                         sx(j1)*occ(m,kpt,iislda)
         fy=fy+Dij(j1,j1,ia,iislda)*dconjg(s(j1))*
     &                         sy(j1)*occ(m,kpt,iislda)
         fz=fz+Dij(j1,j1,ia,iislda)*dconjg(s(j1))*
     &                         sz(j1)*occ(m,kpt,iislda)
         enddo
         endif

40       continue
         

         if(ipsp_type(iitype).eq.2) then
	 fx=0.d0
	 fy=0.d0
	 fz=0.d0
	 do j1=1,nref
         do j2=1,nref
         fx=fx+Dij(j1,j2,ia,iislda)*ssx(j1,j2)
     &        -Qij(j1,j2,iitype)*ssx2(j1,j2)
         fy=fy+Dij(j1,j2,ia,iislda)*ssy(j1,j2)
     &        -Qij(j1,j2,iitype)*ssy2(j1,j2)
         fz=fz+Dij(j1,j2,ia,iislda)*ssz(j1,j2)
     &        -Qij(j1,j2,iitype)*ssz2(j1,j2)
	 enddo
         enddo
         endif


	 fatom(1,ia)=fatom(1,ia)+2*fx*vol**2
	 fatom(2,ia)=fatom(2,ia)+2*fy*vol**2
	 fatom(3,ia)=fatom(3,ia)+2*fz*vol**2
50       continue
*****************************************
*** The factors are right. 
*** The 2 is from the derivative of (\psi*\psi)^2
*****************************************


       return
       end
