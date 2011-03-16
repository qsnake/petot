       subroutine getNLsign(iatom,ntype,is_ref,ip_ref,id_ref)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************


       use fft_data
       use load_data
       use data

       implicit double precision (a-h,o-z)

       include 'param.escan_real'
       integer iatom(matom),iiatom(mtype),icore(mtype),numref(matom)
       integer ityatom(matom)
       integer is_ref(mtype),ip_ref(mtype),id_ref(mtype)
       real*8 occ_t(mtype)
       integer isNL(3,mtype),isNLa(9,matom),ipsp_type(mtype)
       real*8 Dij0(32,32,mtype),Qij(32,32,mtype)

       common /comisNL/isNL
       common /comisNLa/isNLa,Dij0,Qij,ipsp_all,ipsp_type
       common /comNL2/occ_t,iiatom,icore,numref,ityatom
***********************************************

       do ia=1,natom
	do itype=1,ntype
	if(iatom(ia).eq.iiatom(itype)) then
        kk=0
        if(is_ref(itype).eq.1) then
	isNLa(kk+1,ia)=isNL(1,itype)
        kk=kk+1
        endif
        if(ip_ref(itype).eq.1) then
	isNLa(kk+1,ia)=isNL(2,itype)
        isNLa(kk+2,ia)=isNL(2,itype)
        isNLa(kk+3,ia)=isNL(2,itype)
        kk=kk+3
        endif
        if(id_ref(itype).eq.1) then
        isNLa(kk+1,ia)=isNL(3,itype)
	isNLa(kk+2,ia)=isNL(3,itype)
        isNLa(kk+3,ia)=isNL(3,itype)
        isNLa(kk+4,ia)=isNL(3,itype)
        isNLa(kk+5,ia)=isNL(3,itype)
        kk=kk+5
        endif
	goto 11
	endif
	enddo
	write(6,*) "itype not found, stop", iatom(ia),ia
11      continue
      enddo

      return
      end






