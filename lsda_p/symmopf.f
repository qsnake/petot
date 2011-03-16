       subroutine symmopf(smatr,nrot,AL,fatom,
     &       xatom,iatom)
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

       integer smatr(3,3,48)
       real*8 smatrC(3,3,48),tmp(3,3)
       real*8 fatom(3,matom),fatom_t(3,matom)
       real*8 acc(matom)
       real*8 xatom(3,matom)
       integer iatom(matom)
       real*8 AL(3,3)
ccc       real*8 ALI(3,3)      ! ALI is defined in param.escan_real
       real*8 xt(3,6)

       if(nrot.le.1) return

cccccccc generate the sym op for Cartesian Coord.

       do irot=1,nrot 
       do k=1,3
       do j=1,3
       s=0.d0
       do i=1,3
       s=s+AL(k,i)*smatr(j,i,irot)
       enddo
       tmp(k,j)=s
       enddo
       enddo

       do k=1,3
       do k2=1,3
       s=0.d0
       do j=1,3
       s=s+tmp(k,j)*ALI(k2,j)
       enddo
       smatrC(k2,k,irot)=s
       enddo
       enddo

       enddo
*******************************************
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        fatom_t=0.d0
	acc=0.d0
	do 2000 ia=1,natom

	x1=dmod(xatom(1,ia)+3.d0,1.d0)
	x2=dmod(xatom(2,ia)+3.d0,1.d0)
	x3=dmod(xatom(3,ia)+3.d0,1.d0)

	do 1000 irot=1,nrot

	y1=smatr(1,1,irot)*x1+smatr(2,1,irot)*x2+smatr(3,1,irot)*x3
	y2=smatr(1,2,irot)*x1+smatr(2,2,irot)*x2+smatr(3,2,irot)*x3
	y3=smatr(1,3,irot)*x1+smatr(2,3,irot)*x2+smatr(3,3,irot)*x3

	y1=dmod(y1+3.d0,1.d0)
	y2=dmod(y2+3.d0,1.d0)
	y3=dmod(y3+3.d0,1.d0)

	f1=smatrC(1,1,irot)*fatom(1,ia)+smatrC(2,1,irot)*fatom(2,ia)
     &        +smatrC(3,1,irot)*fatom(3,ia)
	f2=smatrC(1,2,irot)*fatom(1,ia)+smatrC(2,2,irot)*fatom(2,ia)
     &        +smatrC(3,2,irot)*fatom(3,ia)
	f3=smatrC(1,3,irot)*fatom(1,ia)+smatrC(2,3,irot)*fatom(2,ia)
     &        +smatrC(3,3,irot)*fatom(3,ia)

	do ia2=1,natom
	if(iatom(ia).eq.iatom(ia2)) then
	z1=dmod(xatom(1,ia2)+3.d0,1.d0)
	z2=dmod(xatom(2,ia2)+3.d0,1.d0)
	z3=dmod(xatom(3,ia2)+3.d0,1.d0)
	if(dabs(z1-y1)+dabs(z2-y2)+dabs(z3-y3).lt.0.0001d0) then
	fatom_t(1,ia2)=fatom_t(1,ia2)+f1
	fatom_t(2,ia2)=fatom_t(2,ia2)+f2
	fatom_t(3,ia2)=fatom_t(3,ia2)+f3
	acc(ia2)=acc(ia2)+1.d0
	goto 1000
	endif
	endif
	enddo

ccccccc should never be here
	write(6,*) "equivalent atom not found under symm op,stop", ia,irot
	stop

1000    continue
2000    continue

	do ia=1,natom
	fatom(1,ia)=fatom_t(1,ia)/acc(ia)
	fatom(2,ia)=fatom_t(2,ia)/acc(ia)
	fatom(3,ia)=fatom_t(3,ia)/acc(ia)
	enddo


       return
       end







