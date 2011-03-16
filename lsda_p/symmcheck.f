       subroutine symmcheck(smatr,AL,nrot,xatom,iatom)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
******************************************


       use fft_data
       use load_data
       use data

       implicit double precision (a-h,o-z)
       include 'mpif.h'
       include 'param.escan_real'

       integer smatr(3,3,48)
       real*8 xatom(3,matom)
       integer iatom(matom)
       real*8 AL(3,3)
       real*8 xt(3,6)
       integer iis_done(48)


       if(nrot.le.1) return

       write(6,*) "nrot=",nrot


cccccccccccc first, check the metric of the supercell.
cccccccccccc  check the six edge of a pyramid do not change.
       xt(1,1)=1
       xt(2,1)=0
       xt(3,1)=0

       xt(1,2)=0
       xt(2,2)=1
       xt(3,2)=0

       xt(1,3)=0
       xt(2,3)=0
       xt(3,3)=1

       xt(1,4)=1
       xt(2,4)=-1
       xt(3,4)=0

       xt(1,5)=1
       xt(2,5)=0
       xt(3,5)=-1

       xt(1,6)=0
       xt(2,6)=1
       xt(3,6)=-1


       do irot=1,nrot
       do ii=1,6
	y1=smatr(1,1,irot)*xt(1,ii)+smatr(2,1,irot)*xt(2,ii)+
     &  	smatr(3,1,irot)*xt(3,ii)
	y2=smatr(1,2,irot)*xt(1,ii)+smatr(2,2,irot)*xt(2,ii)+
     &  	smatr(3,2,irot)*xt(3,ii)
	y3=smatr(1,3,irot)*xt(1,ii)+smatr(2,3,irot)*xt(2,ii)+
     &  	smatr(3,3,irot)*xt(3,ii)

	x1=AL(1,1)*xt(1,ii)+AL(1,2)*xt(2,ii)+AL(1,3)*xt(3,ii)
	x2=AL(2,1)*xt(1,ii)+AL(2,2)*xt(2,ii)+AL(2,3)*xt(3,ii)
	x3=AL(3,1)*xt(1,ii)+AL(3,2)*xt(2,ii)+AL(3,3)*xt(3,ii)
	xx1=dsqrt(x1**2+x2**2+x3**2)

	x1=AL(1,1)*y1+AL(1,2)*y2+AL(1,3)*y3
	x2=AL(2,1)*y1+AL(2,2)*y2+AL(2,3)*y3
	x3=AL(3,1)*y1+AL(3,2)*y2+AL(3,3)*y3
	xx2=dsqrt(x1**2+x2**2+x3**2)

	if(dabs(xx1-xx2).gt.0.0001d0) then
	write(6,*) "the sym. rotation does not preserve the metric, stop", 
     &  	irot
	stop
	endif
       enddo
       enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccc second, check the symmetry of the atoms.
cccccccccc It is a N^2 computation.

	do 2000 ia=1,natom

	x1=dmod(xatom(1,ia)+100.d0,1.d0)
	x2=dmod(xatom(2,ia)+100.d0,1.d0)
	x3=dmod(xatom(3,ia)+100.d0,1.d0)

	do 1000 irot=1,nrot

	y1=smatr(1,1,irot)*x1+smatr(2,1,irot)*x2+smatr(3,1,irot)*x3
	y2=smatr(1,2,irot)*x1+smatr(2,2,irot)*x2+smatr(3,2,irot)*x3
	y3=smatr(1,3,irot)*x1+smatr(2,3,irot)*x2+smatr(3,3,irot)*x3

	y1=dmod(y1+100.d0+0.0001d0,1.d0)
	y2=dmod(y2+100.d0+0.0001d0,1.d0)
	y3=dmod(y3+100.d0+0.0001d0,1.d0)

	do ia2=1,natom
	if(iatom(ia).eq.iatom(ia2)) then
	z1=dmod(xatom(1,ia2)+100.d0+0.0001d0,1.d0)
	z2=dmod(xatom(2,ia2)+100.d0+0.0001d0,1.d0)
	z3=dmod(xatom(3,ia2)+100.d0+0.0001d0,1.d0)
	if(dabs(z1-y1)+dabs(z2-y2)+dabs(z3-y3).lt.0.0001d0) goto 1000
	endif
	enddo

ccccccc should never be here
        if(inode.eq.1) then
	write(6,*) "CHECK: equivalent atom not found under symm op,stop",
     &  ia,irot,nrot
        write(6,*) "y1,y2,y3", y1,y2,y3
        endif
	call mpi_abort(MPI_COMM_WORLD,ierr)

1000    continue
2000    continue


       return
       end







