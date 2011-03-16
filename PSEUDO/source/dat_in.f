      subroutine dat_in(itype,icorr,ispp,lmax,nr,a,b,r,rab,
     1 nameat,norb,ncore,no,lo,so,zo,znuc,zel,cdc,viod,viou,
     2 vid,viu,ev)
c  
c  ***********************************************************
c  *                                                         *
c  *  The routine reads data from file 'datafile.dat'        *
c  *  for latter use                                         *
c  *                                                         *
c  *  Reverse of datout.f Version dated May 1, 1991          *
c  *                                                         *
c  ***********************************************************
c 
      implicit real*8 (a-h,o-z)
c
      dimension r(nr),rab(nr),no(norb),lo(norb),so(norb),zo(norb),
     1 cdc(nr),viod(lmax,nr),viou(lmax,nr),vid(nr),viu(nr),
     2 ev(norb)
c
      character*1 ispp
      character*2 icorr,nameat
c
c  Open and read out data to current file datafile.dat.
c
      open (unit=7,file='datafile.dat',status='unknown',
     1 form='unformatted')
      rewind(7)
      read(7)itype,icorr,ispp,nr,a,b
      read(7)(r(i),i=1,nr) 
      read(7)(rab(i),i=1,nr)
      read(7)lllmax,nameat,norb,ncore 
      if(lllmax .ne. lmax) stop 'lllmax .ne. lmax'
      read(7)(no(i),i=1,norb)
      read(7)(lo(i),i=1,norb) 
      read(7)(so(i),i=1,norb)
      read(7)(zo(i),i=1,norb)     
      read(7)znuc,zel
      read(7)(cdc(i),i=1,nr)
      do 1,j=1,lmax
        read(7)(viod(j,i),i=1,nr)
 1    continue
      do 2,j=1,lmax
        read(7)(viou(j,i),i=1,nr)
 2    continue
      read(7)(vid(i),i=1,nr)
      read(7)(viu(i),i=1,nr)
      read(7)(ev(i),i=1,norb)
      close (unit=7)
c                         
      return
      end
