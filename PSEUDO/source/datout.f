      subroutine datout(itype,icorr,ispp,lmax,nr,a,b,r,rab,
     1 nameat,norb,ncore,no,lo,so,zo,znuc,zel,cdc,viod,viou,
     2 vid,viu,ev)
c  
c  ***********************************************************
c  *                                                         *
c  *  The routine writes needed data to file 'datafile.dat'  *
c  *  for latter use in a minimization program.              *
c  *  Users may want to remove or modify this routine        *
c  *  depending on their needs.                              *
c  *                                                         *
c  *  Version dated May 1, 1991                              *
c  *  njtj                                                   *
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
c  Open and write out data to current file datafile.dat.
c
      open (unit=7,file='datafile.dat',status='unknown',
     1 form='unformatted')
      rewind(7)
      write(7)itype,icorr,ispp,nr,a,b
      write(7)(r(i),i=1,nr) 
      write(7)(rab(i),i=1,nr)
      write(7)lmax,nameat,norb,ncore 
      write(7)(no(i),i=1,norb)
      write(7)(lo(i),i=1,norb) 
      write(7)(so(i),i=1,norb)
      write(7)(zo(i),i=1,norb)     
      write(7)znuc,zel
      write(7)(cdc(i),i=1,nr)
      do 1,j=1,lmax
        write(7)(viod(j,i),i=1,nr)
 1    continue
      do 2,j=1,lmax
        write(7)(viou(j,i),i=1,nr)
 2    continue
      write(7)(vid(i),i=1,nr)
      write(7)(viu(i),i=1,nr)
      write(7)(ev(i),i=1,norb)
      close (unit=7)
c                         
      return
      end
