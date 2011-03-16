c ======================================================================
      subroutine input (twopi,mxatms,iprt,a0,g0,a,b,ainp,binp,natoms,
     *   type,ratom,errflg,ninp,nout,akshift,iz_typ,ntyp)
c
c Read lattice vectors
c
c $Id: input.x,v 1.3 89/05/04 21:26:07 jim Exp $
c
c Input:
      double precision twopi
      integer mxatms,ninp,nout
c Output:
      integer iprt,natoms,type(mxatms)
      integer iz_at,iflag,jj,mtype,iz_typ(ntyp)
      double precision a0,g0,a(3,3),b(3,3),ainp(3,3),binp(3,3),
     *   ratom(3,mxatms)
      logical errflg
c Local:
      integer i,j,l
      double precision shifti(3),shiftc(3),shiftl(3),vol,vola,
     *   volin,volina,xinp(3,3),ratmx(3)
      character*72 title
      character*20 xatom_in
      real*8 akshift(3)
      real*8 tmp1,tmp2,tmp3
      real*8 tmp11,tmp22,tmp33
      real*8 ALL(3,3)
c Constants:
      double precision bohr,zero,one
      parameter (bohr = 0.529177249d0)
      parameter (zero = 0.0d0, one = 1.0d0)
      integer p(3)
      data p / 2,3,1 /
c
c Title
c -----
c

c      write (nout,'('' Input summary'',/,1x,13(''-''))')
c      read (ninp,'(a)') title
c      write (nout,'(/,1x,''Title:'',/,1x,a)') title

       read(ninp,*) xatom_in
       write(6,*) xatom_in
ccccc LWW


c
c Lattice data
c ------------
c
c Scale factor (lattice constant) and print flag:
c
      read (ninp,*) a0
      iprt=2
      g0 = twopi/a0
c
c Read basis for input of other vectors, and generate its
c reciprocal-lattice:
c * ainp(j,i) is the jth Cartesian component of the ith vector.
c * binp(j,i) is the jth Cartesian component of the ith vector.
c
      do 10 i = 1,3
         read (ninp,*) (ainp(j,i),j=1,3),akshift(i)
   10 continue
      call recip (one,ainp,binp,volin)
      volina = volin * a0*a0*a0




c
c Real-space primitive basis vectors:
c * xinp(j,i) is the jth component of the ith vector, in the input
c   basis.
c * a(j,i) is the jth Cartesian component of the ith vector.
c


ccccccccc LWW
      open(14,file=xatom_in)
      rewind(14)
      read(14,*) natoms
      do i=1,3
      read(14,*) tmp1,tmp2,tmp3
      ALL(1,i)=tmp1
      ALL(2,i)=tmp2
      ALL(3,i)=tmp3

      xinp(1,i)=(tmp1*binp(1,1)+tmp2*binp(2,1)+tmp3*binp(3,1))/a0
      xinp(2,i)=(tmp1*binp(1,2)+tmp2*binp(2,2)+tmp3*binp(3,2))/a0
      xinp(3,i)=(tmp1*binp(1,3)+tmp2*binp(2,3)+tmp3*binp(3,3))/a0

         call cnvrt (3,0,ainp,xinp(1,i),a(1,i))
      enddo


ccc      do 20 i = 1,3
ccc         read (ninp,*) (xinp(j,i),j=1,3)
ccc         call cnvrt (3,0,ainp,xinp(1,i),a(1,i))
ccc   20 continue
c
c Generate reciprocal lattice:
c * b(j,i) is the jth Cartesian component of the ith vector,
c   in units of g0 at this point.
c
      call recip (one,a,b,vol)
      vola = vol * a0*a0*a0
c
c Print lattice info:
c
      write (nout,'(//,'' Lattice information'',/,1x,19(''-''))')
      write (nout,1000) a0,g0,volin,volina,vol,vola
      write (nout,1005) ((ainp(i,j),i=1,3),j=1,3)
      write (nout,1006) ((binp(i,j),i=1,3),j=1,3)
      write (nout,1007) ((xinp(i,j),i=1,3),j=1,3)
      write (nout,1010) ((a(i,j),i=1,3),j=1,3)
      write (nout,1020) ((b(i,j),i=1,3),j=1,3)
 1000 format (/,1x,'Lattice constant and volume:',
     *   /,4x,'a0 =',f12.7,
     *   /,4x,'g0 =',f12.7,
     *   /,4x,'v (input-basis cell) =',f12.7,' (a0**3) =',f13.7,
     *   /,4x,'v (primitive cell) =',f12.7,' (a0**3) =',f13.7)
 1005 format (/,' Input-basis vectors (in rows; units of a0):',
     *   /,(4x,3f12.7))
 1006 format (/,' Input-basis reciprocal-lattice vectors (in rows; ',
     *   'units of g0):',
     *   /,(4x,3f12.7))
 1007 format (/,' Lattice vectors (in rows; input-basis coordinates):',
     *   /,(4x,3f12.7))
 1010 format (/,' Lattice vectors (in rows; units of a0):',
     *   /,(4x,3f12.7))
 1020 format (/,' Reciprocal-lattice vectors (in rows; units of g0):',
     *   /,(4x,3f12.7))
c
c Atomic data
c -----------
c
c Read origin shift (input basis), natoms:
c
cccc      read (ninp,*) (shifti(i),i=1,3)
cccc      read (ninp,*) natoms
         shifti(1)=0.d0
         shifti(2)=0.d0
         shifti(3)=0.d0
cccccc LWW
c
c Convert shift to cartesian and lattice bases, and print:
c
      call cnvrt (3,0,ainp,shifti,shiftc)
      call cnvrt (3,1,b,shiftc,shiftl)
      write (nout,'(//,'' Atomic data'',/,1x,11(''-''))')
      write (nout,1004) (shifti(i),i=1,3)
      write (nout,1011) (shiftc(i),i=1,3)
      write (nout,1012) (shiftl(i),i=1,3)
 1004 format (/,1x,'Origin shift:',/,6x,'Input basis:',3(1x,f12.7))
 1011 format (8x,'Cartesian:',3(1x,f12.7))
 1012 format (4x,'Lattice basis:',3(1x,f12.7))
c
c Error check:
c
      errflg = .false.
      if (natoms .gt. mxatms) then
         write (nout,3030) natoms,mxatms
         errflg = .true.
         return
      endif
 3030 format (/,2x,'ERROR (input):',/,
     *   5x,'natoms (',i2,') greater than dimension allows (',i2,')')
c
c Read atomic positions (input basis), print, convert to cartesian:
c
      write (nout,1030)
      mtype=0
      do 60 i = 1,natoms
c         read (ninp,*) type(i),(ratmx(j),j=1,3)
cccc LWW
          read (14,*) iz_at,tmp11,tmp22,tmp33
       iflag=0
       do jj=1,mtype
       if(iz_typ(jj).eq.iz_at) then
       type(i)=jj
       iflag=1
       endif
       enddo
       if(iflag.eq.0) then
       mtype=mtype+1
         if(mtype.gt.ntyp) then
         write(6,*) "mtype.gt.ntyp, stop"
         stop
         endif
       type(i)=mtype
       iz_typ(mtype)=iz_at
       endif


      tmp1=ALL(1,1)*tmp11+ALL(1,2)*tmp22+ALL(1,3)*tmp33
      tmp2=ALL(2,1)*tmp11+ALL(2,2)*tmp22+ALL(2,3)*tmp33
      tmp3=ALL(3,1)*tmp11+ALL(3,2)*tmp22+ALL(3,3)*tmp33

      ratmx(1)=(tmp1*binp(1,1)+tmp2*binp(2,1)+tmp3*binp(3,1))/a0
      ratmx(2)=(tmp1*binp(1,2)+tmp2*binp(2,2)+tmp3*binp(3,2))/a0
      ratmx(3)=(tmp1*binp(1,3)+tmp2*binp(2,3)+tmp3*binp(3,3))/a0

         do 70 l = 1,3
            ratmx(l) = ratmx(l) - shifti(l)
   70    continue
         write (nout,1040) i,type(i),(ratmx(j),j=1,3)
         call cnvrt (3,0,ainp,ratmx,ratom(1,i))
   60 continue
      close(14)
 1030 format (/,1x,'Atomic positions (Input basis):')
 1040 format (4x,i2,' (type ',i2,') at ',3f12.7)
c
c Print cartesian atomic positions:
c
      write (nout,1100) (i,type(i),(ratom(j,i),j=1,3),i=1,natoms)
 1100 format (/,1x,'Atomic positions (Cartesian coordinates):',
     *   /,(4x,i2,' (type ',i2,') at ',3f12.7))
c
c Print atomic positions in lattice coordinates and apply
c scale factor to cartesian atomic positions:
c
      write (nout,1050)
      do 200 i = 1,natoms
         call cnvrt (3,1,b,ratom(1,i),ratmx)
         write (nout,1060) i,type(i),(ratmx(j),j=1,3)
         do 80 l = 1,3
            ratom(l,i) = ratom(l,i) * a0
   80    continue
  200 continue
 1050 format (/,1x,'Atomic positions (lattice coordinates):')
 1060 format (4x,i2,' (type ',i2,') at ',3f12.7)
c
c Apply scale factors to lattice vectors:
c
      vol = vola
      do 25 j = 1,3
         do 26 i = 1,3
            a(i,j) = a(i,j) * a0
            b(i,j) = b(i,j) * g0
   26    continue
   25 continue
c
      return
      end
