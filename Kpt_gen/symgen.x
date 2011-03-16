cccc This program has been modified by L.W. Wang, so it only output
cccc operations which is sommorphic (which is required for EtotMV).
c ======================================================================
      subroutine symgen (natmty,natmi,a,b,atmpos,mtrans,
     *   iprt,twopi,mtran0,mtrx,tnp,icof,namegr,nout)
c
c Determine the symmetry operations of the crystal, and interface to
c the XPW program.  This takes the place of the input of generators
c for the group, but still allows the input of supergroup generators
c in addition to the operations generated here.
c
c @(#) symgen -- 1.2 -- 3/5/89 10:11:17
c
c Parameters:
#ifdef PARAM_H
#include PARAM_H
#else
      integer ntyp,natp
      parameter ( ntyp = ntyp )
      parameter ( natp = natp )
#endif
c Input:
      integer natmty,natmi(natmty),iprt,nout
      double precision a(3,3),b(3,3),atmpos(3,natp),twopi
      character*(*) namegr
c Output:
      integer mtran0,mtrx(3,3,48),icof(48)
      double precision tnp(3,48)
c Input and Output:
      integer mtrans
c Local:
      double precision ratom(3,natp)
      integer natoms,type(natp),op(48),att(48,natp),mt(48,48),hg,nops,
     *   i,j,k,l,class(48),inv(48),first,last
      double precision rc(3,3,48),ai(3,3),tp(3,48),tc(3,48),sum,
     *   r1(3,3),rx(3,natp)
      logical nptflg(48),nsmmrf,hasinv,hex,errflg
c Constants:
      double precision zero,one,eps
      parameter (zero = 0.0d0, one = 1.0d0)
      parameter (eps = 1.0d-4)
c
c Input:
c * natmty is the number of types of atoms.
c * natmi(i) is the number of atoms of type i.
c * iprt is the print flag (0 -> min, 1 -> most, 2 -> all).
c * nout is the unit to use for printing.
c * a(i,j) is the i-th cartesian component of the j-th primitive
c   translation vector of the real-space lattice.
c * b(i,j) is the i-th cartesian component of the j-th primitive
c   translation vector of the reciprocal-space lattice.
c * atmpos(j,i) is the j-th lattice-basis component of the position
c   vector for the i-th atom in the unit cell.
c * twopi is obvious.
c * namegr is some identifier for the group.
c * abs(mtrans) is the number of operations given explicitly in the
c   input file.  If mtrans > 0, symgen does not determine any group
c   operations.  If mtrans <=0, symgen determines the symmetry
c   operations of the crystal.  If mtrans < 0, abs(mtrans) additional
c   operations are to be given explicitly in the input file.
c
c Output:
c * mtran0 returns the number of symmetry operations found here.
c * mtrans returns the total number of operations (mtran0 +
c   abs(mtrans)), or 48, whichever is smaller.
c * mtrx(.,.,i) is the matrix for operation i (in reciprocal-space
c   lattice basis).
c * tnp(.,i) is the fractional translation of operation i (in real-
c   space lattice basis).
c * icof(i) = 1, i=1,mtran0, identifies the returned operations as
c   a (possibly improper) subgroup of the total set of operations.
c
c Local:
c * type(i) is an integer distinguishing atoms of different type.
c * ratom(j,i) is the j-th cartesion component of the position vector
c   for the i-th atom in the unit cell.
c * nops is the order of the group.
c * rc(.,.,i) is the cartesian representation of the matrix of
c   operation i.
c * tc(.,i) is the cartesian representation of the primitive
c   translation of operation i.
c * tp(.,i) is the fractional translation of operation i (in real-
c   space lattice basis).
c * op(i) is an index pointing to operation i in a master list.
c * att(op(i),j) is the index of the atom to which operation op(i)
c   maps atom j.
c * mt(op(i),op(j)) is the index of the operation which is the product
c   of operations op(i) and op(j) (in that order).
c * class(i) is the conjugate class to which operation i belongs.
c * inv(i) is the index of the operation which is the inverse of i.
c * ai is the matrix inverse of a.
c * hg is the holohedral group number.
c * nptflg(i) indicates whether there is a non-primitive translation
c   part of operation i.
c * hex is .true. is the group is hexagonal (hence having non-integral
c   matrix elements).
c * hasinv is .true. if the group contains the inversion.
c * nsmmrf is .true. if the group is non-symmorphic.
c * r1, and rx are scratch arrays.
c
c
c Interface to input from XPW:
c
#ifdef DEBUG
c
c Print out input to check:
c
      write (nout,8000)
      write (nout,8010) ((a(i,j),i=1,3),j=1,3)
      write (nout,8020) ((b(i,j),i=1,3),j=1,3)
      write (nout,8030)
      first = 1
      do 8100 i = 1,natmty
         last = first + natmi(i) - 1
         write (nout,8040) natmi(i),i,((atmpos(j,k),j=1,3),k=first,last)
         first = last + 1
 8100 continue
 8000 format (/,' Summary of input to symgen:',/)
 8010 format (' a matrix:',/,(3x,3(1x,f9.5)))
 8020 format (' b matrix:',/,(3x,3(1x,f9.5)))
 8030 format (' Atomic positions:')
 8040 format (4x,i2,' atoms of type ',i2,' at',/,(3x,3(1x,f9.5)))
#endif
c
c Return silently if mtrans > 0 (user wants to specify all of
c the operations explicitly):
c
      if (mtrans .gt. 0) then
         mtran0 = 0
         return
      endif
c
c Print header:
c
      write (nout,1300)
 1300 format (/,' symmetry information',/,1x,20('*'))
c
c Construct the inverse of a (and check):
c
      do 10 i = 1,3
         do 20 j = 1,3
            ai(j,i) = b(i,j) / twopi
   20    continue
   10 continue
      sum = zero
      do 11 j = 1,3
         do 12 i = 1,3
            if (i .eq. j) then
               sum = sum + abs(one -
     *            (ai(j,1)*a(1,i) + ai(j,2)*a(2,i) + ai(j,3)*a(3,i)))
            else
               sum = sum +
     *            abs(ai(j,1)*a(1,i) + ai(j,2)*a(2,i) + ai(j,3)*a(3,i))
            endif
   12    continue
   11 continue
      if (sum .gt. eps) then
         write (nout,1001)
         write (nout,160) (i,(ai(j,i),j=1,3),(a(j,i),j=1,3),i=1,3)
         stop
      endif
 1001 format (/,2x,'ERROR (symgen):',/,
     *   5x,'Trans(b)/2*pi is not the inverse of a')
  160 format (/,21x,'ainv',33x,'a',
     *   /,(2x,i1,3x,3(1x,f9.4),5x,3(1x,f9.4)))
c
c Convert atomic positions to cartesian coordinates, set
c up type array, and determine total number of atoms:
c
      first = 1
      do 1 i = 1,natmty
         do 2 j = first, first+natmi(i)-1
            type(j) = i
            do 3 k = 1,3
               ratom(k,j) = (a(k,1)*atmpos(1,j) + a(k,2)*atmpos(2,j)
     *            + a(k,3)*atmpos(3,j)) / twopi
    3       continue
    2    continue
         first = first + natmi(i)
    1 continue
      natoms = first - 1
c Error check:
      if (natoms .gt. natp) then
         write (nout,1000) natoms,natp
         stop
      endif
 1000 format (/,2x,'ERROR (symgen):',/,
     *   5x,'natoms (',i2,') greater than dimension allows (',i2,')')
c
c
c Calculate the symmetry operations:
c
c Determine the point group of the lattice and the crystal system:
c
      call pgl (a,ai,rc,nops,hg,op,hex)
c
c Determine the point group of the crystal, the fractional translations,
c and other useful information about the group:
c
      call atftmt (natoms,type,ratom,a,ai,rc,iprt,nops,op,tp,tc,
     *   att,mt,class,inv,nptflg,nsmmrf,hasinv,hex,errflg,rx,nout)
c
c Print standard info (omit cartesian matrices if iprt < 2):
c
      write(nout,1500) nops
 1500 format(1x,' nops=',i5)
      do 1505 i=1,nops
      write(nout,1502) i,op(i)
 1505 continue
 1502 format(1x,'op=',i2,'is the number',i2)
c
      if (iprt .gt. 0) then
         call symprt (natoms,type,nops,hg,op,hasinv,nsmmrf,rc,tc,
     *      mt,att,class,inv,nptflg,hex,iprt-1,nout)
      endif
c
c
c Interface to output for XPW:
c
c Copy nops to mtran0, set icof(i) = 1 (0 < i <= mtran0), convert
c rotation matrices to reciprocal-space lattice coordinates and
c copy to mtrx, and copy non-primitive translations to tnp
c (multiplying them by 2*pi; also check that mtrans + mtran0 <= 48):
c

      mtran0 = nops
      mtrans = abs(mtrans) + mtran0
      if (mtrans .gt. 48) then
c This will cause the input to blow up, but the alternative is
c to risk scribbling on other variables.
         mtrans = 48
      endif
      do 2001 i = 1,nops
         l = op(i)
         icof(i) = 1
         do 2010 j = 1,3
            tnp(j,i) = tp(j,l) * twopi
#ifdef CYBER
cdir$  novector
#endif
            do 2011 k = 1,3
               r1(j,k) = rc(j,1,l)*ai(k,1) + rc(j,2,l)*ai(k,2)
     *            + rc(j,3,l)*ai(k,3)
 2011       continue
#ifdef CYBER
cdir$  vector
#endif
 2010    continue
         do 2020 j = 1,3
            do 2021 k = 1,3
               mtrx(j,k,i) =
     *            nint(a(1,j)*r1(1,k) + a(2,j)*r1(2,k) + a(3,j)*r1(3,k))
 2021       continue
 2020    continue
 2001 continue
c
c Print number of ops generated if detailed output not requested:
c
      if (iprt .lt. 1) then
         write (nout,5600) nops
         return
      endif
 5600 format (/,i8,' symmetry elements generated')
c
c Print rotation matrices and fractional translations in
c xpw input format:
c
      if (iprt .lt. 2) return
      write (nout,440) mtrans,namegr
c      write (nout,698)

      open(18,file='symm.file')
      rewind(18)
      write(18,*) mtrans

      

      do 480 i = 1,mtrans
       write(18,221)  i
221    format(3x, i2,"-----------------")
       do k=1,3
       write(18,222) (mtrx(j,k,i),j=1,3)
       enddo
222    format(5x,3(i5,5x))

         l = op(i)
         if (nptflg(i)) then
            write (nout,460)
     *         ((mtrx(j,k,i),k=1,3),j=1,3),i,(tp(j,l),j=1,3)
C-----------PP input format------
c	       write (nout,601) i
c	       write (nout,602) (mtrx(1,k,l),k=1,3)	
c	       write (nout,602) (mtrx(2,k,l),k=1,3)	
c	       write (nout,602) (mtrx(3,k,l),k=1,3)	
         else
             write (nout,470) ((mtrx(j,k,i),k=1,3),j=1,3),i
C-----------PP input format------
c	       write (nout,601) i
c	       write (nout,602) (mtrx(1,k,l),k=1,3)	
c	       write (nout,602) (mtrx(2,k,l),k=1,3)	
c	       write (nout,602) (mtrx(3,k,l),k=1,3)	
         endif
  480 continue
       write(18,*) "-------------------------------"
       close(18)
  440 format (/,/,' Symmetry operations in xpw input format',
     *   /,/,i5,'  b=2  p=0  s=0',a)
  460 format (3(1x,3i3),1x,'x ',i5,/,3(1x,f9.7))
  470 format (3(1x,3i3),1x,'  ',i5)
c
c  698 format (/,/,1x,'Matrices in PP INPUT FORMAT',/,1x,27('-'),/)
c  601 format (1x,i2,'---------------------------')
c  602 format (3(6x,i3))
c
      return
      end
c ======================================================================
      subroutine atftmt (natoms,type,ratom,a,ai,rc,iprt,nops,op,
     *   tp,tc,att,mt,class,inv,nptflg,nsmmrf,hasinv,hex,errflg,rx,nerr)
c
c Determine the point group of a crystal, the atom transformation
c table, att, the fractional translations, tp, associated with each
c rotation, and the multiplication table, mt, for the point group
c of the crystal.  op now contains operations in the p.g. of the
c crystal, and nops is the order of this group.
c
c @(#) atftmt -- 1.3 -- 3/3/89 17:22:53
c
c Input:
      integer natoms,type(natoms),iprt,nerr
      double precision ratom(3,natoms),a(3,3),ai(3,3),rc(3,3,48)
      logical hex
c Output:
      integer att(48,natoms),mt(48,48),class(48),inv(48)
      double precision tp(3,48),tc(3,48)
      logical nptflg(48),nsmmrf,hasinv,errflg
c Input and Output:
      integer nops,op(48)
c Local:
      integer nopsx,n,l,k,i,j,k1,k2,ks,k3,k4,n1,n2,m,n3,
     *   il,invers
      double precision vr(3),vt(3),xb(3),r1(3,3),dif,da,sum
      logical erri,errc
c Storage:
      double precision rx(3,natoms)
c Numeric constants:
      double precision zero,one,eps1,eps2
      parameter (zero = 0.0d0, one = 1.0d0)
      parameter (eps1 = 1.0d-4, eps2 = 1.0d-3)
c
c Initialization (invers indexes the inversion):
c
      nopsx = 0
      hasinv = .false.
      errflg = .false.
      if (hex) then
         invers = 13
      else
         invers = 25
      endif
c
c Select the crystal symmetry operations, and construct the atom
c transformation table:
c
      do 120 n = 1,nops
         l = op(n)
c Apply rotation to all atoms in the basis:
         do 30 k = 1,natoms
            do 31 i = 1,3
               rx(i,k) = zero
               do 20 j = 1,3
                  rx(i,k) = rx(i,k) + rc(i,j,l)*ratom(j,k)
   20          continue
   31       continue
   30    continue
c Check to see whether the operation is a symmetry op.:
         do 90 k1 = 1,natoms
            do 91 k2 = 1,natoms
               if (type(k1) .ne. type(k2)) goto 91
               do 40 i = 1,3
                  xb(i) = rx(i,k1) - ratom(i,k2)
   40          continue
               il = 0
c rlv removes a direct-lattice vector from xb, leaving the
c remainder in vr.  If a non-zero lattice vector was removed, il is
c made non-zero.  vr stands for v-reference.
               call rlv (ai,xb,vr,il)
               ks = 0
               do 80 k3 = 1,natoms
                  do 70 k4 = 1,natoms
                     if (type(k3) .ne. type(k4)) goto 70
                     do 50 i = 1,3
                        xb(i) = rx(i,k3) - ratom(i,k4)
   50                continue
                     call rlv (ai,xb,vt,il)
c vt stands for v-test.
                     dif = zero
                     do 60 i = 1,3
                        da = abs(vr(i)-vt(i)) + eps1
                        dif = dif + mod(da,one)
   60                continue
                     if (dif .gt. eps2) goto 70
                     att(l,k3) = k4
c att is the function f0 defined in Maradudin and Vosko,
c Rev. Mod. Phys. 40, 1 (1968), by eq.(2.35).  It defines
c the atom transformation table.
                     ks = ks + k4
                     if (ks .eq. natoms*(natoms+1)/2) goto 100
                     goto 80
   70             continue
                  goto 91
   80          continue
   91       continue
   90    continue
         goto 120
  100    nopsx = nopsx+1
         do 110 i = 1,3
            tp(i,l) = vr(i)
  110    continue
c tp(i,l) is the i-th primitive component of the fractional
c translation associated with the rotation rc(l).
         op(nopsx) = l
         if (l .eq. invers) hasinv = .true.
  120 continue
      nops = nopsx
c
c Check to see whether the group is symmorphic:
c

      nsmmrf = .false.

      do 170 n = 1,nops
         l = op(n)
         sum = abs(tp(1,l)) + abs(tp(2,l)) + abs(tp(3,l))
         nptflg(n) = sum .gt. eps1
         nsmmrf = nsmmrf .or. nptflg(n)
  170 continue

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc Modified by L.W. Wang, May, 5, 1999. 
ccccc We need only symmorphic group for our EtotMV program
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      mm=0
      do n=1,nops
      if(.NOT.nptflg(n)) then
      mm=mm+1
      nptflg(mm)= .false.
      op(mm)=op(n)
      endif
      enddo

      nops=mm
      nsmmrf = .false.

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc end of the modification by L.W. Wang
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



c
c Construct the multiplication table:
c
      do 260 n1 = 1,nops
         do 261 n2 = 1,nops
            l = op(n1)
            m = op(n2)
            do 230 i = 1,3
               do 231 j = 1,3
                  r1(i,j) = zero
                  do 220 k = 1,3
                     r1(i,j) = r1(i,j) + rc(i,k,l)*rc(k,j,m)
  220             continue
  231          continue
  230       continue
            do 250 n3 = 1,nops
               n = op(n3)
               sum = zero
               do 240 i = 1,3
                  do 241 j = 1,3
                     sum = sum + abs(r1(i,j)-rc(i,j,n))
  241             continue
  240          continue
               if (sum .le. eps2) then
                  mt(l,m) = n
                  goto 261
               endif
  250       continue
  261    continue
  260 continue
c
c Check multiplication table, and determine class and
c inverse for each element:
c
      call mtanal (nops,mt,op,class,inv,erri,errc)
      if (errc) then
         write (nerr,3010)
         errflg = .true.
      endif
      if (erri) then
         write (nerr,3020)
         errflg = .true.
      endif
 3010 format (' ERROR (atftmt):  Group multiplication table ',
     *   'fails closure test.')
 3020 format (' ERROR (atftmt):  Zero or multiple inverses ',
     *   'found for some group element(s).')
c
c Convert fractional translations to Cartesian coordinates.
c tc(j,l) is in Cartesian coordinates, in the same units as
c the given lattice vectors and atomic positions.
c
      do 1000 i = 1,nops
         l = op(i)
         do 1002 j = 1,3
            tc(j,l) = tp(1,l)*a(j,1) + tp(2,l)*a(j,2) + tp(3,l)*a(j,3)
 1002    continue
 1000 continue
c
      return
      end
c ======================================================================
      subroutine mtanal (nops,mt,op,class,inv,erri,errc)
c
c Determine the inverses and the conjugate classes of the
c operations of a group, given the multiplication table.
c Also, checks the structure of the table for closure and
c the presence of a unique inverse for each operation.
c Only integer arithmetic is needed.
c
c @(#) mtanal -- 1.1 -- 3/5/89 09:45:01
c
c Input:
      integer nops,mt(48,nops),op(nops)
c Output:
      integer class(48),inv(48)
      logical erri,errc
c Local:
      integer curcls,a,u,b,ua,aop,bop,uop,bu
c Notes:
c * nops is the number of operations in the group.
c * op(a) is an index pointing to operation a in a master list.
c * mt(op(a),op(b)) = op(c), where c = a * b.
c * class(a) is an integer identifying the class of operation a.
c * inv(a) is the operation b, such that b * a = 1
c   (or mt(op(b),op(a)) = op(1) := 1).
c * erri returns true if an element does not have an inverse
c   or appears to have more than one.
c * errc returns true if a product is not in the group.
c
c Initialization:
      do 5 a = 1,48
         class(a) = 0
         inv(a) = 0
    5 continue
      erri = .false.
      errc = .false.
      inv(1) = 1
      class(1) = 1
      curcls = 1
c Determine the class structure and inverses:
      do 10 a = 2,nops
         aop = op(a)
         do 20 u = 2,nops
            uop = op(u)
            ua = mt(uop,aop)
            if (ua .eq. 1) then
               if (inv(a) .eq. 0) then
                  inv(a) = u
               else
                  erri = .true.
               endif
            endif
            if (u .ne. a) then
               do 30 b = 2,nops
                  if ((b .ne. a) .and. (b .ne. u)) then
                     bop = op(b)
                     bu = mt(bop,uop)
                     if (bu .eq. ua) then
                        if (class(a) .ne. 0) then
                           class(b) = class(a)
                        elseif (class(b) .ne. 0) then
                           class(a) = class(b)
                        else
                           curcls = curcls + 1
                           class(a) = curcls
                           class(b) = curcls
                        endif
                     endif
                  endif
   30          continue
            endif
   20    continue
         if (inv(a) .eq. 0) then
            erri = .true.
         endif
         if (class(a) .eq. 0) then
            curcls = curcls + 1
            class(a) = curcls
         endif
   10 continue
c
      return
      end
c ======================================================================
      subroutine pgl (a,ai,r,nc,hg,ib,hex)
c
c Determine the point group of the lattice and the crystal system.
c The array ib indexes the locations of the group operations, and
c nc is the order of the group.  hg indicates the crystal system
c (holohedral group number).
c
c @(#) pgl -- 1.1 -- 3/5/89 09:46:42
c
c Input:
      double precision a(3,3),ai(3,3)
c Output:
      double precision r(3,3,48)
      integer nc,hg,ib(48)
      logical hex
c Local:
      double precision vr(3),xa(3),tr
      integer nr,i,j,k,n,lx
c Symbolic constants:
      integer tricl,monocl,rhomb,tetrag,cubic,trigrh,trighx,hexag
      parameter (tricl  = 1)
      parameter (monocl = 2)
      parameter (rhomb  = 3)
      parameter (tetrag = 4)
      parameter (cubic  = 5)
      parameter (trigrh = 6)
      parameter (trighx = 7)
      parameter (hexag  = 8)
c Numeric constants:
      double precision zero,eps
      parameter (zero = 0.0d0)
      parameter (eps = 1.0d-4)
c
c Determine point group of the lattice by applying possible operations
c to each primitive lattice vector and checking to see whether the
c result is a lattice vector:
c
      hex = .true.
      nr = 24
   10 nc = 0
      call rot (nr,r)
      do 50 n = 1,nr
         ib(n) = 0
         tr = zero
         do 40 k = 1,3
            do 20 i = 1,3
               xa(i) = zero
               do 21 j = 1,3
                  xa(i) = xa(i) + r(i,j,n)*a(j,k)
   21          continue
   20       continue
            call rlv (ai,xa,vr,lx)
            do 30 i = 1,3
               tr = tr + abs(vr(i))
   30       continue
            if (tr .gt. eps) goto 50
   40    continue
         nc = nc+1
         ib(nc) = n
   50 continue
      if (hex) then
         if (nc .gt. 12) then
            hg = hexag
         elseif (nc .eq. 12) then
            hg = trighx
         else
             nr = 48
             hex = .false.
             go to 10
         endif
      else
         if (nc .gt. 16) then
            hg = cubic
         elseif (nc .eq. 16) then
            hg = tetrag
         elseif (nc .eq. 12) then
            hg = trigrh
         elseif (nc .eq. 8) then
            hg = rhomb
         elseif (nc .eq. 4) then
            hg = monocl
         elseif (nc .lt. 4) then
            hg = tricl
         endif
      endif
c
      return
      end
c ======================================================================
      subroutine rlv (p,g,y,l)
c
c Remove a lattice vector from g, if possible.  The remainder is
c returned in y in lattice coordinates, and l is non-zero iff a
c vector is removed.
c
c @(#) rlv -- 1.1 -- 3/5/89 09:47:27
c
c Input:
      double precision p(3,3),g(3)
c Output:
      double precision y(3)
      integer l
c Local:
      double precision del(3),ts
      integer i,j
c Constants:
      double precision zero,one,eps,p9
      parameter (zero = 0.0d0, one = 1.0d0)
      parameter (eps = 1.0d-4)
      parameter (p9 = 0.9d0)
c
      l = 0
      ts = zero
      do 10 i = 1,3
         del(i) = zero
         y(i) = zero
         ts = ts + abs(g(i))
   10 continue
      if (ts .le. eps) return
      do 30 i = 1,3
         do 20 j = 1,3
            y(i) = y(i) + p(i,j)*g(j)
   20    continue
         if (y(i) .gt. p9) then
            del(i) = eps
         else if (y(i) .lt. -p9) then
            del(i) = -eps
         endif
c del is added to eliminate roundoff errors in the mod function
         y(i) = y(i) + del(i)
         l = l + abs(y(i))
         y(i) = -mod(y(i),one)
         y(i) = y(i) + del(i)
   30 continue
c
      return
      end
c ======================================================================
      subroutine rot (nr,r)
c
c Set up rotation matrices r for hexagonal group if nr .le. 24, or
c for cubic group if nr .gt. 24.
c
c @(#) rot -- 1.1 -- 3/5/89 09:48:13
c
c Input:
      integer nr
c Output:
      double precision r(3,3,48)
c Local:
      integer i,j,k,n,nv
c Constants:
      double precision zero,half,one,three,rt3o2
      parameter (zero = 0.0d0, half = 0.5d0, one = 1.0d0, three = 3.0d0)
      rt3o2 = sqrt(three) * half
c
c Initialization:
c
      do 10 n = 1,nr
         do 11 i = 1,3
            do 12 j = 1,3
               r(i,j,n) = zero
   12       continue
   11    continue
   10 continue
c
      if (nr .le. 24) then
c
c Define the generators for the rotation matrices -- hexagonal group:
c
         r(1,1,2) = half
         r(1,2,2) = -rt3o2
         r(2,1,2) = rt3o2
         r(2,2,2) = half
         r(1,1,7) = -half
         r(1,2,7) = -rt3o2
         r(2,1,7) = -rt3o2
         r(2,2,7) = half
         do 30 n = 1,6
            r(3,3,n) = one
            r(3,3,n+18) = one
            r(3,3,n+6) = -one
            r(3,3,n+12) = -one
   30    continue
c
c Generate the rest of the rotation matrices -- hexagonal group:
c
         do 40 i = 1,2
            r(i,i,1) = one
            do 41 j = 1,2
               r(i,j,6) = r(j,i,2)
               do 42 k = 1,2
                  r(i,j,3) = r(i,j,3) + r(i,k,2)*r(k,j,2)
                  r(i,j,8) = r(i,j,8) + r(i,k,2)*r(k,j,7)
                  r(i,j,12) = r(i,j,12) + r(i,k,7)*r(k,j,2)
   42          continue
   41       continue
   40    continue
         do 50 i = 1,2
            do 51 j = 1,2
               r(i,j,5) = r(j,i,3)
               do 52 k = 1,2
                  r(i,j,4) = r(i,j,4) + r(i,k,2)*r(k,j,3)
                  r(i,j,9) = r(i,j,9) + r(i,k,2)*r(k,j,8)
                  r(i,j,10) = r(i,j,10) + r(i,k,12)*r(k,j,3)
                  r(i,j,11) = r(i,j,11) + r(i,k,12)*r(k,j,2)
   52          continue
   51       continue
   50    continue
         do 60 n = 1,12
            nv = n+12
            do 61 i = 1,2
               do 62 j = 1,2
                  r(i,j,nv) = -r(i,j,n)
   62          continue
   61       continue
   60    continue
      else
c
c Define the generators for the rotation matrices -- cubic group:
c
         r(1,3,9) = one
         r(2,1,9) = one
         r(3,2,9) = one
         r(1,1,19) = one
         r(2,3,19) = -one
         r(3,2,19) = one
c
c Generate the rest of the rotation matrices -- cubic group:
c
         do 80 i = 1,3
            r(i,i,1) = one
            do 81 j = 1,3
               r(i,j,20) = r(j,i,19)
               r(i,j,5) = r(j,i,9)
               do 82 k = 1,3
                  r(i,j,2) = r(i,j,2) + r(i,k,19)*r(k,j,19)
                  r(i,j,16) = r(i,j,16) + r(i,k,9)*r(k,j,19)
                  r(i,j,23) = r(i,j,23) + r(i,k,19)*r(k,j,9)
   82          continue
   81       continue
   80    continue
         do 90 i = 1,3
            do 91 j = 1,3
               do 92 k = 1,3
                  r(i,j,6) = r(i,j,6) + r(i,k,2)*r(k,j,5)
                  r(i,j,7) = r(i,j,7) + r(i,k,16)*r(k,j,23)
                  r(i,j,8) = r(i,j,8) + r(i,k,5)*r(k,j,2)
                  r(i,j,10) = r(i,j,10) + r(i,k,2)*r(k,j,9)
                  r(i,j,11) = r(i,j,11) + r(i,k,9)*r(k,j,2)
                  r(i,j,12) = r(i,j,12) + r(i,k,23)*r(k,j,16)
                  r(i,j,14) = r(i,j,14) + r(i,k,16)*r(k,j,2)
                  r(i,j,15) = r(i,j,15) + r(i,k,2)*r(k,j,16)
                  r(i,j,22) = r(i,j,22) + r(i,k,23)*r(k,j,2)
                  r(i,j,24) = r(i,j,24) + r(i,k,2)*r(k,j,23)
   92          continue
   91       continue
   90    continue
         do 100 i = 1,3
            do 101 j = 1,3
               do 102 k = 1,3
                  r(i,j,3) = r(i,j,3) + r(i,k,5)*r(k,j,12)
                  r(i,j,4) = r(i,j,4) + r(i,k,5)*r(k,j,10)
                  r(i,j,13) = r(i,j,13) + r(i,k,23)*r(k,j,11)
                  r(i,j,17) = r(i,j,17) + r(i,k,16)*r(k,j,12)
                  r(i,j,18) = r(i,j,18) + r(i,k,16)*r(k,j,10)
                  r(i,j,21) = r(i,j,21) + r(i,k,12)*r(k,j,15)
  102          continue
  101       continue
  100    continue
         do 110 n = 1,24
            nv = n+24
            do 111 i = 1,3
               do 112 j = 1,3
                  r(i,j,nv) = -r(i,j,n)
  112          continue
  111       continue
  110    continue
      endif
c
      return
      end
c ======================================================================
      subroutine symprt (natoms,type,nops,hg,op,hasinv,nsmmrf,rc,tc,
     *   mt,att,class,inv,nptflg,hex,prtlvl,nout)
c
c Print out interesting quantities
c
c @(#) symprt -- 1.1 -- 3/5/89 09:43:04
c
c Input:
      integer hg,nops,natoms,type(natoms),op(48),att(48,natoms),
     *   mt(48,48),class(48),inv(48),prtlvl,nout
      double precision rc(3,3,48),tc(3,48)
      logical hasinv,nsmmrf,nptflg(48),hex
c Local:
      integer opx(48),i,j,k,l,n,il,iu
      character*28 cst(8)
      character*40 fmt
c Symbolic constants:
      integer tricl,monocl,rhomb,tetrag,cubic,trigrh,trighx,hexag
      parameter (tricl = 1)
      parameter (monocl = 2)
      parameter (rhomb = 3)
      parameter (tetrag = 4)
      parameter (cubic = 5)
      parameter (trigrh = 6)
      parameter (trighx = 7)
      parameter (hexag = 8)
c Names:
      data cst (tricl)  /'triclinic'/
      data cst (monocl) /'monoclinic'/
      data cst (rhomb)  /'orthorhombic'/
      data cst (tetrag) /'tetragonal'/
      data cst (cubic)  /'cubic'/
      data cst (trigrh) /'trigonal (rhombohedral axes)'/
      data cst (trighx) /'trigonal (hexagonal axes)'/
      data cst (hexag)  /'hexagonal'/
c
c Description of the group:
c
      if (hg .eq. hexag .and. nops .eq. 24) then
         write (nout,150)
      elseif (hg .eq. cubic .and. nops .eq. 48) then
         write (nout,140)
      else
         write (nout,130) cst(hg)
      endif
      if (hasinv) then
         write (nout,160)
      else
         write (nout,165)
      endif
      if (nsmmrf) then
         write (nout,200) 
      else
         write (nout,180) 
      endif
  130 format (/,' The crystal system is ',a28)
  140 format (/,' The point group of the crystal is the full cubic',
     *   ' group')
  150 format (/,' The point group of the crystal is the full hexagonal',
     *   ' group')
  160 format (' The point group contains the inversion')
  165 format (' The point group does not contain the inversion')
  180 format (' The space group of the crystal is symmorphic')
  200 format (' The space group is non-symmorphic or a non-standard ',
     *   'origin was used')
c
c Atom transformation table:
c
      write (nout,270) 
      il = 1
      if (nops .gt. 24) then
         iu = 24
      else
         iu = nops
      endif
  280 write (nout,290) (op(i),i=il,iu)
      do 340 i = 1,natoms
         do 300 j = 1,nops
            l = op(j)
            opx(j) = att(l,i)
  300    continue
         write (nout,310) type(i),(opx(j),j=il,iu)
  340 continue
      if (nops-iu .gt. 0) then
         il = 25
         iu = nops
         goto 280
      endif
  270 format (/,/,1x,'Atom transformation table',/,1x,25('-'))
  290 format (/,1x,'Type/Op',24i3)
  310 format (1x,i3,4x,24i3)
c
c Multiplication table:
c
      il = 1
      if (nops .gt. 24) then
         iu = 24
      else
         iu = nops
      endif
      write (nout,430) 
  460 write (nout,400)  (op(i),i=il,iu)
      do 510 j = 1,nops
         l = op(j)
         do 470 i = il,iu
            n = op(i)
            opx(i) = mt(l,n)
  470    continue
         write (nout,520)  l,(opx(i),i=il,iu)
  510 continue
      write (fmt,530) (iu-il+1)*3-1
      write (nout,fmt)
      write (nout,540) (class(i),i=il,iu)
      write (nout,550) (op(inv(i)),i=il,iu)
      if (iu .ne. nops) then
         il = 25
         iu = nops
         goto 460
      endif
  400 format (/,8x,24i3)
  430 format (/,/,1x,'Multiplication table',/,1x,20('-'))
  520 format (5x,i2,1x,24i3)
  530 format ('(9x,',i2,'(''-''))')
  540 format (1x,'class:',1x,24i3)
  550 format (1x,'  inv:',1x,24i3)
c
c Rotation matrices and fractional translations (cartesian):
c
      if (prtlvl .lt. 1) return
      if (nsmmrf) then
         write (nout,999) 
      else
         write (nout,1999) 
c         write (nout,1998) 
      endif
      if (hex) then
         do 1010 i = 1,nops
            l = op(i)
            if (nptflg(i)) then
               write (nout,1001)
     *            l,((rc(j,k,l),k=1,3),j=1,3),(tc(j,l),j=1,3)
            else
               write (nout,3004) l,((rc(j,k,l),k=1,3),j=1,3)
            endif
 1010    continue
      else
         do 1020 i = 1,nops
            l = op(i)
            if (nptflg(i)) then
               write (nout,3003)
     *            l,((nint(rc(j,k,l)),k=1,3),j=1,3),(tc(j,l),j=1,3)
C-----------PP input format------SM, 20 Jan 99
c	       write (nout,6001) l
c	       write (nout,6002) (nint(rc(1,k,l)),k=1,3)	
c	       write (nout,6002) (nint(rc(2,k,l)),k=1,3)	
c	       write (nout,6002) (nint(rc(3,k,l)),k=1,3)	
            else
               write (nout,3003) l,((nint(rc(j,k,l)),k=1,3),j=1,3)
C-----------PP input format------
c	       write (nout,6001) l
c	       write (nout,6002) (nint(rc(1,k,l)),k=1,3)	
c	       write (nout,6002) (nint(rc(2,k,l)),k=1,3)	
c	       write (nout,6002) (nint(rc(3,k,l)),k=1,3)	
            endif
 1020    continue
      endif
  999 format (/,/,1x,'Matrices and fractional translations ',
     *   'in Cartesian basis',/,1x,55('-'),/)
c 1998 format (/,/,1x,'Matrices in PP INPUT FORMAT',/,1x,27('-'),/)
 1999 format (/,/,1x,'Matrices in Cartesian basis',/,1x,27('-'),/)
 1001 format (1x,i2,3(1x,3(1x,f8.5)),/,10x,3(1x,f9.5))
 3004 format (1x,i2,3(1x,3(1x,f8.5)))
 3003 format (1x,i2,3(2x,3i3),3x,3(1x,f9.5))
c
c 6001 format (1x,i3,'---------------------------')
c 6002 format (3(6x,i3))
c
      return
      end
