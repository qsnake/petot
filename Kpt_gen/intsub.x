      subroutine intsub(n,xk,xw,imap,nred,nexp,
     +     ntrans,mtrans,icof,mtrx,
     +     nrk,w,rk)
c     
c     subroutine will reduce k points to the irreducible
c     brillouin zone using the symmetry operations
c     $Id: intsub.x,v 1.1 89/03/11 19:09:18 sverre Exp $
c
c     $Log:	intsub.x,v $
c     Revision 1.1  89/03/11  19:09:18  sverre
c     Initial revision
c     
c     Revision 1.1  89/03/11  19:09:18  sverre
c     Initial revision
c     
c     
      implicit DOUBLE PRECISION(a-h,o-z)
#ifdef PARAM_H
#include PARAM_H
#else
      parameter ( nrkp = nrkp )
#endif
      dimension xk(3,n),xw(n),imap(n),
     +     icof(48),mtrx(3,3,48),
     +     w(nrkp),rk(3,nrkp)
c     
c I   n           number of reducible k points
c I   xk(i,j)     i=1,3,j=1,n. irreducible k points
c I   xw(j)       j=1,n. irreducible intergration weights
c I   imap(j)     j=1,n. index mapping reducible to irreducible
c     .           k point
c O   nred        number of distinct reducible k points
c O   nexp        number of expansion k points
c     
c     if the symmetry of the crystal is given as a subgroup
c     of a larger group and the larger group has been
c     properly partitioned into cosets (icof) (by symchk),
c     the larger group is used to symmetrize the k points
c     before reduction by the subgroup
c     
c     the number of fully symmtrized k points is
c     returned in nred and the number of irreducible
c     points that cannot be obtained from the reducible
c     input set without symmetrization is returned in nexp
c     
      dimension jcof(48),yk(3,48),jmap(48),rktr(3,2)
c
c     rcs id string - allows use of ident to identify binaries
c
      character rcsid*50
      rcsid = '$RCSfile: intsub.x,v $$Revision: 1.1 $'
c     
      nrk  = 0
      nred = 0
      nexp = 0
c     
c     outer loop over reducible k points
c     

      do 260 i=1,n
         if (imap(i) .ne. 0) goto 260
c     
c     new irreducible point
c     mark with negative imap
c     
         imap(i) = - (nrk + 1)
         wi = xw(i)
c     
c     initialize jcof = icof
c     
         do 100 j=1,mtrans
            jcof(j) = icof(j)
  100    continue
c     
c     first generate the star of reducible points - yk -
c     related to the current irreducible point by symmetry.
c     jcof is used to keep track of points that are not
c     related by a subgroup symmetry operation.
c     
         nredi = 0
         do 210 j=1,mtrans
c     
c     operate on irreducible xk with the symmetry operation
c     
            do 120 k=1,3
               rktr(k,1) = 0.D0
               do 110 l=1,3
                  rktr(k,1) = rktr(k,1) + DBLE(mtrx(k,l,j))*xk(l,i)
  110          continue
c     
c     translate to interval (0-1) (<1)
c     and find inverse
c     
               kdel = rktr(k,1) + 1.D-6
               if (rktr(k,1) .lt. -1.D-6) kdel = kdel - 1
               rktr(k,1) = rktr(k,1) - DBLE(kdel)
               rktr(k,2) = 1.D0 - rktr(k,1)
               if (rktr(k,1) .lt. 1.D-6) rktr(k,2) = rktr(k,1)
  120       continue
c     
c     loop over transformed xk and its inverse
c     
ccccccccccc note that, k=2 (inversing), can be just considered as 
ccccccccccc another additional symmetry operation (inversion) for
ccccccccccc all the systems (for k-point only). The symm. itself
ccccccccccc mtrx might contain a inversion also. In that case, the inversion
ccccccccccc symm. is operated twice. But that will not screw up the current
ccccccccccc algorithm. The second same operation will not generate any new 
ccccccccccc ireducable kpoints. Note in the original k-point grid, xk(i=1,n), 
ccccccccccc the whole k-point grid generated from kpgen.input are stored (include
ccccccccccc both k and -k (if they both exist). So, only at this step, the -k
ccccccccccc will be reduced due to symm. op (include k=1,2). 
ccccccccccc L.W. Wang

            do 200 k=1,2

ccccccc  L.W. Wang, July 3, 2001 
cccccccc  Since the symmetrized k-point rktr(l,k) might not be within the
cccccccc  original kpoint set xk(l,ii) (in this version, the kpoint set xk(l,ii) 
cccccccc  might not satisfy the symm. op), so, only select those rktr, which 
cccccccc  fall back to the original set xk(l,ii)), otherwise, ignore it. 

cccccc  Test, whether rktr is within xk set

           do ii=1,n
           if((abs(xk(1,ii)-rktr(1,k)).lt.1.D-6).and.
     &     (abs(xk(2,ii)-rktr(2,k)).lt.1.D-6).and.
     &     (abs(xk(3,ii)-rktr(3,k)).lt.1.D-6)) goto 133   ! it is within xk(l,ii)
           enddo
           goto 200       ! it is not within xk, forget about this kpt. 
133        continue
ccccccccccc  L.W. Wang
 
c     
c     check if new reducible point in yk
c     
               if (nredi .eq. 0) goto 160
               do 150 l=1,nredi
                  do 130 m=1,3
                     if (abs(rktr(m,k)-yk(m,l)) .gt. 1.D-6) goto 150
  130             continue
c     
c     point already found - check that jcof agree
c     if not reset jcof for all symmetry elements
c     in the same coset
c     
                  if (jcof(j) .ne. jcof(jmap(l))) then
                     jcofj = jcof(j)
                     do 140 m=1,mtrans
                        if (jcof(m) .eq. jcofj) jcof(m) = jcof(jmap(l))
  140                continue
                  end if
                  goto 200
  150          continue
c     
c     new distinct point - add to array yk
c     
  160          nred  = nred  + 1
               nredi = nredi + 1
               jmap(nredi) = j
               do 170 l=1,3
                  yk(l,nredi) = rktr(l,k)
  170          continue
c     
c     mark (map out) all remaining reducible k points related
c     to current point by symmetry with a positive imap
c     pointing to the current irreducible point
c     
               if (i .eq. n) goto 200
               ip = i+1
               do 190 l=ip,n
                  if (imap(l) .ne. 0) goto 190
                  do 180 m=1,3
                     if (abs(rktr(m,k)-xk(m,l)) .gt. 1.D-6) goto 190
  180             continue
                  wi = wi + xw(l)
                  imap(l) = nrk + 1
  190          continue
  200       continue
  210    continue
c     
c     add inequivalent (with different jcof) points in
c     yk to the irreducible set rk
c     
         nexpi = 0
         do 250 j=1,nredi
            if (jmap(j) .eq. 0) goto 250
c     
c     inequivalent point found
c     
            nexpi = nexpi + 1
            jcofj = jcof(jmap(j))
            jrk = j
c     
c     count (for weights) and remove equivalent points
c     swap points to choose equivalent point with
c     the smaller first components
c     
            neqv = 0
            do 230 k=j,nredi
               if (jmap(k) .eq. 0) goto 230
               if (jcof(jmap(k)) .ne. jcofj) goto 230
               neqv = neqv + 1
               jmap(k) = 0
               do 220 l=1,3
                  dy = yk(l,k) - yk(l,jrk)
                  if (jrk .eq. k .or. dy .gt. 1.D-6) goto 230
                  if (dy .lt. -1.D-6) jrk = k
  220          continue
  230       continue
c     
c     add point to irreducible set rk
c     
            nrk = nrk + 1
            if (nrk .gt. nrkp) call warn(1021,xdum,nrkp)
            do 240 k=1,3
               rk(k,nrk) = yk(k,jrk)
  240       continue

            w(nrk) = wi * DBLE(neqv) / DBLE(nredi)
  250    continue
         nexp = nexp + nexpi - 1
  260 continue
      return
      end
