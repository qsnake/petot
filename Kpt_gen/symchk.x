      subroutine symchk(ierr,
     +     nout,
     +     ntrans,mtrans,icof,mtrx,tnp)
c     
c     check group property of symmetry operations : mtrx, tnp
c     missing operations are added to complete group (if possible)
c     $Id: symchk.x,v 1.1 89/03/11 19:09:05 sverre Exp $
c
c     $Log:	symchk.x,v $
c     Revision 1.1  89/03/11  19:09:05  sverre
c     Initial revision
c     
c     Revision 1.1  89/03/11  19:09:05  sverre
c     Initial revision
c     
c     
      implicit DOUBLE PRECISION(a-h,o-z)
      dimension icof(48),mtrx(3,3,48),tnp(3,48)
c     
c O   ierr                returns status of group
c     .                   ierr = 0  group ok as input
c     .                   ierr < 0  group ok after adding elements
c     .                   ierr > 0  number of suspected operation
c     
      dimension mult(48,48),nerr(48),iinv(48),mtest(3,3)
      common /work3/ mult,nerr,iinv
c
c     rcs id string - allows use of ident to identify binaries
c
      character rcsid*50
      rcsid = '$RCSfile: symchk.x,v $$Revision: 1.1 $'
c     
      pi2 = 8.D0*atan(1.D0)
      ierr = 0
      do 110 i=1,48
         nerr(i) = 0
         iinv(i) = 0
         if (i .gt. mtrans) icof(i) = 0
         do 100 j=1,48
            mult(i,j) = 0
  100    continue
  110 continue
c     
c     test if group is non primitive and
c     check for duplicate operations
c     
      np = 0
      do 180 i=1,mtrans
         do 120 j=1,3
            if (abs(tnp(j,i)) .gt. 1.D-6) np = 1
  120    continue
         im = i-1
         if (im .eq. 0) goto 180
         do 170 j=1,im
            iflg = 0
            do 140 k=1,3
               if (abs(tnp(k,i)-tnp(k,j)) .gt. 1.D-6) iflg = 1
               do 130 l=1,3
                  if (mtrx(k,l,i) .ne. mtrx(k,l,j)) goto 170
  130          continue
  140       continue
            write(nout,150) j,i
            if (iflg .ne. 0) write(nout,160)
  150       format(/,' symmetry operations',i3,' and',i3,' are equal')
  160       format(' apart from the nonprimitive translations')
            ierr = i
            call warn(1011,xdum,ierr)
            return
  170    continue
  180 continue
c     
c     locate identity operation
c     
      do 210 i=1,mtrans
         do 200 j=1,3
            if (mtrx(j,j,i) .ne. 1) goto 210
            do 190 k=1,3
               if (j .ne. k. and. mtrx(j,k,i) .ne. 0) goto 210
  190       continue
  200    continue
         iden = i
         goto 240
  210 continue
c     
c     identity missing - add identity
c     
      if (mtrans .lt. 48) goto 230
      write(nout,220)
  220 format(/,' symmetry operation 1 (identity) is  missing')
      call warn(1012,xdum,0)
      ierr = 1
      return
  230 mtrans = mtrans + 1
      ierr = -1
      iden = mtrans
c     
c     identity must be first operation
c     
  240 icof(iden) = icof(1)
      icof(1) = 1
      do 260 j=1,3
         tnp(j,iden) = tnp(j,1)
         tnp(j,1) = 0.0
         do 250 k=1,3
            mtrx(j,k,iden) = mtrx(j,k,1)
            if (j .ne. k) mtrx(j,k,1) = 0
            if (j .eq. k) mtrx(j,k,1) = 1
  250    continue
  260 continue
      iden = 1
      ntrans = 1
c     
c     construct multiplication table
c     
      nold = 1
  270 nnew = mtrans
      do 370 i=1,nnew
         do 360 j=1,nnew
            if (i .lt. nold .and. j .lt. nold) goto 360
c     
c     mulitiply i and j
c     
            do 300 k=1,3
               do 290 l=1,3
                  mtest(k,l) = 0
                  do 280 m=1,3
                     mtest(k,l) = mtest(k,l) + mtrx(k,m,i)*mtrx(m,l,j)
  280             continue
  290          continue
  300       continue
c     
c     check for match
c     
            do 330 k=1,mtrans
               do 320 l=1,3
                  do 310 m=1,3
                     if (mtest(l,m) .ne. mtrx(l,m,k)) goto 330
  310             continue
  320          continue
c     
c     match - make multiplcation table entry and check for inverse
c     and subgroup property
c     
               mult(i,j) = k
               if (k .eq. iden) iinv(i) = j
               if (icof(i)+icof(j)+icof(k) .ne. 2) goto 360
               icof(i) = 1
               icof(j) = 1
               icof(k) = 1
               goto 360
  330       continue
c     
c     no match - add mtest to operations
c     
            if (mtrans .ge. 48) goto 360
            mtrans = mtrans + 1
            ierr = -1
            mult(i,j) = mtrans
            icof(mtrans) = icof(i)*icof(j)
            do 350 k=1,3
               tnp(k,mtrans) = 0.D0
               do 340 l=1,3
                  mtrx(k,l,mtrans) = mtest(k,l)
  340          continue
  350       continue
  360    continue
  370 continue
c     
c     check nonprimitive translations
c     
      if (np .eq. 0) goto 510
c     
c     compute mtrx(i)*tnp(j) + tnp(i) and compare with tnp(k)
c     
      do 500 i=1,nnew
c     .                                            -1
c     find mtrx(i) for real space (= trans(mtrx(i))  )
c     
         if (iinv(i) .eq. 0) goto 400
c     
c     inverse known
c     
         do 390 j=1,3
            do 380 k=1,3
               mtest(k,j) = mtrx(j,k,iinv(i))
  380       continue
  390    continue
         goto 450
c     
c     inverse unknown - compute
c     
  400    mtest(1,1) = mtrx(2,2,i)*mtrx(3,3,i) - mtrx(2,3,i)*mtrx(3,2,i)
         mtest(1,2) = mtrx(2,3,i)*mtrx(3,1,i) - mtrx(2,1,i)*mtrx(3,3,i)
         mtest(1,3) = mtrx(2,1,i)*mtrx(3,2,i) - mtrx(2,2,i)*mtrx(3,1,i)
         mtest(2,1) = mtrx(3,2,i)*mtrx(1,3,i) - mtrx(3,3,i)*mtrx(1,2,i)
         mtest(2,2) = mtrx(3,3,i)*mtrx(1,1,i) - mtrx(3,1,i)*mtrx(1,3,i)
         mtest(2,3) = mtrx(3,1,i)*mtrx(1,2,i) - mtrx(3,2,i)*mtrx(1,1,i)
         mtest(3,1) = mtrx(1,2,i)*mtrx(2,3,i) - mtrx(1,3,i)*mtrx(2,2,i)
         mtest(3,2) = mtrx(1,3,i)*mtrx(2,1,i) - mtrx(1,1,i)*mtrx(2,3,i)
         mtest(3,3) = mtrx(1,1,i)*mtrx(2,2,i) - mtrx(1,2,i)*mtrx(2,1,i)
         idet = mtrx(1,1,i)*mtest(1,1)
     +        + mtrx(1,2,i)*mtest(1,2)
     +        + mtrx(1,3,i)*mtest(1,3)
         if (iabs(idet) .ne. 1) then
            write(nout,410) i
  410       format(/,' symmetry operation',i3,' is not pure rotation')
            ierr = i
            call warn(1011,xdum,ierr)
            return
         end if
  420    do 440 j=1,3
            do 430 k=1,3
               mtest(j,k) = idet * mtest(j,k)
  430       continue
  440    continue
c     
c     check translation - set mult(i,j) to -1 on error
c     
  450    do 490 j=1,nnew
            if (i .lt. nold .and. j .lt. nold) goto 490
            k = mult(i,j)
            if (k .le. 0) goto 490
            do 480 l=1,3
               ttest = tnp(l,i)
               do 460 m=1,3
                  ttest = ttest + DBLE(mtest(l,m))*tnp(m,j)
  460          continue
               if (k .le. nnew) goto 470
               itest = (ttest + 1.D-4) / pi2
               if (ttest .lt. -1.D-4) itest = itest - 1
               ttest = ttest - pi2 * DBLE(itest)
               tnp(l,k) = ttest
               goto 480
  470          ttest = abs(ttest-tnp(l,k))/pi2
               itest = ttest + 1.D-3
               if (abs(ttest-DBLE(itest)) .lt. 1.D-4) goto 480
               mult(i,j) = -1
               goto 490
  480       continue
  490    continue
  500 continue
c     
c     go back and recompute if operations have been added
c     
  510 nold = nnew + 1
      if (nnew .ne. mtrans) goto 270
c     
c     check multiplication table
c     
      do 530 i=1,mtrans
         do 520 j=1,mtrans
            if (mult(i,j) .gt. 0) goto 520
            nerr(i) = nerr(i) + 1
            nerr(j) = nerr(j) + 1
  520    continue
  530 continue
c     
c     find element with max error
c     
      maxerr = 0
      do 540 i=1,mtrans
         if (nerr(i) .le. maxerr) goto 540
         maxerr = nerr(i)
         ierr = i
  540 continue
      if (ierr .gt. 0) goto 630
c     
c     examine subgroup and generate cosets
c     
      do 560 i=1,mtrans
         do 550 j=1,mtrans
            if (icof(j) .ne. 1) goto 550
            k = mult(j,i)
            if (icof(k) .eq. 0) icof(k) = i
  550    continue
  560 continue
c     
c     sort group so that subgroup appears first
c     
      maxsub = 0
      minoth = mtrans + 1
  570 do 580 i=1,mtrans
         ii = mtrans - i + 1
         if (icof(i)  .eq. 1) maxsub = i
         if (icof(ii) .ne. 1) minoth = ii
  580 continue
      if (maxsub .lt. minoth) goto 620
c     
c     change elements maxsub and minoth
c     
      icof(maxsub) = icof(minoth)
      icof(minoth) = 1
      do 600 i=1,3
         tnpi = tnp(i,maxsub)
         tnp(i,maxsub) = tnp(i,minoth)
         tnp(i,minoth) = tnpi
         do 590 j=1,3
            mtrxij = mtrx(i,j,maxsub)
            mtrx(i,j,maxsub) = mtrx(i,j,minoth)
            mtrx(i,j,minoth) = mtrxij
  590    continue
  600 continue
c     
c     fix indices icof
c     
      do 610 i=minoth,mtrans
         if (icof(i) .eq. minoth) icof(i) = maxsub
  610 continue
      goto 570
c     
c     subgroup sorted - return
c     
  620 ntrans = maxsub
      return
c     
c     multiplication table in error - printout
c     
  630 write(nout,640)
  640 format('1multiplication table',/)
      do 660 i=1,mtrans
         write(nout,650) icof(i),(mult(i,j),j=1,mtrans)
  650    format(1x,i2,2x,48i2)
  660 continue
      call warn(1011,xdum,ierr)
      return
      end
