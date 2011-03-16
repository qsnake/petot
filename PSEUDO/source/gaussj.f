      subroutine gaussj(a,n,np,b,m,mp)
c 
c ****************************************************************
c *                                                              *
c *    Gauss-Jordan routine used by pseud* to find polynominal   *
c *  constants.  Taken from Numerical Recipes, page 28.          *
c *  Version dated May 7, 1991                                   *
c *  njtj                                                        *
c *                                                              *
c ****************************************************************
c
      implicit real*8 (a-h,o-z)
c
      parameter (nmax=50,zero=0.0d0,one=1.0d0)
c
      dimension a(np,np),b(np,mp),ipiv(nmax),indxr(nmax),indxc(nmax)
c
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=zero
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                write(6,100)
                stop 'gaussj one'
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.zero) then
          write(6,100)
          stop 'gaussj two'
        endif
        pivinv=one/a(icol,icol)
        a(icol,icol)=one
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=zero
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
 100  format(//,'Singular matrix, stopped in gaussj')
      end
