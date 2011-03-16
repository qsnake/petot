      subroutine polcoe(x,y,n,cof)
c
c ************************************************
c *  njtj                                        *
c *  Returns the coefficients of a polynominal.  *
c *  Taken from numerical recipes, page 93.      *
c *  njtj                                        *
c ************************************************
c
c  njtj
c  ###  Cray conversions  
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c  
      implicit double precision (a-h,o-z)
c
      parameter (nmax=10,zero=0.D0,one=1.D0)
Cray      parameter (nmax=10,zero=0.0,one=1.0)
c
      dimension x(n),y(n),cof(n),s(nmax)
      do 11 i=1,n
        s(i)=zero
        cof(i)=zero
11    continue
      s(n)=-x(1)
      do 13 i=2,n
        do 12 j=n+1-i,n-1
          s(j)=s(j)-x(i)*s(j+1)
12      continue
        s(n)=s(n)-x(i)
13    continue
      do 16 j=1,n
        phi=n
        do 14 k=n-1,1,-1
          phi=k*s(k+1)+x(j)*phi
14      continue
        ff=y(j)/phi
        b=one
        do 15 k=n,1,-1
          cof(k)=cof(k)+b*ff
          b=s(k)+x(j)*b
15      continue
16    continue
      return
      end
