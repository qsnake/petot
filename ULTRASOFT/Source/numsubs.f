c
c Copyright (C) 2002 Vanderbilt group
c This file is distributed under the terms of the GNU General Public
c License as described in the file 'License' in the current directory.
c
c-------------------------------------------------------------------------
c
c     ************************************************************
      subroutine simeq(ndim,n,b,x,y)
c     ************************************************************
c     solve b * y = x
      implicit double precision (a-h,o-z)
      parameter(idm8=24)
      dimension b(ndim,ndim),bb(idm8,idm8),x(ndim),y(ndim),iw(idm8)
      if (n.gt.ndim.or.ndim.gt.idm8) stop 24
c     -------
c     nag lib
c     -------
c     ifail=0
c     call f04atf(b,ndim,x,n,y,bb,4,w1,w2,ifail)
c     if (ifail.ne.0) stop 16
c     -------
c     num rec
c     -------
      do 22 i=1,n
      y(i)=x(i)
      do 22 j=1,n
   22 bb(i,j)=b(i,j)
      call ludcmp(bb,n,idm8,iw,dwrk)
      call lubksb(bb,n,idm8,iw,y)
c     -------
      return
      end
c
c----------------------------------------------------------------------------
c ===============================================================
c     driver routines first
c ===============================================================
c     ***************************************************************
      subroutine nrinv(a,na,n,b,nb,indx)
c     ***************************************************************
c     inverts matrix using numerical recipes subroutines
c     note matrix a is destroyed
c
c ---------------------------------------------------------------
      implicit double precision (a-h,o-z)
      intrinsic abs,max,min,sqrt,exp,log,sin,cos,tan,asin,acos,
     1   atan,sinh,cosh,tanh
c ---------------------------------------------------------------
c
      dimension a(na,n),b(nb,n),indx(na)
c
c     set up identity
      do 20 i=1,n
      do 10 j=1,n
   10 b(i,j)=0.d0
   20 b(i,i)=1.d0
c
      call ludcmp(a,n,na,indx,d)
c
      do 30 i=1,n
      call lubksb(a,n,na,indx,b(1,i))
   30 continue
c
      return
      end
c ===============================================================
c     now copies of numer recipes - converted to dp - note params
c ===============================================================
c     ***************************************************************
      subroutine ludcmp(a,n,np,indx,d)
c     ***************************************************************
c     -----------------------------------
      implicit double precision (a-h,o-z)
c     -----------------------------------
      parameter (nmax=100,tiny=1.0d-20)
      dimension a(np,np),indx(n),vv(nmax)
      if (n.gt.nmax) stop 33
      d=1.d0
      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
11      continue
        if (aamax.eq.0.d0) then
           write(*,*) 'singular matrix.'
           stop
        endif
        vv(i)=1.d0/aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.d0
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*dabs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.d0)a(j,j)=tiny
        if(j.ne.n)then
          dum=1.d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      end
c     ***************************************************************
      subroutine lubksb(a,n,np,indx,b)
c     ***************************************************************
c     -----------------------------------
      implicit double precision (a-h,o-z)
c     -----------------------------------
      dimension a(np,np),indx(n),b(n)
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      end
c     ***************************************************************
      subroutine spline1(x,y,n,yp1,ypn,y2)
c     ***************************************************************
c     -----------------------------------
      implicit double precision (a-h,o-z)
c     -----------------------------------
      parameter (nmax=1000)
      dimension x(n),y(n),y2(n),u(nmax)
      if (n.gt.nmax) stop 33
      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     *      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      end
c     ***************************************************************
      subroutine splint3(xa,ya,y2a,n,x,y)
c     ***************************************************************
c     ------------------------------------------------
c     my version of splint; produces y', y'', y''' too
c     ------------------------------------------------
c     -----------------------------------
      implicit double precision (a-h,o-z)
c     -----------------------------------
      dimension xa(n),ya(n),y2a(n),y(4)
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.d0) then
         write(*,*) 'bad xa input.'
         stop
      endif
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y(1)=a*ya(klo)+b*ya(khi)+
     *      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
      y(2)=(-ya(klo)+ya(khi))/h+
     *      ((-3*a**2+1)*y2a(klo)+(3*b**2-1)*y2a(khi))*h/6.d0
      y(3)=a*y2a(klo)+b*y2a(khi)
      y(4)=(-y2a(klo)+y2a(khi))/h
      return
      end
c      *****************************************************************
c
      subroutine splift(x,y,yp,ypp,n,w,ierr,isx,a1,b1,an,bn)
      implicit double precision(a-h,o-z)
c
c     sandia mathematical program library
c     applied mathematics division 2613
c     sandia laboratories
c     albuquerque, new mexico  87185
c     control data 6600/7600  version 7.2  may 1978
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                    issued by sandia laboratories                     *
c  *                   a prime contractor to the                       *
c  *                united states department of energy                 *
c  * * * * * * * * * * * * * * * notice  * * * * * * * * * * * * * * * *
c  * this report was prepared as an account of work sponsored by the   *
c  * united states government.  neither the united states nor the      *
c  * united states department of energy nor any of their employees,    *
c  * nor any of their contractors, subcontractors, or their employees  *
c  * makes any warranty, express or implied, or assumes any legal      *
c  * liability or responsibility for the accuracy, completeness or     *
c  * usefulness of any information, apparatus, product or process      *
c  * disclosed, or represents that its use would not infringe          *
c  * owned rights.                                                     *
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  * the primary document for the library of which this routine is     *
c  * part is sand77-1441.                                              *
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     written by rondall e. jones
c
c     abstract
c         splift fits an interpolating cubic spline to the n data points
c         given in x and y and returns the first and second derivatives
c         in yp and ypp.  the resulting spline (defined by x, y, and
c         ypp) and its first and second derivatives may then be
c         evaluated using splint.  the spline may be integrated using
c         spliq.  for a smoothing spline fit see subroutine smoo.
c
c     description of arguments
c         the user must dimension all arrays appearing in the call list,
c         e.g.   x(n), y(n), yp(n), ypp(n), w(3n)
c
c       --input--
c
c         x    - array of abscissas of data (in increasing order)
c         y    - array of ordinates of data
c         n    - the number of data points.  the arrays x, y, yp, and
c                ypp must be dimensioned at least n.  (n .ge. 4)
c         isx  - must be zero on the initial call to splift.
c                if a spline is to be fitted to a second set of data
c                that has the same set of abscissas as a previous set,
c                and if the contents of w have not been changed since
c                that previous fit was computed, then isx may be
c                set to one for faster execution.
c         a1,b1,an,bn - specify the end conditions for the spline which
c                are expressed as constraints on the second derivative
c                of the spline at the end points (see ypp).
c                the end condition constraints are
c                        ypp(1) = a1*ypp(2) + b1
c                and
c                        ypp(n) = an*ypp(n-1) + bn
c                where
c                        abs(a1).lt. 1.0  and  abs(an).lt. 1.0.
c
c                the smoothest spline (i.e., least integral of square
c                of second derivative) is obtained by a1=b1=an=bn=0.
c                in this case there is an inflection at x(1) and x(n).
c                if the data is to be extrapolated (say, by using splint
c                to evaluate the spline outside the range x(1) to x(n)),
c                then taking a1=an=0.5 and b1=bn=0 may yield better
c                results.  in this case there is an inflection
c                at x(1) - (x(2)-x(1)) and at x(n) + (x(n)-x(n-1)).
c                in the more general case of a1=an=a  and b1=bn=0,
c                there is an inflection at x(1) - (x(2)-x(1))*a/(1.0-a)
c                and at x(n) + (x(n)-x(n-1))*a/(1.0-a).
c
c                a spline that has a given first derivative yp1 at x(1)
c                and ypn at y(n) may be defined by using the
c                following conditions.
c
c                a1=-0.5
c
c                b1= 3.0*((y(2)-y(1))/(x(2)-x(1))-yp1)/(x(2)-x(1))
c
c                an=-0.5
c
c                bn=-3.0*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)/(x(n)-x(n-1))
c
c       --output--
c
c         yp   - array of first derivatives of spline (at the x(i))
c         ypp  - array of second derivatives of spline (at the x(i))
c         ierr - a status code
c              --normal code
c                 1 means that the requested spline was computed.
c              --abnormal codes
c                 2 means that n, the number of points, was .lt. 4.
c                 3 means the abscissas were not strictly increasing.
c
c       --work--
c
c         w    - array of working storage dimensioned at least 3n.
      dimension x(n),y(n),yp(n),ypp(n),w(n,3)
c
      if (n.lt.4) go to 200
      nm1  = n-1
      nm2  = n-2
      if (isx.gt.0) go to 40
      do 5 i=2,n
      if (x(i)-x(i-1)) 300,300,5
    5 continue
c
c     define the tridiagonal matrix
c
      w(1,3) = x(2)-x(1)
      do 10 i=2,nm1
      w(i,2) = w(i-1,3)
      w(i,3) = x(i+1)-x(i)
   10 w(i,1) = 2.d0*(w(i,2)+w(i,3))
      w(1,1) = 4.d0
      w(1,3) =-4.d0*a1
      w(n,1) = 4.d0
      w(n,2) =-4.d0*an
c
c     l u decomposition
c
      do 30 i=2,n
      w(i-1,3) = w(i-1,3)/w(i-1,1)
   30 w(i,1)   = w(i,1) - w(i,2)*w(i-1,3)
c
c     define *constant* vector
c
   40 ypp(1) = 4.d0*b1
      dold   = (y(2)-y(1))/w(2,2)
      do 50 i=2,nm2
      dnew   = (y(i+1) - y(i))/w(i+1,2)
      ypp(i) = 6.d0*(dnew - dold)
      yp(i)  = dold
   50 dold   = dnew
      dnew   = (y(n)-y(n-1))/(x(n)-x(n-1))
      ypp(nm1) = 6.d0*(dnew - dold)
      ypp(n) = 4.d0*bn
      yp(nm1)= dold
      yp(n)  = dnew
c
c     forward substitution
c
      ypp(1) = ypp(1)/w(1,1)
      do 60 i=2,n
   60 ypp(i) = (ypp(i) - w(i,2)*ypp(i-1))/w(i,1)
c
c     backward substitution
c
      do 70 j=1,nm1
      i = n-j
   70 ypp(i) = ypp(i) - w(i,3)*ypp(i+1)
c
c     compute first derivatives
c
      yp(1)  = (y(2)-y(1))/(x(2)-x(1)) - (x(2)-x(1))*(2.d0*ypp(1)
     1         + ypp(2))/6.d0
      do 80 i=2,nm1
   80 yp(i)  = yp(i) + w(i,2)*(ypp(i-1) + 2.d0*ypp(i))/6.d0
      yp(n)  = yp(n) + (x(n)-x(nm1))*(ypp(nm1) + 2.d0*ypp(n))/6.d0
c
      ierr = 1
      return
  200 ierr = 2
      write(6,210)
  210 format(47h in splift, there were less than 4 data values.)
      return
  300 ierr = 3
      write(6,310)
  310 format(11h in splift,,
     144h the abscissas were not strictly increasing.)
      return
      end
c
c-------------------------------------------------------------------------
c
      subroutine bndglm(a,np,n,mtp1,nt,emach)
      implicit double precision (a-h,o-z)
      dimension a(np,mtp1),nt(n)
      m=(mtp1-1)/3
      mdp1=m+m+1
      do 2 i=1,n-1
        nrow=min(m+1,n-i+1)
        in=idamax(nrow,a(i+nrow-1,m-nrow+2),np-1)
        lpiv=nrow-in+i
        nt(i)=lpiv
        if (lpiv .ne. i) then
          call dswap1(mdp1,a(i,m+1),np,a(lpiv,i-lpiv+m+1),np)
        endif
        if ( abs(a(i,m+1)) .lt. emach) a(i,m+1)=emach
        do 1 j=1,nrow-1
          if (abs(a(i+j,m+1-j)) .gt. emach ) then
            a(i+j,m+1-j)=a(i+j,m+1-j)/a(i,m+1)
            call daxpy(m+m,-a(i+j,m+1-j),a(i,m+2),np,a(i+j,m+2-j),np)
          endif
1       continue
2     continue
      if (abs(a(n,m+1)) .lt. emach ) a(n,m+1)=emach
      return
      end
      subroutine bndsub(a,np,n,mtp1,nt,x)
      implicit double precision (a-h,o-z)
      dimension a(np,mtp1),nt(n),x(n)
      m=(mtp1-1)/3
      do 1 i=1,n-1
        if ( nt(i) .ne. i ) then
          in=nt(i)
          dum=x(i)
          x(i)=x(in)
          x(in)=dum
        endif
        nrow=min(m+1,n-i+1)
        call daxpy(nrow-1,-x(i),a(i+nrow-1,m-nrow+2),np-1,x(i+1),-1)
1     continue
      x(n)=x(n)/a(n,m+1)
      do 2 i=n-1,1,-1
        lrow=min(n-i,m*2)
        x(i)=x(i)-ddot(lrow,a(i,m+2),np,x(i+1),1)
2       x(i)=x(i)/a(i,m+1)
      return
      end
      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     overwrite double precision dy with double precision da*dx + dy.
c     for i = 0 to n-1, replace       dy(ly+i*incy) with da*dx(lx+i*incx) +
c      dy(ly+i*incy), where lx = 1 if incx .ge. 0, else lx = (-incx)*n,
c      and ly is defined in a similar way using incy.
c
      double precision dx(*),dy(*),da
      if(n.le.0.or.da.eq.0.d0) return
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue
c
c       code for nonequal or nonpositive increments.
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
      dy(iy) = dy(iy) + da*dx(ix)
      ix = ix + incx
      iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop so remaining vector length is a multiple of 4.
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
      dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
      dy(i) = dy(i) + da*dx(i)
      dy(i + 1) = dy(i + 1) + da*dx(i + 1)
      dy(i + 2) = dy(i + 2) + da*dx(i + 2)
      dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
c
c       code for equal, positive, nonunit increments.
c
   60 continue
      ns = n*incx
        do 70 i=1,ns,incx
        dy(i) = da*dx(i) + dy(i)
   70        continue
      return
      end
      double precision function ddot(n,dx,incx,dy,incy)
c
c     returns the dot product of double precision dx and dy.
c     ddot = sum for i = 0 to n-1 of  dx(lx+i*incx) * dy(ly+i*incy)
c     where lx = 1 if incx .ge. 0, else lx = (-incx)*n, and ly is
c     defined in a similar way using incy.
c
      double precision dx(*),dy(*)
      ddot = 0.d0
      if(n.le.0)return
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue
c
c        code for unequal or nonpositive increments.
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
       ddot = ddot + dx(ix)*dy(iy)
      ix = ix + incx
      iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1.
c
c
c       clean-up loop so remaining vector length is a multiple of 5.
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
       ddot = ddot + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
       ddot = ddot + dx(i)*dy(i) + dx(i+1)*dy(i+1) +
     $ dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
      return
c
c        code for positive equal increments .ne.1.
c
   60 continue
      ns = n*incx
        do 70 i=1,ns,incx
        ddot = ddot + dx(i)*dy(i)
   70        continue
      return
      end
      subroutine dswap1(n,dx,incx,dy,incy)
c
c     interchange double precision dx and double precision dy.
c     for i = 0 to n-1, interchange  dx(lx+i*incx) and dy(ly+i*incy),
c     where lx = 1 if incx .ge. 0, else lx = (-incx)*n, and ly is
c     defined in a similar way using incy.
c
      double precision dx(1),dy(1),dtemp1,dtemp2,dtemp3
      if(n.le.0)return
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue
c
c      code for unequal or nonpositive increments.
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
      dtemp1 = dx(ix)
      dx(ix) = dy(iy)
      dy(iy) = dtemp1
      ix = ix + incx
      iy = iy + incy
   10 continue
      return
c
c      code for both increments equal to 1
c
c
c      clean-up loop so remaining vector length is a multiple of 3.
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
      dtemp1 = dx(i)
      dx(i) = dy(i)
      dy(i) = dtemp1
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
      dtemp1 = dx(i)
      dtemp2 = dx(i+1)
      dtemp3 = dx(i+2)
      dx(i) = dy(i)
      dx(i+1) = dy(i+1)
      dx(i+2) = dy(i+2)
      dy(i) = dtemp1
      dy(i+1) = dtemp2
      dy(i+2) = dtemp3
   50 continue
      return
   60 continue
c
c     code for equal, positive, nonunit increments.
c
      ns = n*incx
      do 70 i=1,ns,incx
      dtemp1 = dx(i)
      dx(i) = dy(i)
      dy(i) = dtemp1
   70      continue
      return
      end
      integer function idamax(n,dx,incx)
c
c     find smallest index of maximum magnitude of double precision dx.
c     idamax =      first i, i = 1 to n, to minimize  abs(dx(1-incx+i*incx))
c
      double precision dx(*),dmax,xmag
      idamax = 0
      if(n.le.0) return
      idamax = 1
      if(n.le.1)return
      if(incx.eq.1)goto 20
c
c       code for increments not equal to 1.
c
      dmax = dabs(dx(1))
      ns = n*incx
      ii = 1
        do 10 i = 1,ns,incx
        xmag = dabs(dx(i))
        if(xmag.le.dmax) go to 5
        idamax = ii
        dmax = xmag
    5        ii = ii + 1
   10        continue
      return
c
c       code for increments equal to 1.
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
        xmag = dabs(dx(i))
        if(xmag.le.dmax) go to 30
        idamax = i
        dmax = xmag
   30 continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine simpagn(n,r,f,sum,iflog)
c
c     simpson's rule for integration
c     iflog  eq  1,   log mesh
c     iflog  eq  0,   uniform mesh
c
      implicit double precision (a-h,o-z)
      dimension r(n),f(n),a(4)
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
      if(n.lt.8) stop 'simpagn'
      a(1)=17.0/48.0
      a(2)=59.0/48.0
      a(3)=43.0/48.0
      a(4)=49.0/48.0
c
      if(iflog.eq.0) then
c
        delta=(r(n)-r(1))/dfloat(n-1)
        sum=a(1)*(f(1)+f(n))+a(2)*(f(2)+f(n-1))
     c     +a(3)*(f(3)+f(n-2))+a(4)*(f(4)+f(n-3))
        do 10 i=5,n-4
 10     sum=sum+f(i)
        sum=sum*delta
c
      else
c
c     check if the mesh is logarithmic
        delta=log(r(n)/r(1))/dfloat(n-1)
        absdf=0.0
        do 20 i=2,n
 20     absdf=absdf+abs(log(r(i)/r(i-1))-delta)
        if(absdf.gt.1.0d-3) then
          write(stderr,*) 'absdf = ',absdf
          stop 'absdf'
        endif
c     sum
        sum=a(1)*(f(1)*r(1)+f(n)*r(n))+a(2)*(f(2)*r(2)+f(n-1)*r(n-1))
     c +a(3)*(f(3)*r(3)+f(n-2)*r(n-2))+a(4)*(f(4)*r(4)+f(n-3)*r(n-3))
        do 30 i=5,n-4
 30     sum=sum+f(i)*r(i)
        sum=sum*delta
c
      endif
c
      return
      end
