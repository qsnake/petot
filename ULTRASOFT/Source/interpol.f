c
c Copyright (C) 2002 Vanderbilt group
c This file is distributed under the terms of the GNU General Public
c License as described in the file 'License' in the current directory.
c
      subroutine interpol(xin,yin,nin,xout,yout,nout,y2,u)
      implicit none
      integer nin, nout
      real*8 xin(nin), yin(nin), xout(nout), yout(nout),
     *     y2(nin), u(nin), yp1, ypn
***
***
      yp1 = (yin(2)-yin(1))/(xin(2)-xin(1))
      ypn = 0.0

      call spline(xin,yin,nin,yp1,ypn,y2,u)
      call splint(nin,xin,yin,y2,nout,xout,yout)

      return                                              
      end

      subroutine splint (nspline,xspline,yspline,ysplin2,
     *                   nfit,xfit,yfit)

      implicit real*8 (a-h, o-z)
      dimension xspline(nspline), yspline(nspline), ysplin2(nspline),
     *          xfit(nfit), yfit(nfit)

      klo = 1
      do i=1, nfit
        do k=klo+1, nspline
            if(xspline(k).ge.xfit(i)) then
               if(xspline(k-1).le.xfit(i)) then
                  khi = k
                  klo = k-1
               else
                  if (k-1.eq.1 .and. i.eq.1) then
                     stop '  SPLINT: xfit(1) < xspline(1)'
                  else
                     stop '  SPLINT: xfit not properly ordered'
                  end if
               end if
               h= xspline(khi) - xspline(klo) 
               a= (xspline(khi)-xfit(i))/h
               b= (xfit(i)-xspline(klo))/h

               yfit(i) = a*yspline(klo) + b*yspline(khi)
     *              +( (a**3-a)*ysplin2(klo) +
     *              (b**3-b)*ysplin2(khi) ) *h*h/6
               go to 10
            end if
         end do
*        stop '  SPLINT: out of bounds'
* This is for the unlikely event that rmax exceed r(mesh)
	yfit(i)=0.0
 10     continue
      end do

      return
      end

      subroutine spline(x,y,n,yp1,ypn,y2,u)
      implicit real*8 (a-h, o-z)
      dimension x(n),y(n),y2(n),u(n)

      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     *      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      end

      subroutine deriv(f,df,r,rab,mesh,idim1)
      implicit none
      integer mesh,i,idim1
      real *8 f(idim1),df(idim1),r(idim1),rab(idim1)
      real *8 f1,f2,a,b,c
c.... calculates radial derivative of f

c .. first df/di
      do i=3,mesh-2
         f1= (f(i+1)-f(i-1))/2.
         f2= (f(i+2)-f(i-2))/4.
         df(i) = (4.*f1-f2)/3.
      enddo
       df(3) = (f(4)-f(2))/2.
c      df(2) =(f(3)-f(1))/2.
      f1=f(4)-f(3)
      f2=f(5)-f(4)
      a = (f2-f1)/2.
      b = (f1+f2)/2.-(f2-f1)*4.
      c = f(3)-a*3.*3.-b *3.
      df(2) = 2.*df(3)-df(4)
      df(1) =  2.*df(2)-df(3)
c      df(2) = a*4.+b*2.+c
c      df(1) = a*1.+b*1.+c
      df(mesh-1) = 2.*df(mesh-2)-df(mesh-3)
      df(mesh) =  2.*df(mesh-1)-df(mesh-2)

c now multiply di/dr     
      do i=1,mesh
         df(i) = df(i)/rab(i)
      enddo
      end
