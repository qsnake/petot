      subroutine potran(i,vd,r,nr,zion,a,b,c)
c
c ***********************************************************
c *                                                         *
c *    This is a plotting routine; the user should adjust   *
c *  for their own needs.  The potential is fitted with a   *
c *  second degree polynomial, which is muliplied with the  *
c *  appropriate functions and then integrated by parts     *
c *  to find the fourier transform.  The result is then     *
c *  printed to the current plot.dat file (unit=3) for      *
c *  later plotting.  A marker(marker fn#) is placed at     *
c *  the end of each set of data.                           *
c *                                                         *
c ***********************************************************
c
      implicit real*8 (a-h,o-z)
c
      parameter (zero=0.0d0,one=1.0d0)
c
      dimension vd(nr),r(nr),a(nr),b(nr),c(nr),vql(100)
c
c  The potential times r is fitted to the polynominal
c  a + bx + cx^2 at every other point.
c
      rm=zero
      vm=2*zion
      do 130 k=2,nr,2
        r0=r(k)
        v0=r0*vd(k)+2*zion
        rp=r(k+1)
        vp=rp*vd(k+1)+2*zion                       
        d1=1/((rp-rm)*(r0-rm))
        d2=1/((rp-r0)*(rm-r0))
        d3=1/((r0-rp)*(rm-rp))
        a(k)=vm*d1+v0*d2+vp*d3
        b(k)=-vm*(r0+rp)*d1-v0*(rm+rp)*d2-vp*(rm+r0)*d3
        c(k)=vm*r0*rp*d1+v0*rm*rp*d2+vp*rm*r0*d3
        rm=rp
        vm=vp
 130  continue
c
c  Find the fourier transform q^2/4pi/zion*vql. Everything is 
c  rescaled  by zion.   
c          
      do 150 j=1,94
        q=one/4*j
        q2=q*q
        vql(j)=zero
        rm=zero
        do 140 k=2,nr-1,2
          rp=r(k+1)
          vql(j)=vql(j)+(2*a(k)*rp+b(k))/q*sin(q*rp)
     1     -((a(k)*rp+b(k))*rp+c(k)-2*a(k)/q2)*cos(q*rp)
     2     -(2*a(k)*rm+b(k))/q*sin(q*rm)
     3     +((a(k)*rm+b(k))*rm+c(k)-2*a(k)/q2)*cos(q*rm)
          rm=rp
 140    continue
        vql(j)=vql(j)/2/zion-one
 150  continue        
c
c  Print out the transforms( really q^2/(4pi*zion)*v(q) ) to 
c  the current plot.dat file (unit=3) for latter plotting.
c
      do 170 j=1,48
        write(3,6000)one/4*j,vql(j)
 170  continue
      write(3,6008)i
      return
c
c  format statements
c
 6000 format(1x,f7.4,3x,f10.6)
 6008 format(1x,'marker fn',i1)
c
      end
