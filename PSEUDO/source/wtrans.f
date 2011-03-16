      subroutine wtrans(vd,r,nr,rab,l,ist,b)
c
c **********************************************************
c *    
c *    This is a plotting routine; the user should adjust
c *  for their own needs.  The result
c *  is then printed to the current plot.dat file (unit=3)
c *  for later plotting of the data.  A marker (marker fw#)
c *  is placed at the end of each set of data.
c *   
c **********************************************************
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
      parameter (zero=0.D0,one=1.D0,big=17280.0D0,p5=.05D0)
Cray      parameter (zero=0.0,one=1.0,big=17280.0,p5=.05)
c
      dimension vd(nr),r(nr),rab(nr),b(nr),vql(48),vql2(48),
     1 a(2000),vdpp(2000),r2(2000),v(2000),w(4000)
c
      do 1 i=1,48
        vql(i)=zero
 1    continue
c
c  The wavefuncion(rR) times r times rab.
c
      if (abs(ist) .eq. 2) goto 400
      pi4=16*atan(one)
      do 10 k=2,nr
        if (r(k)-r(k-1).gt. p5) then
          nr2=k
          goto 20
        endif
 10   continue
 20   nr2=7*(nr2/7)+1
      nr3=nr2-7
      do 130 k=2,nr2
        b(k)=vd(k)*r(k)*rab(k)
 130  continue
      do 150 k=nr2,nr
        a(k-nr2+1)=vd(k)*r(k)
 150  continue 
      isx = 0
      a1 = -p5*10
      an = -p5*10
      b1 = zero
      bn = zero
      nrm=nr-nr2+1
      call splift(r(nr2),a,r2,vdpp,nrm,w,ierr,isx,a1,b1,an,bn)
      if(ierr.ne.1) then
        stop
      endif
      nr4=0
      do 155 ak=r(nr2),100.0D0,0.05D0
        nr4=nr4+1
        r2(nr4)=ak
 155  continue
      call splint(r(nr2),a,vdpp,nrm,r2,v,w,w(2000),nr4,kerr) 
c
c  Find the fourier transform-vql. 
c          
      do 140 j=1,48
        q=one/4*j
        vql(j)=zero
        a(1)=zero
        do 135 k=2,nr2
          a(k)=b(k)*sbessj(l,q*r(k))
 135    continue
c
c  Due to the high number of occilations in the intagrand,
c  an eight point Newton-Cotes intagration method is used.
c  See  Abramowitz and Stegun Eq. 25.4.17
c
        do 145 k=1,nr3,7
          vql(j)=vql(j)+751*(a(k)+a(k+7))+3577*(a(k+1)+a(k+6))+
     1     1323*(a(k+2)+a(k+5))+2989*(a(k+3)+a(k+4))
 145    continue
        vql(j)=pi4*7*vql(j)/big
        do 160 k=1,nr4
          a(k)=v(k)*sbessj(l,q*r2(k))
 160    continue
        vql2(j)=zero
        do 165 kk=8,nr4,7
          k=kk-7
          vql2(j)=vql2(j)+751*(a(k)+a(k+7))+3577*(a(k+1)+
     1     a(k+6))+1323*(a(k+2)+a(k+5))+2989*(a(k+3)+a(k+4))
 165    continue 
        vql2(j)=0.35D0*pi4*vql2(j)/big 
        vql(j)=vql(j)+vql2(j)
 140  continue
c
c  Print out the transform vql(q) to the current plot.dat 
c  file (unit=3) for latter plotting.
c
 400  do 170 j=1,48
        write(3,6000)one/4*j,ist*vql(j)
 170  continue
      write(3,6001)l
      return
c
c  format statements
c
 6000 format(1x,f7.4,3x,f10.6)
 6001 format(1x,'marker fw',i1)
      end
