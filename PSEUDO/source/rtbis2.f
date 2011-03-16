      subroutine rtbis2(x1,x2,rc1,rc2,rc3,rc4,rc5,rc6,rc7,
     1 rc8,lp,arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,
     2 gamma,alpha,alpha1,alpha2,alpha3,alpha4,ar)
c
c *************************************************************
c *  njtj
c *  Finds the value of gamma for the v"(0)=0 criteria.
c *  The method used is bisection.  This routine
c *  was taken from Numerical Recipes, page 247.
c *  njtj
c *************************************************************
c
c  njtj
c  ###  Cray conversions  
c  ###    1)Comment out the implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      parameter (jmax=80,pfive=0.5D0,zero=0.D0,xacc=1.D-10) 
Cray      parameter (jmax=80,pfive=0.5,zero=0.0,xacc=1.E-10) 
c
      dimension r(jrc),rab(jrc),ar(jrc)
c
      call gamfn2(rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,lp,
     1 arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,x1,
     2 alpha,alpha1,alpha2,alpha3,alpha4,f,ar)
      call gamfn2(rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,lp,
     1 arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,x2,
     2 alpha,alpha1,alpha2,alpha3,alpha4,fmid,ar)
      if(f*fmid.ge.zero) then
        write(6,4000)
        call ext(840+lp)
      endif
      if(f.lt.zero)then
        gamma=x1
        dx=x2-x1
      else
        gamma=x2
        dx=x1-x2
      endif
      do 11 j=1,jmax
        dx=dx*pfive
        xmid=gamma+dx
        call gamfn2(rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,lp,
     1   arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,
     2   xmid,alpha,alpha1,alpha2,alpha3,alpha4,fmid,ar)
        if(fmid.lt.zero)gamma=xmid
        if(abs(dx).lt.xacc .or. fmid.eq. zero) return
11    continue
      write(6,4001)
      call ext(850+lp)
 4000 format(' error in bisection method(rtbistk)',
     1 ' - root must be bracketed.',
     2 /,'a b o r t i n g   p r o g r a m')
 4001 format(' error in bisection method(rtbistk)',
     1 ' - too many bisections used',
     2 /,'a b o r t i n g   p r o g r a m')
      end
