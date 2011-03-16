      subroutine zrbac2(x1,x2,rc1,rc2,rc3,rc4,rc5,rc6,rc7,
     1 rc8,lp,arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,
     2 gamma,alpha,alpha1,alpha2,alpha3,alpha4,ar)
c
c **********************************************************
c *  njtj
c *    Routine brackets the root of the given function.
c *    Taken from Numerical Recipes page 245.
c *  njtj
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
      parameter (factor=1.6D0,ntry=50) 
Cray      parameter (factor=1.6,ntry=50) 
c
      dimension r(jrc),rab(jrc),ar(jrc)
c
      call gamfn2(rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,lp,
     1 arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,x1,
     2 alpha,alpha1,alpha2,alpha3,alpha4,f1,ar)
      call gamfn2(rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,lp,
     1 arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,x2,
     2 alpha,alpha1,alpha2,alpha3,alpha4,f2,ar)
c
      do 11 j=1,ntry
        if(f1*f2.lt.0.0)return
        if(abs(f1).lt.abs(f2))then
          x1=x1+factor*(x1-x2)
          call gamfn2(rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,lp,
     1     arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,x1,
     2     alpha,alpha1,alpha2,alpha3,alpha4,f1,ar)
        else
          x2=x2+factor*(x2-x1)
          call gamfn2(rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,lp,
     1     arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,x2,
     2     alpha,alpha1,alpha2,alpha3,alpha4,f2,ar)
        endif
11    continue
c
c  failure, abort program
c
      write(6,1000)lp
      call ext(830+lp)
 1000 format(//,'error in zbractk - can not bracket orbital ',i2)
      return
      end
