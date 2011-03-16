      subroutine gamfn2(rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,lp,
     1 arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,gamma,
     2 alpha,alpha1,alpha2,alpha3,alpha4,v0pp,ar)
c
c *********************************************************
c *                                                       *
c *  njtj                                                 *
c *   Retuns the values of delta, alpha, alpha1, alpha2,  *
c *   alpha3, and alpha4 given a fixed value of gamma.    *
c *   Returns V"(0) for the braketing and bisection       *
c *   routines.  Subroutine used in pseudtk routine.      *
c *  njtj                                                 *
c *                                                       *
c *********************************************************
c
      implicit real*8 (a-h,o-z)
c
      dimension r(jrc),rab(jrc),aj(5,5),bj(5),ar(jrc)  
c
      parameter (zero=0.0d0,pfive=0.5d0,one=1.0d0,errmin=1.d-12)
c
      rc9  = rc8*rc1
      rc10 = rc8*rc2
      rc11 = rc8*rc3
      rc12 = rc8*rc4
      delta=zero
      bj(1)=log(arc/rc1**lp)-gamma*rc2
      bj1=bj(1)
      bj(2)=brc-lp/rc1-2*gamma*rc1 
      bj2a=bj(2)+2*gamma*rc1
      bj2=bj(2)
      bj(3)=vrc-ev-2*lp/rc1*bj2a-bj2a**2-2*gamma
      bj3=bj(3)
      bj3a=bj(3)+2*gamma
      bj(4)=vap+2*lp/rc2*bj2a-2*lp/rc1*bj3a-2*bj2a*bj3a
      bj4=bj(4)
      bj(5)=vapp-4*lp/rc3*bj2a+4*lp/rc2*bj3a-2*lp/rc1*bj4-2*bj3a**2
     1 -2*bj2a*bj4
      bj5=bj(5)
      aj(1,1)=rc4
      aj(1,2)=rc6
      aj(1,3)=rc8
      aj(1,4)=rc10 
      aj(1,5)=rc12
      aj(2,1)=4*rc3
      aj(2,2)=6*rc5
      aj(2,3)=8*rc7 
      aj(2,4)=10*rc9
      aj(2,5)=12*rc11
      aj(3,1)=12*rc2
      aj(3,2)=30*rc4
      aj(3,3)=56*rc6
      aj(3,4)=90*rc8
      aj(3,5)=132*rc10
      aj(4,1)=24*rc1
      aj(4,2)=120*rc3
      aj(4,3)=336*rc5
      aj(4,4)=720*rc7
      aj(4,5)=1320*rc9
      aj(5,1)=24*one
      aj(5,2)=360*rc2
      aj(5,3)=1680*rc4
      aj(5,4)=5040*rc6
      aj(5,5)=11880*rc8
      call gaussj(aj,5,5,bj,1,1)
      alpha=bj(1)
      alpha1=bj(2)
      alpha2=bj(3)
      alpha3=bj(4) 
      alpha4=bj(5)
c
c   start iteration loop to find delta(with gamma fixed)
c
      do 550 j=1,100
c
c   generate pseudo wavefunction-note missing factor exp(delta) 
c
        do 560 k=1,jrc
          rp=r(k)
          r2=rp*rp
          polyr = r2*(((((alpha4*r2+alpha3)*r2+alpha2)*r2+
     1     alpha1)*r2+ alpha)*r2+gamma)
          ar(k) = rp**lp * exp(polyr)
 560    continue
c
c   integrate pseudo charge density from r = 0 to rc
c
        ll = 2
        cdps = - ar(jrc) * ar(jrc) * rab(jrc)
        if (jrc .ne. 2*(jrc/2)) then
          do 120 k=jrc,1,-1
            cdps = cdps +  ll * ar(k) * ar(k) * rab(k)
            ll = 6 - ll
 120      continue
        else
          do 121 k=jrc,4,-1
            cdps = cdps +  ll * ar(k) * ar(k) * rab(k)
            ll = 6 - ll
 121      continue
          cdps = cdps - ar(4) * ar(4) * rab(4)
          cdps = cdps + 9 * ( ar(1) * ar(1) * rab(1) +
     1     3 * ar(2) *ar(2) * rab(2) + 
     2     3 * ar(3) *ar(3) * rab(3) +
     3     ar(4) * ar(4) * rab(4))/8
        endif
        cdps = cdps/3
c
c   Calculate new delta(with gamma fixed), uses false position
c
        fdnew = log(cdrc/cdps) - 2*delta
        if (abs(fdnew) .lt. errmin) then
          v0pp=8*((2*one*(lp-one)+5*one)*alpha+gamma**2)
          return
        endif
        if (j .eq. 1) then
          ddelta=-pfive
        else 
          ddelta = - fdnew * ddelta / (fdnew-fdold)
        endif
        delta = delta + ddelta    
        bj(1)=bj1-delta
        bj(2)=bj2
        bj(3)=bj3
        bj(4)=bj4
        bj(5)=bj5
        aj(1,1)=rc4
        aj(1,2)=rc6
        aj(1,3)=rc8
        aj(1,4)=rc10 
        aj(1,5)=rc12
        aj(2,1)=4*rc3
        aj(2,2)=6*rc5
        aj(2,3)=8*rc7 
        aj(2,4)=10*rc9
        aj(2,5)=12*rc11
        aj(3,1)=12*rc2
        aj(3,2)=30*rc4
        aj(3,3)=56*rc6
        aj(3,4)=90*rc8
        aj(3,5)=132*rc10
        aj(4,1)=24*rc1
        aj(4,2)=120*rc3
        aj(4,3)=336*rc5
        aj(4,4)=720*rc7
        aj(4,5)=1320*rc9
        aj(5,1)=24*one
        aj(5,2)=360*rc2
        aj(5,3)=1680*rc4
        aj(5,4)=5040*rc6
        aj(5,5)=11880*rc8
        call gaussj(aj,5,5,bj,1,1)
        alpha=bj(1)
        alpha1=bj(2)
        alpha2=bj(3)
        alpha3=bj(4) 
        alpha4=bj(5)
        fdold = fdnew
 550  continue
      write(6,1000)
 1000 format(//, 'error in gamfind - delta not found')
      call ext(860+lp) 
      end
