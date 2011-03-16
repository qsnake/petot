      subroutine difrel(iter,iorb,v,ar,br,lmax,nr,a,b,r,rab,
     1 norb,no,lo,so,znuc,viod,viou,vid,viu,ev,rabkar,
     2 rabai,fa,fb,evi)
c
c  difrel integrates the relativistic Dirac equation
c  it finds the eigenvalue ev, the major and minor component
c  of the wavefunction, ar and br.  It uses an intial guess
c  for the eigenvalues from dsolv1
c
c  njtj  ***  modifications  ***
c    This routine has major modifications.
c    1)The data needed inside the loops has been calculated
c    outside the main loop(increases speed for non-opt
c    compiliers, i.e. dumb compiliers).  
c    2)The predict/correct values are placed in an array.    
c    Output is unchanged
c  njtj  ***  modifications  ***
c
c  njtj
c  ###  Cray conversions  
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch the 3 double precision parameter
c  ###      to single precision parameter statements. 
c  ###  Cray conversions
c  njtj  
c
c  njtj
c  &&&  Machine dependent Parameter 
c  &&&    The value of expzer is machine dependent.
c  &&&    The user must switch in the correct value for 
c  &&&    the machine in use from the list, or find 
c  &&&    it for their machine.
c  &&&  Machine dependent Parameter
c  njtj
c                
      implicit double precision (a-h,o-z)
c
      parameter (zero=0.D0,pnine=0.9D0,one=1.D0,ai=2*137.0360411D0)
      parameter (etol=-1.D-7)
Cray      parameter (zero=0.0,pnine=0.9,one=1.0,ai=2*137.0360411)
Cray      parameter (etol=-1.E-7)
c
c  Tolernce
c
      parameter (tol = 1.D-10,five=5.0D0)
Cray      parameter (tol = 1.E-10,five=5.0)
c
c  Integration coefficients
c
      parameter(abc1=190.1D0/72,abc2=-138.7D0/36,abc3=10.9D0/3,
     1 abc4=-63.7D0/36,abc5=25.1D0/72,amc0=25.1D0/72,amc1=32.3D0/36,
     2 amc2=-1.1D0/3,amc3=5.3D0/36,amc4=-1.9D0/72)
Cray      parameter(abc1=190.1/72,abc2=-138.7/36,abc3=10.9/3,
Cray     1 abc4=-63.7/36,abc5=25.1/72,amc0=25.1/72,amc1=32.3/36,
Cray     2 amc2=-1.1/3,amc3=5.3/36,amc4=-1.9/72)
c
      dimension v(nr),ar(nr),br(nr),r(nr),rab(nr),
     1 no(norb),lo(norb),so(norb),viod(lmax,nr),viou(lmax,nr),
     2 vid(nr),viu(nr),ev(norb),rabkar(nr),rabai(nr),
     3 fa(nr),fb(nr),evi(norb)
c
      dimension rs(5)
c
c------Machine dependent parameter-
c------Require exp(-2*expzer) to be within the range of the machine
c
      expzer = 3.7D2
cApollo      expzer = 3.7D2
cSun      expzer = 3.7D2
cVax      expzer = 44.D0
Cray      expzer = 2.8E3
c
      itmax = 100
      ai2 = ai * ai
      az = znuc/(2*ai)
      ka = lo(iorb)+1
      if (so(iorb) .lt. 0.1 .and. lo(iorb) .ne. 0) ka=-lo(iorb)
c
c  determine effective charge and vzero for startup of
c  outward integration
c  ar = r**s * (1  + a1 r + a2 r**2 + ... )
c  br = r**s * (b0 + b1 r + b2 r**2 + ... )
c  s = sqrt (ka**2 - az**2)    b0 = - az / (s + ka)
c  an = (az (v0 - e) a(n-1) - (s + n + ka) (v0 - e - ai**2) b(n-1))
c        / (n ai (2 s + n))
c  bn = ((v0 - e) a(n-1) - 2 znuc an ) / ( ai (s + n + ka))
c
      s = sqrt(ka*ka-az*az)
      if (ka .gt. 0) then
        b0 = -az/(s+ka)
      else
        b0 = (s-ka)/az
      endif
      if (so(iorb) .lt. 0.1) then
        vzero=vid(2)
      else
        vzero=viu(2)
      endif
c
c  njtj  ***  start major modification  ***  
c    Loop data calculated only once.
c    Set ar() and br() to zero.
c
      do 1 j=1,nr
        ar(j) = zero
        br(j) = zero
 1    continue
      do 3 j=2,nr
        rabkar(j)=rab(j)*ka/r(j)
 3    continue 
      do 4 j=2,nr
        rabai(j)=rab(j)/ai
 4    continue
      do 5 j=2,5
        rs(j)=r(j)**s
 5    continue
c
c  set the underflow trap, error from Berkeley version,
c  fixed by Troy Barbee, sqrt(expzer) should be expzer/2,
c  4/17/90.
c
      juflow=1
      do 42 j=2,nr
        if (s*abs(log(r(j))) .ge. expzer/2) juflow = j
 42   continue 
c  njtj *** end major modification  ***
c
      emax = zero
      emin = -one*100000
      if (ev(iorb) .gt. emax) ev(iorb) = emax 
 10   if (itmax .lt. 2) write(6,15) iorb,iter,ev(iorb),nodes
 15   format(' iorb =',i3,' iter =',i3,' ev =',e18.10,' nodes =',i2)
      if (itmax .eq. 0) return
      if (ev(iorb) .gt. zero) then
        write(6,1000)iorb
        call ext(620+iorb)
      endif  
 1000 format(//,' error in difrel - ev(',i2,
     1 ') greater then v(infinty)')
c
c  Find practical infinity ninf and classical turning
c  point nctp for orbital.
c
      icount=0
 20   icount=icount+1
      do 22 j=nr,2,-1
        temp = v(j) - ev(iorb)
        if (temp .lt. zero) temp = zero
        if (r(j)*sqrt(temp) .lt. expzer) goto 23
 22   continue
 23   ninf=j
      nctp = ninf - 5 
      do 25 j=2,ninf-5
        if (v(j) .lt. ev(iorb)) nctp = j
 25   continue 
      if (ev(iorb) .ge. etol*100) nctp=ninf-5
      if (ev(iorb) .ge. etol) ev(iorb)=zero
      if (evi(iorb) .ne. zero) then
        ev(iorb)=evi(iorb) 
        do 26 j=2,nr
          if (r(j) .lt. five) nctp=j
 26     continue
      endif               
      if (nctp .le. 6) then
        ev(iorb) = pnine*ev(iorb)
        if (icount .gt. 100) then
          write(6,1010)iorb
          call ext(650+iorb)
        endif
        goto 20
      endif  
 1010 format(//,'error in difrel - cannot find classical',
     1 /,'turning point in orbital ',i2)
c
c  Outward integration from 1 to nctp, startup.
c
      a1 = (az*(vzero-ev(iorb))-(s+1+ka)*(vzero-ev(iorb)-ai2)*b0)
     1   / (ai*(2*s+1))
      b1 = ((vzero-ev(iorb))-2*znuc*a1) / (ai*(s+1+ka))
      a2 = (az*(vzero-ev(iorb))*a1-(s+2+ka)*(vzero-ev(iorb)-ai2)*b1)
     1   / (2*ai*(2*s+2))
      b2 = ((vzero-ev(iorb))*a1-2*znuc*a2) / (ai*(s+2+ka))
      do 35 j=2,5
        ar(j) = rs(j) * (1 +(a1+a2*r(j))*r(j))
        br(j) = rs(j) * (b0+(b1+b2*r(j))*r(j))
 35   continue
      fa(1) = zero
      fb(1) = zero
      fa(2) = rabkar(2)*ar(2)+(ev(iorb)-v(2)+ai2)*br(2)*rabai(2)
      fb(2) = -rabkar(2)*br(2)-(ev(iorb)-v(2))*ar(2)*rabai(2)
      fa(3) = rabkar(3)*ar(3)+(ev(iorb)-v(3)+ai2)*br(3)*rabai(3)
      fb(3) = -rabkar(3)*br(3)-(ev(iorb)-v(3))*ar(3)*rabai(3)
      fa(4) = rabkar(4)*ar(4)+(ev(iorb)-v(4)+ai2)*br(4)*rabai(4)
      fb(4) = -rabkar(4)*br(4)-(ev(iorb)-v(4))*ar(4)*rabai(4)
      fa(5) = rabkar(5)*ar(5)+(ev(iorb)-v(5)+ai2)*br(5)*rabai(5)
      fb(5) = -rabkar(5)*br(5)-(ev(iorb)-v(5))*ar(5)*rabai(5)
c
c  Intergration loop.
c
      nodes = 0
      do 40 j=6,nctp
c
c  Predictor (Adams-Bashforth).
c
        evvai2=ev(iorb)-v(j)+ai2 
        evv=ev(iorb)-v(j)
        arp = ar(j-1) + abc1*fa(j-1)+abc2*fa(j-2)+abc3*fa(j-3)
     1   +abc4*fa(j-4)+abc5*fa(j-5)
        brp = br(j-1) + abc1*fb(j-1)+abc2*fb(j-2)+abc3*fb(j-3)
     1   +abc4*fb(j-4)+abc5*fb(j-5)
        fa(j) = rabkar(j)*arp+evvai2*brp*rabai(j)
        fb(j) = -rabkar(j)*brp-evv*arp*rabai(j)
c
c  Corrector (Adams-Moulton).
c
        arc = ar(j-1) + amc0*fa(j)+amc1*fa(j-1)+amc2*fa(j-2)
     1   +amc3*fa(j-3)+amc4*fa(j-4)
        brc = br(j-1) + amc0*fb(j)+amc1*fb(j-1)+amc2*fb(j-2)
     1   +amc3*fb(j-3)+amc4*fb(j-4)
        faj = rabkar(j)*arc+evvai2*brc*rabai(j)
        fbj = -rabkar(j)*brc-evv*arc*rabai(j)
c
c  Error reduction step.
c
        ar(j) = arc + amc0*(faj-fa(j))
        br(j) = brc + amc0*(fbj-fb(j))
        fa(j) = rabkar(j)*ar(j)+evvai2*br(j)*rabai(j)
        fb(j) = -rabkar(j)*br(j)-evv*ar(j)*rabai(j)
c
c  Count nodes - if no underflow.
c
        if(j.gt.juflow.and.ar(j)*ar(j-1).lt.zero)nodes=nodes+1
 40   continue
       arout = ar(nctp)
       arpout = fa(nctp)
c
c  End outward integration.
c  If number of nodes correct, start inward integration
c  else modify energy stepwise and try again.
c
      if (evi(iorb) .ne. zero) goto 111
      if (nodes .ne. no(iorb)-lo(iorb)-1) then
c
c  too many nodes decrease ev
c
        if (nodes .gt. no(iorb)-lo(iorb)-1) then
          if (ev(iorb) .lt. emax) emax = ev(iorb)
          ev(iorb) = ev(iorb) + ev(iorb)/10
c
c  too few nodes increase ev
c
        else
          if (ev(iorb) .gt. emin) emin = ev(iorb)
          ev(iorb) = ev(iorb) - ev(iorb)/10
        endif
        itmax = itmax-1
        goto 10
      endif
c
c  Inward integration from ninf to nctp startup.
c
      do 70 j=ninf,ninf-4,-1
        alf = v(j) - ev(iorb)
        if (alf .lt. zero) alf = zero
        alf = sqrt(alf)
        ar(j) = exp(-alf*r(j))
        br(j) = ai*(alf+ka/r(j))*ar(j)/(v(j)-ev(iorb)-ai2)
 70   continue
      fa(ninf) = rabkar(ninf)*ar(ninf)+
     1    (ev(iorb)-v(ninf)+ai2)*br(ninf)*rabai(ninf)
      fb(ninf) = -rabkar(ninf)*br(ninf)
     1    -(ev(iorb)-v(ninf))*ar(ninf)*rabai(ninf)
      fa(ninf-1) = rabkar(ninf-1)*ar(ninf-1)+
     1    (ev(iorb)-v(ninf-1)+ai2)*br(ninf-1)*rabai(ninf-1)
      fb(ninf-1) = -rabkar(ninf-1)*br(ninf-1)
     1    -(ev(iorb)-v(ninf-1))*ar(ninf-1)*rabai(ninf-1)
      fa(ninf-2) = rabkar(ninf-2)*ar(ninf-2)
     1    +(ev(iorb)-v(ninf-2)+ai2)*br(ninf-2)*rabai(ninf-2)
      fb(ninf-2) = -rabkar(ninf-2)*br(ninf-2)
     1    -(ev(iorb)-v(ninf-2))*ar(ninf-2)*rabai(ninf-2)
      fa(ninf-3) = rabkar(ninf-3)*ar(ninf-3)
     1    +(ev(iorb)-v(ninf-3)+ai2)*br(ninf-3)*rabai(ninf-3)
      fb(ninf-3) = -rabkar(ninf-3)*br(ninf-3)
     1    -(ev(iorb)-v(ninf-3))*ar(ninf-3)*rabai(ninf-3)
      fa(ninf-4) = rabkar(ninf-4)*ar(ninf-4)
     1    +(ev(iorb)-v(ninf-4)+ai2)*br(ninf-4)*rabai(ninf-4)
      fb(ninf-4) = -rabkar(ninf-4)*br(ninf-4)
     1    -(ev(iorb)-v(ninf-4))*ar(ninf-4)*rabai(ninf-4)
c
c  Integration loop.
c
      istop = ninf-nctp
      if (istop .lt. 5) goto 222
      do 80 j=ninf-5,nctp,-1
c
c  Predictor (Adams-Bashforth).
c
        evvai2=ev(iorb)-v(j)+ai2 
        evv=ev(iorb)-v(j)
        arp = ar(j+1)-(abc1*fa(j+1)+abc2*fa(j+2)+abc3*fa(j+3)
     1   +abc4*fa(j+4)+abc5*fa(j+5))
        brp = br(j+1)-(abc1*fb(j+1)+abc2*fb(j+2)+abc3*fb(j+3)
     1   +abc4*fb(j+4)+abc5*fb(j+5))
        fa(j) = rabkar(j)*arp+evvai2*brp*rabai(j)
        fb(j) = -rabkar(j)*brp-evv*arp*rabai(j)
c
c  Corrector (Adams-Moulton).
c
        arc = ar(j+1)-(amc0*fa(j)+amc1*fa(j+1)+amc2*fa(j+2)
     1   +amc3*fa(j+3)+amc4*fa(j+4))
        brc = br(j+1)-(amc0*fb(j)+amc1*fb(j+1)+amc2*fb(j+2)
     1   +amc3*fb(j+3)+amc4*fb(j+4))
        faj = rabkar(j)*arc+evvai2*brc*rabai(j)
        fbj = -rabkar(j)*brc-evv*arc*rabai(j)
c
c  Error reduction step.
c
        ar(j) = arc + amc0*(faj-fa(j))
        br(j) = brc + amc0*(fbj-fb(j))
        fa(j) = rabkar(j)*ar(j)+evvai2*br(j)*rabai(j)
        fb(j) = -rabkar(j)*br(j)-evv*ar(j)*rabai(j)
 80   continue
 222  arin = ar(nctp)
      arpin = fa(nctp)
c
c  End inward integration
c  Rescale ar and br outside nctp to match ar(nctp) from
c  outward integration.
c
      factor = arout/arin
      do 90 j=nctp,ninf
        ar(j) = factor * ar(j)
        br(j) = factor * br(j)
 90   continue
      arpin = factor * arpin
c
c  Find the normalizing factor.
c
      factor = zero
      ll = 4
      do 100 j=2,ninf
        factor = factor + ll*(ar(j)*ar(j)+br(j)*br(j))*rab(j)
        ll = 6 - ll
 100  continue
      factor = factor / 3
c
c  Modify the eigenvalue ev.
c
      dev = arout * (arpout-arpin) / (factor * rab(nctp))
      if (5*abs(dev) .gt. -ev(iorb)) dev=sign(ev(iorb),dev)/5
      itmax = itmax-1
      evold = ev(iorb)
      ev(iorb) = ev(iorb) + dev
      if (ev(iorb) .gt. emax) then
        ev(iorb) = (evold + emax) / 2
      elseif (ev(iorb) .lt. emin) then
        ev(iorb) = (evold + emin) / 2 
      endif
      if (abs(dev) .gt. tol*(1-ev(iorb))) goto 10
c
c  Normalize the wavefunction.
c
      factor = 1 / sqrt(factor)
      do 110 j=1,ninf
        ar(j) = factor*ar(j)
        br(j) = factor*br(j)
 110  continue
 111  continue
      if (evi(iorb) .ne. zero) then
        factor = zero
        ll = 4
        do 112 j=2,nctp
          factor = factor + ll*(ar(j)*ar(j)+br(j)*br(j))*rab(j)
          ll = 6 - ll
 112    continue
        factor = factor / 3
        factor = 1 / sqrt(factor)
        do 113 j=1,nctp
          ar(j) = factor*ar(j)
          br(j) = factor*br(j)
 113    continue
      endif
      return
      end
