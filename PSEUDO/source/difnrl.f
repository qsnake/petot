      subroutine difnrl(iter,iorb,v,ar,br,lmax,
     1 nr,a,b,r,rab,norb,no,lo,so,znuc,viod,viou,
     2 vid,viu,ev,iflag,rab2,fa,fb,evi)
c
c    difnrl integrates the Schroedinger equation
c    if finds the eigenvalue ev, the wavefunction ar
c    and the derivative br = d(ar)/dr
c
c  njtj  ***  modifications  ***
c    This routine has had major modifications.  Some
c    of the data used inside the main loop has been 
c    calculated outside the main loop to reduce the number
c    of operations(uses extra array space to gain speed)
c    and are passed as work arrays form the main.
c    The predictor-corrector functions have been put 
c    into a array.
c    The iflag variable was added to indicate nonconvergence
c    for other programs.  It has no use in the atom program
c    and can be removed by the user.
c    All output from the routine is compatible to
c    the Berkeley/Sverre Froyen version.
c  njtj  ***  modifications  ***
c
c  njtj
c  ###  Cray conversions  
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch the double precision parameter statements
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
      parameter(zero=0.D0,pnine=0.9D0,two=2.D0,etol=-1.D-7)
Cray      parameter(zero=0.0,pnine=0.9,two=2.0,etol=-1.E-7)
c  
c  Tolerence
c
      parameter(tol=1.D-10,five=5.0D0)
Cray      parameter(tol=1.E-10,five=5.0)
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
      dimension v(nr),ar(nr),br(nr),r(nr),rab(nr),no(norb),
     1 lo(norb),so(norb),viod(lmax,nr),viou(lmax,nr),
     2 vid(nr),viu(nr),ev(norb),evi(norb)
c
c  njtj  *** start modification  ***
c    Arrays added to gain speed.
c
      dimension rabrlo(5),rlp(5),rab2(nr),fa(nr),fb(nr)
c
c  njtj  ***  end modification  ***
c
c------Machine dependent parameter-
c------Require exp(-2*expzer) to be within the range of the machine
c
      expzer = 3.7D2
cApollo      expzer = 3.7D2
cSun      expzer = 3.7D2 
cVax      expzer = 44.D0
Cray      expzer =  2.8E3
c
c  njtj  *** major modification start  ***
c    Loop data calculated outside loop to gain speed.
c
      itmax = 100         
      iflag = 0
      lp = lo(iorb)+1
      ar(1) = zero 
      if (lo(iorb) .eq. 0) then
        br(1) = b*a
      else
        br(1) = zero
      endif
      do 1 j=2,nr
        ar(j) = zero
 1    continue
      do 2 j=2,nr
        br(j) =zero
 2    continue  
      do 4 j=2,5
        rlp(j)=r(j)**lp
 4    continue                       
      do 5 j=2,5
        rabrlo(j)=rab(j)*r(j)**lo(iorb)
 5    continue
      do 6 j=1,nr
        rab2(j)=rab(j)*rab(j)
 6    continue
c
c   set underflow trap, error from Berkeley version,
c   fixed by Troy Barbee sqrt(expzer) should be expzer/2 
c   4/17/90
c
      juflow=1
      do 42 j=2,nr
        if (lp*abs(log(r(j))) .ge. expzer/2) juflow = j
 42   continue
c
c  njtj  *** end major modification  ***
c
c   determine effective charge and vzero for startup of
c   outward integration
c   ar = r**(l+1) * (1 + aa r + bb r**2 + ... )
c   aa = -znuc / lp     bb = (-2 znuc aa + v(0) - e)/(4 l + 6)
c
      zeff = zero
      if (so(iorb) .lt. 0.1 .and. viod(lp,2) .lt. -0.1) zeff=znuc
      if (so(iorb) .gt. 0.1 .and. viou(lp,2) .lt. -0.1) zeff=znuc
      aa = -zeff/lp
      vzero = -2*zeff*aa
      if (zeff .eq. zero) then
        if (so(iorb) .lt. 0.1 ) then
          vzero=vzero+viod(lp,2)/r(2)
        else
          vzero=vzero+viou(lp,2)/r(2)
        endif
      endif
      if (so(iorb) .lt. 0.1) then
        vzero=vzero+vid(2)
      else
        vzero=vzero+viu(2)
      endif
      var0 = zero
      if (lo(iorb) .eq. 0) var0=-2*zeff
      if (lo(iorb) .eq. 1) var0=two
c
      emax = zero
      emin = -two*100000
      if (ev(iorb) .gt. emax) ev(iorb) = emax
 10   if (itmax .lt. 2) write(6,15) iorb,iter,ev(iorb),nodes
 15   format(' iorb =',i3,' iter =',i3,' ev =',e18.10,' nodes =',i2)
      if (itmax .eq. 0) then
        iflag =1
        return
      endif
      if (ev(iorb) .gt. zero) then
        write(6,1000)iorb
        call ext(620+iorb)
      endif  
 1000 format(//,' error in difnrl - ev(',i2,
     1 ') greater then v(infinty)')
c
c   find practical infinity ninf and classical turning
c   point nctp for orbital
c
      icount=0
 20   icount=icount+1
      do 22 j=nr,2,-1
        temp = v(j) -ev(iorb)
        if (temp .lt. zero) temp = zero
        if (r(j)*sqrt(temp) .lt. expzer) goto 23
 22   continue
 23   ninf=j
      nctp = ninf - 5 
      do 25 j=2,ninf-5
        if (v(j) .lt. ev(iorb)) nctp = j
 25   continue                
      if (ev(iorb) .ge. etol*10) nctp=ninf-5
      if (ev(iorb) .ge. etol) ev(iorb)=zero
      if (evi(iorb) .ne. zero) then
        ev(iorb) = evi(iorb)
        do 26 j=1,nr
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
 1010 format(//,'error in difnrl - cannot find the classical '
     1 ,/' turning point for orbital ',i2)
c
c   outward integration from 1 to nctp
c   startup
c
      bb = (vzero-ev(iorb))/(4*lp+2)
      do 35 j=2,5
        ar(j) = rlp(j) * (1+(aa+bb*r(j))*r(j))
        br(j) = rabrlo(j) * (lp+(aa*(lp+1)+bb*(lp+2)*r(j))*r(j))
 35   continue
c
c  njtj  ***  start major modification  ***
c    Predictor-corrector array added.
c
      fa(1) = br(1)
      fb(1) = b*br(1) + rab2(1)*var0
      fa(2) = br(2)
      fb(2) = b*br(2) + rab2(2)*(v(2)-ev(iorb))*ar(2)
      fa(3) = br(3)
      fb(3) = b*br(3) + rab2(3)*(v(3)-ev(iorb))*ar(3)
      fa(4) = br(4)
      fb(4) = b*br(4) + rab2(4)*(v(4)-ev(iorb))*ar(4)
      fa(5) = br(5)
      fb(5) = b*br(5) + rab2(5)*(v(5)-ev(iorb))*ar(5)
c
c   intergration loop
c
      nodes = 0
      do 40 j=6,nctp 
c
c   predictor (Adams-Bashforth)
c                                                             
        j1=j-1
        j2=j-2
        j3=j-3
        j4=j-4
        j5=j-5
        vev=v(j)-ev(iorb)
        arp = ar(j1) + abc1*fa(j1)+abc2*fa(j2)+abc3*fa(j3)+
     1   abc4*fa(j4)+abc5*fa(j5)
        brp = br(j1) + abc1*fb(j1)+abc2*fb(j2)+abc3*fb(j3)+
     1   abc4*fb(j4)+abc5*fb(j5)
        fb1 = b*brp + rab2(j)*vev*arp
c
c   corrector (Adams-Moulton)
c
        arc = ar(j1) + amc0*brp+amc1*fa(j1)+amc2*fa(j2)+
     1   amc3*fa(j3)+amc4*fa(j4)
        brc = br(j1) + amc0*fb1+amc1*fb(j1)+amc2*fb(j2)+
     1   amc3*fb(j3)+amc4*fb(j4)
        fb0 = b*brc + rab2(j)*vev*arc
c
c   error reduction step
c
        ar(j) = arc + amc0*(brc-brp)
        br(j) = brc + amc0*(fb0-fb1)
        fa(j) = br(j)
        fb(j) = b*br(j) + rab2(j)*vev*ar(j)
c
c   count nodes - if no underflow
c
        if(j.gt.juflow.and.ar(j)*ar(j-1).lt.zero)nodes=nodes+1
 40   continue
c
c  njtj  ***  end major modification  ***
c
      arctp = ar(nctp)
      brctp = br(nctp)
c
c   end outward integration
c
c   if number of nodes correct, start inward integration
c   else modify energy stepwise and try again
c
      if (evi(iorb) .ne. zero) goto 111
      if (nodes .ne. no(iorb)-lo(iorb)-1) then
        if (nodes .lt. no(iorb)-lo(iorb)-1) then
c
c  too few nodes; increase ev
c
          if (ev(iorb) .gt. emin) emin = ev(iorb)
          ev(iorb) = ev(iorb) - ev(iorb)/10
        else
c
c  too many nodes; decrease ev
c
          if (ev(iorb) .lt. emax) emax = ev(iorb)
          ev(iorb) = ev(iorb) + ev(iorb)/10
        endif
        itmax = itmax-1
        goto 10
      endif
c
c   inward integration from ninf to nctp
c   startup
c                               
      do 71 j=ninf,ninf-4,-1
        alf = v(j) - ev(iorb)
        if (alf .lt. zero) alf = zero
        alf = sqrt(alf)
        ar(j) = exp(-alf*r(j))
        br(j) = -rab(j)*alf*ar(j)
 71   continue
c
c  njtj  ***  start major modification  ***
c    Array for predictor-corrector added.
c
      fa(ninf) = br(ninf)
      fb(ninf) = b*br(ninf) + rab2(ninf)*
     1 (v(ninf)-ev(iorb))*ar(ninf)
      ninf1 = ninf - 1
      fa(ninf1) = br(ninf1)
      fb(ninf1) = b*br(ninf1) + rab2(ninf1)*
     1       (v(ninf1)-ev(iorb))*ar(ninf1)
      ninf2 = ninf - 2
      fa(ninf2) = br(ninf2)
      fb(ninf2) = b*br(ninf2) + rab2(ninf2)*
     1       (v(ninf2)-ev(iorb))*ar(ninf2)
      ninf3 = ninf - 3
      fa(ninf3) = br(ninf3)
      fb(ninf3) = b*br(ninf3) + rab2(ninf3)*
     1       (v(ninf3)-ev(iorb))*ar(ninf3)
      ninf4 = ninf - 4
      fa(ninf4) = br(ninf4)
      fb(ninf4) = b*br(ninf4) + rab2(ninf4)*
     1       (v(ninf4)-ev(iorb))*ar(ninf4)
c
c   integration loop
c
      istop = ninf - nctp
      if (istop .lt. 5) goto 222
      do 80 j=ninf-5,nctp,-1
c
c   predictor (Adams-Bashforth)
c
        j1 = j + 1
        j2 = j + 2
        j3 = j + 3
        j4 = j + 4
        j5 = j + 5  
        vev = v(j)-ev(iorb)
        arp = ar(j1) - (abc1*fa(j1)+abc2*fa(j2)+abc3*fa(j3)+
     1   abc4*fa(j4)+abc5*fa(j5))
        brp = br(j1) - (abc1*fb(j1)+abc2*fb(j2)+abc3*fb(j3)+
     1   abc4*fb(j4)+abc5*fb(j5))
        fb0 = b*brp + rab2(j)*vev*arp
c
c   corrector (Adams-Moulton)
c
        arc = ar(j1) - (amc0*brp+amc1*fa(j1)+amc2*fa(j2)+
     1   amc3*fa(j3)+amc4*fa(j4))
        brc = br(j1) - (amc0*fb0+amc1*fb(j1)+amc2*fb(j2)+
     1   amc3*fb(j3)+amc4*fb(j4))
c
        fb1 = b*brc + rab2(j)*vev*arc
c
c   error reduction step
c
        ar(j) = arc - amc0*(brc-brp)
        br(j) = brc - amc0*(fb1-fb0)
        fa(j) = br(j)
        fb(j) = b*br(j) + rab2(j)*vev*ar(j)
 80   continue
c
c   end inward integration
c
c  njtj  *** end major modification  ***
c
c   rescale ar and br outside nctp to match ar(nctp) from
c   outward integration
c
  222 factor = arctp/ar(nctp)
      do 90 j=nctp,ninf
        ar(j) = factor * ar(j)
        br(j) = factor * br(j)
 90   continue
c
c   find normalizing factor
c
      factor = zero
      ll = 4
      do 100 j=2,ninf
        factor = factor + ll*ar(j)*ar(j)*rab(j)
        ll = 6 - ll
 100  continue
      factor = factor / 3
c
c   modify eigenvalue ev
c
      dev = arctp * (brctp-br(nctp)) / (factor * rab(nctp))
      if (5*abs(dev) .gt. -ev(iorb)) dev=sign(ev(iorb),dev)/5
      itmax = itmax-1
      evold = ev(iorb)
      ev(iorb) = ev(iorb) + dev
      if (ev(iorb) .gt. emax) ev(iorb) = (evold + emax) / 2
      if (ev(iorb) .lt. emin) ev(iorb) = (evold + emin) / 2
      if (abs(dev) .gt. tol*(1-ev(iorb))) goto 10
c
c   normalize wavefunction and change br from d(ar)/dj to d(ar)/dr
c
      factor = 1 / sqrt(factor)
      do 110 j=1,ninf
        ar(j) = factor*ar(j)
        br(j) = factor*br(j) / rab(j)
 110  continue
 111  continue
      if (evi(iorb) .ne. zero) then
        factor = zero
        ll = 4
        do 112 j=2,nctp
          factor = factor + ll*ar(j)*ar(j)*rab(j)
          ll = 6 - ll
 112    continue
        factor = factor / 3
        factor = 1 / sqrt(factor)
        do 113 j=1,nctp
          ar(j) = factor*ar(j)
          br(j) = factor*br(j) / rab(j)
 113    continue
      endif
      return            
      end
