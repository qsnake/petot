      subroutine pseud2(itype,icorr,ispp,lmax,nr,a,b,r,rab,
     1 nameat,norb,ncore,no,lo,so,zo,znuc,zel,cdd,cdu,cdc,
     2 viod,viou,vid,viu,vod,vou,etot,ev,ek,ep,wk1,wk2,
     3 wk3,wk4,wk5,wk6,wk7,nops,v,ar,br,wkb,evi,nval_orig)
c
c *************************************************************
c *                                                           *
c *     This routine was written by Norman J. Troullier Jr.   *
c *   April 1990, while at the U. of Minnesota, all           *
c *   comments concerning this routine should be directed     *
c *   to him.                                                 *
c *                                                           *
c *     troullie@128.101.224.101                              *
c *     troullie@csfsa.cs.umn.edu                             *
c *     612 625-0392                                          *
c *                                                           *
c *     pseud2 generates a pseudopotential using the          *
c *   improved scheme of N. Troullier and J. L. Martins.      *
c *   The general format of this routine is the same as the   *
c *   pseudo and pseudk routines.  Output/input is            *
c *   compatible.                                             *
c *                                                           *
c *************************************************************
c
c mmga
c mmga
c mmga
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
      parameter (zero=0.D0,one=1.D0,tpfive=2.5D0,ecuts=1.0D-3)  
      parameter (small=1.D-12,pnine=0.9D0,ai=2*137.0360411D0,sml=0.1D0)
Cray      parameter (zero=0.0,one=1.0,tpfive=2.5,ecuts=1.0E-3)  
Cray      parameter (small=1.E-12,pnine=0.9,ai=2*137.0360411,sml=0.1)
c
      character*1 ispp,blank,il(5)
      character*2 icorr,nameat
      character*3 irel
      character*4 nicore
      character*10 iray(6),ititle(7)
c
      dimension r(nr),rab(nr),no(norb),lo(norb),so(norb),zo(norb),
     1 cdd(nr),cdu(nr),cdc(nr),viod(lmax,nr),viou(lmax,nr),
     2 vid(nr),viu(nr),vod(nr),vou(nr),ev(norb),ek(norb),ep(norb),
     3 wk1(nr),wk2(nr),wk3(nr),wk4(nr),wk5(nr),wk6(nr),wk7(nr),
     4 wkb(3*nr),nops(norb),v(nr),ar(nr),br(nr),evi(norb)
c
      dimension indd(5),indu(5),rc(5),rcut(10),
     1 etot(10),aa(7),rr(7),coe(7),aj(5,5),bj(5)
c
      data il/'s','p','d','f','g'/
      if (ncore .eq. norb) return
      ifcore = itype-1
      pi = 4*atan(one)
      do 3 i=1,5
        indd(i)=0
        indu(i)=0
 3    continue
      do 4 i=1,40
        nops(i) = 0
 4    continue
c
c  read rc(s),rc(p),rc(d),rc(f),rc(g),cfac,rcfac
c
c    cfac is used for the pseudocore - the pseudocore stops where
c  the core charge density equals cfac times the renormalized
c  valence charge density (renormalized to make the atom neutral).
c  If cfac is input as negative, the full core charge is used,
c  if cfac is input as zero, it is set equal to one.
c    rcfac is used for the pseudocore cut off radius.  If set
c  to less then or equal to zero cfac is used.  cfac must be 
c  set to greater then zero.
c
cc      read(5,10) (rc(i),i=1,5),cfac,rcfac

ccccccccccccccccccc  changed by Lin-Wang Wang
      cfac=0.d0
      rcfac=0.d0
      do i=1,5
      rc(i)=0.d0
      enddo
      read(5,*) (rc(i),i=1,nval_orig)
ccccccccccccccccccc

 10   format(7f10.5)
      if (cfac .eq. 0.D0) cfac=one
c
c  Reset vod and vou to zero,
c  they are here used to store the pseudo valence charge density.
c
      do 15 i=1,nr
        vod(i) = zero
 15   continue
      do 16 i=1,nr
        vou(i) = zero
 16   continue
c
c  print heading
c
      write(6,20) nameat
 20   format(//,1x,a2,' pseudopotential generation using the ',
     1 'Improved Troullier and Martins method',/,1x,60('-'),//,
     2 ' nl    s    eigenvalue',6x,'rc',10x,'cdrc',7x,'delta',/)
c
c  Start loop over valence orbitals, only one orbital for each 
c  angular momentum and spin can exist.
c
      ncp = ncore+1
      do 190 i=ncp,norb
        lp = lo(i) + 1
        llp = lo(i)*lp
        if (so(i) .lt. 0.1) then
          if (indd(lp) .ne. 0) then
            write(6,1000)lp-1
            call ext(800+lp)
          else
            indd(lp) = i
          endif
        else
          if (indu(lp) .ne. 0) then
            write(6,1010)lp-1
            call ext(810+lp)       
          else
            indu(lp) = i
          endif
        endif              
 1000 format(//,'error in pseud2 - two down spin orbitals of the same ',
     1 /,'angular momentum (',i1,') exist')
 1010 format(//,'error in pseud2 - two up spin orbitals of the same ',
     1 /,'angular momentum (',i1,') exist')
c
c  Find the all electron wave function.
c
        do 29 j=1,nr
          ar(j) = zero
 29     continue
        if (so(i) .lt. 0.1) then
          do 30 j=2,nr
            v(j) = viod(lp,j)/r(j) + vid(j)
 30       continue
        else
          do 31 j=2,nr
            v(j) = viou(lp,j)/r(j) + viu(j)
 31       continue
        endif
        if (ispp .ne. 'r') then
          do 32 j=2,nr
            v(j) = v(j) + llp/r(j)**2
 32       continue
        endif
c
c  The parameter iflag has been added as a nonconvegence
c  indicator for auxillary routines.  Its value does 
c  not change its operation.  iflag is a returned value,
c  set to 1 for none convergence.
c
        if (ispp .ne. 'r') then
          iflag=0
          call difnrl(0,i,v,ar,br,lmax,nr,a,b,
     1     r,rab,norb,no,lo,so,znuc,viod,viou,
     2     vid,viu,ev,iflag,wk1,wk2,wk3,evi)
        else
          call difrel(0,i,v,ar,br,lmax,nr,a,b,r,
     1     rab,norb,no,lo,so,znuc,viod,viou,vid,viu,
     2     ev,wk1,wk2,wk3,wk4,evi)
         endif
c
c  Find last zero and extremum
c
        ka = lo(i)+1
        if (so(i) .lt. 0.1 .and. lo(i) .ne. 0) ka=-lo(i)
        nextr = no(i)-lo(i)
        rzero = zero
        arp = br(2)
c
        if (ispp .eq. 'r') then
          if (so(i) .lt. 0.1) then
            arp = ka*ar(2)/r(2) + (ev(i) - viod(lp,2)/r(2)
     1       - vid(2) + ai*ai) * br(2) / ai
          else
            arp = ka*ar(2)/r(2) + (ev(i) - viou(lp,2)/r(2)
     1       - viu(2) + ai*ai) * br(2) / ai
          endif
        endif
c
        do 40 j=3,nr-7
          if (nextr .eq. 0) goto 50
          if (ar(j-1)*ar(j) .le. zero .and. evi(i) .eq. zero)
     1     rzero = (ar(j)*r(j-1)-ar(j-1)*r(j)) / (ar(j)-ar(j-1))
          arpm = arp
          arp = br(j)
c
          if (ispp .eq. 'r') then
            if(so(i) .lt. 0.1) then
              arp = ka*ar(j)/r(j) + (ev(i) - viod(lp,j)/r(j)
     1         - vid(j) + ai*ai) * br(j) / ai
            else
              arp = ka*ar(j)/r(j) + (ev(i) - viou(lp,j)/r(j)
     1         - viu(j) + ai*ai) * br(j) / ai
            endif
          endif
c
          if (arp*arpm .gt. zero) goto 40
          rextr = (arp*r(j-1)-arpm*r(j)) / (arp-arpm)
          nextr = nextr - 1
 40     continue
 50     if (rzero .lt. r(2)) rzero = r(2)
c
c  Check rc if inside rzero,
c  reset to .9 between rmax and rzero if inside
c  if rc(lp) is negative, rc(lp) is percent of way 
c  betweeen rzero and rmax.
c
        if (rc(lp) .gt. rzero) then
        elseif(rc(lp) .ge. zero) then
          rc(lp) = rzero + pnine*(rextr-rzero)
        else
          rc(lp) = rzero - rc(lp)*(rextr-rzero)
        endif
c
c  Find the index for odd grid point closest to rc.
c
        do 70 j=1,nr
          if (r(j) .gt. rc(lp)) goto 80
 70     continue
 80     jrc=j-1
        rc(lp)=r(jrc)
c      
c  njtj  ***  plotting routines ***
c  potrw is called to save an usefull number of points
c  of the wave function to make a plot.  The info is
c  written to the current plot.dat file.
c   
        ist=1
        if (ar(jrc) .lt. zero) ist=-1
        call potrw(ar,r,nr-85,lo(i),1,ist)
        do 41 j=1,nr
          ar(j)=ar(j)*ist
          br(j)=br(j)*ist
 41     continue
c
c  njtj  ***  user should adjust for their needs  ***
c
c
c  Reset n quantum numbers.
c
        nops(i) = lp
c
c  Find the integrated charge inside rc(1-charge outside).
c
        ll = 2
        if (ispp .eq. 'r') then
          cdrc = -(ar(jrc)*ar(jrc)+br(jrc)*br(jrc))*rab(jrc)
          if (jrc .ne. 2*(jrc/2)) then
            do 102 k=jrc,1,-1
              cdrc = cdrc+ll*(ar(k)*ar(k)+br(k)*br(k))*rab(k)
              ll = 6 - ll
 102        continue
          else
            do 103 k=jrc,4,-1
              cdrc = cdrc+ll*(ar(k)*ar(k)+br(k)*br(k))*rab(k)
              ll = 6 - ll
 103        continue
            cdrc = cdrc-(ar(4)*ar(4)+br(4)*br(4))*rab(4)
            cdrc = cdrc+9*((ar(1)*ar(1)+br(1)*br(1))*rab(1)+
     1       3*(ar(2)*ar(2)+br(2)*br(2))*rab(2)+ 
     2       3*(ar(3)*ar(3)+br(3)*br(3))*rab(3)+
     3       (ar(4)*ar(4)+br(4)*br(4))*rab(4))/8
          endif
          cdrc = cdrc/3
        else
          cdrc = - ar(jrc) * ar(jrc) * rab(jrc)
          if (jrc .ne. 2*(jrc/2)) then
            do 100 k=jrc,1,-1
              cdrc = cdrc +  ll * ar(k) * ar(k) * rab(k)
              ll = 6 - ll
 100        continue
          else
            do 101 k=jrc,4,-1
              cdrc = cdrc +  ll * ar(k) * ar(k) * rab(k)
              ll = 6 - ll
 101        continue
            cdrc = cdrc - ar(4) * ar(4) * rab(4)
            cdrc = cdrc + 9 * ( ar(1) * ar(1) * rab(1) +
     1       3 * ar(2) *ar(2) * rab(2) + 
     2       3 * ar(3) *ar(3) * rab(3) +
     3       ar(4) * ar(4) * rab(4))/8
          endif
          cdrc = cdrc/3
        endif
c
c  Find the values for wave(arc), d(wave)/dr(arp), potential(vrc),
c  d(potential)/dr(vrp), and d2(potential)/dr2(vrpp)
c
        rc1 = r(jrc)
        rc2 = rc1 * rc1
        rc3 = rc2 * rc1
        rc4 = rc2 * rc2 
        rc5 = rc4 * rc1
        rc6 = rc4 * rc2
        rc7 = rc4 * rc3
        rc8 = rc4 * rc4   
        rc9 = rc4 * rc5
        rc10= rc4 * rc6
        arc = ar(jrc)
        arp = br(jrc)
        if (ispp .eq. 'r') then
          if (so(i) .lt. 0.1) then
            arp=ka*ar(jrc)/r(jrc) + (ev(i) - viod(lp,jrc)/r(jrc)
     1       - vid(jrc) + ai*ai) * br(jrc)/ai
          else
            arp=ka*ar(jrc)/r(jrc) + (ev(i) - viou(lp,jrc)/r(jrc)
     1       - viu(jrc) + ai*ai) * br(jrc)/ai
          endif
        endif
        arp =arp 
        brc = arp / arc
c 
        if (so(i) .lt. 0.1) then 
          vrc = viod(lp,jrc)/r(jrc) + vid(jrc)
          aa(1)=viod(lp,jrc-3)/r(jrc-3) + vid(jrc-3)
          aa(2)=viod(lp,jrc-2)/r(jrc-2) + vid(jrc-2)
          aa(3)=viod(lp,jrc-1)/r(jrc-1) + vid(jrc-1)
          aa(4)=vrc
          aa(5)=viod(lp,jrc+1)/r(jrc+1) + vid(jrc+1)
          aa(6)=viod(lp,jrc+2)/r(jrc+2) + vid(jrc+2)
          aa(7)=viod(lp,jrc+3)/r(jrc+3) + vid(jrc+3)
       else
          vrc = viou(lp,jrc)/r(jrc) + viu(jrc)
          aa(1)=viou(lp,jrc-3)/r(jrc-3) + viu(jrc-3)
          aa(2)=viou(lp,jrc-2)/r(jrc-2) + viu(jrc-2)
          aa(3)=viou(lp,jrc-1)/r(jrc-1) + viu(jrc-1)
          aa(4)=vrc
          aa(5)=viou(lp,jrc+1)/r(jrc+1) + viu(jrc+1)
          aa(6)=viou(lp,jrc+2)/r(jrc+2) + viu(jrc+2)
          aa(7)=viou(lp,jrc+3)/r(jrc+3) + viu(jrc+3)
        endif
        rr(1)=r(jrc-3)-r(jrc)
        rr(2)=r(jrc-2)-r(jrc)
        rr(3)=r(jrc-1)-r(jrc)
        rr(4)=zero
        rr(5)=r(jrc+1)-r(jrc)
        rr(6)=r(jrc+2)-r(jrc)
        rr(7)=r(jrc+3)-r(jrc)
        call polcoe(rr,aa,7,coe)
        vap   = coe(2)
        vapp  = 2*coe(3)
c
c   Set up matrix without the d2(potential(0)/dr2=0 condition
c   to find an intial guess for gamma.
c
        delta=zero
        bj(1)=log(arc/rc1**lp)
        bj1=bj(1)
        bj(2)=brc-lp/rc1
        bj2=bj(2)
        bj(3)=vrc-ev(i)-2*lp/rc1*bj2-bj2**2
        bj3=bj(3)
        bj(4)=vap+2*lp/rc2*bj2-2*lp/rc1*bj3-2*bj2*bj3
        bj4=bj(4)
        bj(5)=vapp-4*lp/rc3*bj2+4*lp/rc2*bj3-2*lp/rc1*bj4-2*bj3**2
     1   -2*bj2*bj4
        bj5=bj(5)
        aj(1,1)=rc2
        aj(1,2)=rc4
        aj(1,3)=rc6
        aj(1,4)=rc8 
        aj(1,5)=rc10
        aj(2,1)=2*rc1
        aj(2,2)=4*rc3
        aj(2,3)=6*rc5
        aj(2,4)=8*rc7 
        aj(2,5)=10*rc9
        aj(3,1)=2*one
        aj(3,2)=12*rc2
        aj(3,3)=30*rc4
        aj(3,4)=56*rc6
        aj(3,5)=90*rc8
        aj(4,1)=zero
        aj(4,2)=24*rc1
        aj(4,3)=120*rc3
        aj(4,4)=336*rc5
        aj(4,5)=720*rc7
        aj(5,1)=zero
        aj(5,2)=24*one
        aj(5,3)=360*rc2
        aj(5,4)=1680*rc4
        aj(5,5)=5040*rc6
        call gaussj(aj,5,5,bj,1,1)
        gamma=bj(1)
        alpha=bj(2)
        alpha1=bj(3)
        alpha2=bj(4)
        alpha3=bj(5) 
c
c  Start iteration loop to find delta, uses false postion.
c
        do 150 j=1,50
c
c  Generate pseudo wavefunction-note missing factor exp(delta).
c
          do 110 k=1,jrc
            rp=r(k)
            r2=rp*rp
            polyr = ((((alpha3*r2+alpha2)*r2+
     1       alpha1)*r2+ alpha)*r2+gamma)*r2
            ar(k) = rp**lp * exp(polyr)
 110      continue
c
c  Integrate pseudo charge density from r = 0 to rc.
c
          ll = 2
          cdps = - ar(jrc) * ar(jrc) * rab(jrc)
          if (jrc .ne. 2*(jrc/2)) then
            do 120 k=jrc,1,-1
              cdps = cdps +  ll * ar(k) * ar(k) * rab(k)
              ll = 6 - ll
 120        continue
          else
            do 121 k=jrc,4,-1
              cdps = cdps +  ll * ar(k) * ar(k) * rab(k)
              ll = 6 - ll
 121        continue
            cdps = cdps - ar(4) * ar(4) * rab(4)
            cdps = cdps + 9 * ( ar(1) * ar(1) * rab(1) +
     1       3 * ar(2) *ar(2) * rab(2) + 
     2       3 * ar(3) *ar(3) * rab(3) +
     3       ar(4) * ar(4) * rab(4))/8
          endif
          cdps = cdps/3
c
c   Calculate new delta
c
          fdnew = log(cdrc/cdps) - 2*delta
          if (abs(fdnew) .lt. small) goto 160
          if (j .eq. 1) then
            ddelta=-0.5
          else
            ddelta = - fdnew * ddelta / (fdnew-fdold)
          endif
          delta = delta + ddelta
          bj(1)=bj1-delta
          bj(2)=bj2
          bj(3)=bj3
          bj(4)=bj4
          bj(5)=bj5
          aj(1,1)=rc2
          aj(1,2)=rc4
          aj(1,3)=rc6
          aj(1,4)=rc8 
          aj(1,5)=rc10
          aj(2,1)=2*rc1
          aj(2,2)=4*rc3
          aj(2,3)=6*rc5
          aj(2,4)=8*rc7 
          aj(2,5)=10*rc9
          aj(3,1)=2*one
          aj(3,2)=12*rc2
          aj(3,3)=30*rc4
          aj(3,4)=56*rc6
          aj(3,5)=90*rc8
          aj(4,1)=zero
          aj(4,2)=24*rc1
          aj(4,3)=120*rc3
          aj(4,4)=336*rc5
          aj(4,5)=720*rc7
          aj(5,1)=zero
          aj(5,2)=24*one
          aj(5,3)=360*rc2
          aj(5,4)=1680*rc4
          aj(5,5)=5040*rc6
          call gaussj(aj,5,5,bj,1,1)
          gamma=bj(1)
          alpha=bj(2)
          alpha1=bj(3)
          alpha2=bj(4)
          alpha3=bj(5) 
          fdold = fdnew
 150    continue
c
c  End iteration loop for delta.
c
        write(6,1020)lp-1
        call ext(820+lp)
 1020 format(//,'error in pseud2 - nonconvergence in finding',
     1 /,' starting delta for angular momentum ',i1)
c
c  Bracket the correct gamma, use gamma and -gamma
c  from above as intial brackets, expands brackets
c  until a root is found..
c
 160    alpha4=zero
        x1=gamma
        x2=-gamma  
c
        call zrbac2(x1,x2,rc1,rc2,rc3,rc4,rc5,rc6,rc7,
     1   rc8,lp,arc,brc,vrc,vap,vapp,ev(i),cdrc,r,rab,
     2   jrc,delta,gamma,alpha,alpha1,alpha2,alpha3,
     3   alpha4,ar)
c
c  Iteration loop to find correct gamma, uses
c  bisection to find gamma.
c
        call rtbis2(x1,x2,rc1,rc2,rc3,rc4,rc5,rc6,rc7,
     1   rc8,lp,arc,brc,vrc,vap,vapp,ev(i),cdrc,r,rab,jrc,delta,
     2   gamma,alpha,alpha1,alpha2,alpha3,alpha4,ar)
c    
c  Augment charge density and invert schroedinger equation
c  to find new potential.
c
 645    expd = exp(delta)
        if (so(i) .lt. 0.1) then
          do 169 j=1,jrc 
            r2=r(j)*r(j)
            poly=r2*(((((alpha4*r2+alpha3)*r2+alpha2)*r2+alpha1)*
     1       r2+alpha)*r2+gamma) 
            ar(j) = r(j)**lp * expd * exp(poly)
            vod(j) = vod(j) + zo(i)*ar(j)*ar(j)
            xlamda=((((12*alpha4*r2+10*alpha3)*r2+8*alpha2)*r2+
     1       6*alpha1)*r2+4*alpha)*r2+2*gamma
            vj = ev(i) + xlamda * (2 * lp + xlamda * r2)
     1       +((((132*alpha4*r2+90*alpha3)*r2+56*alpha2)*r2+30*alpha1)*
     2       r2+12*alpha)*r2+2*gamma
            viod(lp,j) = (vj-vid(j)) * r(j)
 169      continue 
          do 168 j=jrc+1,nr
            vod(j) = vod(j) + zo(i)*ar(j)*ar(j)
 168      continue         
        else                                                   
          do 170 j=1,jrc
            r2=r(j)*r(j)
            poly=r2*(((((alpha4*r2+alpha3)*r2+alpha2)*r2+alpha1)*
     1       r2+alpha)*r2+gamma) 
            ar(j) = r(j)**lp * expd * exp(poly)
            vou(j) = vou(j) + zo(i)*ar(j)*ar(j)
            xlamda=((((12*alpha4*r2+10*alpha3)*r2+8*alpha2)*r2+
     1       6*alpha1)*r2+4*alpha)*r2+2*gamma
            vj = ev(i) + xlamda * (2 * lp + xlamda * r(j)**2)
     1       +((((132*alpha4*r2+90*alpha3)*r2+56*alpha2)*r2+30*alpha1)*
     2       r2+12*alpha)*r2+2*gamma
            viou(lp,j) = (vj-viu(j)) * r(j)
 170      continue
          do 171 j=jrc+1,nr
            vou(j) = vou(j) + zo(i)*ar(j)*ar(j)
 171      continue         
        endif 
c      
c  njtj  ***  plotting routines ***
c  potrw is called to save a usefull number of points
c  of the pseudowave function to make a plot.  The 
c  info is written to the current plot.dat file.
c  wtrans is called to fourier transform the the pseudo
c  wave function and save it to the current plot.dat file.
c  
        ist=1
        call potrw(ar,r,nr-85,lo(i),0,ist)
        if (ev(i) .eq. zero .or. evi(i) .ne. zero) ist=2
        call wtrans(ar,r,nr,rab,lo(i),ist,wk1)
c
c  njtj  ***  user should adjust for their needs  ***
c
        write(6,180) nops(i),il(lp),so(i),ev(i),rc(lp),cdrc,delta
 180  format(1x,i1,a1,f6.1,5f12.6)
 190  continue
c
c  End loop over valence orbitals.
c
c  Reset the n quantum numbers to include all valence orbitals.
c  Compute the ratio between the valence charge present and the
c  valence charge of a neutral atom.
c  Transfer pseudo valence charge to charge array
c
      zval = zero
      zratio = zero
      do 200 i=ncp,norb
        nops(i) = lo(i) + 1
        zval = zval + zo(i)
 200  continue
      zion = zval+znuc-zel
      if (zval .ne. zero) zratio=zion/zval
      do 210 i=1,nr
        cdd(i) = vod(i)
 210  continue
      do 211 i=1,nr
        cdu(i) = vou(i)
 211  continue

cccccc
c mmga
c
c  If a core correction is indicated construct pseudo core charge
c  cdc(r) = r^2*exp(ac+bc*r^2+cc*r^4) inside r(icore)
c  if cfac < 0 or the valence charge is zero the full core is used

      if (ifcore .ne. 0) then
        AC = ZERO
        BC = ZERO
        CC = ZERO
        icore = 1
        if (cfac .le. zero .or. zratio .eq. zero) then
          write(6,280) r(icore),AC,BC,CC
        else
          if (rcfac .le. zero) then
            do 220 i=nr,2,-1
              if (cdc(i) .gt. cfac*zratio*(cdd(i)+cdu(i))) goto 230
 220        continue
          else
            do 221 i=nr,2,-1
              if (r(i) .le. rcfac ) goto 230
 221        continue
          endif
 230      icore = i

          if(icorr.eq.'pb') then        ! changed by Lin-Wang Wang
cc          CALL PCC_EXP(NR,ICORE,AC,BC,CC,R,CDC)
          CALL PCC_POLY(NR,ICORE,AC,BC,CC,R,CDC)
          else
          CALL PCC_SIN(NR,ICORE,AC,BC,CC,R,CDC)
          endif

          WRITE(6,280) r(icore),AC,BC,CC

        endif
      endif
 280  format(//,' core correction used',/,
     . ' pseudo core inside r =',f7.3,/,
     . ' ac =',f7.3,' bc =',f7.3,' cc =',f7.3,/)
c 1030 format(//,' error in pseud2 - noncovergence in finding ',
c     1 /,'pseudo-core values')
c
c mmga
cccccc
c
c  End the pseudo core charge.
c  Compute the potential due to pseudo valence charge.
c
c  njtj  ***  NOTE  ***
c  Spin-polarized potentails should be unscreend with
c  spin-polarized valence charge.  This was not
c  done in pseudo and pseudok in earlier versions
c  of this program.       
c  njtj  ***  NOTE  ***
c
 290  if (ispp .eq. 's') then
        blank='s'
      else
        blank=' '
      endif
      zval2=zval
      call velect(0,1,icorr,blank,ifcore,nr,r,rab,zval,
     1 cdd,cdu,cdc,vod,vou,etot,wk1,wk2,wk3,wk4,wk5,wkb)
      if (ifcore .eq. 2) zion=zion+zval-zval2
c
c  Construct the ionic pseudopotential and find the cutoff,
c  ecut should be adjusted to give a reassonable ionic cutoff
c  radius, but should not alter the pseudopotential, ie.,
c  the ionic cutoff radius should not be inside the pseudopotential
c  cutoff radius
c
      ecut=ecuts
      do 315 i=ncp,norb
        lp = lo(i)+1
        if (so(i) .lt. 0.1) then
          do 300 j=2,nr
            viod(lp,j)=viod(lp,j) + (vid(j)-vod(j))*r(j)
            vp2z = viod(lp,j) + 2*zion
            if (abs(vp2z) .gt. ecut) jcut = j
 300      continue
          rcut(i-ncore) = r(jcut)
          do 310 j=jcut,nr
            fcut = exp(-5*(r(j)-r(jcut)))
            viod(lp,j) = - 2*zion + fcut * (viod(lp,j)+2*zion)
 310      continue 
          do 311 j=2,nr
            v(j) = viod(lp,j)/r(j)
 311      continue
c
c  njtj  ***  plotting routines ***
c
          call potran(lo(i)+1,v,r,nr,zion,wk1,wk2,wk3)
          call potrv(v,r,nr-120,lo(i))
c
c  njtj  ***  user should adjust for their needs  ***
c
        else
          do 312 j=2,nr
            viou(lp,j)=viou(lp,j)+ (viu(j)-vou(j))*r(j)
            vp2z = viou(lp,j) + 2*zion
            if (abs(vp2z) .gt. ecut) jcut = j
 312      continue
          rcut(i-ncore) = r(jcut)
          do 313 j=jcut,nr
            fcut = exp(-5*(r(j)-r(jcut)))
            viou(lp,j) = - 2*zion + fcut * (viou(lp,j)+2*zion)
 313      continue
          do 314 j=2,nr
            v(j) = viou(lp,j)/r(j)
 314      continue
c
c  njtj  ***  plotting routines ***
c
          call potran(lo(i)+1,v,r,nr,zion,wk1,wk2,wk3)
          call potrv(v,r,nr-110,lo(i))
c
c  njtj  ***  user should adjust for their needs  ***
c
        endif
 315  continue
c
c  njtj  ***  plotting routines ***
c   The calls to 1)potran take the fourier transform of
c   the potential and saves it in the current plot.dat file,
c   2)potrv saves the potential in the current plot.dat file
c   3)zion is saved to the current plot.dat file wtih a 
c   marker 'zio' for latter plotting
c
      write(3,4559)
      write(3,4560) zion
 4559 format(1x,'marker zio')
 4560 format(2x,f5.2)
c
c  njtj  ***  user should adjust for their needs  ***
c
c   Convert spin-polarized potentials back to nonspin-polarized
c   by occupation weight(zo).  Assumes core polarization is
c   zero, ie. polarization is only a valence effect. 
c
      if (ispp .eq. 's' ) then
        do 500 i=ncp,norb,2
          lp = lo(i)+1
          zot=zo(i)+zo(i+1)
          if (zot .ne. zero) then
            do 505 j=2,nr
              viod(lp,j)=(viod(lp,j)*zo(i)+viou(lp,j)
     1         *zo(i+1))/zot
              viou(lp,j)=viod(lp,j)
 505        continue 
          else
            do 506 j=2,nr
              viod(lp,j)=viod(lp,j)/2+viou(lp,j)/2
              viou(lp,j)=viod(lp,j)
 506        continue 
          endif
 500    continue
      endif
c
      do 320 i=2,nr
        vid(i) = vod(i)
        viu(i) = vou(i)
 320  continue
c 
c   Test the pseudopotential self consistency.  Spin-polarized
c   is tested as spin-polarized(since up/down potentials are
c   now the same)
c
      call dsolv2(0,1,blank,ifcore,lmax,nr,a,b,r,rab,
     1 norb-ncore,0,nops(ncp),lo(ncp),so(ncp),zo(ncp),
     2 znuc,cdd,cdu,cdc,viod,viou,vid,viu,ev(ncp),ek(ncp),
     3 ep(ncp),wk1,wk2,wk3,wk4,wk5,wk6,wk7,evi(ncp))
c
c  Printout the pseudo eigenvalues after cutoff.
c
      write(6,325) (il(lo(i)+1),rcut(i-ncore),i=ncp,norb)
      write(6,326) (ev(i),i=ncp,norb)
 325  format(//,' test of eigenvalues',//,' rcut =',8(2x,a1,f7.2))
 326  format(' eval =',8(2x,f8.5))
c
c  Printout the data for potentials.
c
      write(6,330)
 330  format(///,' l    vps(0)    vpsmin      at r',/)
      do 370 i=1,lmax
        if (indd(i)+indu(i) .eq. 0) goto 370
        if (indd(i) .ne. 0) then
          vpsdm = zero
          do 350 j=2,nr
            if (r(j) .lt. .00001) goto 350
            vps = viod(i,j)/r(j)
            if (vps .lt. vpsdm) then
              vpsdm = vps
              rmind = r(j)
            endif
 350      continue
          write(6,360) il(i),viod(i,2)/r(2),vpsdm,rmind
        endif
        if (indu(i) .ne. 0) then
          vpsum = zero
          do 351 j=2,nr
            if (r(j) .lt. .00001) goto 351
            vps = viou(i,j)/r(j)
            if (vps .lt. vpsum) then
              vpsum = vps
              rminu = r(j)
            endif
 351      continue
          write(6,360) il(i),viou(i,2)/r(2),vpsum,rminu
        endif
 360  format(1x,a1,3f10.3)
 370  continue
c
c   Print out the energies from etotal.
c
      call etotal(itype,one,nameat,norb-ncore,
     1 nops(ncp),lo(ncp),so(ncp),zo(ncp),
     2 etot,ev(ncp),ek(ncp),ep(ncp))
c
c  Find the jobname and date, date is a machine 
c  dependent routine and must be chosen/written/
c  comment in/out in the zedate section.
c
      iray(1) = 'atom-lda  '
      call zedate(iray(2))
      iray(3) = '  Improved'
      iray(4) = ' Troullier'
      iray(5) = ' - Martins'
      iray(6) = ' potential'
c
c  Encode the title array.
c
      do 390 i=1,7
        ititle(i) = '          '
 390  continue
      do 420 i=1,lmax
        if (indd(i) .eq. 0 .and. indu(i) .eq. 0) goto 420
        zelu = zero
        zeld = zero
        if (indd(i) .ne. 0) then
          noi = no(indd(i))
          zeld = zo(indd(i))
        endif
        if (indu(i) .ne. 0) then
          noi = no(indu(i))
          zelu = zo(indu(i))
        endif
        zelt = zeld + zelu
       if (ispp .ne. 's') then
         write(ititle(2*i-1),400) noi,il(i),zelt
         write(ititle(2*i),401)ispp,rc(i)
 400     format(i1,a1,'(',f6.2,')')
 401     format(a1,' rc=',f5.2)
       else
         write(ititle(2*i-1),410) noi,il(i),zeld
         write(ititle(2*i),411)zelu,ispp,rc(i) 
 410     format(i1,a1,'  (',f4.2,',')
 411     format(f4.2,')',a1,f4.2)
        endif
 420  continue
c
c  Construct relativistic sum and difference potentials.
c
      if (ispp .eq. 'r') then
        if (indu(1) .eq. 0) goto 429
        indd(1)=indu(1)
        indu(1)=0
        do 428 j=2,nr
          viod(1,j) = viou(1,j)
          viou(1,j) = zero
 428    continue
 429    do 431 i=2,lmax
          if (indd(i) .eq. 0 .or. indu(i) .eq. 0) goto 431
          do 430 j=2,nr
            viodj = viod(i,j)
            viouj = viou(i,j)
            viod(i,j) = ((i-1)*viodj + i*viouj) / (2*i-1)
            viou(i,j) = 2 * (viouj - viodj) / (2*i-1)
 430      continue
 431    continue
      endif
c
c  Determine the number of  potentials.  Coded them as 
c  two digits, where the first digit is the number
c  of down or sum potentials and the second the number of
c  up or difference potentials.
c
      npotd = 0
      npotu = 0
      do 450 i=1,lmax
        if (indd(i) .ne. 0) npotd=npotd+1
        if (indu(i) .ne. 0) npotu=npotu+1
 450  continue
c
c  Write the heading to the current pseudo.dat 
c  file (unit=1).
c
      ifull = 0
      if (cfac .le. zero .or. zratio .eq. zero) ifull = 1
      if (ifcore .eq. 1) then
        if (ifull .eq. 0) then
          nicore = 'pcec'
        else
          nicore = 'fcec'
        endif
      elseif (ifcore .eq. 2) then
        if (ifull .eq. 0) then
          nicore = 'pche'
        else
          nicore = 'fche'
        endif
      else
        nicore = 'nc  '
      endif
      if (ispp .eq. 's') then
        irel='isp'
      elseif (ispp .eq. 'r') then
        irel='rel'
      else
        irel = 'nrl'
      endif
      rewind 1
      write(1) nameat,icorr,irel,nicore,(iray(i),i=1,6),
     1 (ititle(i),i=1,7),npotd,npotu,nr-1,a,b,zion
      write(1) (r(i),i=2,nr)
c
c  Write the potentials to the current pseudo.dat 
c  file (unit=1).
c
      do 460 i=1,lmax
        if (indd(i) .eq. 0) goto 460
        write(1) i-1,(viod(i,j),j=2,nr)
 460  continue
      do 465 i=1,lmax
        if (indu(i) .eq. 0) goto 465
        write(1) i-1,(viou(i,j),j=2,nr)
 465  continue
c
c  Write the charge densities to the current pseudo.dat 
c  file (unit=1).
c
      if (ifcore .eq. 0) then
        write(1) (zero,i=2,nr)
      else
        write(1) (cdc(i),i=2,nr)
      endif
      write(1) (zratio*(cdd(i)+cdu(i)),i=2,nr)
c
      return
      end
