      subroutine pseudk(itype,icorr,ispp,lmax,nr,a,b,
     1 r,rab,nameat,norb,ncore,no,lo,so,zo,znuc,zel,
     2 cdd,cdu,cdc,viod,viou,vid,viu,vod,vou,etot,ev,
     3 ek,ep,wk1,wk2,wk3,wk4,wk5,wk6,wk7,nops,v,ar,br,
     4 wkb,evi,nval_orig)
c
c *************************************************************
c *                                                           *
c *     pseudk generates the pseudo potential using the       *
c *   scheme of G. P. Kerker, J. Phys. C13, L189 (1980).      *
c *                                                           *
c *************************************************************
c
c mmga
c mmga
c mmga
c
c  njtj  ***  modifications  ***
c    The only major modification is in the spin-polarization
c    treatment of the unscreening of the pseudopotential.
c    Spin-polarized potentails should be unscreend with
c    spin-polarized valence charge.  This was not
c    done in pseudo and pseudk in earlier Berkeley/Froyen
c    versions of this program.       
c  njtj  ***  modifications  ***
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
      parameter (zero=0.D0,tpfive=2.5D0,pfive=0.5D0,smtol=1.D-12)
      parameter (one=1.D0,ai=2*137.0360411D0,ecuts=1.0D-3)
Cray       parameter (zero=0.0,tpfive=2.5,pfive=0.5,smtol=1.E-12)
Cray       parameter (one=1.0,ai=2*137.0360411,ecuts=1.0D-3)
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
      dimension indd(5),indu(5),rc(5),rcut(10),etot(10)
c
      data il/'s','p','d','f','g'/
      do 3 i=1,5
        indd(i)=0
        indu(i)=0
 3    continue
      if (ncore .eq. norb) return
      if (itype .lt. 1 .or. itype .gt. 3) return
      ifcore = itype-1
      pi = 4*atan(one)
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
      if (cfac .eq. zero) cfac=one
c
c    Reset vod and vou to zero, they are used to store the pseudo 
c  valence charge density.
c
      do 15 i=1,nr
        vod(i) = zero
        vou(i) = zero
 15   continue
c
c  Print the heading.
c
      write(6,20) nameat
 20   format(//,a2,' Pseudopotential generation using the method',
     1 ' of Kerker',/,1x,60('-'),//,
     2 ' nl    s    eigenvalue',6x,'rc',4x,6x,'cdrc',7x,'delta',
     3 7x,/)
c
c      start loop over valence orbitals
c
      ncp = ncore+1
      do 190 i=ncp,norb
        lp = lo(i) + 1
        llp = lo(i)*lp
        if (so(i) .lt. 0.1 .and. indd(lp) .ne. 0) call ext(800+lp)
        if (so(i) .gt. 0.1 .and. indu(lp) .ne. 0) call ext(810+lp)
        if (so(i) .lt. 0.1) then
          indd(lp) = i
        else
          indu(lp) = i
        endif
c
c  Find the all-electron wave function.
c
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
        if (ispp .ne. 'r') then
          call difnrl(0,i,v,ar,br,lmax,nr,a,b,r,rab,norb,no,lo,so,
     1     znuc,viod,viou,vid,viu,ev,iflag,wk1,wk2,wk3,evi)
        else
          call difrel(0,i,v,ar,br,lmax,nr,a,b,r,rab,norb,no,lo,so,
     1     znuc,viod,viou,vid,viu,ev,wk1,wk2,wk3,wk4,evi)
        endif
c      
c  njtj  ***  plotting routines ***
c  potrw is called to save an usefull number of points
c  of the wave function to make a plot.  The info is
c  written to the current plot.dat file.
c   
        ist=1
        if (ar(nr-85) .lt. zero) ist=-1
        call potrw(ar,r,nr-85,lo(i),1,ist) 
c
c  njtj  ***  user should adjust for their needs  ***
c
c  Find the last zero and extremum point.
c
        ka = lo(i)+1
        if (so(i) .lt. 0.1 .and. lo(i) .ne. 0) ka=-lo(i)
        nextr = no(i)-lo(i)
        rzero = zero
        arp = br(2)
c
        if (ispp .eq. 'r') then
          if (so(i) .lt. 0.1) then
            arp=ka*ar(2)/r(2)+(ev(i)-viod(lp,2)/r(2)-vid(2)+ai*ai) *
     1         br(2)/ai
          else
            arp=ka*ar(2)/r(2)+(ev(i)-viou(lp,2)/r(2)-viu(2)+ai*ai) *
     1         br(2)/ai
          endif
        endif
c
        do 40 j=3,nr-7
          if (nextr .eq. 0) goto 50
          if (ar(j-1)*ar(j) .le. zero)
     1     rzero = (ar(j)*r(j-1)-ar(j-1)*r(j)) / (ar(j)-ar(j-1))
          arpm = arp
          arp = br(j)
c
          if (ispp .eq. 'r') then
            if(so(i) .lt. 0.1) then
              arp=ka*ar(j)/r(j)+(ev(i)-viod(lp,j)/r(j)-
     1         vid(j)+ai*ai)*br(j)/ai
            else
              arp=ka*ar(j)/r(j)+(ev(i)-viou(lp,j)/r(j)-
     1         viu(j)+ai*ai)*br(j)/ai
            endif
          endif
c
          if (arp*arpm .gt. zero) goto 40
          rextr = (arp*r(j-1)-arpm*r(j)) / (arp-arpm)
          nextr = nextr - 1
 40     continue
c
c   Check rc, if outside bounds reset.
c
 50     if (rzero .lt. r(2)) rzero = r(2)
        if (rc(lp) .gt. rzero .and. rc(lp) .lt. rextr) goto 60
        if (rc(lp) .lt. zero) rc(lp) = rzero - rc(lp)*(rextr-rzero)
c
c   Find index for grid point closest to rc.
c
 60     do 70 j=1,nr
          if (r(j) .gt. rc(lp)) goto 80
          jrc = j
 70     continue
c
c   Reset the n quantum numbers.
c
 80     rc(lp)=r(jrc)
        do 90 j=1,norb
          nops(j) = 0
 90     continue
        nops(i) = lp
c
c  Find the integrated charge inside rc.
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
c   The initial values for alpha, beta, gamma and delta.
c
         rc2 = r(jrc) * r(jrc)
         rc3 = r(jrc) * rc2
         rc4 = r(jrc) * rc3
         iswtch = 1
         if (ar(jrc) .lt. zero) iswtch = -1
         arc = iswtch * ar(jrc)
         arp = br(jrc)
c
         if (ispp .eq. 'r') then
           if(so(i) .lt. 0.1) then
             arp=ka*ar(jrc)/r(jrc)+(ev(i)-viod(lp,jrc)/r(jrc)-
     1        vid(jrc) + ai*ai) * br(jrc)/ai
           else
             arp=ka*ar(jrc)/r(jrc)+(ev(i)-viou(lp,jrc)/r(jrc)-
     1        viu(jrc) + ai*ai) * br(jrc)/ai
           endif
         endif
c
         brc = arp / ar(jrc)
         if (so(i) .lt. 0.1) then
           vrc = viod(lp,jrc)/r(jrc) + vid(jrc)
         else
           vrc = viou(lp,jrc)/r(jrc) + viu(jrc)
         endif
         alpha = ( 3*log(arc/r(jrc)**lp) - 2*(r(jrc)*brc-lp)
     1    + (rc2*vrc+lp*lp-rc2*(ev(i)+brc*brc))/2 ) / rc4
         beta  = (-8*log(arc/r(jrc)**lp) + 5*(r(jrc)*brc-lp)
     1    - (rc2*vrc+lp*lp-rc2*(ev(i)+brc*brc))   ) / rc3
         gamma = ( 6*log(arc/r(jrc)**lp) - 3*(r(jrc)*brc-lp)
     1    + (rc2*vrc+lp*lp-rc2*(ev(i)+brc*brc))/2 ) / rc2
         delta = zero
c
c  Start the iteration loop to find delta.
c                                         
         do 150 j=1,50
c
c  Generate the pseudo-wavefunction (note missing factor exp(delta)).
c
           do 110 k=1,jrc
             polyr=r(k)*r(k)*((alpha*r(k)+beta)*r(k)+gamma)
             ar(k) = iswtch * r(k)**lp * exp(polyr)
 110       continue
c
c  Integrate  the pseudo charge density from r = 0 to rc.
c
           ll = 2
           cdps = - ar(jrc) * ar(jrc) * rab(jrc)
           if (jrc .ne. 2*(jrc/2)) then
             do 120 k=jrc,1,-1
               cdps = cdps +  ll * ar(k) * ar(k) * rab(k)
               ll = 6 - ll
 120         continue
           else
             do 121 k=jrc,4,-1
               cdps = cdps +  ll * ar(k) * ar(k) * rab(k)
               ll = 6 - ll
 121         continue
             cdps = cdps - ar(4) * ar(4) * rab(4)
             cdps = cdps + 9 * ( ar(1) * ar(1) * rab(1) +
     1        3 * ar(2) *ar(2) * rab(2) + 
     2        3 * ar(3) *ar(3) * rab(3) +
     3        ar(4) * ar(4) * rab(4))/8
           endif
           cdps = cdps/3
c
c  Find the new delta.
c
           fdnew = log(cdrc/cdps) - 2*delta
           if (abs(fdnew) .lt. smtol) goto 160
           if (j .eq. 1) then
             ddelta = pfive
           else
             ddelta = - fdnew * ddelta / (fdnew-fdold)
           endif
           alpha = alpha - 3*ddelta/rc4
           beta  = beta  + 8*ddelta/rc3
           gamma = gamma - 6*ddelta/rc2
           delta = delta + ddelta
           fdold = fdnew
 150     continue
c
c  End the iteration loop for delta.
c
         call ext(820+lp)
c
c    Augment the charge density and invert schroedinger equation
c  to find new potential.
c
 160     expd = exp(delta)
         if (so(i) .lt. 0.1) then
           do 170 j=1,jrc
             ar(j) = expd * ar(j)
             xlamda=(4*alpha*r(j)+3*beta)*r(j)+2*gamma
             vj = ev(i) + xlamda * (2 * lp + xlamda * r(j)**2)
     1        + (12 * alpha * r(j) + 6 * beta) * r(j) + 2 * gamma
             viod(lp,j) = (vj - vid(j)) * r(j)
             vod(j) = vod(j) + zo(i)*ar(j)*ar(j)
 170       continue
           do 171 j=jrc+1,nr
             vod(j) = vod(j) + zo(i)*ar(j)*ar(j)
 171       continue
         else
           do 175 j=1,jrc
             ar(j) = expd * ar(j)
             xlamda=(4*alpha*r(j)+3*beta)*r(j)+2*gamma
             vj = ev(i) + xlamda * (2 * lp + xlamda * r(j)**2)
     1        + (12 * alpha * r(j) + 6 * beta) * r(j) + 2 * gamma
             viou(lp,j) = (vj - viu(j)) * r(j)
             vou(j) = vou(j) + zo(i)*ar(j)*ar(j)
 175       continue
           do 176 j=jrc+1,nr
             vou(j) = vou(j) + zo(i)*ar(j)*ar(j)
 176       continue
         endif
         write(6,180) nops(i),il(lp),so(i),ev(i),rc(lp),cdrc,delta
 180     format(1x,i1,a1,f6.1,5f12.6)
c      
c  njtj  ***  plotting routines ***
c  potrw is called to save a usefull number of points
c  of the pseudowave function to make a plot.  The 
c  info is written to the current plot.dat file.
c  wtrans is called to fourier transform the the pseudo
c  wave function and save it to the current plot.dat file.
c   
         ist=1
         if (ar(nr-85) .lt. zero) ist=-1
         call potrw(ar,r,nr-85,lo(i),0,ist)
         if (ev(i) .eq. zero .or. evi(i) .ne. zero) ist=2
         call wtrans(ar,r,nr,rab,lo(i),ist,wk1)
c
c  njtj  ***  user should adjust for their needs  ***
c

 190   continue
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

          if(icorr.eq.'pb') then
c          CALL PCC_EXP(NR,ICORE,AC,BC,CC,R,CDC)
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
      call velect(0,1,icorr,blank,ifcore,nr,r,rab,zval,
     1 cdd,cdu,cdc,vod,vou,etot,wk1,wk2,wk3,wk4,wk5,wkb)
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
          do 500 j=2,nr
            viod(lp,j)=viod(lp,j) + (vid(j)-vod(j))*r(j)
            vp2z = viod(lp,j) + 2*zion
            if (abs(vp2z) .gt. ecut) jcut = j
 500      continue
          rcut(i-ncore) = r(jcut)
          do 510 j=jcut,nr
            fcut = exp(-5*(r(j)-r(jcut)))
            viod(lp,j) = - 2*zion + fcut * (viod(lp,j)+2*zion)
 510      continue 
          do 511 j=2,nr
            v(j) = viod(lp,j)/r(j)
 511      continue
c
c  njtj  ***  plotting routines ***
c
          call potran(lo(i)+1,v,r,nr,zion,wk1,wk2,wk3)
          call potrv(v,r,nr-120,lo(i))
c
c  njtj  ***  user should adjust for their needs  ***
c
        else
          do 512 j=2,nr
            viou(lp,j)=viou(lp,j)+ (viu(j)-vou(j))*r(j)
            vp2z = viou(lp,j) + 2*zion
            if (abs(vp2z) .gt. ecut) jcut = j
 512      continue
          rcut(i-ncore) = r(jcut)
          do 513 j=jcut,nr
            fcut = exp(-5*(r(j)-r(jcut)))
            viou(lp,j) = - 2*zion + fcut * (viou(lp,j)+2*zion)
 513      continue
          do 514 j=2,nr
            v(j) = viou(lp,j)/r(j)
 514      continue
c
c  njtj  ***  plotting routines ***
c
          call potran(lo(i)+1,v,r,nr,zion,wk1,wk2,wk3)
          call potrv(v,r,nr-120,lo(i))
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
        do 700 i=ncp,norb,2
          lp = lo(i)+1
          zot=zo(i)+zo(i+1)
          if (zot .ne. zero) then
            do 705 j=2,nr
              viod(lp,j)=(viod(lp,j)*zo(i)+viou(lp,j)
     1         *zo(i+1))/zot
              viou(lp,j)=viod(lp,j)
 705        continue 
          else
            do 706 j=2,nr
              viod(lp,j)=viod(lp,j)/2+viou(lp,j)/2
              viou(lp,j)=viod(lp,j)
 706        continue 
          endif
 700    continue
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
      iray(3) = '   Kerker-'
      iray(4) = 'potential '
      do 380 i=5,6
        iray(i) = '          '
 380  continue
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
