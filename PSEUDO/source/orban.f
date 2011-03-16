      subroutine orban(ispp,iorb,ar,br,lmax,nr,a,b,r,rab,
     1 norb,no,lo,zo,so,viod,viou,vid,viu,ev,ek,ep)
c
c  orban is used to analyze and printout data 
c  about the orbital.
c
      implicit real*8 (a-h,o-z)
c
      parameter (ai=2*137.0360411d0,zero=0.0d0)
c
      character*1 ispp
      character*10 name
c
      dimension ar(nr),br(nr),r(nr),rab(nr),no(norb),
     1 lo(norb),zo(norb),so(norb),viod(lmax,nr),viou(lmax,nr),
     2 vid(nr),viu(nr),ev(norb),ek(norb),ep(norb)
c
      dimension rzero(10),rextr(10),aextr(10),bextr(10)
c
c     dimension wk1(1000),wk2(1000),wk3(1000),v(1000)
c
      ka = lo(iorb)+1
      lp = ka
      if (so(iorb) .lt. 0.1 .and. lo(iorb) .ne. 0) ka=-lo(iorb)
c
c      compute zeroes and extrema
c
      nzero = 0
      nextr = 0
      rzero(1) = zero
      arp = br(2)
      if (ispp .eq. 'r') then
        if (so(iorb) .lt. 0.1) then
          arp = ka*ar(2)/r(2) + (ev(iorb) - viod(lp,2)/r(2)
     1     - vid(2) + ai*ai) * br(2) / ai
        else
          arp = ka*ar(2)/r(2) + (ev(iorb) - viou(lp,2)/r(2)
     1     - viu(2) + ai*ai) * br(2) / ai
        endif
      endif
      do 20 i=3,nr
        if (nextr .ge. no(iorb)-lo(iorb)) goto 30
        if (ar(i)*ar(i-1) .gt. zero) goto 10
c
c   zero
c
        nzero = nzero + 1
        rzero(nzero) = (ar(i)*r(i-1)-ar(i-1)*r(i)) / (ar(i)-ar(i-1))
 10     arpm = arp
        arp = br(i)
        if (ispp .eq. 'r') then
          if ( so(iorb) .lt. 0.1) then
            arp = ka*ar(i)/r(i) + (ev(iorb) - viod(lp,i)/r(i)
     1       - vid(i) + ai*ai) * br(i) / ai
          else
            arp = ka*ar(i)/r(i) + (ev(iorb) - viou(lp,i)/r(i)
     1       - viu(i) + ai*ai) * br(i) / ai
          endif
        endif
        if (arp*arpm .gt. zero) goto 20
c
c   extremum
c
        nextr = nextr + 1
        rextr(nextr) = (arp*r(i-1)-arpm*r(i)) / (arp-arpm)
        aextr(nextr) = (ar(i)+ar(i-1))/2
     1   - (arp**2+arpm**2) * (r(i)-r(i-1)) / (4*(arp-arpm))
        bextr(nextr) = br(i)
 20   continue
c
c   find orbital kinetic and potential energy
c   the potential part includes only the interaction with
c   the nuclear part
c
 30   ek(iorb) = br(1)*br(1)*rab(1)
      ep(iorb) = zero
      sa2 = zero
      lp = lo(iorb)+1
      llp = lo(iorb)*lp
      ll = 2
      if (2*(nr/2) .eq. nr) ll=4
      i90=nr
      i99=nr
      do 40 i=nr,2,-1
        ar2 = ar(i)*ar(i)
        br2 = br(i)*br(i)
        deni = ar2
        if (ispp .eq. 'r') deni=deni+br2
        ek(iorb) = ek(iorb) + ll * (br2 + ar2*llp/r(i)**2)*rab(i)
        if (so(iorb) .lt. 0.1) then
          ep(iorb) = ep(iorb) + ll * deni*viod(lp,i)*rab(i)/r(i)
        else
          ep(iorb) = ep(iorb) + ll * deni*viou(lp,i)*rab(i)/r(i)
        endif
        ll = 6 - ll
        if (sa2 .gt. 0.1) goto 40
        sa2 = sa2 + deni*rab(i)
        if (sa2 .le. 0.01) i99 = i
        i90 = i
 40   continue
      ek(iorb) = ek(iorb) / 3
      ep(iorb) = ep(iorb) / 3
      if (ispp .eq. 'r') ek(iorb) = zero
c
c   printout
c
      write(6,80) no(iorb),lo(iorb),so(iorb)
 80   format(/,' n =',i2,'  l =',i2,'  s =',f4.1)
      name = 'a extr    '
      write(6,100) name,(aextr(i),i=1,nextr)
      name = 'b extr    '
      if (ispp .eq. 'r') write(6,100) name,(bextr(i),i=1,nextr)
      name = 'r extr    '
      write(6,100) name,(rextr(i),i=1,nextr)
      name = 'r zero    '
      write(6,100) name,(rzero(i),i=1,nzero)
      name = 'r 90/99 % '
      write(6,100) name,r(i90),r(i99)
      if (ev(iorb) .eq. zero) then
        if (zo(iorb) .ne. zero) then
          write(6,110)zo(iorb)
        else
          write(6,120)
        endif
      endif
 100  format(8x,a10,2x,8f8.3)
 110  format(8x,'WARNING: This orbital is not bound',
     1 ' and contains ',f6.4,' electrons!!')
 120  format(8x,'WARNING:  This orbital is not bound!')
c
c  njtj  ***  plotting routines  ***
c  jlm   uncomment dimension declaration for wk1,etc
c
c    Save plotting information to current plot.dat file
c  (unit = 3),  User must specify what orbital 
c   is to be saved(or all).
c
c      ist=1
c      if (ar(nr-80) .lt. 0.0) ist=-1
c      call potrw(ar,r,nr-85,lo(iorb),1,ist)
c      call wtrans(ar,r,nr,rab,lo(iorb),ist,wk1)
c
c       do i=2,nr
c         v(i)=viod(lo(iorb)+1,i)/r(i)+vid(i)
c       enddo
c       zion=4
c       call potran(lo(iorb)+1,v,r,nr,zion,wk1,wk2,wk3)
c       call potrv(v,r,nr,lo(iorb))
c
c  njtj  ***  user should adjust for their needs  ***
c
       return
       end
