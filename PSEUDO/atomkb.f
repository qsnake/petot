      program atomkb
c            
c  *************************************************************
c  *     Program for atomic calculations                       *
c  *     with a Kleinman-Bylander operator pseudo-potential.   *
c  *     By Norman J. Troullier Jr., Feb 1990                  *
c  *     while at U of Minn, NEC Research                      *
c  *                                                           *
c  *     Send comments/suggestions/bug reports to:             *
c  *     troullie@128.101.224.101,norm@research.nec.com        *
c  *     troullie@csfsa.cs.umn.edu                             *
c  *     612 625-0392                                          *
c  *                                                           *
c  *     Version 1.10, Dated Oct. 14, 1991                     *
c  *                                                           *
c mmga                                                         *
c mmga                                                         *
c mmga                                                         *
c  *                                                           *
c  *************************************************************
c
c    Some parameters are set inside the program,
c  the most important ones are:
c  1)the tolerance for selfconsistency in the screening 
c    potential (set in the main-atm program-named tol),
c  2)the accuracy in the eigenvalue for a given potential
c    (set in difgrd-named tol),
c  3)the dimensions of the work space used: nrmax,
c    norbmx, lmax(needs to be set only in main-atm),
c  4)the machine precision - MACHEP, for use in the 
c    eispack group of fuctions: tinvit, and tridib.
c    (The current value is ok for this application.)
c  5)the machine precision exp(-2*expzer), set in difgrd
c    for the predictor-corrector methods, named expzer.
c  
c    For traditional applications the default values 
c  should be enough for 1-4.  For 5, expzer should be
c  found for the system in use.
c    NOTE: that for most VAX systems the value of expzer
c  can be very very small !!
c
c    The subroutine orban is called once for each orbital
c  after selfconsistency has been achieved.
c  Any additional analysis of the orbitals should therefore
c  be placed in orban.
c
c  tolerance for self-consistency
c
c
      implicit real*8 (a-h,o-z)
      parameter (tol=1.D-8)
c
      parameter (zero=0.0D0,one=1.0D0)
c
      parameter(lmax=5,nrmax=2000,norbmx=40)

c
      dimension r(nrmax),rab(nrmax),no(norbmx),lo(norbmx),
     1 so(norbmx),zo(norbmx),cdd(nrmax),cdu(nrmax),cdc(nrmax),
     2 viod(lmax,nrmax),viou(lmax,nrmax),vid(nrmax),viu(nrmax),
     3 vod(nrmax),vou(nrmax),vn1d(nrmax),vn1u(nrmax),
     4 vn11d(nrmax),vn11u(nrmax),vn2d(nrmax),vn2u(nrmax),
     5 vn22d(nrmax),vn22u(nrmax),ev(norbmx),ek(norbmx),ep(norbmx),
     6 wk1(nrmax),wk2(nrmax),wk3(nrmax),wk4(nrmax),wk5(nrmax),
     7 wk6(nrmax),wk7(nrmax),wk8(nrmax),wk9(nrmax),wkb(6*nrmax),
     8 inorm(norbmx),anorm(norbmx),arw(nrmax,norbmx),vql(nrmax),
     9 evi(norbmx)
c
      dimension econf(100),etot(10)
c
      character*1 ispp
      character*2 icorr,nameat
      character*10 plotfile
c
      do 5 j=1,nrmax
        do 6 i=1,norbmx
          arw(j,i)= zero
 6      continue
 5    continue
c
c    Call to machine dependent cpu time routine.
c    User may want to comment out timing calls -
c    here and at exit of main - atm
c
      call zesec(t1)
c 
      nconf = 0
c
c      open files
c
      open(unit=1,file='pseudo.dat',form='unformatted',
     1  status='unknown')
c      open(unit=3,file='plot.dat',status='new',form='formatted')
      open(unit=5,file='atom.input',status='old',form='formatted')
      open(unit=6,file='atomkb.out',status='unknown',form='formatted')
      rewind(5)
      rewind(6)
c
c   Start of loop over configuration.
c   Read the input data.

c
 20   nr = nrmax
      norb = norbmx
      call input(itype,ikerk,icorr,ispp,zsh,rsh,
     1 nr,a,b,r,rab,nameat,norb,ncore,no,lo,so,zo,
     2 znuc,zel,evi,nval_orig)
c

c  Stop - no more data in input file,
c  Find time taken for total calculation.
c  zesec - machine dependent routine
c
      if (itype .lt. 4) then
        call zesec(t2)
        write(6,2000)t2-t1
 2000 format(//,' The total time for the K-B Atomic',
     1 '  calculation is ',f12.5,' seconds')
        stop
      endif
c
c  Jump charge density 'set up' and ionic data input if
c  configuration test.
c
      itsm=znuc/9+3
c
c  set up ionic potentials
c
        call vionic(ispp,itype,icorr,ifcore,0.0,0.0,
     1   lmax,nr,a,b,r,rab,nameat,ncore,znuc,
     2   cdd,cdu,cdc,viod,viou)
c
c  Read in normalization data for operators.
c
        if (nconf .eq. 0) call input2(nr,lmax,vql,inorm,anorm)
c
c   Set up the electronic potential.
c
      call velect(0,0,icorr,ispp,ifcore,
     1 nr,r,rab,zel,cdd,cdu,cdc,vod,vou,etot,wk1,wk2,
     2 wk3,wk4,wk5,wkb)
c
      do 50 i=1,nr
        vid(i) = vod(i)
        viu(i) = vou(i)
 50   continue
c
c   Start the iteration loop for electronic convergence.
c
      iconv = 0
      icon2 = 0
      maxit = 200
c
c    The screening potential mixing parameter is
c    an empirical function of the nuclear charge.
c    Larger atoms require a slower convergence
c    then smaller atoms.
c
      xmixo = one/log(znuc+7*one)
c
      do 100 iter=1,maxit
        if (iter .eq. maxit) iconv=1
c
c  compute orbitals
c
        if (icon2 .lt. 1) then
          call dsolk1(lmax,nr,a,b,r,rab,norb,ncore,
     1     no,lo,so,zo,cdd,cdu,viod,viou,vid,viu,
     2     ev,vql,inorm,anorm,arw,wk1,wk2,wk3,wk4,
     3     wk5,wk6,wk7,wk8,wk9,wkb)
        endif
        call dsolk2(iter,iconv,ispp,ifcore,lmax,nr,
     1   a,b,r,rab,norb,ncore,no,lo,so,zo,znuc,cdd,
     2   cdu,cdc,viod,viou,vid,viu,ev,ek,ep,vql,
     3   inorm,arw,wk1,wk2,wk3,wk4,wk5,wk6,wk7,evi)
c
c  set up output electronic potential from charge density
c
        if (icon2 .gt. 0) then
        call velect(iter,iconv,icorr,ispp,ifcore,
     1   nr,r,rab,zel,cdd,cdu,cdc,vod,vou,etot,wk1,wk2,
     2   wk3,wk4,wk5,wkb)
c
c  check for convergence
c
        if (iconv .gt. 0) goto 120
        dvmax = zero
        do 60 i=1,nr
          dv = (vod(i)-vid(i))/(1.0+vod(i)+vou(i))
          if (abs(dv) .gt. dvmax) dvmax=abs(dv)
          dv = (vou(i)-viu(i))/(1.0+vou(i)+vod(i))
          if (abs(dv) .gt. dvmax) dvmax=abs(dv)
 60     continue
        if (dvmax .le. tol) iconv=1
        endif
        icon2 = icon2 + 1
c
c    Mix the input and output electronic potentials.
c    The screening potential is initially mixed with a
c    percentage of old and new for itsm iterations.
c    This brings the potential to a stable region
c    after which an Anderson's extrapolation scheme
c    is used.
c
        if (iter .lt. itsm) then
          iiter=2
        else
          iiter=iter-itsm+3
        endif
        call dmixp(vod,vid,xmixo,iiter,3,nr,wk1,wk2,
     1   vn1d,vn11d,vn2d,vn22d)
        call dmixp(vou,viu,xmixo,iiter,3,nr,wk1,wk2,
     1   vn1u,vn11u,vn2u,vn22u)
 100  continue
c
c   End of iteration of electronic convergence loop.
c
      write(6,110) dvmax,xmixo
 110  format(/,' potential not converged - dvmax =',e10.4,
     1 ' xmixo =',f5.3)
       stop
c
c
c  Find the total energy.
c
 120  write(6,121)icon2
 121  format(/,'Total number of iterations needed for',
     1 ' electron screening potential is ',i2,/)
      call etotal(itype,0.0,nameat,norb,
     1 no,lo,so,zo,etot,ev,ek,ep)
c
c   Replace the valence charge density.
c
      if (itype .eq. 5) call vionic(ispp,6,icorr,
     1 ifcore,0.0,0.0,lmax,nr,a,b,r,rab,nameat,
     2 ncore,znuc,cdd,cdu,cdc,viod,viou)
c
      nconf = nconf + 1
c
c   End loop of configuration.
c
      goto 20
      end
C
C
C
       subroutine dsolk1(lmax,nr,a,b,r,rab,norb,ncore,
     1  no,lo,so,zo,cdd,cdu,viod,viou,vid,viu,
     3  ev,vql,inorm,anorm,arw,dk,d,sd,sd2,rv1,rv2,
     4  rv3,rv4,rv5,z)
c
c   dsolk1 finds the non - relativistic wave function
c   using finite differences and matrix diagonalization.
c   An initial guess for the eigenvalues need not be supplied.
c
      implicit real*8 (a-h,o-z)
c
      parameter (zero=0.0,one=1.0,pone=0.1,opf=1.5)
c
      dimension r(nr),rab(nr),no(norb),lo(norb),so(norb),
     1 zo(norb),cdd(nr),cdu(nr),viod(lmax,nr),viou(lmax,nr),
     2 vid(nr),viu(nr),ev(norb),dk(nr),d(nr),sd(nr),sd2(nr),
     3 z(6*nr),rv1(nr),rv2(nr),rv3(nr),rv4(nr),rv5(nr),
cccccc
c mmga
c
     . inorm(lmax),anorm(lmax),arw(nr,norb),vql(nr)
c     4 inorm(lmax),anorm(lmax),arw(1000,5),vql(nr)
c
c mmga
cccccc
c
      dimension nmax(2,5),e(10),ind(10)
c
c   Initialize the charge density arrays.
c
       do 10 i=1,nr
         cdd(i) = zero
         cdu(i) = zero
 10    continue
c
c   Find the max n given l and s.
c   Zero spin is treated as down.
c
      do 20 i=1,2
        do 20 j=1,lmax
          nmax(i,j) = 0
          do 20 k=1,norb
            if (no(k) .le. 0) goto 20
            if (lo(k) .ne. j-1) goto 20
            if ((so(k)-pone)*(i-opf) .lt. zero) goto 20
            nmax(i,j)=no(k)
 20   continue
c
c   Set up hamiltonian matrix for kinetic energy.
c   Only the diagonal depends on the potential.
c
      c2 = -one/b**2
      c1 = -2*one*c2 + one/4
      dk(1)  = c1 / (r(2)+a)**2
      sd(1)  = zero
      sd2(1) = zero
      do 30 i=3,nr
        dk(i-1)  = c1 / (r(i)+a)**2
        sd(i-1)  = c2 / ((r(i)+a)*(r(i-1)+a))
        sd2(i-1) = sd(i-1)**2
 30   continue
c
c   Start loop over spin down=1 and spin up=2.
c
      nrm = nr - 1
      do 80 i=1,2
c
c   Start loop over s p d... states.
c
        do 80 j=1,lmax
          if (nmax(i,j) .eq. 0) goto 80
          llp = j*(j-1)
          prowav = zero
          viod(j,1) = zero
          do 34 k=nr+1,nr+4
            viod(j,k) = zero
 34       continue
          do 35  k=1,nr,4
            prowav = prowav + 7*(viod(j,k)*rab(k)*r(k)**2
     1       +viod(j,k+4)*rab(k+4)*r(k+4)**2)+
     2       32*(viod(j,k+1)*rab(k+1)*r(k+1)**2+
     3       viod(j,k+3)*rab(k+3)*r(k+3)**2)+
     4       12*viod(j,k+2)*rab(k+2)*r(k+2)**2
 35       continue
          prowav = inorm(j)*2*prowav/45
          do 40 k=2,nr
            if (i .eq. 1) then
              d(k-1)=dk(k-1)+vql(k)+llp/r(k)**2+vid(k)+
     1         prowav*viod(j,k)
            else
              d(k-1)=dk(k-1)+vql(k)+llp/r(k)**2+viu(k)+
     1         prowav*viod(j,k)
            endif
 40       continue
c
c   Diagonalize the matrix.
c
          eps = -one
          call tridib(nrm,eps,d,sd,sd2,bl,bu,1,
     1     nmax(i,j),e,ind,ierr,rv4,rv5)
          if (ierr .ne. 0) write(6,50) ierr
 50   format(/,' error in tridib ****** ierr =',i3,/)
          call tinvit(nrm,nrm,d,sd,sd2,nmax(i,j),e,ind,z,ierr,
     1     rv1,rv2,rv3,rv4,rv5)
          if (ierr .ne. 0) write(6,55) ierr
 55   format(/,' error in tinvit ****** ierr =',i3,/) 
c
c   Save the energy levels and add to charge density.
c
          ki = 1
          kn = 0
        do 70 k=1,norb
          if (no(k) .le. 0) goto 70
          if (lo(k) .ne. j-1) goto 70
          if ((so(k)-pone)*(i-opf) .lt. zero) goto 70
          ev(k) = e(ki)
          do 60 l=2,nr
            arw(l,k) = z(kn+l-1) + 1.D-30
            denr = zo(k) * z(kn+l-1)**2 / rab(l)
            if (i .eq. 1) then
              cdd(l) = cdd(l) + denr
            else
              cdu(l) = cdu(l) + denr
            endif
 60       continue
          ki = ki + 1
          kn = kn + nrm
 70     continue
 80   continue
c
c   End loop over s p and d states.
c
      return
      end
C
C
C
      subroutine dsolk2(iter,iconv,ispp,ifcore,lmax,nr,a,b,r,
     1 rab,norb,ncore,no,lo,so,zo,znuc,cdd,cdu,cdc,viod,
     2 viou,vid,viu,ev,ek,ep,vql,inorm,arw,wk1,wk2,wk3,
     3 wk4,v,ar,br,evi)
c
c  dsolk2 finds the non - relativistic wave function using
c  difnrl to intgrate the Scroedinger equation.
c  The energy level from the previous iteration is used
c  as initial guess, and it must therefore be reasonable
c  accurate.
c
      implicit real*8 (a-h,o-z)
c
      character*1 ispp
c
      parameter (zero=0.0,smev=1.E-6)
c
      dimension r(nr),rab(nr),no(norb),lo(norb),so(norb),zo(norb),
     1 cdd(nr),cdu(nr),cdc(nr),viod(lmax,nr),viou(lmax,nr),
     2 vid(nr),viu(nr),ev(norb),ek(norb),ep(norb),
     3 wk1(nr),wk2(nr),wk3(nr),wk4(nr),v(nr),ar(nr),br(nr),
cccccc
c mmga
c
     . inorm(lmax),arw(nr,norb),vql(nr),evi(norb)
c     4 inorm(lmax),arw(1000,5),vql(nr),evi(norb)
c
c mmga
cccccc
c
c  Initialize arrays for charge density.
c
      do 5 i=1,nr
        cdd(i) = zero
 5    continue
      do 10 i=1,nr
        cdu(i) = zero
 10   continue
      if (ifcore .ne. 1) then
        do 15 i=1,nr
          cdc(i)= zero
 15     continue
      endif
c
c  Start the loop over orbitals.
c  Note that spin zero is treated as down.
c
      do 50 i=1,norb
        if (no(i) .le. 0) goto 50
        iback = 0
        prowao = 0.D0
 102    continue
        if (ev(i) .ge. 0.0) ev(i)=-smev
c
c  Set up the potential, set the wave functionc array to zero-ar.
c  Find projector operator factor prowav.
c   
        lp = lo(i)+1
        llp = lp * (lp-1) 
        prowav = zero
        do 75 j=2,nr
          wk3(j) = rab(j)*r(j)*arw(j,i)*viod(lp,j)
 75     continue
        wk3(1) = zero
        do 76 j=nr+1,nr+4
          wk3(j) = zero
 76     continue
        do 70 j=1,nr,4
          prowav = prowav + 7*(wk3(j)+wk3(j+4)) + 
     1     32*(wk3(j+1)+wk3(j+3)) + 12*wk3(j+2)
 70     continue
        prowav = inorm(lp) * 2 * prowav / 45
        if (iback .ge. 1) prowav=(2*prowav+prowao)/3
        print*,'prowav',lo(i),prowav
        do 20 j=2,nr
          v(j) = vql(j)+llp/r(j)**2
     1     +prowav*viod(lp,j)*r(j)/arw(j,i)
 20     continue
c
        do 17 j=1,nr
          ar(j)=zero
 17     continue
        store = viod(lp,2)
        viod(lp,2) = vql(2)+prowav*store*r(j)/arw(2,i) 
        if (so(i) .lt. 0.1) then
          do 18 j=2,nr
            v(j) = v(j) + vid(j)
 18       continue
        else
          do 19 j=2,nr
            v(j) = v(j) + viu(j)
 19       continue
        endif
c
c  Call the integration routine.
c
        iflag=0
        call difnrl(iter,i,v,ar,br,lmax,nr,a,b,r,
     1   rab,norb,no,lo,so,znuc,viod,viou,vid,viu,
     2   ev,iflag,wk1,wk2,wk3,evi)
c
        viod(lp,2) = store
c
c  Store wavefunction
c
         do 60 j=1,nr
           arw(j,i)=(ar(j)+2*arw(j,i))/3+1.D-20
 60      continue
         iback =iback +1
         if (iback .gt. 100) stop 100
         if (iback .le. 1) then
           prowao=prowav
           goto 102
         endif
         if (prowav .ne. zero) then
           if (prowao/prowav .lt. 0.95 .or. prowao/prowav
     1      .gt. 1.05 ) then
             prowao=prowav
             goto 102
           endif
         endif
c
c  Add to the charge density.
c
         if (so(i) .lt. 0.1) then
           do 32 j=1,nr
             denr = zo(i) * ar(j) * ar(j)
             cdd(j) = cdd(j) + denr
 32        continue
         else
           do 33 j=1,nr
             denr = zo(i) * ar(j) * ar(j)
             cdu(j) = cdu(j) + denr
 33        continue
         endif
       if (ifcore .ne. 1 .and. i .le. ncore) then
         do 34 j=1,nr
           denr = zo(i) * ar(j) * ar(j)
           cdc(j)=cdc(j)+denr
 34      continue
       endif
c
c  Compute various quantitities if last iteration.
c
        if (iconv .eq. 1) call orban(ispp,i,ar,br,
     1   lmax,nr,a,b,r,rab,norb,no,lo,zo,so,viod,viou,
     2   vid,viu,ev,ek,ep)
 50   continue
c
c  End loop over orbitals.
c
      return
      end
      subroutine input2(nr,lmax,vql,inorm,anorm)
c
c  input2 reads in the normalization data for the 
c  Kleinman-Bylander operator pseudopotentials
c
      implicit real*8 (a-h,o-z)
c
      dimension inorm(lmax),anorm(lmax),vql(nr)
c
      read(1) (vql(j),j=2,nr)
      read(1) norb
      do 10 i=1,norb
        read(1) inorm(i),anorm(i)
10    continue                        
c
      do 20 i=norb+1,lmax
        inorm(i) = inorm(norb)
        anorm(i) = anorm(norb)
 20   continue
c
      return
      end
