       subroutine dsolv1(lmax,nr,a,b,r,rab,norb,ncore,
     1  no,lo,so,zo,cdd,cdu,viod,viou,vid,viu,
     3  ev,dk,d,sd,sd2,rv1,rv2,rv3,rv4,rv5,z)
c
c   dsolv1 finds the (non)-relativistic wave function
c   using finite differences and matrix diagonalization.
c   An initial guess for the eigenvalues need not be supplied.
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
      parameter (zero=0.D0,one=1.D0,pone=0.1D0,opf=1.5D0)
Cray      parameter (zero=0.0,one=1.0,pone=0.1,opf=1.5)
c
      dimension r(nr),rab(nr),no(norb),lo(norb),so(norb),
     1 zo(norb),cdd(nr),cdu(nr),viod(lmax,nr),viou(lmax,nr),
     2 vid(nr),viu(nr),ev(norb),dk(nr),d(nr),sd(nr),sd2(nr),
     2 z(6*nr),rv1(nr),rv2(nr),rv3(nr),rv4(nr),rv5(nr)
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
          do 40 k=2,nr
            if (i .eq. 1) then
              d(k-1)=dk(k-1)+(viod(j,k)+llp/r(k))/r(k)+vid(k)
            else
              d(k-1)=dk(k-1)+(viou(j,k)+llp/r(k))/r(k)+viu(k)
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
