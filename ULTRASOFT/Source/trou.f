c
c Copyright (C) 2002 Vanderbilt group
c This file is distributed under the terms of the GNU General Public
c License as described in the file 'License' in the current directory.
c
      subroutine pswfnormc(psi,rho,vpot,rcut,l,energy,r,
     $     rab,dx,con,mesh,ifprt,idim1,idim2)
c     -----------------------------------
c
c     recherche des coefficients du polynome
      implicit none
      integer l,mesh,ifprt,ndm,idim1,idim2,info
      real*8 psi(idim1), rho(idim1), vpot(idim1), rcut,energy, 
     +       r(idim1), rab(idim1),dx, con(idim2)
      parameter(ndm=1000)
      real*8 vpot1(ndm)
      real*8 r_,rab_,c,c2,amat,y,rc,aenorm
      integer lam,ipvt,ik
      common/radmes/r_(ndm),rab_(ndm)
      common/fundat/ c(6),c2,lam
     +      /funex / amat(6,6),ipvt(6),y(6),rc,aenorm,ik
      integer stderr,input,iout,ioae,iplot,iologd,iops
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops


      real*8 c2l,c2h,precision
      real*8 funz,rtflsp,chip2,r2
      integer i,i1,i2,i3,i4
      external funz,rtflsp,chip2
c
      do i=1,idim1
         r_(i) = r(i)
         rab_(i)= rab(i)
      enddo
      do i=2,idim1
         vpot1(i) = vpot(i)/r(i)
      enddo
      vpot1(1) = 2.*vpot1(2)-vpot1(3)
      do i = 1,6
         c(i) = 0.0
      end do
c     compute the weight under psi
      rho(1) = 0.0d0
      do 10 i = 2,mesh
        rho(i) = psi(i)**2
   10 continue
      call radin(mesh,rho,rab,aenorm,idim1)
c
      if (ifprt.gt.-2) write (iout,20) aenorm
   20 format (' integral under unpseudised psi is ',f15.10)
c
      do i2 = mesh,1,-2
        if ( r(i2) .le. rcut ) goto 110
      enddo
c
  110 continue
c
c     check that the value of i2 is reasonable
      if ( i2 .lt. 3 .or. i2 .gt. mesh-6 ) then
        write(iout,*) '***error in subroutine pswfnormc'
        write(iout,*) 'rcut =',rcut,' generated i2 =',i2,
     +  ' which is illegal with mesh =',mesh
        call exit(0)
      endif
c
      rc = r(i2)
      i1 = i2 - 1
      i3 = i2 + 1
      i4 = i2 + 2
c
      rcut = rc
      if (ifprt.gt.0) then
        write (iout,*) ' i1, i2, i3, i4, kkbeta =',i1,i2,i3,i4,mesh
        write (iout,120) ' r1, r2, r3, r4 =',(r(i),i=i1,i4)
        write (iout,120) 'ps1,ps2,ps3,ps4 =',(psi(i),i=i1,i4)
  120   format(a18,4f8.4)
        write (iout,130) rcut
  130   format(' rcut =',f8.4)
      endif

      ik = i2
      lam = l
c aenorm is the norm out to the cutoff
      aenorm = rho(ik)
      call fill_matrix(amat,rc,lam)
* calcul de mat = lu
      call dgefa(amat,6,6,ipvt,info)
* calcul de y pour ax = y
      call eval_coeff(r,psi,ik,lam,energy,dx,vpot1,y)

      c2h =  5.0          ! borne superieure pour c2
      c2l = -5.0          ! borne inferieure pour c2
      precision = 1e-12
      c2 = rtflsp(funz,c2l,c2h,precision)
* cherche la valeur de  c2 qui verifie l'equ 29a
      con(1) = c(1)
      con(2) = c2
      do i=2,6
         con(i+1)=c(i)
      enddo
c
      psi(1) = 0.0d0
      do  i = 2,ik
         r2 = r(i)*r(i)
         psi(i) = r(i)**(lam+1) * exp(((((((c(6)*r2+c(5))*r2+c(4))*
     $        r2+c(3))*r2+c(2))*r2+c2)*r2)+c(1))
      enddo
      do i=1,mesh
         rho(i) = psi(i)**2
      enddo
      call radin(mesh,rho,rab,aenorm,idim1)
      if (ifprt.gt.-2) then
         write (iout,200) aenorm
 200     format (' integral under pseudised psi is ',f15.10)
      endif
c test for rhodium
      return
      end

c
c ---------------------------------------------------------------
      subroutine fill_matrix(amat,rc,l)
c     --------------------------------
      implicit none
      real*8 amat(6,6),rc
      integer l
c      
c     this routine fills the matrix withn the coefficients taken
c     from  p,p',p'',p''', p(iv), where p is
c     p(r)= c0 + c4 r^4 + c6 r^6 + c8 r^8 + ...
      integer pr1(6),cr1(6),pr(6),cr(6),i,j
c matrices representing the coefficients (cr) and the powers of r ( pr)
      data pr1/0,4,6,8,10,12/
      data cr1/1,1,1,1,1,1/

      do i=1,6
         pr(i) = pr1(i)
         cr(i) = cr1(i)
      end do
      do i = 1,5
c     fill matrix row by row
         do j = 1,6
            amat(i,j) = cr(j) * rc**(pr(j)*1.0)
         end do
c     derivate polynomial expression
         do j = 1,6 
            cr(j) = cr(j) * pr(j)
            pr(j) = max(0,pr(j)-1)
         end do
      end do
      do j = 1,6
         amat(6,j) = 0.0
      end do
      amat(6,2) = 2.0*l +5.0

      return
      end

c ----------------------------------------------------------------
      subroutine eval_coeff (r,psi,ik,l,energy,dx,vpot,y)
c ----------------------------------------------------------------
c    calcule les coefficients dependant de la fct d'onde calculee
c    avec tous les electrons. ces coefficients ervent a la resolution
c    du systeme lineaire.
c    en entree : ik,nx,vpot comme dans le programme principal
c    en sortie : une matrice colonne y contenant les coefficients
c
      implicit none
      integer ik,l,lp1
      real*8 r(ik+3),psi(ik+3),vpot(ik+3),energy,dx,y(6)
      real*8 p,dp,d,vae,dvae,ddvae,rc,rc2,rc3
      real*8   deriv_7pts,deriv2_7pts
      external deriv_7pts,deriv2_7pts,polyn2
      real*8 val,d1,d2
c
c   evaluer p et p' ( dp )
      rc = r(ik)
c      p  = psi(ik)
c      dp = deriv_7pts(psi,ik,rc,dx)
      call polyn2(rc,val,r(ik-6),psi(ik-6),-1)
      call polyn2(rc,val,r(ik-6),psi(ik-6), 0)
      call polyn2(rc, d1,r(ik-6),psi(ik-6),+1)
c      write(*,*) 'dpsi',p,dp,val,d1
c      write(*,*) 'psi',(psi(i),i=ik-1,ik+6)
      p=val
      dp = d1
      if (p.lt.0.0) then
          p  = -p
          dp = -dp
      endif
      d = dp /p
c   evaluer vae, dvae, ddvae
c      vae   = vpot(ik)
c      dvae  =  deriv_7pts(vpot,ik,rc,dx)
c      ddvae = deriv2_7pts(vpot,ik,rc,dx)
      call polyn2(rc,val,r(ik-6),vpot(ik-6),-1)
      call polyn2(rc,val,r(ik-6),vpot(ik-6), 0)
      call polyn2(rc, d1,r(ik-6),vpot(ik-6),+1)
      call polyn2(rc, d2,r(ik-6),vpot(ik-6),+2)
      write(*,*) 'vpd',vae,val,dvae,d1,ddvae,d2
      vae = val
      dvae= d1
      ddvae= d2
c   calcul de parametres intervenant dans les calculs successifs
      lp1 = l + 1
      rc2 = rc * rc
      rc3 = rc2* rc
c
      y(1) = log ( p / rc**lp1 )
      y(2) = dp/p - (lp1 / rc)
      y(3) = vae - energy + (lp1*lp1)/rc2 - d*d
      y(4) = dvae - 2*(vae - energy + l*lp1/rc2)*d - 2*(lp1*lp1)/rc3
     &            + 2*(d*d*d)
      y(5) = ddvae - 2*(dvae - 2*l*lp1/rc3)*d + 6*lp1*lp1/(rc3*rc)
     &       - 2*(vae - energy + l*lp1/rc2 - 3*d*d)*
     &           (vae - energy + l*lp1/rc2 - d*d)

      return
      end
c ----------------------------------------------------------------
       function deriv2_7pts(f,ik,rc,h)
c ----------------------------------------------------------------
c
c      evaluates the second derivative of function f, the function
c      is given numerically on logarithmic mesh r.
c      nm = dimension of mesh
c      ik = integer : position of the point in which the derivative
c      will be evaluated.
c      h is distance between x(i) and x(i+1) where
c      r(i) = exp(x(i))/zed & r(j) = exp(x(j))/zed
c
       implicit none
       integer a(7),i,ik
       real*8 f(ik+3),rc,h,sum,sum1,deriv_7pts,deriv2_7pts
c       integer a(7),i,ik
c      coefficients for the formula in abramowitz & stegun p.914
c      the formula is used for 7 points.
       data a/4,-54,540,-980,540,-54,4/   ! these are coefficients

c formula for linear mesh 
       sum = 0.0
       do i=1,7
          sum = sum + a(i)*f(i-4+ik)
       end do
       sum = 2.0*sum/(720.0*h**2)
c transform to logarithmic mesh
       sum1 = deriv_7pts(f,ik,rc,h) 
       deriv2_7pts = sum/(rc*rc) - sum1 /rc

       return
       end
c ---------------------------------------------------------------     
       function deriv_7pts(f,ik,rc,h)
c ---------------------------------------------------------------
c      evaluates the first derivative of function f, the function
c      is given numerically on logarithmic mesh r.
c      nm = dimension of mesh
c      ik = integer : position of the point in which the derivative
c      will be evaluated.
c      h is distance between x(i) and x(i+1) where
c      r(i) = exp(x(i))/zed & r(j) = exp(x(j))/zed
c
       implicit none
       integer a(7),ik,i
       real*8 f(ik+3),rc,h,sum,deriv_7pts
c       integer a(7),ik,i
c      coefficients for the formula in abramowitz & stegun p.914
       data a/-12,108,-540,0,540,-108,12/

c formula for linear mesh 
       sum = 0
       do i=1,7
          sum = sum + a(i)*f(i-4+ik)
       end do
       deriv_7pts = sum/(720.0*h)
c transform to logarithmic mesh
       deriv_7pts = deriv_7pts /rc

       return
       end
c
c ---------------------------------------------------------------
       real*8 function rtflsp(func,x1,x2,xacc)
c ---------------------------------------------------------------
c       using the false position method, find the root of a 
c      function func, known to lie between x1 and x2.
c       the root, returned as rtflsp, is refined until its accuracy
c      is +- xacc.
c      nb.: routine taken from numerical recipes, fortran edition
c           p 248.
c 
       implicit none
       real*8 func,x1,x2,xacc,fl,fh,zzz,zh,swap,f,del,xl,xh,dx
       integer maxit,ll,j
       external func
       parameter ( maxit = 1000 )


       fl = func(x1)
       fh = func(x2)
       do while (fl*fh.gt.0)       
          print '(/,''   root not bracketed'')'
          print '(''   actual value of xmin xmax   '',2f10.3)',x1,x2
          zzz = x1
          zh = abs(x2 -x1)/10
          print '(/,3x,''     x          f(x)'')'
          do ll = 0,9
             zzz = zzz + zh
             fh = func(zzz)
             print '(3x,2f10.3)',zzz,fh
             if (fl*fh .lt. .0d0) then
                x1 = zzz-zh
                fl = func(x1)
                x2 = zzz
                goto  100
             endif
          end do
          zzz = (x2-x1)
          x2 = x2+zzz
          x1 = x1 - zzz
          fl = func(x1)
          fh = func(x2)         
       end do
 100   if (fl.lt.0) then
            xl = x1
            xh = x2
       else
            xl = x2
            xh = x1
            swap = fl
            fl = fh
            fh = swap
       endif
       dx = xh - xl
       do j = 1,maxit
           rtflsp = xl + dx*fl/(fl-fh)
           f = func(rtflsp)
           if (f.lt.0) then
                del = xl - rtflsp
                xl = rtflsp
                fl = f
           else
                del = xh - rtflsp
                xh = rtflsp
                fh = f
           endif
           dx = xh - xl
           if (abs(del).lt.xacc.or.f.eq.0) return
       end do

       stop ' rtflsp max iter execeeded'
       end

c --------------------------------------------------------
      function funz(x)
c --------------------------------------------------------
c    cette fonction est la fonction de c2 qu'il faut annuler
c    pour trouver une valeur de c2 qui verifie la conservation 
c    de la charge de coeur.
c    cette fonction calcule de plus a chaque fois les ci qui 
c    verifient les equations lineaires avec un c2 donne.
c    en entree : une valeur de c2 donnee dans x
c    en sortie, la valeur de la fonction pour cette valeur de x
c    : c'est la fonction qui correspond a l'equation integrale
c
      implicit none
      integer ndm,lam,ipvt,ik,i
      real*8 funz, x, c, c2, amat, y, rc, aenorm, r,rab,
     $        chip2, psnorm
      parameter(ndm=1000)
      common/radmes/ r(ndm),rab(ndm)
      common /fundat/c(6),c2,lam
     +       /funex/ amat(6,6),ipvt(6),y(6),rc,aenorm,ik
      external chip2,radin
      real*8 fint(ndm)
c    resolution du systeme lineaire pour cette valeur de x (=c2)
      c(1) = y(1) - x*rc**2  ! ajoute les termes en c2
      c(2) = y(2) - 2*x*rc
      c(3) = y(3) - 2*x
      c(4) = y(4)       ! pas de coeff en c2
      c(5) = y(5)
      c(6) = -x*x
      call dgesl(amat,6,6,ipvt,c,0)    ! resoud le systeme
c    calcul de la norme de la pseudo-fonction d'onde
      psnorm = 0.0
      fint(1) = 0
      do i=2,ik
         fint(i) = chip2(c,x,lam,r(i))
      enddo
      call radin(ik,fint,rab,psnorm,ndm)
      funz = log( psnorm / aenorm )
c
      return
      end

c --------------------------------------------------------
      function chip2(c,c2,l,r)
c --------------------------------------------------------
      implicit none
      real*8 chip2,c(6),c2,r,r2
      integer l

      r2 = r**2
      chip2 = r2**(l+1) * exp(2.0*
     * (((((((c(6)*r2+c(5))*r2+c(4))*r2+c(3))*r2+c(2))*r2+c2)*r2)+c(1)))

      return
      end
c
c -------------------------------------------------------------------

c ----------------------------------------------------------------
       function der3an(l,c,c2,rc)
c ----------------------------------------------------------------
c
c 3rd derivative of r^(l+1) e^p(r)
c
       implicit none
       integer l
       real*8  c(6), c2, rc,pr,dexpr,d2expr,d3expr, der3an
       external pr,dexpr,d2expr,d3expr

       der3an= (l-1)*l*(l+1)*rc**(l-2) *exp(pr(c,c2,rc)) +
     $        3 *    l*(l+1)*rc**(l-1) * dexpr(c,c2,rc)  +
     $        3 *      (l+1)*rc**(l  ) *d2expr(c,c2,rc)  +
     $                       rc**(l+1) *d3expr(c,c2,rc) 

       return
       end

c ----------------------------------------------------------------
       function der4an(l,c,c2,rc)
c ----------------------------------------------------------------
c
c 4rd derivative of r^(l+1) e^p(r)
c
       implicit none
       integer l
       real*8  c(6), c2, rc,pr,dexpr,d2expr,d3expr,d4expr,der4an
       external pr,dexpr,d2expr,d3expr,d4expr

       der4an = (l-2)*(l-1)*l*(l+1)*rc**(l-3) * exp(pr(c,c2,rc)) +
     $         4 *    (l-1)*l*(l+1)*rc**(l-2) *  dexpr(c,c2,rc)  +
     $         6 *          l*(l+1)*rc**(l-1) * d2expr(c,c2,rc)  +
     $         4 *            (l+1)*rc**(l  ) * d3expr(c,c2,rc)  +
     $                              rc**(l+1) * d4expr(c,c2,rc)

       return
       end
c ----------------------------------------------------------------
       function dexpr(c,c2,rc)
c ----------------------------------------------------------------
c
c 1st derivative of e^p(r)
c
       implicit none
       real*8  c(6), c2, rc, pr, dpr, dexpr
       external pr, dpr

       dexpr= exp(pr(c,c2,rc)) * dpr(c,c2,rc)

       return
       end
c ----------------------------------------------------------------
       function d2expr(c,c2,rc)
c ----------------------------------------------------------------
c
c 2nd derivative of e^p(r)
c
       implicit none
       real*8  c(6), c2, rc, pr, dpr, d2pr, d2expr
       external pr, dpr, d2pr

       d2expr= exp(pr(c,c2,rc)) * ( dpr(c,c2,rc)**2+d2pr(c,c2,rc) )

       return
       end
c ----------------------------------------------------------------
       function d3expr(c,c2,rc)
c ----------------------------------------------------------------
c
c 3nd derivative of e^p(r)
c
       implicit none
       real*8  c(6), c2, rc, pr, dpr, d2pr, d3pr, d3expr
       external pr, dpr, d2pr, d3pr

       d3expr= exp(pr(c,c2,rc)) * (   dpr(c,c2,rc)**3
     $                            +  3*dpr(c,c2,rc)*d2pr(c,c2,rc)
     $                            +   d3pr(c,c2,rc)               )

       return
       end
c ----------------------------------------------------------------
       function d4expr(c,c2,rc)
c ----------------------------------------------------------------
c
c 4th derivative of e^p(r)
c
       implicit none
       real*8  c(6), c2, rc, pr,  dpr, d2pr, d3pr, d4pr, d4expr
       external pr, dpr, d2pr, d3pr, d4pr

       d4expr= exp(pr(c,c2,rc)) * (     dpr(c,c2,rc)**4
     $                            +  6* dpr(c,c2,rc)**2*d2pr(c,c2,rc)
     $                            +  3*d2pr(c,c2,rc)**2
     $                            +  4* dpr(c,c2,rc)*d3pr(c,c2,rc)
     $                            +    d4pr(c,c2,rc)                  )

       return
       end
c
c --------------------------------------------------------
      function pr(c,c2,x)
c --------------------------------------------------------
c     cette fonction evalue le polynome dont les coefficients
c     sont dans c d'apres la forme suivante
c     p(x) =  c(2)*x^4 + c(3) x^6 + c(4) x^8 + c(5) x^10
c            +c(6) x^12 + c2*x^2 + c(1)
      implicit none
      real*8 c(6), c2, x, y, pr

      y = x*x
      pr = (((((c(6)*y+c(5))*y+c(4))*y+c(3))*y+c(2))*y+c2)*y
     +          + c(1)

      return
      end
c
      function dpr(c,c2,x)
c --------------------------------------------------------
      implicit none
      real*8 c(6),c2,x,dpr

      dpr  = 2*c2*x + 4*c(2)*x**3 + 6*c(3)*x**5 + 8*c(4)*x**7
     &     + 10*c(5)*x**9 + 12*c(6)*x**11

      return
      end
c 
      function d2pr(c,c2,x)
c --------------------------------------------------------
      implicit none
      real*8 c(6),c2,x,d2pr

      d2pr = 2*c2 + 12*c(2)*x**2 + 30*c(3)*x**4 + 56*c(4)*x**6
     &     + 90*c(5)*x**8 + 132*c(6)*x**10

      return
      end
c
      function d3pr(c,c2,x)
c --------------------------------------------------------
      implicit none
      real*8 c(6),c2,x,d3pr

      d3pr  = 24*c(2)*x + 120*c(3)*x**3 + 336*c(4)*x**5
     &     + 720*c(5)*x**7 + 1320*c(6)*x**9

      return
      end
c 
      function d4pr(c,c2,x)
c --------------------------------------------------------
      implicit none
      real*8 c(6),c2,x,d4pr

      d4pr  = 24*c(2) + 360*c(3)*x**2 + 1680*c(4)*x**4
     &     + 5040*c(5)*x**6 + 11880*c(6)*x**8

      return
      end

      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
      INTEGER LDA,N,IPVT(1),INFO
      DOUBLE PRECISION A(LDA,1)
C
C     DGEFA FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION.
C
C     DGEFA IS USUALLY CALLED BY DGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR DGEFA) .
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT DGESL OR DGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN DGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DSCAL,IDAMAX
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION T
      INTEGER IDAMAX,J,K,KP1,L,NM1
C
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      END

      SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)
      INTEGER LDA,N,IPVT(1),JOB
      DOUBLE PRECISION A(LDA,1),B(1)
C
C     DGESL SOLVES THE DOUBLE PRECISION SYSTEM
C     A * X = B  OR  TRANS(A) * X = B
C     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE OUTPUT FROM DGECO OR DGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM DGECO OR DGEFA.
C
C        B       DOUBLE PRECISION(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE
C                            TRANS(A)  IS THE TRANSPOSE.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF DGECO HAS SET RCOND .GT. 0.0
C        OR DGEFA HAS SET INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DDOT
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,NM1
C
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE  L*Y = B
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
         DO 60 K = 1, N
            T = DDOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
C
C        NOW SOLVE TRANS(L)*X = Y
C
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END

* ======================================================================
* NIST Guide to Available Math Software.
* Fullsource for module DSCAL from package BLAS1.
* Retrieved from NETLIB on Fri Jun 21 17:32:50 1996.
* ======================================================================
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
