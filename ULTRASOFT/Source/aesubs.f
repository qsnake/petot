c
c Copyright (C) 2002 Vanderbilt group
c This file is distributed under the terms of the GNU General Public
c License as described in the file 'License' in the current directory.
c
c-------------------------------------------------------------------------
c
      subroutine schsol(ruae,r,rab,sqr,a,b,v,z,ee,
     +  nnlz,wwnl,nkkk,snl,snlo,ncspvs,mesh,thresh,ifprt,idim1,idim2)
c
c     routine to compute non relativistic solutions to
c     schroedinger equation
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1),sqr(idim1)
c.....potentials for the all electron calculation
      dimension ruae(idim1),v(idim1)
c.....wavefunctions
      dimension snl(idim1,idim2),snlo(idim1)
c.....shell labels, energies, occupancies and real space cutoff
      dimension nnlz(idim2),ee(idim2),wwnl(idim2),nkkk(idim2)
c.....scratch
      dimension yval(-1:1)
c
c.....variable for file = 0
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c--------------------------------------------------------------------------
c
c               r o u t i n e  i n i t i a l i s a t i o n
c
c     set the maximum number of iterations for improving wavefunctions
      itmax = 50
c
      do 10 j = 1,ncspvs
      do 10 i = 1,mesh
        snl(i,j) = 0.0d0
  10  continue
c
c--------------------------------------------------------------------------
c
c               m a i n  l o o p  o v e r  o r b i t a l s
c
      do 100 iorb = 1,ncspvs
c
c       check that this orbital is required
        if ( wwnl(iorb) .lt. 0.0d0 ) then
          ee(iorb) = 0.0d0
          nkkk(iorb) = 1
          snl(1,iorb) = 0.0d0
          goto 100
        endif
c
c       compute the quantum numbers
        ncur = nnlz(iorb) / 100
        lcur = nnlz(iorb) / 10 - ncur * 10
        llp = lcur * ( lcur + 1 )
c
c       set the starting energies for the iterative loop
        ecur = ee(iorb)
        emax = + 1.0d0
        emin = - 1.0d+10
c
c       try uncommenting the "write(42..." statements if you
c       encounter convergence problems in the following loop
c
c       write(42,'(/a,i5/)') 'nnlz',nnlz(iorb)
c
c-------------------------------------------------------------------------
c
c            i t e r a t i o n  o n  e i g e n f u n c t i o n s
c
        do 200 iter = 1,itmax
c
c         compute the sums of potential and centrifugal terms
c         required for the radial integrator
c
          call setv(r,rab,ruae,b,ecur,llp,v,mesh,idim1)
c
c----------------------------------------------------------------------------
c
c         c l a s s i c a l  t u r n i n g  a n d  i n f i n i t y
c
c         first attempt to locate the classical turning point
          do 210 i = mesh,2,-1
            if ( v(i) .lt. 0.0d0 ) goto 220
            if ( i .eq. 2 ) then
c             if the energy is reasonable this point should not be reached
              write(iout,*) '***error in subroutine schsol'
              write(iout,500) nnlz(iorb),iter,ecur,mesh
  500         format(' no tuning point fot state nnlz =',i4,
     +        ' on iteration',i3,/,' with ecur = ',f10.5,
     +        ' and v(mesh) = ',f10.5)
              call exit(1)
            endif
  210     continue
c
  220     continue
c
          if ( i .gt. mesh - 20 ) then
c           state unbound or close to being unbound
            nctp = mesh - 20
            ninf = mesh
          else
c           state well bound so i is nctp
            nctp = i
c           find the pracitcal infinity
            tolinf = 1.0d-18
            tolinf = dlog(tolinf) ** 2
            do 230 i = nctp,mesh
              alpha2 = v(i) * (dble(i-nctp)**2)
              if ( alpha2 .gt. tolinf ) goto 240
  230       continue
c
  240       continue
            ninf = i
          endif
c
          if ( ninf .eq. mesh ) then
            write(iout,*) '***warning in subroutine schsol'
            write(iout,540) nnlz(iorb),ecur,iter
  540       format(' state ',i4,' with energy',f9.5,
     +      ' either unbound or requires larger',/,
     +      ' radial mesh with energy on iteration',i4)
          endif
c
c         if ( ifprt .eq. 1 ) write(iout,520) nnlz(iorb),
c    +      iter,nctp,ninf
  520     format(' state nnlz=',i4,' on iteration',i3,' required',
     +    ' nctp =',i5,' and ninf =',i5)
c
c---------------------------------------------------------------------------
c
c               o u t w a r d  i n t e g r a t i o n 
c
c         compute the starting value of psi using taylor expansion 
          zz = 0.0d0
          vzero = ruae(2) / r(2)
          call orgsnl(r,snlo,1,3,lcur,zz,vzero,ecur,idim1)
c
c         convert to the phi function
          snlo(1) = snlo(1) / sqr(1)
          snlo(2) = snlo(2) / sqr(2)
          snlo(3) = snlo(3) / sqr(3)
c
c         call the integration routine
          call difsol(v,snlo,4,nctp,idim1)
c
          snlctp = snlo(nctp)
c
c--------------------------------------------------------------------------
c
c               i n w a r d  i n t e g r a t i o n
c
          if ( ninf .lt. mesh ) then
c           exponential start up of the wavefunction
            alpha = sqrt(v(ninf))
            snlo(ninf) = exp( - alpha * dble(ninf - nctp ) )
            alpha = sqrt(v(ninf-1))
            snlo(ninf-1) = exp ( - alpha * dble(ninf - nctp - 1) )
          else
c           linear start up
            snlo(ninf) = 0.0d0
            snlo(ninf-1) = r(ninf) - r(ninf-1)
          endif
c
          call difsol(v,snlo,ninf-2,nctp,idim1)
c
c--------------------------------------------------------------------------
c
c         c h e c k  t h a t  t h e  n u m b e r  o f  n o d e s  o k
c
c         renormalise snlo so continuous at nctp
          factor = snlctp / snlo(nctp)
          do 350 i = nctp,ninf
            snlo(i) = snlo(i) * factor
  350     continue
c
c         count the number of nodes
          call nodeno(snlo,2,ninf,nodes,idim1)
c
c--------------
c
          if ( nodes .lt. ncur - lcur - 1 ) then
c           energy is too low
            emin = ecur
c           write(42,'(i5,3f12.5,2i5)')
c    +         iter,emin,ecur,emax,nodes,ncur-lcur-1
            if ( ecur * 0.9d0 .gt. emax ) then
              ecur = 0.5d0 * ecur + 0.5d0 * emax 
            else
              ecur = 0.9d0 * ecur
            endif
            goto 200
          endif
c
          if ( nodes .gt. ncur - lcur - 1 ) then
c           energy is too high
            emax = ecur
c           write(42,'(i5,3f12.5,2i5)')
c    +         iter,emin,ecur,emax,nodes,ncur-lcur-1
            if ( ecur * 1.1d0 .lt. emin ) then
              ecur = 0.5d0 * ecur + 0.5d0 * emin
            else
              ecur = 1.1d0 * ecur
            endif
            goto 200
          endif
c
c------------------------------------------------------------------------------
c
c         v a r i a t i o n a l  i m p ro v e m e n t  o f  e c u r
c
c         find normalisation of psi using simpson's rule
          factor = 0.0d0
          xw = 4.0d0
          do 360 i = 2,ninf
            factor = factor + xw * rab(i) * (snlo(i)*sqr(i))**2
            xw = 6.0d0 - xw
  360     continue
          factor = factor / 3.0d0
c
c         compute the noumerov discontinuity in the phi function
          do 370 i = -1,1
            yval(i) = snlo(nctp+i) * ( 1.0d0 - v(nctp+i) / 12.0d0 )
  370     continue
          d2p = v(nctp) * snlctp
          disc = - yval(-1) + 2.0d0*yval(0) - yval(1) + d2p
c
c         variational estimate for the change in ecur
          decur = snlctp*disc*sqr(nctp)**2 / (factor*rab(nctp))
c
c         to prevent convergence problems:
c         do not allow decur to exceed 20% of | ecur |
c         do not allow decur to exceed 70% of distance to emin or emax
          if (decur.gt.0.) then
            emin=ecur
            decurp=min(decur,-0.2*ecur,0.7*(emax-ecur))
          else
            emax=ecur
            decurp=-min(-decur,-0.2*ecur,0.7*(ecur-emin))
          endif
c
c         write(42,'(i5,3f12.5,1p2e12.4)')
c    +         iter,emin,ecur,emax,decur,decurp
c
c         test to see whether eigenvalue converged
          if ( abs(decur) .lt. thresh ) goto 400
c
          ecur = ecur + decurp
c
          if ( iter .eq. itmax ) then
c           eigenfunction has not converged in allowed number of iterations
            write(iout,*) '***error in subroutine schsol'
            write(iout,530) nnlz(iorb),ee(iorb),ecur
  530       format(' state',i4,' could not be converged.',/,
     +      ' starting energy for calculation was',f15.5,
     +      ' and end value =',f15.5)
            call exit(1)
          endif
c
c       close iteration loop
  200   continue
c
c-----------------------------------------------------------------------------
c
  400   continue
c
c       update the energy
        ee(iorb) = ecur
c       update the real-space cutoff array
        nkkk(iorb) = ninf
c
c       update the wavefunction array
        factor = 1.0d0 / sqrt( factor )
        do 410 i = 1,ninf
          snl(i,iorb) = snlo(i) * sqr(i) * factor
  410   continue
c
c     close loop over bands
  100 continue
c
      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine difsol(g,p,jj1,jj2,idim1)
c
c---------------------------------------------------------------------------
c     routine for solving second order differential equation of the form
c
c         2
c        d p(x)
c        ------  = g(x) p(x)
c           2
c         dx
c
c     using the method due to noumerov.
c     method taken from g.w. pratt, phys. rev. 88 p.1217 (1952).
c
c     routine integrates from jj1 to jj2 and can cope with both cases
c     jj1 < jj2 and jj1 > jj2.  first two starting values of p must be
c     provided by the calling program.
c---------------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
c
      dimension g(idim1),p(idim1)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c     i n i t i a l i s a t i o n
c
      r12 = 1.0d0 / 12.0d0
c
c     decide whether integrating from:
c            left to right ---> isgn = + 1
c        or  right to left ---> isgn = - 1
c
      isgn = ( jj2 - jj1 ) / iabs( jj2 - jj1 )
c
c     run some test just to be conservative
      if ( isgn .eq. + 1 ) then
        if ( jj1 .le. 2 .or. jj2 .gt. idim1 ) then
          write(iout,10) isgn,jj1,jj2,idim1
          call exit(1)
        endif
      elseif ( isgn .eq. - 1 ) then
        if ( jj1 .ge. ( idim1 - 1 ) .or. jj2 .lt. 1 ) then
          write(iout,10) isgn,jj1,jj2,idim1
          call exit(1)
        endif
      else
        write(iout,10) isgn,jj1,jj2,idim1
      endif
c
  10  format(' ***error in subroutine difsol',/,
     +' isgn =',i2,' jj1 =',i5,' jj2 =',i5,' idim1 =',i5,
     +' are not allowed')
c
c     initialise current second derivative of p (d2p), ycur and yold
c
      i = jj1 - isgn
      d2p = g(i) * p(i)
      ycur = ( 1.0d0 - r12 * g(i) ) * p(i)
      i = i - isgn
      yold = ( 1.0d0 - r12 * g(i) ) * p(i)
c
c
c     b e g i n  i n t e g r a t i o n  l o o p
c
      do 100 i = jj1,jj2,isgn
        ynew = 2.0d0 * ycur - yold + d2p
        p(i) = ynew / ( 1.0d0 - r12 * g(i) )
        d2p  = g(i) * p(i)
        yold = ycur
        ycur = ynew
100   continue
c
      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine setv(r,rab,ruae,b,ecur,llp,v,mesh,idim1)
c
c     compute the sums of potential and centrifugal terms
c     required for the radial integrator
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1)
c.....potentials for the all electron calculation
      dimension ruae(idim1),v(idim1)
c
c
      v(1) = 0.0d0
      do 100 i = 2,mesh
        v(i) = b * b / 4.0d0 + rab(i) * rab(i) *
     +    ( ruae(i) / r(i) + dble(llp) / ( r(i) * r(i) ) - ecur )
  100 continue
c
      return
      end
c
c-------------------------------------------------------------------------
c
      subroutine vhxc(ipass,ifpcor,ruexch,ruhar,reexch,rscore,
     +  rsvale,rsatom,rspsco,rstot,exfact,r,rab,
     +  xi,xj,fu,ehar,emvxc,mesh,idim1)
c
c     calculate hartree and exchange-correlation potentials
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1)
c.....hartree, exchange and correlation potentials
      dimension ruhar(idim1),reexch(idim1),ruexch(idim1)
c.....charge densities
      dimension rscore(idim1),rsvale(idim1),rsatom(idim1),rspsco(idim1)
c.....scratch arrays
      dimension xi(idim1),xj(idim1),fu(idim1),rstot(idim1)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
      data pi    /3.141592654/

      nneg=0
      do 100 i=2,mesh
c       charge depends on whether this is ae or pseudo calculation
        if ( ipass .eq. 1 ) then
          rsatom(i) = rscore(i) + rsvale(i)
        elseif ( ipass .eq. 2 ) then
          rsatom(i) = rsvale(i)
        endif
c
c       section of code added 12/19/91 to prevent problems of
c       pseudo-ion having incorrect charge
c
        rstot(i)=rsatom(i)
        if (ipass.eq.2.and.ifpcor.gt.0) rstot(i)=rsatom(i)+rspsco(i)
c
c       take measures to ensure that rho is greater than zero
        if ( rstot(i) .lt. 1.0d-30 ) then
c         check if charge density is negative
          if ( rstot(i) .lt. -1.0d-30 ) then
            write(iout,*) 'rstot =',rstot(i),' at r =',r(i)
            nneg=nneg+1
          endif
        endif
        fu(i) = rstot(i)/(4.*pi*r(i)**2)
 100  continue
      if (nneg.ne.0) then
        write (iout,*) ' vhxc: negative density at',nneg,' r-values'
        write (stderr,*) '*** stopping in vhxc: negative densities ***'
        call exit(1)
      endif
      fu(1) = fu(3) - (fu(4)-fu(2))/2./rab(3)*(r(3)-r(1))
c      write(*,*)  'test',fu(1),fu(2),fu(3),fu(4)
c specify cut-off density         
c specify cutoff density below which the derivatives are exponentially damped
      rhoc1 = 1.d-10
c specify cutoff density below which the derivatives are zero
      rhoc2 = 1.d-20

      do i=2,mesh
         if (fu(i) .gt. rhoc1) iexp1 = i
         if (fu(i) .gt. rhoc2) iexp2 = i
      enddo
      do i=1,mesh
         xj(i) = rstot(i)
      enddo
      call radin(mesh,xj,rab,asum,idim1)
c     fu contains the real density, now calculate the derivatives
      call deriv(fu,xi,r,rab,mesh,idim1)
c   xi the 1. derivative
c   xj the 2. derivative, calculate later
c   first fit exponential tail to the density
       if(iexp1 .lt. mesh-2) then   
         bexp = xi(iexp1)/fu(iexp1) +  xi(iexp1-1)/fu(iexp1-1) +
     +      xi(iexp1+1)/fu(iexp1+1)
         bexp = bexp /3.
         x = -bexp*r(iexp1)
         aa=-1.*(xj(mesh)-xj(iexp1))/((x*x+2.*x+2.)*
     +      exp(-x)*4.*pi)*bexp**3
      endif
      call deriv(xi,xj,r,rab,mesh,idim1)
      if(iexp1 .lt. mesh-2) then
         do i=iexp1,iexp2
            fu(i) =  aa * exp(bexp *r(i))
            xi(i) = bexp * aa* exp(bexp * r(i))
            xj(i) = bexp * xi(i)
         enddo
      endif

c     rhoc2 = fu(iexp2)
      do i=iexp2,mesh
         xi(i) = 0.0d0
         xj(i) = 0.0d0
      enddo
      
      call xctype(exfact,rhoc2,0.0d0,0.0d0,r(iexp2),uxcc,excc)
      coef1 = ( 3.d0 * excc - uxcc ) / rhoc2
      coef2 = (-2.d0 * excc + uxcc ) / rhoc2**2
c
c     calculate total charge density and exchange-correlation
      ruexch(1)=0.0d0
      reexch(1)=0.0d0
      do i=2,mesh
         rho=fu(i)
        if (rho.gt.rhoc2) then
           call xctype(exfact,rho,xi(i),xj(i),r(i),uxc,exc)
        else
          exc =     coef1 * rho +     coef2 * rho**2
          uxc = 2 * coef1 * rho + 3 * coef2 * rho**2
        endif
c
        ruexch(i)=r(i)*uxc
        reexch(i)=r(i)*exc
      enddo

c     ruexch contains exchange-correlation potential
c     reexch contains exchange-correlation energy
c
c     calculate atomic coulomb potential
      xi(1)=0.0d0
      xj(1)=0.0d0
      do 110 i=2,mesh
        xi(i)=rsatom(i)
        xj(i)=rsatom(i)/r(i)
  110 continue
      call radin(mesh,xi,rab,xim,idim1)
      call radin(mesh,xj,rab,xjmesh,idim1)
c
      do 120 i=1,mesh
        rehar=xi(i)+r(i)*(xjmesh-xj(i))
        ruhar(i)=rehar*2.0
  120 continue
c
      fu(1)=0.
      do 130 i=2,mesh
        fu(i)=0.5*ruhar(i)/r(i)*rsatom(i)
  130 continue
      call radin(mesh,fu,rab,ehar,idim1)
      fu(1)=0.
      do 140 i=2,mesh
        fu(i)=( reexch(i)*rstot(i) - ruexch(i)*rsatom(i) )/r(i)
  140 continue
      call radin(mesh,fu,rab,emvxc,idim1)
c
      return
      end
c
c     ****************************************************************
      subroutine xctype(alphax,rho,rho1,rho2,r,uxc,exc)
c     ****************************************************************
c
c ---------------------------------------------------------------
      implicit double precision (a-h,o-z)
      real*8  LAPrho
c ---------------------------------------------------------------
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c ---------------------------------------------------------------
c
c     subroutine to calculate exchange-correlation potential and energy
c
c     alphax : key to type of exchange-correlation
c        alphax =  0. -->  perdew-zunger param of ceperley-alder
c        alphax = -1. -->  wigner interpolation
c        alphax = -2. -->  hedin-lundquist
c        alphax = -3. -->  gunnarson-lundquist
c        alphax = -4. -->  no exchange or correlation at all
c                  1  LDA (CA) + Becje88 exchange +  LYP correlation
c                  2  LDA (CA) + Becke88 exchange 
c                  3  LDA (CA) + Becke88 exchange + Perdew86 correlation
c                  4  LDA (CA) + PW(91)
c                  5  PBE
c

c
c     rho = density in electrons per au**3 (input)
c
c     uxc = exchange-correlation potential in ry (output)
c     exc = exchange-correlation energy    in ry (output)
c
      data pi    /3.141592654/
      data third /0.33333333333333333/
      data qq    /0.620350491/
      data cc    /0.916330586/
c
c     qq = (3./(4.*pi))**(1./3.)
c     cc=pi*(1.5/pi)**(5./3.)
c
c     perdew-zunger exchange-correlation
      data beta1 /1.0529/
      data beta2 /0.3334/
      data gamma /-0.2846/
      data pza /0.0622/
      data pzb /-0.096/
      data pzc /0.0040/
      data pzd /-0.0232/
c     wigner exchange-correlation
      data aa /0.88/
      data bb /7.8/
c     hedin-lundquist exchange-correlation
      data ah /21./
      data ch /0.045/
c     gunnarson-lundquist exchange-correlation
      data ag /11.4/
      data cg /0.0666/
c
      rs=qq/rho**third
      ex=-cc/rs
      ux=4.0*ex/3.0
c     ex and ux are exchange energy and potential respectively
c
      if (alphax.ge. 0.) then
c
c       perdew-zunger correlation
        if (rs.gt.1.) then
          dd=1.0+beta1*sqrt(rs)+beta2*rs
          dn=1.0+(7./6.)*beta1*sqrt(rs)+(4./3.)*beta2*rs
          ec=gamma/dd
          ecrs=-1.*gamma/dd/dd*(.5*beta1/sqrt(rs)+beta2)
          uc=ec*dn/dd
        else
          rl=log(rs)
          ec= pza*rl + pzb + pzc*rs*rl + pzd*rs
          ecrs=pza/rs+pzc*(rl+1.)+pzd
          uc= pza*rl + (pzb-pza/3.) + (2.*pzc/3.)*rs*rl
     1       + ((2*pzd-pzc)/3.)*rs
        endif
        uxc=ux+uc
        exc=ex+ec
c
      else if (alphax.eq.-1.) then
c
c       wigner correlation
        ec=-aa/(rs+bb)
        ecmuc=rs*aa/(3.0*(rs+bb)**2)
        uxc=ux+ec-ecmuc
        exc=ex+ec
c
      else if (alphax.eq.-2.) then
c
c       hedin-lundquist correlation
        xx=rs/ah
        al=log(1.+1./xx)
        uc=-ch*al
        ecmuc=-ch*(al*xx**3+xx/2.-xx*xx-1./3.)
        uxc=ux+uc
        exc=ex+uc+ecmuc
c
      else if (alphax.eq.-3.) then
c
c       gunnarson-lundquist correlation
        xx=rs/ag
        al=log(1.+1./xx)
        uc=-cg*al
        ecmuc=-cg*(al*xx**3+xx/2.-xx*xx-1./3.)
        uxc=ux+uc
        exc=ex+uc+ecmuc
c
      else if (alphax.eq.-4.) then
c
c       no exchange or correlation at all (for debugging only)
        uxc=0.
        exc=0.
c
c      else if (alphax.gt. 0.) then
c
c       slater exchange-correlation
c       alphax = 2/3 corresponds to pure exchange
c        uxc=ux * 1.5 * alphax
c        exc=ex * 1.5 * alphax
c
      else 
c
        write(iout,*) '***error in subroutine xctype'
        write(iout,*) 'alphax = ',alphax,' is out of range'
        call exit(1)
c
      endif
      rmin = 1.d-20
c      write(6,*)rmin, rho1,rho2
c      if (( rho1 .le. 1.d0-40) .and. (rho1 .ge. -1.d0-40)) rho1=0.d0
c      if (abs(rho2) .le. 1.0d-40) rho2=0.d0
      if (alphax.gt. 0. .and. (abs(rho1) + abs(rho2)).gt.rmin) then
c add gradient corrections
         ecg=0.d0
         exg=0.d0
         ucg=0.d0
         uxg=0.d0


         LAPrho=rho2 + 2.d0*rho1/R
         rk = rho1
         rk2=rk*rk
         rnorm=abs(rk)

C
C 
C GRADIENT CORRECTIONS TO THE CORRELATION POTENTIAL  --  Lee-Yang-Parr
C   
      if(alphax.eq.1.) then
         A=0.04918d0 
         B=0.132d0 
         C=0.2533d0 
         D=0.349d0 
         CF=2.871233d0          !  3/10 (3pi^2)^(2/3)
        T1 = rho**third
        T2 = 1/T1
        T4 = 1+D*T2
        T5 = 1/T4
        T6 = T1**2
        T7 = 1/T6
        T8 = T6**2
        T9 = T8*T1
        T11 = rho1**2
        T16 = CF*T9-17.D0/72.D0*T11/rho+7.D0/24.D0*LAPrho
        T19 = EXP(-C*T2)
        T20 = T16*T19
        T23 = rho+B*T7*T20
        T27 = T4**2
        T28 = 1/T27
        T35 = 1/T9
        T41 = rho**2
        T42 = 1/T41
        T49 = C*T19
        T59 = rho1*T19
        T63 = A*T5*B*T35*T59
        T71 = 1/T41/rho
        T73 = T11*T19*D
        T77 = A*T28*B*T71*T73
        T79 = T8**2
        T81 = T11*T19
        T85 = A*T5*B/T79*T81
        T87 = rho2*T19
        T91 = A*T5*B*T35*T87
        T97 = A*T5*B*T71*T11*T49
        T118 = 1/T79/T6
        T119 = D**2
        T150 = C**2
        ecg = -A*T5*T23
        ecg = ecg/rho 

        DNF = -A*T28*T23*D/T8/3
     $      -A*T5*(1-2.D0/3.D0*B*T35*T20+B*T7*(5.D0/3.D0*CF*T6
     $      +17.D0/72.D0*T11*T42)*T19+B*T42*T16*T49/3)
        DN1F = 17.D0/36.D0*T63
        DNLF = -7.D0/24.D0*A*T5*B*T7*T19
        DRN1F = 17.D0/108.D0*T77-85.D0/108.D0*T85+17.D0/36.D0*T91
     $      +17.D0/108.D0*T97
        DRNLF = -7.D0/72.D0*A*T28*B*T42*T19*D*rho1+7.D0/36.D0*T63
     $     -7.D0/72.D0*A*T5*B*T42*C*T59
        DR2NLF = -7.D0/108.D0*A/T27/T4*B*T118*T19*T119*T11
     $     +7.D0/27.D0*T77-7.D0/108.D0*A*T28*B*T118*C*T73
     $     -7.D0/72.D0*A*T28*B*T42*T19*D*rho2-35.D0/108.D0*T85
     $     +7.D0/36.D0*T91+7.D0/27.D0*T97
     $     -7.D0/72.D0*A*T5*B*T42*C*T87
     $     -7.D0/216.D0*A*T5*B*T118*T150*T81
        ucg = DNF-(2/R*DN1F+DRN1F)+(2*DRNLF+DR2NLF)
c    due to the LYP  e_c
        uc = 0.0
        ec = 0.0
      endif
C GRADIENT CORRECTIONS TO THE CORRELATION POTENTIAL  --  PERDEW (1986)
c    PRB33, 8822(1986)
        if(alphax.eq.3. ) then
        rs2=rs*rs
        ccn=1.d0+8.723d0*rs+0.472d0*rs2+7.389d-02*rs2*rs
        ccn=1.d0/ccn*(0.002568d0+0.023266d0*rs+7.389d-06*rs2)
        ccn=ccn+0.001667d0
        ccinf=0.004235d0
        rapp=RNORM/RHO**(1.166666666d0)
        cfi=0.19195d0*ccinf/ccn*rapp
        ecg=exp(-cfi)*ccn*rapp*rapp
        
        rrho=RHO
        drho=RNORM
        rapp=drho/rrho
        rapp2=rapp*rapp
        rlarho=LAPrho

        dcdrs=((0.023266d0+2.d0*7.389d-06*rs)*(1.d0+8.723d0*rs+
     $    0.472d0*rs*rs+7.389d-02*rs*rs*rs)-
     $   (0.002568d0+rs*(0.023266d0+7.389d-06*rs))*(8.723d0+
     $    2.d0*0.472d0*rs+3.d0*7.389d-02*rs*rs))/
     $   (1.d0+rs*(8.723d0+rs*(0.472d0+7.389d-02*rs)))**2

        dcdn=-dcdrs/3.d0*rs/rrho
        equal=dcdn*exp(-cfi)*rnorm*rnorm/rrho**(4.d0/3.d0)*
     $   (cfi*cfi-cfi-1.d0)
        fa1=4.d0/3.d0-11.d0/3.d0*cfi+7.d0/6.d0*cfi*cfi
        t1=(2.d0-cfi)*rlarho/rrho
        t2=-fa1*(drho/rrho)**2
        auxb=rho2/rrho
        t3=cfi*(cfi-3.d0)*auxb
        ucg=-exp(-cfi)*ccn*(t1+t2+t3)/rrho**(1.d0/3.d0)+equal
c        write(*,*) 'pw86',r,ec,ecg,uc,ucg
      endif
C
C GRADIENT CORRECTIONS TO THE Correlation POTENTIAL  ---  Perdew Wang (1991)
C   PRB 46, 6671 (1992)
      if (alphax .eq. 4.) then
         zet = 0.0
         g=1.d0
c test
c         zet=1.0
c         g=1./2.**third
c test end

         rkf= (2.25*pi)**third/rs
         rks= sqrt(4.*rkf/pi)
         t= rnorm/(rho*2.*rks*g)
         uu=rnorm*rho2/(rho**2*(2.*rks*g)**3)
         vv=laprho/(rho*(2.*rks*g)**2)
         
         ww=0.0
         eczet=.0
c         write(*,*) 'CA',ec,ecrs,eczet,uc
         call CLSDPW91(rs,zet,ec,uc1,uc2,ecrs,eczet,alfc)
         uc=uc1+uc2
         call cggapw91(rs,rkf,zet,g,t,uu,vv,ww,ec,
     +        ecrs,eczet,ecg,ucg1,ucg2)
C conversion from hartree to rydberg         
         ec = 2.*ec
         ucg = .5*(ucg1+ucg2)

c         uc = 2*uc1
c         ucg = ucg1
c         write(*,*) 'PW',ec,ecrs,eczet,uc
c        write(*,*) 'pw92',r,ec,ecg,uc,ucg

      endif
C
C GRADIENT CORRECTIONS TO THE EXCHANGE POTENTIAL  ---  Perdew Wang (1991)
C   PRB 46, 6671 (1992)
      if (alphax .eq. 4.) then
c    H test begin
c          RS= rs /  2. ** third
c        RK2=RK2*4.d0
c        RNORM=SQRT(RK2)
c        rho= 2.*rho
c        laprho=2.*laprho
c        rho2=2.*rho2
c    H test end


        rkf= (2.25*pi)**third/rs
         s=rnorm/(2.*rho*rkf)
         u=rnorm*rho2/(rho**2*(2.*rkf)**3)
         v=laprho/(rho*(2.*rkf)**2)
         call exchpw91(rho,s,u,v,exg,uxg)
c test
c          exg=.5*exg
c          uxg = .5*uxg
c    H test end

         ex=0.0
         ux=0.0
      endif

c
c*********** GGA - PBE approximations of XC-energy and potential 
c (written by Keith Ng 2/9/97 ) ?
      if (alphax .eq. 5.) then
c     do correlation first, unpolarized
         zet = 0.0
         g=1.d0
         rkf= (2.25*pi)**third/rs
         rks= sqrt(4.*rkf/pi)
         t= rnorm/(rho*2.*rks*g)
         uu=rnorm*rho2/(rho**2*(2.*rks*g)**3)
         vv=laprho/(rho*(2.*rks*g)**2)
         ww=0.0
         eczet=.0
c
         lgga=1
         lpot=1
c     flags to do GGA and return V_XC
         call CORPBE(RS,ZET,T,UU,VV,WW,lgga,lpot,ec,uc1,uc2,
     1                  ecg,ucg1,ucg2)
C conversion from hartree to rydberg         
         ec = 2.*ec
         uc=uc1+uc2
         ucg = .5*(ucg1+ucg2)
c     now do exchange, initialize ex & ux, use PW version.
         ex=0.0
         ux=0.0
c  
         rkf= (2.25*pi)**third/rs
         s=rnorm/(2.*rho*rkf)
         u=rnorm*rho2/(rho**2*(2.*rkf)**3)
         v=laprho/(rho*(2.*rkf)**2)
         call EXCHPBE(rho,S,U,V,lgga,lpot,exg,uxg)
         endif



*********************** GRADIENT CORRECTIONS ***********************
C
C GRADIENT CORRECTIONS TO THE EXCHANGE POTENTIAL  ---  A.D.BECKE (1988)
C   PRA 38, 3098 (1988)
       if(alphax.eq.1. .or. alphax.eq.2 .or.alphax.eq.3 ) then
         bet=0.0042d0
        RS= rs *  2. ** third
        RK2=RK2/4.d0
        RNORM=SQRT(RK2)
        RHO=0.5D0*RHO
        LAPrho=0.5D0*LAPrho
        rrho=RHO
        y03=RHO**third
        varx=rnorm/(RHO*y03)
        varx2=varx*varx
        sq=dsqrt(1.d0+varx2)
        shm1=log(varx+sq)
        deno=1.d0+6.d0*bet*varx*shm1
        exg = -bet*y03*varx2/deno

        caf=1.d0/sq
        auxb=rho2/2.d0/(rrho*rrho*rrho)*rnorm
        auxb=6.d0*bet*y03*(auxb-4.d0*varx**3/3.d0)
        auxb1=3.d0*shm1*(1.d0+2.d0*bet*varx*shm1)
        auxb1=(auxb1+4.d0*varx*caf*(1.d0-3.d0*bet*varx*varx*caf))
        auxb1=auxb1/deno+varx*caf*caf*caf
        auxb=auxb*auxb1
        uxg=auxb+4.d0/3.d0*y03*varx2*deno
        uxg=uxg-2.d0*LAPrho/rrho/y03*
     +       (1.d0+3.d0*bet*varx*(shm1-varx*caf))
        uxg=-bet*uxg/(deno*deno)
        RS= rs /  2. ** third
        RHO=RHO+RHO
        LAPrho=LAPrho+LAPrho
      endif

c factor of 2, conversion from hartree to rydberg
      uxc=ux+ uc+2.d0* (uxg+ucg)
      exc=ex+ec+2.d0*(ecg + exg)

c
c      write(21,*) r,rho,rho1,rho2,uxc
* the factor 2 is added because of unity reasons (Ha --> Ryd)
      endif
c
c
      return
      end

      SUBROUTINE CLSDPW91(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
C  UNIFORM-GAS CORRELATION OF PERDEW AND WANG 1991
C  INPUT: SEITZ RADIUS (RS), RELATIVE SPIN POLARIZATION (ZET)
C  OUTPUT: CORRELATION ENERGY PER ELECTRON (EC), UP- AND DOWN-SPIN
C     POTENTIALS (VCUP,VCDN), DERIVATIVES OF EC WRT RS (ECRS) & ZET (ECZET)
C  OUTPUT: CORRELATION CONTRIBUTION (ALFC) TO THE SPIN STIFFNESS
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA GAM,FZZ/0.5198421D0,1.709921D0/
      DATA THRD,THRD4/0.333333333333D0,1.333333333333D0/
      F = ((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
      CALL GCPW91(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,
     1    0.49294D0,1.00D0,RS,EU,EURS)
      CALL GCPW91(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,
     1    0.62517D0,1.00D0,RS,EP,EPRS)
      CALL GCPW91(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0,
     1    0.49671D0,1.00D0,RS,ALFM,ALFRSM)
C  ALFM IS MINUS THE SPIN STIFFNESS ALFC
      ALFC = -ALFM
      Z4 = ZET**4
      EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
C  ENERGY DONE. NOW THE POTENTIAL:
      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
      FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
      ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU
     1        -(1.D0-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3.D0-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
      RETURN
      END
      SUBROUTINE GCPW91(A,A1,B1,B2,B3,B4,P,RS,GG,GGRS)
C  CALLED BY SUBROUTINE CLSDpw91
      IMPLICIT REAL*8 (A-H,O-Z)
      P1 = P + 1.D0
      Q0 = -2.D0*A*(1.D0+A1*RS)
      RS12 = DSQRT(RS)
      RS32 = RS12**3
      RSP = RS**P
      Q1 = 2.D0*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RSP)
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/RS12+2.D0*B2+3.D0*B3*RS12+2.D0*B4*P1*RSP)
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1**2+Q1)
      RETURN
      END

      SUBROUTINE CGGAPW91(RS,FK,ZET,G,T,UU,VV,WW,EC,ECRS,ECZET,
     +     H,DVCUP,DVCDN)
C  GGA91 CORRELATION
C  INPUT RS: SEITZ RADIUS
C  INPUT FK: Fermi wave wector
C  INPUT ZET: RELATIVE SPIN POLARIZATION
C  INPUT G  : ((1+zet)^(2/3)+(1-zet)^(2/3))/2
C  INPUT T: ABS(GRAD D)/(D*2.*KS*G)
C  INPUT UU: (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KS*G)**3)
C  INPUT VV: (LAPLACIAN D)/(D * (2*KS*G)**2)
C  INPUT WW:  (GRAD D)*(GRAD ZET)/(D * (2*KS*G)**2
C  INPUT EC: CORRELATION ENERGY PER ELECTRON 
C  INPUT ECRS: DERIVATIVES OF EC WRT RS
C  INPUT ECZET: DERIVATIVES OF EC WRT ZET
C  OUTPUT H: NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
C  OUTPUT DVCUP,DVCDN:  NONLOCAL PARTS OF CORRELATION POTENTIALS

      IMPLICIT REAL*8 (A-H,O-Z)
      DATA XNU,CC0,CX,ALF/15.75592D0,0.004235D0,-0.001667212D0,0.09D0/
      DATA C1,C2,C3,C4/0.002568D0,0.023266D0,7.389D-6,8.723D0/
      DATA C5,C6,A4/0.472D0,7.389D-2,100.D0/
      DATA THRDM,THRD2/-0.333333333333D0,0.666666666667D0/
      data pi    /3.141592654/
      SK=DSQRT(4.d0*FK/PI)
      BET = XNU*CC0
      DELT = 2.D0*ALF/BET
      G3 = G**3
      G4 = G3*G
      PON = -DELT*EC/(G3*BET)
      B = DELT/(DEXP(PON)-1.D0)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      T6 = T4*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.D0+B*T2
      Q5 = 1.D0+B*T2+B2*T4
      Q6 = C1+C2*RS+C3*RS2
      Q7 = 1.D0+C4*RS+C5*RS2+C6*RS3
      CC = -CX + Q6/Q7
      R0 = (SK/FK)**2
      R1 = A4*R0*G4
      COEFF = CC-CC0-3.D0*CX/7.D0
      R2 = XNU*COEFF*G3
      R3 = DEXP(-R1*T2)
      H0 = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T2/Q5)
      H1 = R3*R2*T2
      H = H0+H1
C  LOCAL CORRELATION OPTION:
C     H = 0.0D0
C  ENERGY DONE. NOW THE POTENTIAL:
      CCRS = (C2+2.*C3*RS)/Q7 - Q6*(C4+2.*C5*RS+3.*C6*RS2)/Q7**2
      RSTHRD = RS/3.D0
      R4 = RSTHRD*CCRS/COEFF
      GZ = ((1.D0+ZET)**THRDM - (1.D0-ZET)**THRDM)/3.D0
      FAC = DELT/B+1.D0
      BG = -3.D0*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1.D0+2.D0*B*T2
      H0B = -BET*G3*B*T6*(2.D0+B*T2)/Q8
      H0RS = -RSTHRD*H0B*BEC*ECRS
      FACT0 = 2.D0*DELT-6.D0*B
      FACT1 = Q5*Q9+Q4*Q9*Q9
      H0BT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
      H0RST = RSTHRD*T2*H0BT*BEC*ECRS
      H0Z = 3.D0*GZ*H0/G + H0B*(BG*GZ+BEC*ECZET)
      H0T = 2.*BET*G3*Q9/Q8
      H0ZT = 3.D0*GZ*H0T/G+H0BT*(BG*GZ+BEC*ECZET)
      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
      H0TT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
      H1RS = R3*R2*T2*(-R4+R1*T2/3.D0)
      FACT4 = 2.D0-R1*T2
      H1RST = R3*R2*T2*(2.D0*R4*(1.D0-R1*T2)-THRD2*R1*T2*FACT4)
      H1Z = GZ*R3*R2*T2*(3.D0-4.D0*R1*T2)/G
      H1T = 2.D0*R3*R2*(1.D0-R1*T2)
      H1ZT = 2.D0*GZ*R3*R2*(3.D0-11.D0*R1*T2+4.D0*R1*R1*T4)/G
      H1TT = 4.D0*R3*R2*R1*T*(-2.D0+R1*T2)
      HRS = H0RS+H1RS
      HRST = H0RST+H1RST
      HT = H0T+H1T
      HTT = H0TT+H1TT
      HZ = H0Z+H1Z
      HZT = H0ZT+H1ZT
      COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
      PREF = HZ-GZ*T2*HT/G
      FACT5 = GZ*(2.D0*HT+T*HTT)/G
      COMM = COMM-PREF*ZET-UU*HTT-VV*HT-WW*(HZT-FACT5)
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
C  LOCAL CORRELATION OPTION:
C     DVCUP = 0.0D0
C     DVCDN = 0.0D0
      RETURN
      END

      SUBROUTINE EXCHPW91(D,S,U,V,EX,VX)
C  GGA91 EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
C  INPUT D : DENSITY
C  INPUT S:  ABS(GRAD D)/(2*KF*D)
C  INPUT U:  (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KF)**3)
C  INPUT V: (LAPLACIAN D)/(D*(2*KF)**2)
C  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA A1,A2,A3,A4/0.19645D0,0.27430D0,0.15084D0,100.D0/
      DATA AX,A,B1/-0.7385588D0,7.7956D0,0.004D0/
      DATA THRD,THRD4/0.333333333333D0,1.33333333333D0/
      FAC = AX*D**THRD
      S2 = S*S
      S3 = S2*S
      S4 = S2*S2
      P0 = 1.D0/DSQRT(1.D0+A*A*S2)
      P1 = DLOG(A*S+1.D0/P0)
      P2 = DEXP(-A4*S2)
      P3 = 1.D0/(1.D0+A1*S*P1+B1*S4)
      P4 = 1.D0+A1*S*P1+(A2-A3*P2)*S2
      F = P3*P4
      EX = FAC*F
C  LOCAL EXCHANGE OPTION
C     EX = FAC
C  ENERGY DONE. NOW THE POTENTIAL:
      P5 = B1*S2-(A2-A3*P2)
      P6 = A1*S*(P1+A*S*P0)
      P7 = 2.D0*(A2-A3*P2)+2.D0*A3*A4*S2*P2-4.D0*B1*S2*F
      FS = P3*(P3*P5*P6+P7)
      P8 = 2.D0*S*(B1-A3*A4*P2)
      P9 = A1*P1+A*A1*S*P0*(3.D0-A*A*S2*P0*P0)
      P10 = 4.D0*A3*A4*S*P2*(2.D0-A4*S2)-8.D0*B1*S*F-4.D0*B1*S3*FS
      P11 = -P3*P3*(A1*P1+A*A1*S*P0+4.D0*B1*S3)
      FSS = P3*P3*(P5*P9+P6*P8)+2.D0*P3*P5*P6*P11+P3*P10+P7*P11
      VX = FAC*(THRD4*F-(U-THRD4*S3)*FSS-V*FS)
C  LOCAL EXCHANGE OPTION:
C     VX = FAC*THRD4
      RETURN
      END

c
c----------------------------------------------------------------------------
c
      subroutine radin(mesh,func,rab,asum,idim1)
c
c     simpson's rule integrator for function stored on the
c     logarithmic mesh
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c
c.....logarithmic radial mesh information
      dimension rab(idim1)
c.....function to be integrated
      dimension func(idim1)
c
c.....variable for file = 0
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c     routine assumes that mesh is an odd number so run check
      if ( mesh - ( mesh / 2 ) * 2 .ne. 1 ) then
        write(iout,*) '***error in subroutine radin'
        write(iout,*) 'routine assumes mesh is odd but mesh =',mesh
        call exit(1)
      endif
c
c     should really be no case where func(1) .ne. 0
      if ( abs(func(1)) .gt. 1.0d-10 ) then
        write(iout,*) '***error in subroutine radin'
        write(iout,*) 'func(1) is non-zero, func(1) =',func(1)
        call exit(1)
      endif
c
      asum = 0.0d0
      r12 = 1.0d0 / 12.0d0
      f3  = func(1) * rab(1) * r12
      func(1) = 0.0d0
c
      do 100 i = 2,mesh-1,2
        f1 = f3
        f2 = func(i) * rab(i) * r12
        f3 = func(i+1) * rab(i+1) * r12
        asum = asum + 5.0d0*f1 + 8.0d0*f2 - 1.0d0*f3
        func(i) = asum
        asum = asum - 1.0d0*f1 + 8.0d0*f2 + 5.0d0*f3
        func(i+1) = asum
100   continue
c
      return
      end
c
c-------------------------------------------------------------------------
c
      subroutine rsae(snl,wwnl,nkkk,nnlz,rscore,rsvale,ncores,ncspvs,
     +  mesh,alpha,nbl,nbl0,rab,ipass,keyps,kkbeta,qfunc,dum,
     +  ifprt,idim1,idim2,idim3,idim5)
c
c     compute the charge density for the all electron calculation
c     and compute valence corrections for generalized eigenvalue case
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c.....logarithmic radial mesh information
      dimension rab(idim1)
c.....charge densities
      dimension rscore(idim1),rsvale(idim1)
c.....wavefunctions
      dimension snl(idim1,idim2)
c.....shell labels, energies, occupancies and real-space cutoff index
      dimension nnlz(idim2),wwnl(idim2),nkkk(idim2)
c.....q function
      dimension qfunc(idim1,idim3,idim3)
c.....work space for solution to inhomogeneous equation and overlaps
      dimension alpha(idim3,idim2)
c.....indexing for the beta functions
      dimension nbl0(idim5),nbl(idim5)
c.....scratch
      dimension dum(idim1)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c--------------------------------------------------------------------------
c
c     i n i t i a l i s a t i o n
c
      if ( ifprt .ge. 1) write (iout,'(/a)')
     +      ' subroutine rsae: construct charge density'
c
c     zero the core and valence charge density arrays
      do 10 i = 1,mesh
        rscore(i) = 0.0d0
        rsvale(i) = 0.0d0
   10 continue
c
c--------------------------------------------------------------------------
c
c     c o r e  c h a r g e  d e n s i t y  (  o n l y  i f  a e  )
c
      if ( ipass .eq. 1 ) then
c
        do 100 j = 1,ncores
        do 100 i = 1,nkkk(j)
          rscore(i) = rscore(i) + wwnl(j) * snl(i,j)**2
  100   continue
c
      endif
c
c--------------------------------------------------------------------------
c
c            v a l e n c e  c h a r g e  d e n s i t y
c
      do 110 j = ncores+1,ncspvs
c
        do 120 i = 1,nkkk(j)
          rsvale(i) = rsvale(i) + wwnl(j) * snl(i,j)**2
  120   continue
c
c       possible extra contribution from vanderbilt scheme
        if ( ipass .eq. 2 .and. keyps .eq. 3 ) then
c
          ncur = nnlz(j) / 100
          lcur = nnlz(j) / 10 - 10 * ncur
          lp   = lcur + 1
c
          do 210 i = 1,mesh
            dum(i) = rsvale(i)
  210     continue
          call radin(mesh,dum,rab,asum,idim1)
          write (iout,'(3x,a,i2,f12.8)')
     +      'intermed rsvale for lp=',lp,asum
          write (iout,'(5x,a,2f12.8)')
     +      'alpha=',(alpha(nbl0(lp)+ib,j),ib=1,nbl(lp))
c
          do 220 ib = 1,nbl(lp)
          do 220 jb = 1,nbl(lp)
            ibeta = nbl0(lp) + ib
            jbeta = nbl0(lp) + jb
            fac = wwnl(j) * alpha(ibeta,j) * alpha(jbeta,j)
            do 230 i = 1,kkbeta
              rsvale(i) = rsvale(i) + fac * qfunc(i,ibeta,jbeta)
  230       continue
  220     continue
c
          do 240 i = 1,mesh
            dum(i) = rsvale(i)
  240     continue
          call radin(mesh,dum,rab,asum,idim1)
          write (iout,'(5x,a,i2,f12.8)')
     +      'integrated rsvale after lp',lp,asum
c
        endif
c
c     close loop over valence orbitals
  110 continue
c
c-----------------------------------------------------------------------
c
c         b a l a n c e  v a l e n c e  c h a r g e
c
      if ( ipass .eq. 2 .and. keyps .eq. 3 ) then
c
        wwsum = 0.0d0
        do 300 j = ncores+1,ncspvs
          wwsum = wwsum + wwnl(j)
  300   continue
c
        do 310 i = 1,mesh
          dum(i) = rsvale(i)
  310   continue
c
        call radin(mesh,dum,rab,asum,idim1)
c
        write (iout,'(3x,a,f12.8,a,f12.8)')
     +      'balancing valence density, asum =',asum,' wwsum =',wwsum
c
        if ( abs( asum - wwsum ) .gt. 8.0d-05 ) then
          write(iout,*) '***error in subroutine rsae'
          write(iout,*) 'asum and wwsum are not consistent'
          call exit(1)
        endif
c
        if ( abs(asum) .lt. 1.0d-05 ) return
c
        fac = wwsum / asum
        do 320 i = 1,mesh
          rsvale(i) = fac * rsvale(i)
  320   continue
c
      endif
c
      if ( ifprt .ge. 1) write (iout,'(a)') ' leaving rsae'
c
      return
      end
c
c------------------------------------------------------------------------
c
      subroutine mixer(ruae,runuc,ruhar,ruexch,rudif,damp,
     +  tol,delta,iselfc,mesh,idim1)
c
c     routine for mixing old and new charge potentials.
c     if the potentials are found to be self-consistent iselfc = +1
c                                             otherwise iselfc = -1
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c.....potentials for the all electron calculation
      dimension ruae(idim1),runuc(idim1),rudif(idim1)
c.....hartree, exchange and correlation potentials
      dimension ruhar(idim1),ruexch(idim1)
c
c     initialisation
      iselfc = - 1
      delta = 0.0d0
c
c     compute the difference between new and old all electron potentials
      do 100 i = 2,mesh
        rudif(i) = runuc(i) + ruhar(i) + ruexch(i) - ruae(i)
        if ( delta .lt. abs( rudif(i) ) ) delta = abs( rudif(i) )
  100 continue
c
c     check whether the potential is self consistent
      if ( delta .lt. tol ) then
        iselfc = + 1
c       on last iteration, update potential without damping
        do 190 i = 2,mesh
          ruae(i) = ruae(i) + rudif(i)
  190   continue
      else
        do 200 i = 2,mesh
          ruae(i) = ruae(i) + damp * rudif(i)
  200   continue
      endif
c
      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine orgsnl(r,snlo,j1,j2,lcur,zz,vzero,ecur,idim1)
c
c     generate an analytic starting guess for the wavefunction
c     at the origin using taylor series expansion.
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c
c.....logarithmic radial mesh information
      dimension r(idim1)
c.....wavefunctions
      dimension snlo(idim1)
c
c
c     compute the starting value of psi according to
c
c     psi(r) = r**(l+1) * ( 1 + aa * r + bb * r**2 + ... )
c
c     where:
c           aa = - z / (l+1) 
c           bb = ( 2 * z * z / (l+1) + v(0) - e ) / (4*l+6)
c
c
      lp = lcur + 1
      aa = - zz / dble(lp)
      bb = ( 2 * zz * zz / dble(lp) + vzero - ecur ) / dble(4*lcur+6)
c
      do 100 i = j1,j2
        snlo(i) = (r(i)**lp) * ( 1.0d0 + aa*r(i) + bb*r(i)*r(i) )
100   continue
c
      return
      end
c
c----------------------------------------------------------------------------
c
c     ***************************************************************
      subroutine lderiv(ifprt,ipass,keyps,ilogd,klogd,irel,emax,emin,
     +  nnt,zz,vzero,b,r,rab,sqr,ruae,v,snlo,yy,flname,mesh,
     +  vloc,beta,nbl,nbl0,gamma,alpha,eta,kkbeta,ddd,qqq,
     +  bndmat,idim1,idim3,idim5,idim7)
c     ***************************************************************
c     compute logarithmic derivatives over an energy range and
c       print comparisons between ae and pseudo case
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
      parameter ( idim10 = 81 , idim11 = 10 )
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1),sqr(idim1)
c.....potentials for the all electron calculation
      dimension ruae(idim1),v(idim1)
c.....wavefunctions
      dimension snlo(idim1),yy(idim1,2)
c.....beta and q functions and q pseudization coefficients
      dimension beta(idim1,idim3)
c.....angular momenta, reference energy, qij and dij of vanderbilt scheme
      dimension qqq(idim3,idim3),ddd(idim3,idim3)
c.....work space for solution to inhomogeneous equation and overlaps
      dimension eta(idim1,idim7),alpha(idim3)
c.....indexing for the beta functions
      dimension nbl0(idim5),nbl(idim5)
c.....local component of the pseudopotential
      dimension vloc(idim1)
c.....scratch arrays
      dimension gamma(idim1),bndmat(idim1,2,2)
c.....logarithmic derivative arrays
      dimension dlsav(idim10,idim11),ltab(idim11)
c.....symbols for the line printer plot out
      dimension isymbl(9)
c.....file name array
      character*40 flname(6)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c     initialise the symbol array
      data  isymbl / 1h!,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8 /
c
      save isymbl, dlsav, ltab
c
c----------------------------------------------------------------------------
c
c                    i n i t i a l i s a t i o n
c
      write (iout,'(/a,i3)') ' subroutine lderiv: pass ',ipass
c
      nnt1 = nnt + 1
c
      if ( ipass .lt. 1 .or. ipass .gt. 2 ) then
        write(iout,*) '***error in subroutine lderiv'
        write(iout,*) 'ipass =',ipass,' is out of range'
        call exit(1)
      endif
c
      if ( ipass .eq. 1 ) then
c
c       check that idim10 is sufficient
        if ( nnt1 .gt. idim10 ) then
          write(iout,*) '***error in subroutine lderiv'
          write(iout,*) 'nnt+1=',nnt1,' but idim10=',idim10
          call exit(1)
        endif
c
c       check that idim11 is sufficient
        if ( 2 * ilogd + 2 .gt. idim11 ) then
          write(iout,*) '***error in subroutine lderiv'
          write(iout,*) '2*ilogd+2=',2*ilogd+2,' but idim11',idim11
          call exit(1)
        endif
c
c       open the file for storing the logarithmic derivative information
        open( unit = iologd , file = flname(5) , status = 'unknown' ,
     +    form = 'unformatted' )
c
c       zeroing
        do 20 j = 1,2*ilogd+2
        do 20 i = 1,nnt1
          dlsav(i,j) = 0.0d0
   20   continue
c
      endif
c
c---------------------------------------------------------------------------
c
c          c o m p u t e  t h e  l o g  d e r i v a t i v e s
c
      dele = ( emax - emin ) / dble(nnt)
c
c     loop over the angular momenta
      do 100 lcur = 0,ilogd-1
c
        index = 3 + ilogd * ( ipass - 1 ) + lcur
        ltab(index) = lcur
c
        do 110 i = 1,nnt1
c
          ecur = emin + dele * dble(i-1)
c
          if ( ipass .eq. 1 ) then
c
            if ( irel .eq. 0 ) then
c
c             schroedinger equation is required
              call phase(zz,vzero,ecur,lcur,b,r,rab,sqr,
     +          ruae,v,snlo,mesh,klogd,dlwf,idim1)
c
            elseif ( irel .eq. 2 ) then
c
c             koelling and harmon scalar relativistic equation
              call phakoe(ruae,r,rab,bndmat,ecur,lcur,
     +          yy,v,klogd,dlwf,idim1)
c
            else
c
              write(iout,*) '***error in subroutine lderiv'
              write(iout,*) 'irel =',irel,' is not programmed'
              call exit(1)
c
            endif
c
c
          elseif ( ipass .eq. 2 .and. keyps .eq. 3 ) then
c
            call schgef(ecur,lcur,b,r,rab,sqr,vloc,v,
     +        snlo,mesh,klogd,dlwf,beta,nbl,nbl0,gamma,alpha,
     +        eta,kkbeta,ddd,qqq,idim1,idim3,idim5,idim7)
c
          else
            write(iout,*) '***error in subroutine lderiv'
            write(iout,*) 'ipass=',ipass,'and keyps=',keyps,'illegal'
            call exit(1)
          endif
c
          dlsav(i,index) = dlwf
c
  110   continue
c
        write (iologd) ilogd,nnt1,emin,dele,(dlsav(i,index),i=1,nnt1)
c
  100 continue
c
      if ( ipass .eq. 1 .or. ifprt .lt. 1) then
        write (iout,'(/1x,a)') 'leaving lderiv'
        return
      endif
c
c---------------------------------------------------------------------------
c
c     p r i n t  o u t  o f  l o g  d e r i v a t i v e s
c
      ndl = 2 + ilogd * ipass
      write (iout,'(/5x,a)') 'logarithmic derivatives'
      if (ndl-2.le.2) then
        write (iout,202) (ltab(m),m=3,ndl)
      else if (ndl-2.le.4) then
        write (iout,204) (ltab(m),m=3,ndl)
      else if (ndl-2.le.6) then
        write (iout,206) (ltab(m),m=3,ndl)
      else
        write (iout,208) (ltab(m),m=3,ndl)
      endif
      write (iout,*)
      do i=1,nnt1
        e=emin+dele*dble(i-1)
        dlsav(i,1)=e
        dlsav(i,2)=0.
        if (ndl-2.le.6) then
          write (iout,230) e,(dlsav(i,m),m=3,ndl)
        else
          write (iout,238) e,(dlsav(i,m),m=3,ndl)
        endif
      end do
c
  202 format (/6x,'e',2x,2(8x,2hl=,i1,5x))
  204 format (/6x,'e',2x,4(8x,2hl=,i1,5x))
  206 format (/6x,'e',2x,6(8x,2hl=,i1,5x))
  208 format (/6x,'e',2x,8(6x,2hl=,i1,4x))
  230 format (3x,f6.3,6(3x,e13.5))
  238 format (3x,f6.3,8(1x,e12.4))
c
      write (iout,'(/3x,a)') 'chart: log derivatives'
c
c     line printer plot out of results
      call lpplot(8,dlsav,nnt1,ndl,nnt1,0,idim10,idim11,isymbl,-4.,4.)
c
      write (iout,'(/1x,a)') 'leaving lderiv'
c
      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine phase(zz,vzero,ecur,lcur,b,r,rab,sqr,
     +  ruae,v,snlo,mesh,klogd,dlwf,idim1)
c
c     routine for computing logarithmic derivatives
c
c     output in dlwf is d ln ( r*psi(r) ) / dr
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1),sqr(idim1)
c.....potentials for the all electron calculation
      dimension ruae(idim1)
c.....wavefunctions
      dimension snlo(idim1)
c.....scratch arrays
      dimension v(idim1)
c
c     compute the potential for this angular momentum and energy
      llp = lcur * ( lcur + 1 )
      call setv(r,rab,ruae,b,ecur,llp,v,mesh,idim1)
c
c     analytic starting value for snlo using taylor expansion
      call orgsnl(r,snlo,1,3,lcur,zz,vzero,ecur,idim1)
c
c     convert to the phi function
      snlo(1) = snlo(1) / sqr(1)
      snlo(2) = snlo(2) / sqr(2)
      snlo(3) = snlo(3) / sqr(3)
c
c     integrate the schroedinger equation outwards
      call difsol(v,snlo,4,klogd,idim1)
c
c     compute the logarithmic derivative
c
      snlom = snlo(klogd-1) * sqr(klogd-1)
      snlop = snlo(klogd)   * sqr(klogd)
c
      dwf =   snlop - snlom
      wf  = ( snlop + snlom ) * 0.5d0
      dr  =   r(klogd) -    r(klogd-1)
c
      dlwf = dwf / ( dr * wf )
c
      return
      end
c
c----------------------------------------------------------------------------
c
c     *******************************************
      subroutine nodeno(snlo,jj1,jj2,nodes,idim1)
c     *******************************************
c
c     routine counts the number of nodes of the wavefunction snlo
c     between the points jj1 and jj2
c
      implicit double precision (a-h,o-z)
c
c.....wavefunction array
      dimension snlo(idim1)
c
      nodes = 0
c
      do 100 i = jj1+1,jj2
        if ( snlo(i-1) * snlo(i) .lt. 0.0d0 ) nodes = nodes + 1
  100 continue
c
      return
      end
c
c---------------------------------------------------------------------------
