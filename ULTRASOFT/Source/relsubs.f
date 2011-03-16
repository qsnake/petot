c
c Copyright (C) 2002 Vanderbilt group
c This file is distributed under the terms of the GNU General Public
c License as described in the file 'License' in the current directory.
c
c-------------------------------------------------------------------------
c
      subroutine dirsol(ruae,r,rab,zz,ee,nnlz,wwnl,nkkk,krel,
     +  snl,yy,ncspvs,mesh,thresh,ifprt,idim1,idim2)
c
c     subroutine to compute solutions to the full dirac equation
c
c     the dirac equation in rydberg units reads:
c
c     df(r)     k                     alpha
c     ----- = + - f(r) - ( e - v(r) ) ----- g(r)
c      dr       r                       2
c
c     dg(r)     k                          4      alpha
c     ----- = - - g(r) + ( e - v(r) +  -------- ) ----- f(r)
c      dr       r                      alpha**2     2
c
c     where 
c            alpha is the fine structure constant
c            f(r) is r*minor component
c            g(r) is r*major component
c            k is quantum number 
c               k = - (l+1)    if  j = l+0.5
c               k = + l        if  j = l-0.5
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1)
c.....potentials for the all electron calculation
      dimension ruae(idim1),zz(idim1,2,2)
c.....wavefunctions
      dimension snl(idim1,idim2),yy(idim1,2)
c.....shell labels, energies, occupancies and real space cutoff, 
c     relativistic quantum numbers
      dimension nnlz(idim2),ee(idim2),wwnl(idim2),nkkk(idim2),
     +  krel(idim2)
c     logic to determine whether intialization of quantum nos needed
      logical linit
c
c.....variable for file = 0
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
      save linit
c
      data linit /.false./
c
c--------------------------------------------------------------------------
c
c               r o u t i n e  i n i t i a l i s a t i o n
c
c     set the maximum number of iterations for improving wavefunctions
      itmax = 50
c
c     set ( 2 / fine structure constant )
      tbya = 2.0d0 * 137.04d0
c     set ( fine structure constant / 2 )
      abyt = 1.0d0 / tbya
c
c
      if ( .not. linit ) then
c
c       split the orbitals into l-0.5 and l+0.5 pairs
c
c       first of all determine the total number of orbitals
        ip = 0
        do 10 iorb = 1,ncspvs
c
          ncur = nnlz(iorb) / 100
          lcur = nnlz(iorb) / 10  - 10 * ncur
c
          if ( lcur .eq. 0 ) then
            ip = ip + 1
          else
            ip = ip + 2
          endif
c
   10   continue
c
c       check that there will be sufficient space to hold expanded 
c       orbital lists
        if ( ip .gt. idim2 ) then
          write(iout,*) '***error in subroutine dirsol'
          write(iout,*) 'idim2 is too small for the dirac case'
          call exit(1)
        endif
c
c       update ncspvs, nnlz, wwnl, ee and initialize krel
        jp = ip + 1
        do 20 iorb = ncspvs,1,-1
c
          ncur = nnlz(iorb) / 100
          lcur = nnlz(iorb) / 10 - 10 * ncur
c
c         deal with case jcur = lcur + 0.5 first
          jp = jp - 1
          nnlz(jp) = nnlz(iorb)
          ee(jp)   = ee(iorb)
          wwnl(jp) = wwnl(iorb) * dble(lcur+1) / dble(2*lcur+1)
          krel(jp) = - ( lcur + 1 )
c
c         case jcur = lcur - 0.5 (only if lcur .ne. 0)
          if ( lcur .ne. 0 ) then
            jp = jp - 1
            nnlz(jp) = nnlz(iorb)
            ee(jp)   = ee(iorb)
            wwnl(jp) = wwnl(iorb) * dble(lcur) / dble(2*lcur+1)
            krel(jp) = lcur
          endif
c
   20   continue
c
        ncspvs = ip
c
        linit = .true.
c
      endif
c
c     zero the wavefunction arrays
      do 30 iorb = 1,ncspvs
      do 30 ir   = 1,mesh
        snl(ir,iorb) = 0.0d0
   30 continue
c
c-----------------------------------------------------------------------
c
c               m a i n  l o o p  o v e r  o r b i t a l s
c
      do 100 iorb = 1,ncspvs
c
c       check that this state is required
        if ( wwnl(iorb) .lt. 0.0d0 ) then
          ee(iorb) = 0.0d0
          nkkk(iorb) = 1
          snl(1,iorb) = 0.0d0
          goto 100
        endif
c
c       set the current quantum numbers
        kcur = krel(iorb)
        ecur = ee(iorb)
        ncur = nnlz(iorb) / 100
        lcur = nnlz(iorb) / 10 - 10 * ncur
c
c       set initial upper and lower bounds for the eigen value
        emin = - 1.0d10
        emax = 1.0d0
c
c       try uncommenting the "write(43..." statements if you
c       encounter convergence problems in the following loop
c
c       write(43,'(/a,i5/)') 'nnlz',nnlz(iorb)
c
c       ======================================================
c       i t e r a t i v e  s o l u t i o n  f o r  e n e r g y
c       ======================================================
c
        do 200 iter = 1,itmax
c
c         ===================
c         define the zz array
c         ===================
c
          if ( iter .eq. 1 ) then
            do 210 ir = 2,mesh
              zz(ir,1,1) = rab(ir) * dble(kcur) / r(ir)
              zz(ir,2,2) = - zz(ir,1,1)
  210       continue
          endif
c
          do 220 ir = 2,mesh
            zz(ir,1,2) = - rab(ir) * ( ecur - ruae(ir) / r(ir) ) * abyt
            zz(ir,2,1) = - zz(ir,1,2) + rab(ir) * tbya
  220     continue
c
c         ==============================================
c         classical turning point and practical infinity
c         ==============================================
c
          do 230 nctp = mesh,10,-1
            if ( zz(nctp,1,2) .lt. 0.0d0 ) goto 240
            if ( ir .eq. 10 ) then
              write(iout,*) '***error in dirsol'
              write(iout,*) 'no classical turning point found'
              call exit(1)
            endif
  230     continue
c
c         jump point out of classical turning point loop
  240     continue
c
          if ( nctp .gt. mesh - 10 ) then
            write(iout,*) '***error in dirsol'
            write(iout,*) 'classical turning point too close to mesh'
            call exit(1)
          endif
c
          tolinf = dlog(thresh) ** 2
          do 250 ninf = nctp+10,mesh
            alpha2 = (ruae(ninf)/r(ninf)-ecur) * (r(ninf) - r(nctp))**2
            if ( alpha2 .gt. tolinf ) goto 260
  250     continue
c
c         jump point out of practical infinity loop
  260     continue
c
          if ( ninf .ge. mesh ) then
            write(iout,*) '***error in dirsol'
            write(iout,*) 'ninf=',ninf,' too big for this mesh',mesh
            call exit(1)
          endif
c
c         ===========================================================
c         analytic start up of minor and major components from origin
c         ===========================================================
c
c         with finite nucleus so potential constant at origin we have
c
c         f(r) = sum_n f_n r ** ( ig + n )
c         g(r) = sum_n g_n r ** ( ig + n )
c
c         with
c
c         f_n+1 = - (ecur-v(0)) * abyt * g_n / ( ig - kcur + n + 1 )
c         g_n+1 = (ecur-v(0)+tbya**2 ) * abyt * f_n / ( ig + kcur + n + 1)
c
c         if kcur > 0  ig = + kcur , f_0 = 1 , g_0 = 0
c         if kcur < 0  ig = - kcur , f_0 = 0 , g_1 = 1
c
          vzero = ruae(2) / r(2)
c
c         set f0 and g0
          if ( kcur .lt. 0 ) then
            ig = - kcur
            f0 = 0
            g0 = 1
          else
            ig = kcur
            f0 = 1
            g0 = 0
          endif
c
          f1 = - (ecur-vzero) * abyt * g0 / dble( ig - kcur + 1 )
          g1 = (ecur-vzero+tbya**2) * abyt * f0 / dble( ig + kcur + 1 )
          f2 = - (ecur-vzero) * abyt * g1 / dble( ig - kcur + 2 )
          g2 = (ecur-vzero+tbya**2) * abyt * f1 / dble( ig + kcur + 2 )
c
          yy(1,1) = 0.0d0
          yy(1,2) = 0.0d0
c
          do 300 ir = 2,6
            yy(ir,1) = r(ir)**ig * ( f0 + r(ir) * ( f1 + r(ir) * f2 ) )
            yy(ir,2) = r(ir)**ig * ( g0 + r(ir) * ( g1 + r(ir) * g2 ) )
  300     continue

c         ===========================
c         outward integration to nctp
c         ===========================
c
c         fifth order predictor corrector integration routine
          call cfdsol(zz,yy,7,nctp,idim1)
c
c         save major component and its gradient at nctp
          gout = yy(nctp,2)
          gpout = zz(nctp,2,1)*yy(nctp,1) + zz(nctp,2,2)*yy(nctp,2)
          gpout = gpout / rab(nctp)
c
c         ==============================================
c         start up of wavefunction at practical infinity
c         ==============================================
c
          do 310 ir = ninf,ninf-4,-1
            alpha = sqrt( ruae(ir) / r(ir) - ecur )
            yy(ir,2) = exp ( - alpha * ( r(ir) - r(nctp) ) )
            yy(ir,1) = ( dble(kcur)/r(ir) - alpha ) * yy(ir,2)*tbya /
     +        ( ecur - ruae(ir)/r(ir) + tbya ** 2 )
  310     continue
c
c         ==========================
c         inward integration to nctp
c         ==========================
c
c         fifth order predictor corrector integration routine
          call cfdsol(zz,yy,ninf-5,nctp,idim1)
c
c         save major component and its gradient at nctp
          gin = yy(nctp,2)
          gpin = zz(nctp,2,1)*yy(nctp,1) + zz(nctp,2,2)*yy(nctp,2)
          gpin = gpin / rab(nctp)
c
c
c         ===============================================
c         rescale tail to make major component continuous
c         ===============================================
c
          factor = gout / gin
          do 320 ir = nctp,ninf
            yy(ir,1) = factor * yy(ir,1)
            yy(ir,2) = factor * yy(ir,2)
  320     continue
c
          gpin = gpin * factor
c
c
c         =================================
c         check that the number of nodes ok
c         =================================
c
c         count the number of nodes in major component
          call nodeno(yy(1,2),2,ninf,nodes,idim1)
c
          if ( nodes .lt. ncur - lcur - 1 ) then
c           energy is too low
            emin = ecur
c           write(43,'(i5,3f12.5,2i5)')
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
c           write(43,'(i5,3f12.5,2i5)')
c    +         iter,emin,ecur,emax,nodes,ncur-lcur-1
            if ( ecur * 1.1d0 .lt. emin ) then
              ecur = 0.5d0 * ecur + 0.5d0 * emin
            else
              ecur = 1.1d0 * ecur
            endif
            goto 200
          endif
c
c
c         =======================================================
c         find normalisation of wavefunction using simpson's rule
c         =======================================================
c
          factor = 0.0d0
          xw = 4.0d0
          do 360 ir = 2,ninf
            factor = factor + xw * rab(ir) * (yy(ir,1)**2 + yy(ir,2)**2)
            xw = 6.0d0 - xw
  360     continue
          factor = factor / 3.0d0
c
c
c         =========================================
c         variational improvement of the eigenvalue
c         =========================================
c
          decur = gout * ( gpout - gpin ) / factor
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
c         write(43,'(i5,3f12.5,1p2e12.4)')
c    +         iter,emin,ecur,emax,decur,decurp
c
c         test to see whether eigenvalue converged
          if ( abs(decur) .lt. thresh ) goto 400
c
          ecur = ecur + decurp
c
c
c         =======================================================
c         check that the iterative loop is not about to terminate
c         =======================================================
c
          if ( iter .eq. itmax ) then
c           eigenfunction has not converged in allowed number of iterations
            write(iout,*) '***error in subroutine dirsol'
            write(iout,999) nnlz(iorb),ee(iorb),ecur
  999       format(' state',i4,' could not be converged.',/,
     +      ' starting energy for calculation was',f10.5,
     +      ' and end value =',f10.5)
            call exit(1)
          endif
c
c       close iterative loop
  200   continue
c
c       jump point on successful convergence of eigenvalue
  400   continue
c
c       update the energy
        ee(iorb) = ecur
c       update the real-space cutoff array
        nkkk(iorb) = ninf
c
c       update the wavefunction array
        factor = 1.0d0 / sqrt( factor )
        do 410 ir = 1,ninf
          snl(ir,iorb) = sqrt( yy(ir,1)**2 + yy(ir,2)**2 ) * factor
  410   continue
c
c     close loop over bands
  100 continue
c
c     stop
      return
      end
c
c-------------------------------------------------------------------------
c
      subroutine koesol(ruae,r,rab,zz,ee,nnlz,wwnl,nkkk,
     +  snl,yy,vme,ncspvs,mesh,thresh,ifprt,idim1,idim2)
c
c     subroutine to compute solutions of koelling and harmon's
c     scalar relativistic equation (j. phys. c 10, 3107 (1977))
c
c     in rydberg units the koelling harmon equations are:
c
c     df(r)     f(r)    l(l+1)
c     ----- = - ---- + -------- g(r) + ( v(r) - e ) g(r)
c      dr        r     m(r)*r*r
c
c     dg(r)                 g(r)
c     ----- =  m(r) f(r) +  ---
c      dr                    r
c
c     m(r) = 1 - alpha**2 * ( v(r) - e ) / 4
c
c     where g(r) is r*(spin-orbit averaged major component) and f(r)
c     is an analogue of r*(minor component) in the full dirac
c     equation.  the charge density is evaluated by taking the
c     square of the spin-orbit averaged major component.
c
c     to transform the above equations into form given in paper
c     must make substitutions phi(r) = f(r) / r , g(r) = g(r) / r
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1)
c.....potentials for the all electron calculation
      dimension ruae(idim1),zz(idim1,2,2)
c.....wavefunctions
      dimension snl(idim1,idim2),yy(idim1,2)
c.....shell labels, energies, occupancies and real space cutoff
      dimension nnlz(idim2),ee(idim2),wwnl(idim2),nkkk(idim2)
c.....scratch
      dimension vme(idim1)
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
c     set ( 2 / fine structure constant )
      tbya = 2.0d0 * 137.04d0
c
c     *******test*****test******test******test****test
c     tbya = 1.0d6
c
c     set ( fine structure constant / 2 )
      abyt = 1.0d0 / tbya
c
c     zero the wavefunction arrays
      do 30 iorb = 1,ncspvs
      do 30 ir   = 1,mesh
        snl(ir,iorb) = 0.0d0
   30 continue
c
c     diagonal elements of zz are constant
      do 40 ir = 2,mesh
        zz(ir,1,1) = - rab(ir) / r(ir)
        zz(ir,2,2) = + rab(ir) / r(ir)
   40 continue
c
c-----------------------------------------------------------------------
c
c               m a i n  l o o p  o v e r  o r b i t a l s
c
      do 100 iorb = 1,ncspvs
c
c       check that this state is required
        if ( wwnl(iorb) .lt. 0.0d0 ) then
          ee(iorb) = 0.0d0
          nkkk(iorb) = 1
          snl(1,iorb) = 0.0d0
          goto 100
        endif
c
c       set the current quantum numbers
        ecur = ee(iorb)
        ncur = nnlz(iorb) / 100
        lcur = nnlz(iorb) / 10 - 10 * ncur
        llp  = lcur * ( lcur + 1 )
c
c       set initial upper and lower bounds for the eigen value
        emin = - 1.0d10
        emax = 1.0d0
c
c       try uncommenting the "write(44..." statements if you
c       encounter convergence problems in the following loop
c
c       write(44,'(/a,i5/)') 'nnlz',nnlz(iorb)
c
c       ======================================================
c       i t e r a t i v e  s o l u t i o n  f o r  e n e r g y
c       ======================================================
c
        do 200 iter = 1,itmax
c
c         ===================
c         define the zz array
c         ===================
c
          do 220 ir = 2,mesh
c
            vme(ir) = ( ruae(ir) / r(ir) - ecur )
            xm  = ( 1.0d0 - vme(ir) * abyt**2 ) 
c
            zz(ir,1,2) = rab(ir) * 
     +        ( dble(llp) / ( xm * r(ir)**2 ) + vme(ir) )
            zz(ir,2,1) = rab(ir) * xm
c
  220     continue
c
c         ==============================================
c         classical turning point and practical infinity
c         ==============================================
c
          do 230 nctp = mesh,10,-1
            if ( vme(nctp) .lt. 0.0d0 ) goto 240
            if ( ir .eq. 10 ) then
              write(iout,*) '***error in koesol'
              write(iout,*) 'no classical turning point found'
              call exit(1)
            endif
  230     continue
c
c         jump point out of classical turning point loop
  240     continue
c
          if ( nctp .gt. mesh - 10 ) then
            write(iout,*) '***error in koesol'
            write(iout,*) 'classical turning point too close to mesh'
            call exit(1)
          endif
c
          tolinf = dlog(thresh) ** 2
          do 250 ninf = nctp+10,mesh
            alpha2 = vme(ninf) * (r(ninf) - r(nctp))**2
            if ( alpha2 .gt. tolinf ) goto 260
  250     continue
c
c         jump point out of practical infinity loop
  260     continue
c
          if ( ninf .ge. mesh ) then
            write(iout,*) '***error in koesol'
            write(iout,*) 'ninf',ninf,' too big for this mesh',mesh
            call exit(1)
          endif
c
c         ===========================================================
c         analytic start up of minor and major components from origin
c         ===========================================================
c
c         with finite nucleus so potential constant at origin we have
c
c         f(r) = sum_n f_n r ** ( ig + n - 1 )
c         g(r) = sum_n g_n r ** ( ig + n )
c
c         with
c
c         g_n+1 = m(0) (v(0)-e) g_n-1 / ( (n+1) * (2l+n+2) )
c         f_n   = ( ig + n - 1 ) g_n / m(0)
c
c         ig = l+1 , g_0 = 1
c
          vmezer = vme(2)
          xmzer  = 1.0d0 - abyt**2 * vmezer
c
          ig = lcur + 1
c
          g0 = 1.0d0
          g1 = 0.0d0
          g2 = xmzer * vmezer * g0 / dble ( 4*lcur + 6 )
          f0 = dble( lcur ) * g0 / xmzer
          f1 = 0.0d0
          f2 = dble( lcur + 2 ) * g2 / xmzer
c
          yy(1,1) = 0.0d0
          yy(1,2) = 0.0d0
c
          do 300 ir = 2,6
            yy(ir,1) = r(ir)**(ig-1)*( f0 + r(ir)*( f1 + r(ir) * f2 ) )
            yy(ir,2) = r(ir)**ig * ( g0 + r(ir) * ( g1 + r(ir) * g2 ) )
  300     continue

c         ===========================
c         outward integration to nctp
c         ===========================
c
c         fifth order predictor corrector integration routine
          call cfdsol(zz,yy,7,nctp,idim1)
c
c         save major component and its gradient at nctp
          gout = yy(nctp,2)
          gpout = zz(nctp,2,1)*yy(nctp,1) + zz(nctp,2,2)*yy(nctp,2)
          gpout = gpout / rab(nctp)
c
c         ==============================================
c         start up of wavefunction at practical infinity
c         ==============================================
c
          do 310 ir = ninf,ninf-4,-1
            alpha = sqrt( vme(ir) )
            yy(ir,2) = exp ( - alpha * ( r(ir) - r(nctp) ) )
            yy(ir,1) = - alpha * yy(ir,2)
  310     continue
c
c         ==========================
c         inward integration to nctp
c         ==========================
c
c         fifth order predictor corrector integration routine
          call cfdsol(zz,yy,ninf-5,nctp,idim1)
c
c         save major component and its gradient at nctp
          gin = yy(nctp,2)
          gpin = zz(nctp,2,1)*yy(nctp,1) + zz(nctp,2,2)*yy(nctp,2)
          gpin = gpin / rab(nctp)
c
c
c         ===============================================
c         rescale tail to make major component continuous
c         ===============================================
c
          factor = gout / gin
          do 320 ir = nctp,ninf
            yy(ir,1) = factor * yy(ir,1)
            yy(ir,2) = factor * yy(ir,2)
  320     continue
c
          gpin = gpin * factor
c
c
c         =================================
c         check that the number of nodes ok
c         =================================
c
c         count the number of nodes in major component
          call nodeno(yy(1,2),2,ninf,nodes,idim1)
c
c         write(iout,*) 'number of nodes =',nodes,ncur-lcur-1
c
          if ( nodes .lt. ncur - lcur - 1 ) then
c           energy is too low
            emin = ecur
c           write(44,'(i5,3f12.5,2i5)')
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
c           write(44,'(i5,3f12.5,2i5)')
c    +         iter,emin,ecur,emax,nodes,ncur-lcur-1
            if ( ecur * 1.1d0 .lt. emin ) then
              ecur = 0.5d0 * ecur + 0.5d0 * emin
            else
              ecur = 1.1d0 * ecur
            endif
c           write(iout,*) 'classical turning point =',r(nctp)
c           write(iout,'(4(f10.5,e15.5))') (r(ir),yy(ir,2),ir=1,ninf)
            goto 200
          endif
c
c
c         =======================================================
c         find normalisation of wavefunction using simpson's rule
c         =======================================================
c
          factor = 0.0d0
          xw = 4.0d0
          do 360 ir = 2,ninf
            factor = factor + xw * rab(ir) * yy(ir,2)**2
            xw = 6.0d0 - xw
  360     continue
          factor = factor / 3.0d0
c
c
c         =========================================
c         variational improvement of the eigenvalue
c         =========================================
c
          decur = gout * ( gpout - gpin ) / factor
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
c         write(44,'(i5,3f12.5,1p2e12.4)')
c    +         iter,emin,ecur,emax,decur,decurp
c
c         test to see whether eigenvalue converged
          if ( abs(decur) .lt. thresh ) goto 400
c
          ecur = ecur + decurp
c
c
c         =======================================================
c         check that the iterative loop is not about to terminate
c         =======================================================
c
          if ( iter .eq. itmax ) then
c           eigenfunction has not converged in allowed number of iterations
            write(iout,*) '***error in subroutine koesol'
            write(iout,999) nnlz(iorb),ee(iorb),ecur
  999       format(' state',i4,' could not be converged.',/,
     +      ' starting energy for calculation was',f10.5,
     +      ' and end value =',f10.5)
            call exit(1)
          endif
c
c       close iterative loop
  200   continue
c
c       jump point on successful convergence of eigenvalue
  400   continue
c
c       update the energy
        ee(iorb) = ecur
c       update the real-space cutoff array
        nkkk(iorb) = ninf
c
c       update the wavefunction array
        factor = 1.0d0 / sqrt( factor )
        do 410 ir = 1,ninf
          snl(ir,iorb) = yy(ir,2) * factor
  410   continue
c
c
c     close loop over bands
  100 continue
c
c     stop
      return
      end
c
c-------------------------------------------------------------------------
c
      subroutine phakoe(ruae,r,rab,zz,ecur,lcur,
     +  yy,vme,klogd,dlwf,idim1)
c
c     subroutine to compute log derivatives of koelling and harmon's
c     scalar relativistic equation (j. phys. c 10, 3107 (1977))
c
c     in rydberg units the koelling harmon equations are:
c
c     df(r)     f(r)    l(l+1)
c     ----- = - ---- + -------- g(r) + ( v(r) - e ) g(r)
c      dr        r     m(r)*r*r
c
c     dg(r)                 g(r)
c     ----- =  m(r) f(r) +  ---
c      dr                    r
c
c     m(r) = 1 - alpha**2 * ( v(r) - e ) / 4
c
c     where g(r) is r*(spin-orbit averaged major component) and f(r)
c     is an analogue of r*(minor component) in the full dirac
c     equation.
c
c     to transform the above equations into form given in paper
c     must make substitutions phi(r) = f(r) / r , g(r) = g(r) / r
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1)
c.....potentials for the all electron calculation
      dimension ruae(idim1),zz(idim1,2,2)
c.....wavefunctions
      dimension yy(idim1,2)
c.....scratch
      dimension vme(idim1)
c
c.....variable for file = 0
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c-----------------------------------------------------------------------
c
c               r o u t i n e  i n i t i a l i s a t i o n
c
c     check that klogd is reasonable
      if ( klogd .lt. 10 ) then
        write(iout,*) '***error in subroutine phakoe'
        write(iout,*) 'klogd =',klogd,' is too small'
        call exit(1)
      endif
c
c     set ( 2 / fine structure constant )
      tbya = 2.0d0 * 137.04d0
c
c     *******test*****test******test******test****test
c     tbya = 1.0d6
c
c     set ( fine structure constant / 2 )
      abyt = 1.0d0 / tbya
c
      llp  = lcur * ( lcur + 1 )
c
c-----------------------------------------------------------------------
c
c     ===================
c     define the zz array
c     ===================
c
c     set the diagonal elements of zz
      do 100 ir = 2,klogd
        zz(ir,1,1) = - rab(ir) / r(ir)
        zz(ir,2,2) = + rab(ir) / r(ir)
  100 continue
c
c     set the off diagonal elements of zz
      do 110 ir = 2,klogd
c
        vme(ir) = ( ruae(ir) / r(ir) - ecur )
        xm  = ( 1.0d0 - vme(ir) * abyt**2 ) 
c
        zz(ir,1,2) = rab(ir) * 
     +    ( dble(llp) / ( xm * r(ir)**2 ) + vme(ir) )
        zz(ir,2,1) = rab(ir) * xm
c
  110 continue
c
c
c     ===========================================================
c     analytic start up of minor and major components from origin
c     ===========================================================
c
c     with finite nucleus so potential constant at origin we have
c
c     f(r) = sum_n f_n r ** ( ig + n - 1 )
c     g(r) = sum_n g_n r ** ( ig + n )
c
c     with
c
c     g_n+1 = m(0) (v(0)-e) g_n-1 / ( (n+1) * (2l+n+2) )
c     f_n   = ( ig + n - 1 ) g_n / m(0)
c
c     ig = l+1 , g_0 = 1
c
      vmezer = vme(2)
      xmzer  = 1.0d0 - abyt**2 * vmezer
c
      ig = lcur + 1
c
      g0 = 1.0d0
      g1 = 0.0d0
      g2 = xmzer * vmezer * g0 / dble ( 4*lcur + 6 )
      f0 = dble( lcur ) * g0 / xmzer
      f1 = 0.0d0
      f2 = dble( lcur + 2 ) * g2 / xmzer
c
      yy(1,1) = 0.0d0
      yy(1,2) = 0.0d0
c
      do 300 ir = 2,6
        yy(ir,1) = r(ir)**(ig-1)*( f0 + r(ir)*( f1 + r(ir) * f2 ) )
        yy(ir,2) = r(ir)**ig * ( g0 + r(ir) * ( g1 + r(ir) * g2 ) )
  300 continue

c     ============================
c     outward integration to klogd
c     ============================
c
c     fifth order predictor corrector integration routine
      call cfdsol(zz,yy,7,klogd,idim1)
c
c     compute the logarithmic derivative as major component wavefunction
c
      snlom = yy(klogd-1,2)
      snlop = yy(klogd,2)
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
c-----------------------------------------------------------------------
c
      subroutine cfdsol(zz,yy,jj1,jj2,idim1)
c
c-----------------------------------------------------------------------
c
c     routine for solving coupled first order differential equations
c
c      d yy(x,1)
c      ---------   =  zz(x,1,1) * yy(x,1) + zz(x,1,2) * yy(2,1)
c         dx
c
c      d yy(x,2)
c      ---------   =  zz(x,2,1) * yy(x,1) + zz(x,2,2) * yy(2,1)
c         dx
c     
c
c     using fifth order predictor corrector algorithm
c
c     routine integrates from jj1 to jj2 and can cope with both cases
c     jj1 < jj2 and jj1 > jj2.  first five starting values of yy must 
c     be provided by the calling program.
c
c-----------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
c
      dimension zz(idim1,2,2),yy(idim1,2)
      dimension fa(0:5),fb(0:5),abp(1:5),amc(0:4)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c-----------------------------------------------------------------------
c
c                   i n i t i a l i s a t i o n
c
c     decide whether integrating from:
c            left to right ---> isgn = + 1
c        or  right to left ---> isgn = - 1
c
      isgn = ( jj2 - jj1 ) / iabs( jj2 - jj1 )
c
c     run some test just to be conservative
      if ( isgn .eq. + 1 ) then
        if ( jj1 .le. 5 .or. jj2 .gt. idim1 ) then
          write(iout,10) isgn,jj1,jj2,idim1
          call exit(1)
        endif
      elseif ( isgn .eq. - 1 ) then
        if ( jj1 .ge. ( idim1 - 4 ) .or. jj2 .lt. 1 ) then
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
c     integration coefficients
c
      abp(1) = 1901.d0 / 720.d0
      abp(2) = -1387.d0 / 360.d0
      abp(3) = 109.d0 / 30.d0
      abp(4) = -637.d0 / 360.d0
      abp(5) = 251.d0 / 720.d0
      amc(0) = 251.d0 / 720.d0
      amc(1) = 323.d0 / 360.d0
      amc(2) = -11.d0 / 30.d0
      amc(3) = 53.d0 / 360.d0
      amc(4) = -19.d0 / 720.d0
c
c     set up the arrays of derivatives
      do 20 j = 1,5
        ip = jj1 - isgn * j
        fa(j) = zz(ip,1,1) * yy(ip,1) + zz(ip,1,2) * yy(ip,2)
        fb(j) = zz(ip,2,1) * yy(ip,1) + zz(ip,2,2) * yy(ip,2)
   20 continue
c
c-----------------------------------------------------------------------
c
c                i n t e g r a t i o n  l o o p
c
      do 100 j = jj1,jj2,isgn
c
c       predictor (adams-bashforth)
c
        arp = yy(j-isgn,1)
        brp = yy(j-isgn,2)
        do 110 i = 1,5
          arp = arp + dble(isgn) * abp(i) * fa(i)
          brp = brp + dble(isgn) * abp(i) * fb(i)
  110   continue
c
        fa(0) = zz(j,1,1) * arp + zz(j,1,2) * brp
        fb(0) = zz(j,2,1) * arp + zz(j,2,2) * brp
c
c       corrector (adams-moulton)
c
        yy(j,1) = yy(j-isgn,1)
        yy(j,2) = yy(j-isgn,2)
        do 120 i = 0,4,1
          yy(j,1) = yy(j,1) + dble(isgn) * amc(i) * fa(i)
          yy(j,2) = yy(j,2) + dble(isgn) * amc(i) * fb(i)
  120   continue
c
c       book keeping
c
        do 130 i = 5,2,-1
          fa(i) = fa(i-1)
          fb(i) = fb(i-1)
  130   continue
        fa(1) = zz(j,1,1) * yy(j,1) + zz(j,1,2) * yy(j,2)
        fb(1) = zz(j,2,1) * yy(j,1) + zz(j,2,2) * yy(j,2)
c
  100 continue
c
      return
      end
