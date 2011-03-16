c
c
      subroutine exccor( rho, pol, exc, vxc, vpol, ixctyp, imvex )
c
c     input
c
c     rho    the number density in inverse bohr radii cubed
c     pol    the difference of density (majority less minority)
c            divided by the sum majority plus minority
c     ixctyp the parameterization of the exchange correlation energy
c
c     ixctyp= 1 von Barth - Hedin polarization applied to
c               Wigner interpolation
c               Phys. Rev. 46, 1002 (1934)
c     ixctyp= 2 Hedin - Lundqvist
c               J. Phys. C, 4, 2064 (1971)
c     ixctyp= 3 Perdew - Zunger
c               Phys. Rev. B 23, 5048 (1981)
c     ixctyp= 4 Vosko - Wilk - Nusair
c               Can. J.  Phys. 58, 1200 (1980)
c               Phys. Rev. B 22, 3812 (1980)
c     ixctyp= 5 Perdew-Yang 92
c               Phys. Rev. B 45 13244 (1992)
c
c     imvex  inclusion of the Macdonald-Vosko correction for
c            the high density limit in the exchange
c            J. Phys. C, 12, 2977 (1979)
c
c            Not implemented for wigner exchange correlation
c
c     imvex= non-zero include it
c     imvex= 0 do not include
c
c     output
c
c     exc    the exchange correlation energy density in hartrees
c     vxc    the exchange correlation potential in hartrees
c            average over spins for polarized calculations
c     vpol   the difference in exchange correlation potentials
c
c (Note: majority spins sit in the potential vxc+vpol;
c        minority spins sit in the potential vxc-vpol.)
c
c     p refers to the paramagnetic or unpolarized limit
c     f refers to the ferromagnetic or polarized limit
c
      implicit double precision (a-h,o-z)
      parameter (zero= 0.0d0)
      parameter (one= 1.0d0)
      parameter (two= 2.0d0)
      parameter (three= 3.0d0)
      parameter (four= 4.0d0)
      parameter (half= 0.5d0)
      parameter (third= one/three)
      parameter (fthird= four*third)
      parameter (two4o3= 2.51984209978975d0)
      parameter (conrs= .6203504908993999d0)

      parameter (tiny = 1.d-12)
c
      if( rho .eq. zero ) then
         exc= zero
         vxc= zero
         vpol= zero
         return
      else if (rho .lt. zero) then
         print *, 'exccor:  rho ', rho, ' must be non-negative'
         stop 999
      endif
c
      if ( dabs(pol) .gt. one ) then
c Allow for pol to be very slightly greater than one, e.g., rounding problem.
         if (dabs(pol) .gt. one+tiny) then
             print *, 'exccor:  pol ', pol, ' must be <= 1 in abs value'
             print *, 'exccor:  dabs(pol)-1 ', dabs(pol)-one
             stop 999
         else if (pol .gt. zero) then
             pol = one
         else
             pol = -one
         endif
      endif
c
c     ------------------------
      if( ixctyp .eq. 1 ) then
c
c     wigner interpolation with rpa scaling for polarization
c     as in von barth - hedin
c
c     paramagentic limit
c
         x= rho**third
         call wigcor( x , ep, vp )
c
         if( pol .eq. zero ) then
c
c     unpolarized (paramagnetic)
c
            exc= ep
            vxc= vp
            vpol= zero
c
         else
c
c     ferromagentic limit
c
            xf= x * two4o3
            call wigcor( xf, ef, vf )
            ef= half * ef
            vf= half * vf
c
c     mix the paramagnetic and ferromagnetic limits
c
            call vbhmix(pol,ep,ef,vp,vf,exc,vxc,vpol)
c
         endif
c
c     -----------------------------
      else if( ixctyp .eq. 2 ) then
c
c     hedin - lundquist parameterization
c
         rs= conrs / rho**third
         rsin= one / rs
c
         if( pol .eq. zero ) then
c
c     unpolarized (paramagnetic) exchange
c
            call rpaexp(rsin,exp,vxp,imvex)
c
c     correlation in paramagnetic limit
c
            call hlcor( rsin, ecp, vcp )
c
            exc= exp + ecp
            vxc= vxp + vcp
            vpol= zero
c
         else
c
c     find the exchange energy and potential
c
            call rpaexf(rsin,pol,ex,vx,vxpol,imvex)
c
c     correlation in paramagnetic limit
c
            call hlcor( rsin, ecp, vcp )
c
c     correlation in ferromagnetic limit
c
            rsinf= rsin * two4o3
            call hlcor( rsinf, ecf, vcf )
c
c     scale by a factor of a half
c
            ecf= half * ecf
            vcf= half * vcf
c
c     mix the paramagnetic and ferromagnetic limits
c
            call vbhmix(pol,ecp,ecf,vcp,vcf,ec,vc,vcpol)
c
c     exchange and correlation
c
            exc= ex + ec
            vxc= vx + vc
            vpol= vxpol + vcpol
c
         endif
c
c     -----------------------------
      else if( ixctyp .eq. 3 ) then
c
c     perdew-zunger exchange and correlation
c
         rs= conrs / rho**third
         rsin= one / rs
c
         if( pol .eq. zero ) then
c
c     unpolarized (paramagnetic) exchange
c
            call rpaexp(rsin,exp,vxp,imvex)
c
c     correlation in paramagnetic limit
c
            call pzcop(rs,ecp,vcp)
c
            exc= exp + ecp
            vxc= vxp + vcp
            vpol= zero
c
         else
c
c     find the exchange energy and potential
c
            call rpaexf(rsin,pol,ex,vx,vxpol,imvex)
c
c     correlation in paramagnetic limit
c
            call pzcop(rs,ecp,vcp)
c
c     find the ferromagentic limit of the correlation
c
            call pzcof(rs,ecf,vcf)
c
c     interpolation between polarized and unpolarized
c
            call vbhmix(pol,ecp,ecf,vcp,vcf,ec,vc,vcpol)
c
c     exchange and correlation
c
            exc= ex + ec
            vxc= vx + vc
            vpol= vxpol + vcpol
c
         endif
c
c     -----------------------------
      else if( ixctyp .eq. 4 ) then
c
c     vosko - wilk - nusair correlation plus exchange
c
         rs= conrs / rho**third
         rsin= one / rs
         x= dsqrt( rs )
c
         if( pol .eq. zero ) then
c
c     unpolarized (paramagnetic) exchange
c
            call rpaexp(rsin,exp,vxp,imvex)
c
c     correlation in paramagnetic limit
c
            call vwncop(x,ecp,vcp)
c
            exc= exp + ecp
            vxc= vxp + vcp
            vpol= zero
c
         else
c
c     find the exchange energy and potential
c
            call rpaexf(rsin,pol,ex,vx,vxpol,imvex)
c
c     correlation in paramagnetic limit
c
            call vwncop(x,ecp,vcp)
c
c     correlation in ferromagnetic limit
c
            call vwncof(x,ecf,vcf)
            call vwncoa(x,eca,vca)
c
c     interpolation between the paramagnetic and ferromagnetic limits
c     for the correlation - not simple exchange scaling
c
            call vwnmix(pol,ecp,ecf,eca,vcp,vcf,vca,ec,vc,vcpol)
c
c     exchange and correlation
c
            exc= ex + ec
            vxc= vx + vc
            vpol= vxpol + vcpol
c
         endif

c     -----------------------------
      else if( ixctyp .eq. 5 ) then
         rs= conrs / rho**third
         rsin= one / rs
         
c     find the exchange energy and potential
         if( pol .eq. zero ) then
            call rpaexp(rsin,ex,vx,imvex)
            vxpol = 0.d0
         else
            call rpaexf(rsin,pol,ex,vx,vxpol,imvex)
         endif
c
c     Perdew-Wang 92 correlation
c
         call pwlsd (rs, pol, ec, vc, vcpol)

         exc = ex + ec
         vxc = vx + vc
         vpol= vxpol + vcpol
             
c     -----------------------------
      else
	  print *, 'exccor:  ixctyp ', ixctyp, ' out of range (1-5)'
          stop 999
      endif
c
      return
      end
c
c
c
      subroutine rpaexf(rsin,pol,ex,vx,vxpol,imvex)
c
c     rpa exchange energy for a spin polarized (ferromagnetic) electron density
c
      implicit double precision (a-h,o-z)
      parameter (zero= 0.0d0)
      parameter (one= 1.0d0)
      parameter (two= 2.0d0)
      parameter (three= 3.0d0)
      parameter (third= one/three)
      parameter (half= 0.5d0)
      parameter (tthird= two*third)
      parameter (cexrs= 0.4581652932831428d0)
      parameter (twotrd= 1.25992104989487d0)
      parameter (cmacv= 1.400477567739614E-02)
      parameter (thhalf= 1.5d0)
c
      fup= half * ( one + pol )
      fdn= half * ( one - pol )
      tmpup= fup**third
      tmpdn= fdn**third
c
c     exchange energy density is up plus down exchange energies
c     divided by the total charge density
c
      exup= -          twotrd * cexrs * fup * tmpup * rsin
      exdn= -          twotrd * cexrs * fdn * tmpdn * rsin
      vxup= - tthird * twotrd * cexrs *       tmpup * rsin
      vxdn= - tthird * twotrd * cexrs *       tmpdn * rsin
c
      if( imvex .ne. 0 ) then
c
c     macdonald-vosko correction for spin up
c
         beta= cmacv * twotrd * tmpup * rsin
         eta= dsqrt( one + beta**2 )
         xi= dlog( beta + eta )
         if( beta .gt. zero ) then
            scex= one - thhalf * ( ( beta*eta - xi ) / beta**2 )**2
            scvx= ( three * xi / ( two * beta * eta ) ) - half
         else
            scex= one
            scvx= one
         endif
         exup= exup * scex
         vxup= vxup * scvx
c
c     macdonald-vosko correction for spin down
c
         beta= cmacv * twotrd * tmpdn * rsin
         eta= dsqrt( one + beta**2 )
         xi= dlog( beta + eta )
         if( beta .gt. zero ) then
            scex= one - thhalf * ( ( beta*eta - xi ) / beta**2 )**2
            scvx= ( three * xi / ( two * beta * eta ) ) - half
         else
            scex= one
            scvx= one
         endif
         exdn= exdn * scex
         vxdn= vxdn * scvx
c
      endif
c
c     combined
c
      ex=    exup + exdn
      vx=    vxup + vxdn
      vxpol= vxup - vxdn
c
      return
      end
c
c
c
      subroutine rpaexp(rsin,ex,vx,imvex)
c
c     rpa exchange energy for a unpolarized (paramagnetic) electron density
c
      implicit double precision (a-h,o-z)
      parameter (one= 1.0d0)
      parameter (two= 2.0d0)
      parameter (three= 3.0d0)
      parameter (third= one/three)
      parameter (half= 0.5d0)
      parameter (tthird= two*third)
      parameter (cexrs= 0.4581652932831428d0)
      parameter (twotrd= 1.25992104989487d0)
      parameter (hlftrd= 0.7937005259840997d0)
      parameter (cmacv= 1.400477567739614E-02)
      parameter (thhalf= 1.5d0)
c
c     exchange energy density is up plus down exchange energies
c     divided by the total charge density
c
      ex= -                cexrs * rsin
      vx= - two * tthird * cexrs * rsin
c
      if( imvex .ne. 0 ) then
c
c     macdonald-vosko correction
c
         beta= cmacv * rsin
         eta= dsqrt( one + beta**2 )
         xi= dlog( beta + eta )
         scex= one - thhalf * ( ( beta*eta - xi ) / beta**2 )**2
         scvx= ( three * xi / ( two * beta * eta ) ) - half
         ex= ex * scex
         vx= vx * scvx
c
      endif
c
      return
      end
c
c
c
      subroutine vwncop(x,ec,vc)
c
c     Vosko - Wilk - Nusair parameterization of Ceperly - Alder
c     correlation energy for the paramagnetic limit
c     Can. J.  Phys. 58, 1200 (1980)
c     Phys. Rev. B 22, 3812 (1980)
c     on input x is the square root of rs
c
      implicit double precision (a-h,o-z)
      parameter (one= 1.0d0)
      parameter (two= 2.0d0)
      parameter (three= 3.0d0)
      parameter (third= one/three)
      parameter (a= 0.0310907d0)
      parameter (b= 3.72744d0)
      parameter (c= 12.9352d0)
      parameter (q= 6.15199081975908d0)
      parameter (x0= -0.10498d0)
      parameter (c1= two * b / q )
      parameter (c2= two * ( b + two * x0 ) / q )
      parameter (c3= b * x0 / ( c + b * x0 + x0**2 ) )
c
      x2= x*x
      xox= x2 + b * x + c
      taninq= datan( q / ( two * x + b ) )
      xxb= ( x2 + b*x ) / xox
      ec= dlog( x2 / xox )
     $     + c1 * taninq
     $     - c3 * ( dlog( (x-x0)**2 / xox )
     $            + c2 * taninq )
      ec= a * ec
      vc= one - xxb - c3 * ( x / (x-x0) - xxb - x*x0 / xox )
      vc= ec -third * a * vc
c
      return
      end
c
c
c
      subroutine vwncof(x,ec,vc)
c
c     Vosko - Wilk - Nusair parameterization of Ceperly - Alder
c     correlation energy for the ferromagnetic limit
c     Can. J.  Phys. 58, 1200 (1980)
c     Phys. Rev. B 22, 3812 (1980)
c     on input x is the square root of rs
c
      implicit double precision (a-h,o-z)
      parameter (one= 1.0d0)
      parameter (two= 2.0d0)
      parameter (three= 3.0d0)
      parameter (third= one/three)
      parameter (a= .01554535d0)
      parameter (b= 7.06042d0)
      parameter (c= 18.0578d0)
      parameter (q= 4.73092690956011d0)
      parameter (x0= -.325d0)
      parameter (c1= two * b / q )
      parameter (c2= two * ( b + two * x0 ) / q )
      parameter (c3= b * x0 / ( c + b * x0 + x0**2 ) )
c
      x2= x*x
      xox= x2 + b * x + c
      taninq= datan( q / ( two * x + b ) )
      xxb= ( x2 + b*x ) / xox
      ec= dlog( x2 / xox )
     $     + c1 * taninq
     $     - c3 * ( dlog( (x-x0)**2 / xox )
     $            + c2 * taninq )
      ec= a * ec
      vc= one - xxb - c3 * ( x / (x-x0) - xxb - x*x0 / xox )
      vc= ec -third * a * vc
c
      return
      end
c
c
c
      subroutine vwncoa(x,ec,vc)
c
c     Vosko - Wilk - Nusair parameterization of Ceperly - Alder
c     correlation contribution to the spin stiffness in the paramagnetic limit
c     Can. J.  Phys. 58, 1200 (1980)
c     Phys. Rev. B 22, 3812 (1980)
c     on input x is the square root of rs
c
      implicit double precision (a-h,o-z)
      parameter (one= 1.0d0)
      parameter (two= 2.0d0)
      parameter (three= 3.0d0)
      parameter (third= one/three)
      parameter (a= -.01688685d0)
      parameter (b= 1.13107d0)
      parameter (c= 13.0045d0)
      parameter (q= 7.12310891781812d0)
      parameter (x0= -.0047584d0)
      parameter (c1= two * b / q )
      parameter (c2= two * ( b + two * x0 ) / q )
      parameter (c3= b * x0 / ( c + b * x0 + x0**2 ) )
c
      x2= x*x
      xox= x2 + b * x + c
      taninq= datan( q / ( two * x + b ) )
      xxb= ( x2 + b*x ) / xox
      ec= dlog( x2 / xox )
     $     + c1 * taninq
     $     - c3 * ( dlog( (x-x0)**2 / xox )
     $            + c2 * taninq )
      ec= a * ec
      vc= one - xxb - c3 * ( x / (x-x0) - xxb - x*x0 / xox )
      vc= ec -third * a * vc
c
      return
      end
c
c
c
      subroutine pzcop(rs,ec,vc)
c
c     Perdew - Zunger  parameterization of Ceperly - Alder
c     correlation energy for the paramagnetic (unpolarized) limit
c     Phys. Rev. B 23, 5048 (1981)
c
      implicit double precision (a-h,o-z)
      parameter (one= 1.0d0)
      parameter ( b= 0.0480d0 )
      parameter ( a= 0.0311d0 )
      parameter ( d= 0.0116d0 )
      parameter ( c= 0.0020d0 )
      parameter ( ao3= a / 3.0d0 )
      parameter ( do3= d / 3.0d0 )
      parameter ( co3= c / 3.0d0 )
      parameter ( gamma= 0.1423d0 )
      parameter ( beta1= 1.0529d0 )
      parameter ( beta2= 0.3334d0 )
      parameter ( gb1o6= gamma * beta1 / 6.0d0 )
      parameter ( gb2o3= gamma * beta2 / 3.0d0 )
c
      if( rs .lt. one ) then
c
c     high density correlation
c
         alrs= dlog( rs )
         ec=  - b
     $        + a * alrs
     $        - d * rs
     $        + c * rs * alrs
         vc=  - rs * ( ao3 / rs
     $        - do3
     $        + co3 * ( one + alrs ) )
c
      else
c
c     low density correlation
c
         srrs= dsqrt( rs )
         denom= one + beta1 * srrs + beta2 * rs
         ec= - gamma / denom
         vc=  - rs * ( ( gb1o6 / srrs
     $                 + gb2o3 ) / denom**2 )
c
      endif
c
      vc= vc + ec
c
      return
      end
c
c
c
      subroutine pzcof(rs,ec,vc)
c
c     Perdew - Zunger  parameterization of Ceperly - Alder
c     correlation energy for the ferromagnetic (polarized) limit
c     Phys. Rev. B 23, 5048 (1981)
c
      implicit double precision (a-h,o-z)
      parameter (one= 1.0d0)
      parameter ( b= 0.0269d0 )
      parameter ( a= 0.01555d0 )
      parameter ( d= 0.0048d0 )
      parameter ( c= 0.0007d0 )
      parameter ( ao3= a / 3.0d0 )
      parameter ( do3= d / 3.0d0 )
      parameter ( co3= c / 3.0d0 )
      parameter ( gamma= 0.0843d0 )
      parameter ( beta1= 1.3981d0 )
      parameter ( beta2= 0.2611d0 )
      parameter ( gb1o6= gamma * beta1 / 6.0d0 )
      parameter ( gb2o3= gamma * beta2 / 3.0d0 )
c
      if( rs .lt. one ) then
c
c     high density correlation
c
         alrs= dlog( rs )
         ec=  - b
     $        + a * alrs
     $        - d * rs
     $        + c * rs * alrs
         vc=  - rs * ( ao3 / rs
     $        - do3
     $        + co3 * ( one + alrs ) )
c
      else
c
c     low density correlation
c
         srrs= dsqrt( rs )
         denom= one + beta1 * srrs + beta2 * rs
         ec= - gamma / denom
         vc=  - rs * ( ( gb1o6 / srrs
     $                 + gb2o3 ) / denom**2 )
c
      endif
c
      vc= vc + ec
c
      return
      end
c
c
c
      subroutine vbhmix(pol,ep,ef,vp,vf,e,v,vpol)
c
c     mixing between paramagnetic limits and ferromagnetic limits
c     based on rpa scaling for exchange
c     von Barth - Hedin
c     J. Phys. C 5, 1629 (1972)
c
      implicit double precision (a-h,o-z)
      parameter (one= 1.0d0)
      parameter (two= 2.0d0)
      parameter (three= 3.0d0)
      parameter (four= 4.0d0)
      parameter (third= one/three)
      parameter (fthird= four*third)
      parameter (cmix1= 1.92366105093154d0)
c
      fup= one + pol
      fdn= one - pol
      fupth= fup**third
      fdnth= fdn**third
      fmix= ( fup*fupth + fdn*fdnth - two ) * cmix1
      dfmix= fthird * ( fupth - fdnth ) * cmix1
c
      e=                ep * ( one - fmix ) + ef * fmix
      vpol= dfmix * ( ef - ep )
      v= - pol * vpol + vp * ( one - fmix ) + vf * fmix
c
      return
      end
c
c
c
      subroutine wigcor( x, exc, vxc )
c
c     wigner interpolation
c     Phys. Rev. 46, 1002 (1934)
c
      implicit double precision (a-h,o-z)
      parameter (one= 1.0d0)
c
      exc= -0.738d0 * x * ( one + 0.959d0 / ( one + 12.57d0 * x ) )
      vxc= -x * ( 0.984d0 + ( 0.943656d0 + 8.8963d0 * x ) /
     $     ( one + 12.57d0 * x)**2 )
c
      return
      end
c
c
c
      subroutine hlcor( rsin, ec, vc )
c
c     Hedin - Lundqvist parameterization of Singwi et al correlation
c     J. Phys. C, 4, 2064 (1971)
c
      implicit double precision (a-h,o-z)
      parameter (rscut= 0.01d0)
      parameter (rp= 21.0d0)
      parameter (cp= 0.0225d0)
      parameter (fiotw= 2.5d0)
      parameter (tquart= 0.75d0)
      parameter (one= 1.0d0)
      parameter (three= 3.0d0)
      parameter (four= 4.0d0)
      parameter (half= 0.5d0)
      parameter (third= one/three)
c
      x= rp * rsin
      x2= x * x
      x3= x2 * x
c
      temp= dlog( one + x )
      vc= -cp * temp
c
      if( x .lt. rscut ) then
         ec= tquart * x * ( one
     $                     - x  * ( four  - x3 ) / 10.d0
     $                     + x2 * ( three - x3 ) / 13.5d0
     $                     - x3 * ( fiotw - x3 ) / 17.5d0 )
      else
         ec= temp + ( temp
     $               - ( x - half * x2 + third * x3 ) ) / x3
      endif
c
      ec= -cp * ec
c
      return
      end
c
c
c
      subroutine vwnmix(pol,ecp,ecf,eca,vcp,vcf,vca,ec,vc,vcpol)
c
c     mixing between paramagnetic limits and ferromagnetic limits
c     Vosko - Wilk - Nusair
c     Can. J.  Phys. 58, 1200 (1980)
c     Phys. Rev. B 22, 3812 (1980)
c
      implicit double precision (a-h,o-z)
      parameter (one= 1.0d0)
      parameter (two= 2.0d0)
      parameter (three= 3.0d0)
      parameter (four= 4.0d0)
      parameter (third= one/three)
      parameter (fthird= four*third)
      parameter (cmix1= 1.92366105093154d0)
      parameter (cfppin= 0.5848223622634647d0)
c
      fup= one + pol
      fdn= one - pol
      fupth= fup**third
      fdnth= fdn**third
      fmix= ( fup*fupth + fdn*fdnth - two ) * cmix1
      dfmix= fthird * ( fupth - fdnth ) * cmix1
c
      pol3= pol**3
      pol4= pol3 * pol
      fmpol4= fmix * pol4
      dfmp4= four * pol3 * fmix + pol4 * dfmix
c
      ec=    ecp * (  one - fmpol4 )
     $     + ecf *          fmpol4
     $     + eca * ( fmix - fmpol4 ) * cfppin
      vcpol=
     $     - ecp * dfmp4
     $     + ecf * dfmp4
     $     + eca * ( dfmix - dfmp4 ) * cfppin
      vc=  - pol * vcpol
     $     + vcp * (  one - fmpol4 )
     $     + vcf *          fmpol4
     $     + vca * ( fmix - fmpol4 ) * cfppin
c
      return
      end

c Implementation of Perdew-Wang Phys. Rev. B 45, 13244 (1992) correlation
c simplified version of the subroutine CORLSD distributed by J. P. Perdew
c This version runs about twice as fast and gives almost exactly the same
c results.

c
c input
c  rs: electron gas parameter
c  zet: spin polarization   (rhoUp-rhoDn)/(rhoUp+rhoDn)
c
c output
c  ec: correlation energy per electron in hartree
c  vc: average correlation potential (called COMM in CORLSD)
c  eczet: difference potential, also d ec / d zet
c      vcUp = vc + ecZet
c      vcDn = vc - ecZet
c
c internal: (reported by CORLSD)
c  ALFC = -ALFM is the correlation contribution to the spin stiffness
c  ecrs: d ec / d rs

      subroutine pwlsd (rs, zet, ec, vc, eczet)
      double precision alfm, alfrsm, ec, ecrs, eczet, ep, eprs,
     $                 eu, eurs, f, fz, fzCon, fzzInv, gamInv,
     $                 rs, rs12, thrd, vc,
     $                 z3, z4, zet, zet1m, zet3m, zet1p, zet3p
      double precision a, a1, b1, b2, b3, b4, q0, q1, q2, q3

      parameter (thrd=1.d0/3.d0)
      parameter (gamInv=1.d0/0.5198421D0, fzzInv=1.d0/1.709921D0)
      parameter (fzCon = 4.d0*thrd*gamInv)

      zet1p = 1.d0 + zet
      zet1m = 1.d0 - zet
      zet3p = zet1p**thrd
      zet3m = zet1m**thrd
      f = ( zet1p*zet3p + zet1m*zet3m - 2.d0 ) * gamInv
      rs12 = dsqrt(rs)

      a  = 0.0310907D0
      a1 = 0.21370D0
      b1 = 7.5957D0
      b2 = 3.5876D0
      b3 = 1.6382D0
      b4 = 0.49294D0

      Q0 = -2.D0*A*(1.D0+A1*rs)
      q1 = 2.d0*a* rs12* (b1+ rs12* (b2+ rs12*(b3 + rs12* b4)))
      Q2 = DLOG(1.D0+1.D0/Q1)
      eu = Q0*Q2
      Q3 = A*(B1/RS12+2.D0*B2+ rs12* (3.d0*b3 + rs12* 4.d0*b4))
      eurs = -2.D0*A*A1*Q2-Q0*Q3/(Q1**2+Q1)
  
      a  =  0.01554535D0
      a1 =  0.20548D0
      b1 = 14.1189D0
      b2 =  6.1977D0
      b3 =  3.3662D0
      b4 =  0.62517D0

      Q0 = -2.D0*A*(1.D0+A1*rs)
      q1 = 2.d0*a* rs12* (b1+ rs12* (b2+ rs12*(b3 + rs12* b4)))
      Q2 = DLOG(1.D0+1.D0/Q1)
      ep = Q0*Q2
      Q3 = A*(B1/RS12+2.D0*B2+ rs12* (3.d0*b3 + rs12* 4.d0*b4))
      eprs = -2.D0*A*A1*Q2-Q0*Q3/(Q1**2+Q1)

      a  =  0.0168869D0
      a1 =  0.11125D0
      b1 = 10.357D0
      b2 =  3.6231D0
      b3 =  0.88026D0
      b4 =  0.49671D0

      Q0 = -2.D0*A*(1.D0+A1*rs)
      q1 = 2.d0*a* rs12* (b1+ rs12* (b2+ rs12*(b3 + rs12* b4)))
      Q2 = DLOG(1.D0+1.D0/Q1)
      alfm = Q0*Q2
      Q3 = A*(B1/RS12+2.D0*B2+ rs12* (3.d0*b3 + rs12* 4.d0*b4))
      alfrsm = -2.D0*A*A1*Q2-Q0*Q3/(Q1**2+Q1)

      z3 = zet * zet * zet
      z4 = z3 * zet
      EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4) * fzzInv

C  ENERGY DONE. NOW THE POTENTIAL:
      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4) * fzzInv
      FZ = fzCon * ( zet3p - zet3m )
      ecZet = 4.D0*z3*F*(EP-EU+ALFM*fzzInv)+FZ*(Z4*(EP-EU)
     $        -(1.D0-Z4)*ALFM*fzzInv)
      vc = ec - thrd * rs * ecRs - zet * ecZet

      return
      end


