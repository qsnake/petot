       subroutine velect(iter,iconv,icorr,ispp,ifcore,
     1 nr,r,rab,zel,cdd,cdu,cdc,vod,vou,etot,y,yp,
     2 ypp,s1,s2,w)
c
c    velect generates the electronic output potential from
c    the electron charge density.  The ionic part is 
c    added in dsolv1/dsolv2.
c
c  This version is a  revision of the exchange-correlation part
c  to introduce the interface ATOMXC which allows GGA calculations
c  with Perdew'96 and others xc functionals. 
c  J.M. SOLER and L.C. BALBAS, january 1997.
c
c   version 5.62
c                              
       implicit double precision (a-h,o-z)
c
       character*1 ispp
       character*2 icorr
c
c  njtj  *** modification start  ***
       parameter (zero=0.D0,one=1.D0,pfive=.5D0,opf=1.5D0,pnn=.99D0)
Cray       parameter (zero=0.0,one=1.0,pfive=0.5,opf=1.5,pnn=0.99)
c  njtj  *** modification end  ***
c
      dimension r(nr),rab(nr),cdd(nr),cdu(nr),cdc(nr),
     1 vod(nr),vou(nr),etot(10),y(nr),yp(nr),ypp(nr),
     2 s1(nr),s2(nr),w(3*nr)

      parameter ( maxr = 1500 )
      dimension dens(maxr,2),vxc(maxr,2)

      pi=4*atan(one)
c
c------Machine dependent parameter-
c------Require exp(-2*expzer) to be within the range of the machine
c
c      expzer = 3.7D2
cApollo      expzer = 3.7D2
cSun      expzer = 3.7D2
cVax      expzer = 44.D0
Cray      expzer = 2.8E3
c
c      fit cd/r by splines
c
       y(1) = zero
       do 10 i=2,nr
         y(i) = (cdd(i)+cdu(i))/r(i)
 10    continue
       if (ifcore .eq. 2) then
         do 11 i=2,nr
           y(i) = y(i) + cdc(i)/r(i)
 11      continue
       endif
       isx = 0
       a1 = zero
       an = zero
       b1 = zero
       bn = zero
       nrm=nr
       call splift(r,y,yp,ypp,nrm,w,ierr,isx,a1,b1,an,bn)
       if(ierr.ne.1) then
         write(6,20000)ierr
         call ext(420+ierr)
       endif
20000  format(1x,'****** Error in splift ierr =',i2)
c
c      compute the integrals of cd/r and cd from
c      r(1)=0 to r(i)
c
       xlo = zero
       call spliq(r,y,yp,ypp,nrm,xlo,r,nrm,s2,ierr)
       if(ierr.ne.1) then
         write(6,20001)ierr
         call ext(440+ierr)
       endif
20001  format(1x,'****** Error in spliq ierr =',i2)
       do 20 i=1,nr
         ypp(i) = r(i)*ypp(i) + 2*yp(i)
         yp(i)  = r(i)*yp(i)  + y(i)
         y(i)   = r(i)*y(i)
 20    continue
       call spliq(r,y,yp,ypp,nrm,xlo,r,nrm,s1,ierr)
       if(ierr.ne.1) then
         write(6,20002)ierr
         call ext(460+ierr)
       endif
20002  format(1x,'****** Error in spliq ierr =',i2)
c
c      check normalization
c
       xnorm = zero
       if (ifcore .eq. 2 .and. iter .eq. 0 ) zel=s1(nr)
       if (zel .ne. zero) xnorm = zel/s1(nr)
       if (iter .gt. 3 .and. abs(zel-s1(nr)) .gt. 0.01) then
         if (zel .lt. s1(nr)+1.0 ) then
           write(6,24) iter,xnorm
 24    format(/,' warning *** charge density rescaled in',
     1 ' velect',/,' iteration number',i4,3x,
     2 'scaling factor =',f6.3,/)
         else
           xnorm=pnn*xnorm
           write(6,25) iter,xnorm
 25    format(/,' warning *** charge density partially rescaled in',
     1 ' velect',/,' iteration number',i4,3x,
     2 'scaling factor =',f6.3,/)
         endif
       endif
c
c      compute new hartree potential
c      renormalize the charge density
c
       do 30 i=2,nr
         vod(i) = 2 * xnorm*(s1(i)/r(i) + s2(nr) - s2(i))
         vou(i) = vod(i)
         cdd(i) = xnorm*cdd(i)
         cdu(i) = xnorm*cdu(i)
 30    continue
c
c      compute hartree contribution to total energy
c
       if (iconv .eq. 1) then
         ehart = zero
         ll = 4
         do 40 i=2,nr
           ehart = ehart+ll*(cdd(i)+cdu(i))*vod(i)*rab(i)
           ll = 6 - ll
 40      continue
         ehart = ehart / 6 
       endif       
c
c     Add exchange and correlation
c     This part is totally new. It is written to use
c     the package XC of J.M. Soler.
c     J.M. Soler and L.C. Balbas. January 1997  
c
c     Compute dens(i,nspin) = density up, density down
c
      if (nr.gt.maxr) then
        write(6,*) 'velect: maxr too small. nr, maxr =', nr, maxr
        stop
      endif

      do 90 i=2,nr
        if (ispp .eq. 's') then
	  dens(i,1) = cdu(i)/(4.d0*pi*r(i)**2)
	  dens(i,2) = cdd(i)/(4.d0*pi*r(i)**2)
        else
          dens(i,1) = 0.5d0*(cdu(i) + cdd(i))/(4.d0*pi*r(i)**2)
          dens(i,2) = dens(i,1)
        endif
        if (ifcore .ge. 1) then
              dens(i,1) = dens(i,1) + 0.5d0 * cdc(i)/(4.d0*pi*r(i)**2)
              dens(i,2) = dens(i,2) + 0.5d0 * cdc(i)/(4.d0*pi*r(i)**2)
        endif
  90  continue

c
c     Extrapolate the density at r=0  
c
      dens(1,1) = dens(2,1) - (dens(3,1)-dens(2,1))*r(2)/(r(3)-r(2))
      dens(1,2) = dens(2,2) - (dens(3,2)-dens(2,2))*r(2)/(r(3)-r(2))
c
c     Define 'irel' and 'nspin' for the interface ATOMXC
c
      if (ispp .eq. 'r') irel = 1
      if (ispp .ne. 'r') irel = 0
      nspin = 2
c
      r(1) = 0.0d0
           if (icorr .eq. 'ca') then
 	      call atomxc('LDA','ca',irel,nr,maxr,r,nspin,dens,
     .                     ex,ec,dx,dc,vxc)      
	    elseif(icorr .eq. 'pw') then
              call atomxc('LDA','pw92',irel,nr,maxr,r,nspin,dens,
     .                    ex,ec,dx,dc,vxc)
            elseif(icorr .eq. 'pb') then
              call atomxc('GGA','pbe',irel,nr,maxr,r,nspin,dens,   
     .                     ex,ec,dx,dc,vxc)
           endif

c
c   Add vxc to total potential and energies
c   
      do 150 i=2,nr
        vou(i) = vou(i) + vxc(i,1)
        vod(i) = vod(i) + vxc(i,2)
  150 continue

clcb
      xccor=0.0d0
      ll = 4
      do 160 i=2,nr
        xccor = xccor + ll * rab(i)*
     .             (vxc(i,1)*cdu(i) + vxc(i,2)*cdd(i))
        ll = 6 - ll
  160 continue

      etot(4) = ehart
      etot(5) = xccor / 3
      etot(6) = xccor - 4*(ex + ec)
      etot(7) = ex + ec 

clcb
c
c     Obtain total potential at r = 0
c
      vod(1) = vod(2) - (vod(3)-vod(2))*r(2)/(r(3)-r(2))
      vou(1) = vou(2) - (vou(3)-vou(2))*r(2)/(r(3)-r(2))

c     *** lcb-jms modification end ********************
      return
      end

