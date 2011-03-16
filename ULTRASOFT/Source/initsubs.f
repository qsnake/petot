c
c Copyright (C) 2002 Vanderbilt group
c This file is distributed under the terms of the GNU General Public
c License as described in the file 'License' in the current directory.
c
c----------------------------------------------------------------------------
c
      subroutine rinit(r,rab,sqr,rmax,aasf,bbsf,a,b,z,rlogd,klogd,
     +  mesh,ifprt,idim1)
c
c     routine to generate logarithmic radial mesh
c     after algorithm due to sverre froyen
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c
c.....logarithmic radial mesh information
      dimension r(idim1),rab(idim1),sqr(idim1)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c----------------------------------------------------------------------------
c
c           g e n e r a t e  t h e  l o g a r i t h m i c  m e s h
c
      a = exp( - aasf ) / z
      b = 1.0d0 / bbsf
c
c     compute the mesh parameter
      mesh = 2 + int( dlog( rmax / a + 1.0d0 ) / b )
c
c     ensure that mesh is odd for simpson's rule integration
      if ( 2 * ( mesh / 2 ) .eq. mesh ) mesh = mesh + 1
c
      if ( mesh .gt. idim1 ) then
        write(iout,*) '***error in subroutine rinit'
        write(iout,*) 'mesh =',mesh,' gt idim1 =',idim1
        call exit(1)
        endif
c
      do 100 i = 1,mesh
        r(i) = a * ( exp( b*dble(i-1) ) - 1.0d0 )
        rab(i) = b * ( r(i) + a )
        sqr(i) = exp ( b * dble(i) / 2.0d0 )
  100 continue
c
c     further check for self-consistency for development purposes
      if ( r(mesh) .lt. rmax ) then
        write(iout,*) '***error in subroutine rinit'
        write(iout,*) 'r(mesh)=',r(mesh),' .lt. rmax =',rmax
        call exit(1)
        endif
c
      write(iout,110) mesh
  110 format(/' value of mesh generated in rinit is',i6)
c
c----------------------------------------------------------------------------
c
c     p o i n t  t o  e v a l u a t e  l o g  d e r i v a t i v e s
c
c     compute the index for evaluating the logarithmic derivatives
      klogd = 2 + int( dlog( rlogd / a + 1.0d0 ) / b )
c
c     check that i seems reasonable
      if ( klogd .lt. 10 .or. klogd .eq. mesh ) then
        write(iout,*) '***error in subroutine rinit'
        write(iout,*) 'klogd =',klogd,' not allowed with',
     +    ' mesh =',mesh
        call exit(1)
      endif
c
c----------------------------------------------------------------------------
c
c        p r i n t  o u t  s u m m a r y  o f  r e s u l t s
c
      if ( ifprt .ge. 1 ) write(iout,170) rlogd,klogd,r(klogd)
c
  170 format(' rlogd =',f9.5,' generated klogd =',i5,' with',
     +' r(klogd) =',f9.5)
c
      if ( ifprt .ge. 1 ) then
        write(iout,210)
        do 220 i = 1,10
  220     write(iout,230) i,r(i),rab(i)
        do 222 i = 1,3
  222     write(iout,240)
        do 224 i = mesh-10,mesh
  224     write(iout,230) i,r(i),rab(i)
      endif
c
  210 format(/10x,'logarimic radial mesh values:',/,/,
     +5x,'index',12x,'r',18x,'rab'/)
  230 format(4x,i4,6x,f13.8,7x,f13.8)
  240 format(7x,'.',10x,'.',19x,'.')
c
      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine startv(ipass,z,xion,mesh,r,runuc,ruae,
     +  wwnl,wwnlps,nvales,nvalps,ncores,zv,vloc,vloc0,
     +  idim1,idim2,idim4)
c
c     routine establishes an analytic starting potential
c
c             ipass = 1    all electron calculation
c             ipass = 2    vanderbilt pseudopotential
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c.....logarithmic radial mesh information
      dimension r(idim1)
c.....potentials for the all electron calculation
      dimension runuc(idim1),ruae(idim1)
c.....ion occupancies
      dimension wwnl(idim2)
c.....pseudo ion occupancies
      dimension wwnlps(idim4)
c.....local components of the vanderbilt potential
      dimension vloc(idim1),vloc0(idim1)
c
      if ( ipass .eq. 1 ) then
c
        drng=0.675
        zeff= z-(xion+0.5)
        yy=zeff**.4
c
        do 100 i=1,mesh
          tt=dmin1(1.d+02,(r(i)/drng))
          tt=2.0*zeff*(1.-1.0/(drng*yy*(exp(tt)-1.)+1.))
          ruae(i)=tt + runuc(i)
  100   continue
c
      elseif ( ipass .eq. 2 ) then
c
        wsumps = 0.0d0
        do 200 i = 1,nvalps
          wsumps = wsumps + wwnlps(i)
  200   continue
c
        wsum = 0.0d0
        do 210 i = ncores+1,ncores+nvales
          wsum = wsum + wwnl(i)
  210   continue
c
c       avoid neutral starting potential to insure states bound
        wsum = dmin1(wsum,zv-0.25d0)
        wratio = wsum/wsumps
c
c       screening corrected so that right number of electrons are there
        do 220 i = 1,mesh
          vscr = vloc(i) - vloc0(i)
          vloc(i) = vloc0(i) + vscr * wratio
  220   continue
c
      endif
c
      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine setnuc(z,r,runuc,sigma,mesh,idim1,ifprt)
c
c     routine to establish a nuclear potential with gausian
c     cuttoff in the core.  sigma is the effective nuclear
c     radius.
c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
c
c.....logarithmic radial mesh information
      dimension r(idim1)
c.....potentials for the all electron calculation
      dimension runuc(idim1)
c
      integer stderr
      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c     set the value of sigma to the fifth mesh point
c
      sigma = r(30)
      rsigma = 1.0d0 / ( sigma )
c
      do 100 i = 1,mesh
        arg = rsigma * r(i)
        if ( arg .lt. 20.0d0 ) then
          tt = exp ( - arg )
        else
          tt = 0.0d0
        endif
        runuc(i) = - 2.0d0 * z * ( 1.0d0 - tt )
 100  continue
c
      if ( ifprt .eq. 1 ) then
        write(iout,200) sigma
c       write(iout,210)
c       do 150 i = 1,mesh
c         write(iout,220) i,r(i),runuc(i),-2.0d0*z
c 150   continue
      endif
c
  200 format(' subroutine setnuc - nuclear core radius =',g15.5)
  210 format(/,4x,'i',12x,'r',16x,'runuc',10x,' - 2.0d0 * z ',/)
  220 format(i5,5x,f10.5,10x,f10.5,10x,f10.5)
c
      return
      end
c
c-------------------------------------------------------------------------
