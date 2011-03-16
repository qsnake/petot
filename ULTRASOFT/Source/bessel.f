c
c Copyright (C) 2002 Vanderbilt group
c This file is distributed under the terms of the GNU General Public
c License as described in the file 'License' in the current directory.
c
c---------------------------------------------------------------
      subroutine besdiag(ifprt,iout,rm,ecut,lmx,mesh,r,rab,vloc,
     *     kkbeta,idim1,idim3,nbl0,nbl,ddd,qqq,beta,ne)
c---------------------------------------------------------------
c
c-   nonselfconsistent program for eigenvalues in a
c-   spherical bessel functions basis sets
c     units: atomic rydberg units
c
c---------------------------------------------------------------
c
      implicit none
      integer ifprt,                 !-3 no prints, -1 more or less standard,1  everything
     *        iout                   ! fileunit for output
      integer  lmx,mesh,ne,kkbeta,idim1,idim3,nbl0(lmx+1),nbl(lmx+1)
      real*8    ecut,                ! kinetic-energy cutoff
     *          rm,                  ! radius of the box
     *          r(mesh),rab(mesh),   !logarithmic mesh
     *          vloc(mesh),
     *          ddd(idim3,idim3),qqq(idim3,idim3),
     *          beta(idim1,idim3)
      
      integer   nswx,         ! maximum number of spherical waves
     *          lmax,          ! maximum angular momentum
     *          nptx,nbt,mshm

      parameter(nswx=200, lmax=3, nptx=4*nswx)
      parameter (nbt=1501,mshm = 2000)
      integer   nsw(0:lmax),     ! actual number of spherical waves
     *          i,l,
     *           iret         ! number of points for radial integration

      real*8    pi,           !
     *          q(nswx,0:lmax),! quantized momenta in the spherical box
     *          s(nswx,0:lmax),! normalization constants for the bessel's
     *          h(nswx,nswx), chi(nswx,nswx), enl(nswx) , qm(nswx,nswx),
     *          bessel,             ! spherical bessel functions
     *          fint(nbt), work1(2*nswx),
     *          rbt(nbt),vbt(nbt),y2(mshm),u(mshm),dr

      external bessel
      pi = 2.0d0*asin(1.0d0)

31    if(lmx.gt.lmax) call error('besdiag','increase lmax',1)

c init r grid
      dr = rm/(nbt-1)
      do i=1,nbt
         rbt(i)=(i-1)*dr
      end do

c init q grid

      do i=1, nswx
         q(i,0) = i*pi
      end do

      do l = 1,lmx
         do i=1, nswx-l
            call find_root(bessel,l,q(i,l-1),q(i+1,l-1),q(i,l),
     *           1.0d-10,iret)
            call error('besdiag','error in FIND_ROOT',iret)
         end do
      enddo

      do l=0,lmx
         do i=1, nswx-l
            q(i,l) = q(i,l) / rm
         end do
      end do

      do l=0,lmx
         do i=1,nswx-l
            if (q(i,l)**2.gt.ecut) then
	       nsw(l)=i
	       goto 111
            endif
         enddo
         call error('nswx','nswx is too small',nswx)
111      continue
      end do

      if ( ifprt .ge. 1) then
         write(iout,9000 ) ecut,rm,nbt
 9000    format(/3x,'ecut=',f6.1,'  rmax=',f6.1,'   nbt=',i5)
         do l=0,lmx
            write(iout,'(3x,''nsw(l='',i1,'')='', i5)') l, nsw(l)
         end do
      endif
c
      call s_fill(nbt,rbt,dr,fint,nsw,nswx,lmx,q,s)
c
      write(iout,'(/3x,a)') 'eigenvalues of local potential '
      do l = 0, lmx
         call beseigen (mesh,r,vloc,nbt,rbt,vbt,y2,u,dr,fint,kkbeta,
     *       idim1,idim3,nbl0(l+1),nbl(l+1),ddd,qqq,beta,nsw(l),nswx,
     *       l,q(1,l), s(1,l),h,qm,chi,enl,work1,0)
         write(iout,9010) l,(enl(i),i=1,2)
 9010     format(3x,'l=',i2,4x,2f11.6)
      enddo
c
      write(iout,'(/3x,a)') 'eigenvalues of US hamiltonian '
      do l = 0, lmx
         call beseigen (mesh,r,vloc,nbt,rbt,vbt,y2,u,dr,fint,kkbeta,
     *       idim1,idim3,nbl0(l+1),nbl(l+1),ddd,qqq,beta,nsw(l),nswx,
     *       l,q(1,l), s(1,l),h,qm,chi,enl,work1,1)
         write(iout,9020) l,(enl(i),i=1,ne)
 9020     format(3x,'l=',i2,4x,5f11.6)
      enddo
c
      return
      end


      subroutine beseigen (npt,r,vloc,nbt,rbt,vbt,y2,u,dr,fint,kkbeta,
     *   idim1,idim3,nbl0,nbl,ddd,qqq,beta,nsw,nswx,l,q,s,
     *     h,qm,chi,enl,work1,ihnl)
*     ===================
*-----------------------------------------------------------------------
c    solves the radial schrodinger equation for an l channel
c
      implicit none
      integer npt,nbt,nsw,nswx,l,kkbeta,idim1,idim3,nbl0,nbl,ihnl,
     *  ifail
      real*8 r(npt),vloc(npt),rbt(nbt),vbt(nbt),dr,fint(nbt),
     *     ddd(idim3,idim3),y2(npt),u(npt),
     *  qqq(idim3,idim3),beta(idim1,idim3),s(nswx),q(nswx),h(nswx,nswx),
     *  qm(nswx,nswx),chi(nswx,nswx),enl(nswx),work1(2*nswx)
c first find the eigenstates of the local part
      call h_fill(npt,r,vloc,nbt,rbt,vbt,y2,u,dr,fint,nsw,nswx,l,q,s,h)
      call h_us(npt,r,vloc,nbt,rbt,vbt,y2,u,dr,fint,nsw,nswx,l,
     *    kkbeta,idim1,idim3,
     *     nbl0,nbl,ddd,qqq,beta,q,s,h,qm,chi,ihnl)
c
c solve ( h - enl qm) chi
c
c     --- nag library subroutine ---
c      call f02aef(h,nswx,qm,nswx,nsw,enl,chi,nswx,work1,work1(nswx),
c    *     ifail)
c
c     --- essl library subroutine ---
c      call dsygv(0,h,nswx,qm,nswx,enl,chi,nswx,nsw,work1,2*nswx)
c
c     --- eispack library subroutine ---
       call rsg(nswx,nsw,h,qm,enl,1,chi,work1,work1(nsw+1),ifail)
c
      return
      end

*-----------------------------------------------------------------------
      subroutine h_us(npt,r,vloc,nbt,rbt,vbt,y2,u,dr,fint,nsw,nswx,l,
     *     kkbeta,idim1,idim3,nbl0,nbl,ddd,qqq,beta,q,s,h,qm,work,ihnl)
*
* write the ultra soft hamiltonian in a basis of besselfunctions
*-----------------------------------------------------------------------
      implicit none
      integer  npt,nbt,nsw,nswx,l,kkbeta,idim1,idim3,nbl0,nbl,ihnl

      real*8 r(npt), vloc(npt),rbt(nbt),vbt(nbt),dr,fint(nbt),
     *   ddd(idim3,idim3),qqq(idim3,idim3),beta(idim1,idim3),    
     *   s(nswx), q(nswx), h(nswx,nswx), qm(nswx,nswx),
     *   work(nswx,nbl),y2(npt),u(npt)
      
      real *8 asum,bessel
      integer ib,jb,ib0,jb0,n,m,i
      external bessel,radlin,interpol
c   add the nonlocal part to the hamiltonian
c-
c-    bessel transform the projectors
c-
      do ib =1,nbl
         ib0 = ib+nbl0
         call interpol(r,beta(1,ib0),npt,rbt,vbt,nbt,y2,u)
         do n=1,nsw
            do i=1,nbt
               fint(i) = vbt(i)*rbt(i)*bessel(q(n)*rbt(i),l)
            end do
            call radlin(nbt,fint,dr,asum)
            work(n,ib) = asum/s(n)
         end do
      enddo
calculate the hamiltonian and the qnm matrix
      do n = 1,nsw
         do m=1,nsw
            qm(n,m) = 0.d0
         enddo
         qm(n,n) = 1.d0
      enddo
      do n=1,nsw
         do m=1,n
            do ib =1,nbl
               ib0 = ib+nbl0
               do jb =1,nbl
                  jb0=jb+nbl0
                  if (ihnl .gt. 0) h(n,m) = 
     $                 h(n,m)+work(n,ib)*work(m,jb)*ddd(ib0,jb0)
                  qm(n,m) = qm(n,m)+work(n,ib)*work(m,jb)*qqq(ib0,jb0)
               end do
            end do
            h(m,n)=h(n,m)
            qm(m,n)=qm(n,m)
         enddo
      enddo
      return
      end

*-----------------------------------------------------------------------
      subroutine s_fill   (npt,r,dr,fint,nsw,nswx,lmx,q,s)
*     =================
*-----------------------------------------------------------------------
      implicit none
      integer n, lmx, nsw(0:lmx), nswx, l, i, npt
      real*8 s(nswx,0:lmx), q(nswx,0:lmx), r(npt),fint(npt)
      real*8 bessel,asum,dr
      external bessel,radlin

      do l=0,lmx
         do n=1, nsw(l)
            do i=1, npt
               fint(i) = (r(i)*bessel(q(n,l)*r(i),l))**2
            end do
            call radlin(npt,fint,dr,asum)
            s(n,l) = sqrt( asum)
         end do
      end do

      return
      end

*-----------------------------------------------------------------------
      subroutine h_fill(npt,r,vloc,nbt,rbt,vbt,y2,u,dr,fint,nsw,nswx,
     *                l,q,s,h)
*     =================
*-----------------------------------------------------------------------
      implicit none
      integer n, m, nsw, nswx, l, i, npt,nbt
      real*8 r(npt), fint(nbt),asum,rbt(nbt),vbt(nbt),dr,y2(npt),u(npt)
      real*8 s(nswx), q(nswx), h(nswx,nswx), vloc(npt), bessel
      external bessel,radlin,interpol
c v(r) is really rv(r)
      call interpol(r,vloc,npt,rbt,vbt,nbt,y2,u)
      do n=1, nsw
         do m=1,n
            h(n,m) = 0.0
            do i=1, nbt
              fint(i) =  bessel(q(n)*rbt(i),l)*rbt(i)*
     *              bessel(q(m)*rbt(i),l)*vbt(i)
            end do
            call radlin(nbt,fint,dr,asum)
            h(n,m) = asum / s(n)/s(m)
         end do
         h(n,n) = h(n,n) + q(n)**2
      end do

      return
      end

*-----------------------------------------------------------------------
      subroutine find_root   (f,l,xt1,xt2,x0,eps,iret)
*     ====================
*-----------------------------------------------------------------------
      implicit none
      integer iret,l
      real*8 f, xt1, xt2, x1, x2, x0, eps, f1, f2, f0
      external f

      x1 = xt1
      x2 = xt2

      if(x1.gt.x2) then
         x0 = x1
         x1 = x2
         x2 = x0
      end if

      f1 = f(x1,l)
      f2 = f(x2,l)

      iret = 0

      if(sign(f1,f2).eq.f1) then
         iret = 1
         return
      end if

      do while(abs(x2-x1).gt.eps)
         x0 = 0.5*(x1+x2)
         f0 = f(x0,l)
         if(sign(f0,f1).eq.f0) then
            x1 = x0
            f1 = f0
         else
            x2 = x0
            f2 = f0
         end if
      end do

      return
      end

      subroutine error(routin,messag,ierr)
c-----------------------------------------------------------------------
c
      character*(*) routin, messag
      if(ierr.eq.0) return
      write(6,*) ' '
      write(6,'(1x,79(''!''))')
      write(6,'(5x,''from '',a,'' : error #'',i10)') routin,ierr
      write(6,'(5x,a)') messag
      write(6,'(1x,79(''!''))')
      if(ierr.gt.0) then
         write(6,'(''     stopping ...'')')
         stop
      else
         write(6,*) ' '
         return
      end if
      end


      subroutine radlin(mesh,func,dr,asum)
c-----------------------------------------------------------------------
c
c     simpson's rule integrator for function stored on the
c     radial logarithmic mesh
c
c.....logarithmic radial mesh information
c.....function to be integrated
      real*8 func(mesh),dr,asum
      integer mesh,i
c     
c.....variable for file = 0
c
c     routine assumes that mesh is an odd number so run check
      if ( mesh - ( mesh / 2 ) * 2 .ne. 1 ) then
        write(*,*) '***error in subroutine radlin'
        write(*,*) 'routine assumes mesh is odd but mesh =',mesh
        stop
      endif

      asum = func(1)+func(mesh)
      do  i = 2,mesh-1,2
         asum = asum + 4.0d0*func(i)+2.0d0*func(i+1)
      enddo
      asum = asum*dr/3.0d0
      return
      end

