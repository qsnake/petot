       subroutine w_line_usp(vwr_atom,zatom,iiatom,
     &  Ealpha,occ_t,isNL,Dij0,Qij,qijrad,
     &  icore,is_ref,
     &  ip_ref,id_ref,qi,wq,qi2,vq,
     &  rhoq,rhocq,vqT,ri,amr,lll,nbeta,rcut_q1t,
     &  rcut_q2t,qfuncLM0,r_at,a_r,b_r)
c----------------------------------------------------------------------------
*************************************************************************
*** Written by Lin-Wang Wang, 2001
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

      implicit double precision(a-h,o-z)
      include "param.escan_real"

      real*8 rhoc(2000),vlocT(2000),vloc(2000)
      real*8 ch_at(2000),vw(2000),qfuncLM0(12,2000)
      real*8 r_at(2000)

 
      real*8 qi2(mnq),vq(mnq),rhoq(mnq)
      real*8 rhocq(mnq)
      real*8 vqT(mnq)

cccccccc wq has 1/rmask for real space implementation
cccccccc    and has no 1/rmask for q space implementation

      real*8 qi(mnq),wq(mnq,8)
      real*8 ri(201),amr(201)
      real*8 Dij0(32,32),Qij(32,32),qijrad(0:6,32,32)

      integer  iiatom

      real*8  zatom,Ealpha
      real*8  occ_s,occ_p,occ_d,occ_t
      integer isNL(3),icore
      integer is_ref,ip_ref,id_ref
      integer is_TB,ip_TB,id_TB


      character*20  vwr_atom
c----------------------------------------------------------------------------
c----------------------------------------------------------------------------
c

cccc idim3 must be 8, otherwise, many of the other dimensions in other subroutines
cccc need to be changed !
      parameter( idim1 = 1000      , idim2 = 26         )
      parameter( idim3 = 8        , idim4 = 5          )
      parameter( idim5 = 4         , idim6 = 2*idim5-1  )
      parameter( idim7 = 4         , idim8 = 20         )

c     idim1  .ge.  no. of points in the radial mesh
c     idim2  .ge.  no. of shells in the atom
c     idim3  .ge.  no. of beta functions in vanderbilt potential
c     idim4  .ge.  no. of valence states for pseudization
c     idim5  .ge.  lmax + 1, lmax is maximum l value of potential
c     idim7  .ge.  no. of reference states for each angular momentum value
c     idim8  .ge.  no. of terms in taylor expansion of pseudized q_i,j



c
c.....radial mesh information
      dimension r(idim1),rab(idim1),rinner(idim6)
c.....charge densities
      dimension rsatom(idim1),rspsco(idim1)
c.....shell labels, energies, occupancies and real-space cutoff index
      dimension nnlz(idim2),ee(idim2),wwnl(idim2)
c.....pseudo quantum numbers, energies occupancies and cutoff radii
      dimension nnlzps(idim4),eeps(idim4),wwnlps(idim4),rc(idim5)
c.....beta and q functions and q pseudization coefficients
      dimension beta(idim1,idim3),qfunc(idim1,idim3,idim3),
     +qfcoef(idim8,idim6,idim3,idim3),vlocSR(idim1),vloc0(idim1)
c.....indexing for the beta functions
      dimension nbl0(idim5),nbl(idim5)
c.....angular momenta, reference energy, qij and dij of vanderbilt scheme
      dimension lll(idim3),eee(idim3),iptype(idim3),qqq(idim3,idim3),
     +ddd0(idim3,idim3),ddd(idim3,idim3)
c.....version number and date
      dimension iver(3),idmy(3)
c.....logic for logarithmic mesh types
      logical tlog
c.....wave functions
      dimension snl(idim1,idim2)
c
c.....title
      character*20 title,xctype
c.....file name array
      character*40 flname(6)
c
      integer stderr

      real*8 dum(idim1),qrad(0:6,idim3,idim3)
c      common /files/ stderr,input,iout,ioae,iplot,iologd,iops
c
c
c--------------------------------------------------------------------------
c
c       r e a d  p s e u d o p o t e n t i a l  f r o m  f i l e
c

        iops=10
        iout=6
        stderr=6
      

c       save data so compatibility checks can be run later
c        mold  = mesh
c        r2    = r(2)
c        rmesh = r(mesh)
c
        open( unit = iops , file = vwr_atom , status = 'old' ,
     +    form = 'unformatted' )
c
        read (iops) (iver(i),i=1,3),(idmy(i),i=1,3)
c
        if ( iver(1) .gt. 99 .or. iver(1) .lt. 1 .or. 
     +       iver(2) .gt. 9  .or. iver(2) .lt. 0 .or.
     +       iver(3) .gt. 9  .or. iver(3) .lt. 0       ) then
c
c         handle errors caused by pre-version number files
c
          write(stderr,*) '***error reading pseudopotential file'
          write(stderr,*) 'file does not contain version number'
          write(stderr,*) 'taking corrective action - will assume:'
c
          iver(1) = 2
          iver(2) = 1
          iver(3) = 0
c
          idmy(1) = 0
          idmy(2) = 0
          idmy(3) = 0
c
          write(stderr,90) (iver(i),i=1,3),idmy(2),idmy(1),idmy(3)
   90     format('pseudopotential program version',i3,'.',
     +      i1,'.',i1,'   date:',i3,' - ',i2,' - ',i4)
c
          rewind (iops)
c
        endif
c
c       set the mesh type
        if ( iver(1) .eq. 1 ) then
          tlog = .false.
        else
          tlog = .true.
        endif
c
c
        read (iops) title,zps,zvps,exftps,nvalps,mesh,etot
        read (iops) (nnlzps(i),wwnlps(i),eeps(i),i=1,nvalps)
        read (iops) keyps,ifpcor,rinner(1)
c
        if ( keyps .ne. 3 ) then
          write(iout,*) '***error in subroutine rwps'
          write(iout,*) 'keyps =',keyps,' out of programed range'
          call exit(1)
        endif
c
        if ( iver(1) .ge. 3 ) then
          read(iops) nang,lloc,eloc,ifqopt,nqf,qtryc
        endif
c
        if (10*iver(1)+iver(2).ge.51) then
          read (iops) (rinner(i),i=1,nang*2-1)
        else
          if (nang.gt.1) then
            do i=2,2*nang-1
              rinner(i)=rinner(1)
            end do
          endif
        endif
c
        if ( iver(1) .ge. 4 ) then
          read(iops) irelps
        endif
c
c
c       set the number of angular momentum terms in q_ij to read in
        if ( iver(1) .eq. 1 ) then
c         no distinction between nang and nvalps
          nang = nvalps
c         no optimisation of q_ij so 3 term taylor series
          nqf = 3
          nlc = 5
        elseif ( iver(1) .eq. 2 ) then
c         no distinction between nang and nvalps
          nang = nvalps
c         no optimisation of q_ij so 3 term taylor series
          nqf = 3
          nlc = 2 * nvalps - 1
        else
          nlc = 2 * nang - 1
        endif
c
c
        read (iops) (rc(i),i=1,nang)
        read (iops) nbeta,kkbeta
        do 100 j=1,nbeta
          read (iops) lll(j),eee(j),(beta(i,j),i=1,kkbeta)
          do 100 k=j,nbeta
            read (iops) ddd0(j,k),ddd(j,k),qqq(j,k),
     +      (qfunc(i,j,k),i=1,kkbeta),
     +      ((qfcoef(i,lp,j,k),i=1,nqf),lp=1,nlc)
            ddd0(k,j)=ddd0(j,k)
            ddd (k,j)=ddd (j,k)
            qqq (k,j)=qqq (j,k)
            do 100 i=1,kkbeta
              qfunc(i,k,j)=qfunc(i,j,k)
  100   continue
        do 110 i=1,nqf
        do 110 lp = 1,nlc
          qfcoef(i,lp,k,j)=qfcoef(i,lp,j,k)
  110   continue
c
        if (10*iver(1)+iver(2).ge.72) then
          read (iops) (iptype(j),j=1,nbeta),npf,ptryc
        endif
c
        read (iops) rcloc,(vloc0(i),i=1,mesh)
c
c       set index arrays nbl and nbl0
        do 150 lp = 1,nang
          nbl(lp) = 0
  150   continue
        do 160 ib = 1,nbeta
          lp = lll(ib) + 1
          nbl(lp) = nbl(lp) + 1
  160   continue
        lmaxps = lll(nbeta)
        nbl0(1) = 0
        do 170 i = 1,lmaxps
          nbl0(i+1) = nbl0(i) + nbl(i)
  170   continue
c
c       possible readin of the frozen core correction
        if (ifpcor.gt.0) then
           read(iops) rpcor
           read (iops) (rspsco(i),i=1,mesh)
        endif
c
        read (iops) (vlocSR(i),i=1,mesh)
        read (iops) (rsatom(i),i=1,mesh)
c
c       with the logarithmic mesh potentials grid information included
        if ( tlog ) then
          read (iops) (r(i),i=1,mesh)
          read (iops) (rab(i),i=1,mesh)
        endif
        if (iver(1) .ge. 6) then
c           nchi = nvales        ! warning, used before defined ! use nvalps ?
           nchi = nvalps        ! warning, used before defined ! use nvalps ?
           if (iver(1) .ge. 7)  read (iops) nchi
           read (iops) ((snl(i,j), i=1,mesh),j=ncores+1,ncores+nchi)
        endif
c
        close (iops)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   END OF THE READ IN
cccccccccccccccccccccccccccccccccccccccccccccc
ccc  definitions of the read-in functions:
ccc    do 100 i=1,mesh
ccc    r(i)=a*(exp(b*(i-1))-1)
ccc    rab(i)=b*(r(i)+a)
ccc    sqr(i)=exp(b*i/2)
ccc    enddo
ccc
ccc    zps: total atomic charge (atomic number)
ccc    zvps:  valence electron charge
ccc    nvalps: num of valence ele state (not counted spin, not the reference state)
ccc    nbeta:  num of reference states (not counted spin,i.e, just s,p,d,f, etc)
ccc    etot:   total atomic energy
ccc    mesh:   real space radial mesh r(i),i=1,mesh
ccc    kkbeta: end point for ref. state real space r: r(i),i=1,kkbeta
ccc    
ccc    nnlzps(i),i=1,nvalps:  (quantum num)(angular mom)(mz ? 0), eg, 310 for 3p, 400 for 4s
ccc    wwnlps(i),i=1,nvalps:  occupation for the valence state i for neutral atom
ccc    eeps(i), i=1,nvalps:   the eigen energy for valence state i
ccc    lll(j), j=1,nbeta:  the angular momentum for ref. state j (0,1,2,3 for s,p,d,f)
ccc    eee(j), j=1,nbeta:  the energy for ref. state j
ccc
ccc    local ionic potential:     vloc0(i)/r(i)          ! in Ryd
ccc    total screened potential:  vlocSR(i)/r(i)           ! in Ryd
ccc
ccc    atom charge:               rsatom(i)/r(i)**2/(4*pi)
ccc                            normaliz: \int rsatom(i)*dr = zvps
ccc                            This is the real charge of the valence electron
ccc    valence el. wave func:     snl(i,j)/r(i)/sqr(4*pi), j=1,nvalps
ccc                            normaliz: \int snl(i,j)**2 dr < 1.
ccc                            these wavefunctions are not norm conserved.
ccc                            sum_j snl(i,j)**2*wwwnlps(j) =< rsatom(i)
ccc    projector wave func:       beta(i,j)/r(i)/sqr(4*pi), j=1,nbeta, i=1,kkbeta (real space)
ccc
ccc    The nolocal potential is:
ccc    V_NL=sum_L,m sum_i,j  ddd(i,j) |beta(:,i)/r Y_Lm> <beta(:,j)/r Y_Lm|
ccc    note the overlap between L and i,j. Also, ddd(i,j) is zero if L(i).ne.L(j)
ccc    The overlap operator is:
ccc    S=1+sum_L,m sum_i,j qqq(i,j) |beta(:,i)/r Y_Lm> <beta(:.j)/r Y_Lm|
ccc    The orthnormal condition is:
ccc    <psi_i1| S | psi_i2> = delta_i1,i2
ccc    The Schroedinger's equationi s:
ccc    H psi_i1 = E_i1 S psi_i1
ccc    The charge density is:
ccc    rho(r) = sum_i1 |psi_i1(r)|^2+ sum_i,j  qfunc(r,i,j) X 
ccc        sum_i1 <beta(:,i)/r Y_Lm| psi_i1>  < psi_i1|beta(:,j)/r Y_Lm>  
ccc    ddd(j,k)=ddd0(j,k)+\int vlocSR(r)/r(i)*qfunc(i,j,k)/r(i)**2 * r^2 dr
ccc    qqq(j,k)=delta_j,k(same l)*\int qfunc(i,j,k)/r(i)**2 * r^2 dr
ccc    qijrad(L,j,k)=\int qfunc(i,j,k)/r(i)**2* r^(2+L) dr
ccc    qqq(j,k)=delta_j,k(same l)*qijrad(0,j,k)
cccccccccccccccccccccccccccccccccccc
         zatom=zvps
cccccccccccccccccccccccccccccccccccc
         rcut_q1t=r(kkbeta)*0.8
cccccccccccccccccccccccccccccccccccc
         qfuncLM0=0.d0
         r_at=0.d0
         LM0=0
         do j=1,nbeta
         do k=1,j
          if(lll(j).eq.lll(k)) then
          LM0=LM0+1
          if(LM0.gt.12) then
          write(6,*) "LM0.gt.12, stop, w_line"
          stop
          endif
          do i=2,kkbeta
          qfuncLM0(LM0,i)=qfunc(i,j,k)/r(i)**2
          enddo
          qfuncLM0(LM0,1)=qfuncLM0(LM0,2)
          endif
         enddo
         enddo
ccccccccccccccccccccccccccccccccccccccccccccc
cccc add a mask function to qfuncLM0. Note, the 
cccc magnetude is not important for qfuncLM0, it will not normalized, 
cccc but the shape is important. The original amplitude is recorded in qrad(ll,j,k)
ccccccccccccccccccccccccccccccccccccccccccccc
          LM0N=LM0
          b_r=(rab(300)-rab(100))/(r(300)-r(100))
          a_r=r(300)/(dexp(b_r*(300-1))-1.d0)
          call qfuncLM0_maskr()
ccccccccccccccccccccccccccccccccccccccccccccc


          do i=1,kkbetaL+1
          r_at(i)=r(i)
          enddo
cccccccccccccccccccccccccccccccccccc
         do j=1,nbeta
         do k=1,nbeta
         do ll=0,6
          dum=0.d0
          do i=2,kkbeta
          dum(i) = qfunc(i,j,k)*r(i)**ll
          enddo
          dum(1)=dum(2)
         sum=0.d0
         do i=1,kkbeta-1
         sum=sum+(dum(i)+dum(i+1))/2*(r(i+1)-r(i))
         enddo
         qrad(ll,j,k)=sum
         enddo
         enddo
         enddo
cccccccccccccccccccccccccccccccccccc

      pi=4*datan(1.d0)
      nrr=mesh
      icore=0
      if(ifpcor.gt.0) icore=1
      if(icore.eq.1) then
      do i=2,nrr
      rhoc(i)=rspsco(i)/r(i)**2/(4*pi)     ! warning, normalization
      enddo
      rhoc(1)=rhoc(2)
      else
      do i=1,nrr
      rhoc(i)=0.d0
      enddo
      endif
c-------------------------------------------
       do i=2,nrr
       vloc(i)=vloc0(i)/r(i)*0.5d0    ! change to Hartree unit
       enddo
       vloc(1)=vloc(2)

c-------------------------------------------
      pi=4*datan(1.d0)
      ch=zvps
      s=0.d0
      do i=2,nrr-1
      if(r(i).lt.15.d0) then
      s=s+(ch*r(i)+vloc(i)*r(i)**2)*(r(i+1)-r(i-1))/2
      endif
      enddo
      Ealpha=s*4*pi
c-------------------------------------------
c   vlocT is for the use of Thomas procedure
c   Since we don't have vs,vp,vd, etc, we just need to take vloc as vlocT
c-------------------------------------------
      do i=1,nrr
      vlocT(i)=vloc(i)
      enddo
c-------------------------------------------
cc   atomic charge. Note, this is the total charge
cc   if we need the non-norm conserved charge for the wavefunction occupation, 
cc   we need to use the (snl(i,j)/r(i))**2*wwnlps(i)
cc   ch_at-> rhoq is used in getVrhoL as the total charge. It is correct here. 
      sum_test=0.d0
      do i=2,nrr
      ch_at(i)=rsatom(i)/r(i)**2/(4*pi)

      if(i.le.nrr-1) then
      sum_test=sum_test+ch_at(i)*(r(i+1)-r(i-1))/2*
     &   r(i)**2*4*pi
      endif

      enddo
      ch_at(1)=ch_at(2)

      if(inode.eq.1) then
      write(6,*) "input atom charge for iiatom=",iiatom, sum_test
      endif

c-------------------------------------------
c      open(10,file="graph.vr")
c      rewind(10)
c      do i=2,kkbeta
c      vwr1=fac*beta(i,1)/r(i) 
c      vwr2=fac*beta(i,2)/r(i) 
c      vwr3=fac*beta(i,3)/r(i) 
c      vwr4=fac*beta(i,4)/r(i) 
c      vwr5=fac*beta(i,5)/r(i) 
c      vwr6=fac*beta(i,6)/r(i) 
c      write(10,1300) r(i),vwr1,vwr2,vwr3,vwr4,vwr5,vwr6
c      enddo
c1300  format(7(E12.5,1x))
c      close(10)
ccccc test



      qmx=dsqrt(2*Ecut2L)*1.2d0  
      
      occ_t=0.d0
      do j=1,nvalps
      occ_t=occ_t+wwnlps(j)
      enddo

      z=zvps
      rst=10.d0
      do 101 iq=1,mnq
      g1=(iq-1)*qmx/(mnq-1.d0)
      if(g1.lt.1.D-3) g1=1.D-3

      s=0.d0
      s1=0.d0
      s2=0.d0
      s3=0.d0

      do i=1,nrr-1
      if(r(i).gt.rst) then
      rst1=r(i)
      goto 97
      endif
      r1=r(i)
      r2=r(i+1)
      x1=g1*r1
      x2=g1*r2

      b1=vloc(i)*r1
      b2=(vloc(i+1)*r2-b1)/(x2-x1)

      c1=(dcos(x1)-dcos(x2))*b1+(x1-x2)*b2*dcos(x2)+
     &  b2*(dsin(x2)-dsin(x1))
      s=s+c1

      b11=ch_at(i)*r1
      b22=(ch_at(i+1)*r2-b11)/(x2-x1)

      c11=(dcos(x1)-dcos(x2))*b11+(x1-x2)*b22*dcos(x2)+
     &  b22*(dsin(x2)-dsin(x1))
      s1=s1+c11

      b1=vlocT(i)*r1
      b2=(vlocT(i+1)*r2-b1)/(x2-x1)

      c1=(dcos(x1)-dcos(x2))*b1+(x1-x2)*b2*dcos(x2)+
     &  b2*(dsin(x2)-dsin(x1))
      s2=s2+c1

      b1=rhoc(i)*r1
      b2=(rhoc(i+1)*r2-b1)/(x2-x1)
      c1=(dcos(x1)-dcos(x2))*b1+(x1-x2)*b2*dcos(x2)+
     &  b2*(dsin(x2)-dsin(x1))
      s3=s3+c1

      enddo
97    continue

      s2=s2-s

      s=s/g1**2
      s=s*4*pi
      s=s-z*4*pi*dcos(g1*rst1)/g1**2

      s1=s1/g1**2
      s1=s1*4*pi

      s2=s2/g1**2
      s2=s2*4*pi

      s3=s3/g1**2
      s3=s3*4*pi

      qi2(iq)=g1
      vq(iq)=s
      rhoq(iq)=s1
      vqT(iq)=s2
      rhocq(iq)=s3
101   continue

c-------------------------------------------
c-------------------------------------------
ccccccc the projector vw(i) is normalized to 1 for all ibeta 

      rst=rcut/1.03
      fac=1.d0/dsqrt(4*pi)
      do 1000 ibeta=1,nbeta


      do i=2,kkbeta
      if(r(i).gt.rst) then
      amrI=0.d0
      else
      r2=r(i)/rcut
      ir=1+r2*200.d0
      f1=(ri(ir+1)-r2)/(ri(ir+1)-ri(ir))
      f2=(r2-ri(ir))/(ri(ir+1)-ri(ir))
      amrI=amr(ir)*f1+amr(ir+1)*f2
      amrI=1.d0/amrI
      endif
      vw(i)=fac*beta(i,ibeta)/r(i)    

      vw(i)=vw(i)*amrI
      enddo
      vw(1)=vw(2)


ccccccccccccccccccc   test
c       sum=0.d0
c       do i=2,kkbeta-1
c       sum=sum+vw(i)**2*(r(i+1)-r(i-1))/2*r(i)**2
c       enddo
c       sum=sum*4*pi
c       write(6,*) "ibeta,sum=",ibeta, sum
ccccccccccccccccccc   test
      

      if(lll(ibeta).eq.0) then
      do 202 iq=1,mnq
      g1=(iq-1)*qmx/(mnq-1.d0)

      if(g1.lt.1.D-3) g1=1.D-3

      qi(iq)=g1

      s=0.d0
      do i=2,kkbeta-1
      if(r(i).gt.rst) goto 96
      x=r(i)*g1
      if(x.gt.0.01d0) then
      aj0=dsin(x)/x
      else
      aj0=1.d0-x**2/6.d0+x**4/120.d0-x**6/5040.d0
      endif
      s=s+aj0*r(i)**2*vw(i)*(r(i+1)-r(i-1))/2
      enddo
96    continue
      s=s*4*pi
      qi(iq)=g1
      wq(iq,ibeta)=s
202   enddo
      endif     ! lll(ibeta).eq.0

c--------------------------------------------

      if(lll(ibeta).eq.1) then
      do 200 iq=1,mnq
      g1=(iq-1)*qmx/(mnq-1.d0)

      if(g1.lt.1.D-3) g1=1.D-3

      qi(iq)=g1

      s=0.d0
      do i=2,kkbeta-1
      if(r(i).gt.rst) goto 99
      x=r(i)*g1
      if(x.gt.0.1) then
      aj1=dsin(x)/x**2-dcos(x)/x
      else
      aj1=x/3.d0-x**3/30.d0+x**5/840.d0-x**7/45360.d0
      endif
      s=s+aj1*r(i)**2*vw(i)*(r(i+1)-r(i-1))/2
      enddo
99    continue
      s=s*4*pi
      qi(iq)=g1
      wq(iq,ibeta)=s
200   enddo
      endif    ! lll(ibeta).eq.1

c--------------------------------------------
      if(lll(ibeta).eq.2) then
      do 201 iq=1,mnq
      g1=(iq-1)*qmx/(mnq-1.d0)

      if(g1.lt.1.D-3) g1=1.D-3

      qi(iq)=g1

      s=0.d0
      do i=2,kkbeta-1
      if(r(i).gt.rst) goto 98
      x=r(i)*g1

      if(x.gt.0.2d0) then
      aj2=(3/x**3-1/x)*dsin(x)-3*dcos(x)/x**2
      else
      aj2=x**2/15.d0-x**4/210.d0+x**6/7560.d0-x**8/498960.d0
      endif

      s=s+aj2*r(i)**2*vw(i)*(r(i+1)-r(i-1))/2
      enddo
98    continue
      s=s*4*pi
      qi(iq)=g1
      wq(iq,ibeta)=s
201   enddo
      endif     ! lll(ibeta).eq.2

      if(lll(ibeta).gt.2) then
      write(6,*) "lll(ibeta).gt.2, not programed yet, stop",lll(ibeta)
      stop
      endif

1000  continue

c      open(10,file="graph.vq")
c      rewind(10)
c      do i=1,mnq
c      write(10,1300) qi(i),(wq(i,ibeta),ibeta=1,nbeta)
c      enddo
c1300  format(7(E12.5,1x))
c       do ibeta=1,nbeta
c       sum=0.d0
c       do i=2,mnq-1
c       sum=sum+wq(i,ibeta)**2*(qi(i+1)-qi(i-1))/2*qi(i)**2
c       enddo
c       sum=sum*4*pi/(2*pi)**3
c       write(6,*) "ibeta, sumwq=",ibeta,sum
c       enddo
c      close(10)
ccccc test
ccccccccccccccccccccccccccccccccccccccccccccccc
ccc  ddd(ibeta1,ibeta2), ibeta only counts l, not m
ccc  Dij0(iref1,iref2),   iref counts both l and m
ccc  same for qqq and Qij

      Dij0=0.d0
      Qij=0.d0
      iref1=0
      do ibeta1=1,nbeta

      iref2=0
      do ibeta2=1,nbeta
      if(lll(ibeta1).eq.lll(ibeta2)) then
        do ll=1,2*lll(ibeta1)+1
        Dij0(iref1+ll,iref2+ll)=0.5d0*ddd0(ibeta1,ibeta2)    ! change into Hartree unit
cc         Dij0(iref1+ll,iref2+ll)=0.5d0*(ddd(ibeta1,ibeta2)
cc     &       -ddd0(ibeta1,ibeta2))    ! change into Hartree unit
cccc        Qij(iref1+ll,iref2+ll)=qqq(ibeta1,ibeta2)
        Qij(iref1+ll,iref2+ll)=qrad(0,ibeta1,ibeta2)     ! num consistent with rho_nL
        enddo

      endif
      iref2=iref2+2*lll(ibeta2)+1
      enddo
      iref1=iref1+2*lll(ibeta1)+1
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccc

      iref2=0
      do ibeta2=1,nbeta
      do ll2=1,2*lll(ibeta2)+1
      iref2=iref2+1
      iref1=0
      do ibeta1=1,nbeta
      do ll1=1,2*lll(ibeta1)+1
      iref1=iref1+1
      do jj=0,6
      qijrad(jj,iref1,iref2)=qrad(jj,ibeta1,ibeta2)
      enddo
      enddo
      enddo
      enddo
      enddo

cccccccc  test
c      if(inode.eq.1) then
c      write(6,*) "******* Qij *********"
c      do j=1,10
c      write(6,333) (Qij(i,j),i=1,10)
c      enddo
c      write(6,*) "******* Qij *********"

c      write(6,*) "******* Dij_change *********"
c      do j=1,10
c      write(6,333) (Dij0(i,j),i=1,10)
c      enddo
c      write(6,*) "******* Dij_change *********"
c333   format(10(E11.5,1x))
c      endif


      return
c
c----------------------------------------------------------------------------

      contains
      

       subroutine qfuncLM0_maskr()
ccccccccccccccccccccccccccccccccccccccccccccc
cccc add a mask function to qfuncLM0. Note, the 
cccc magnetude is not important for qfuncLM0, it will not normalized, 
cccc but the shape is important. The original amplitude is recorded in qrad(ll,j,k)
ccccccccccccccccccccccccccccccccccccccccccccc
       implicit double precision (a-h,o-z)
       real*8 vr_tmp(2000),vq_tmp(mnq),qi_tmp(mnq)

          pi=4*atan(1.d0)
          rcut_q2t=r(kkbeta)*1.2
ccccccc  This cut off is important. If one uses 1.5, we can get better results
ccccccc  But that might be too costly. 1.2 is a compromise
ccccccc

          kkbetaL=log(rcut_q2t/a_r+1.d0)/b_r+1.00001d0
          rcut_q2t=r(kkbetaL)
          LM0N=LM0
          qmx_tmp=dsqrt(2*Ecut2L)

cccccccccccccccccccccccccc
          if(inode.eq.1) then
          open(10, file="qfunLM0.orig")
          rewind(10)
          write(6,*) "LM0N=",LM0N
          do i=1,kkbetaL
          write(10,400) r(i),qfuncLM0(1,i),qfuncLM0(2,i),
     &  qfuncLM0(3,i),qfuncLM0(4,i),qfuncLM0(5,i)
          enddo
          close(10)
400       format(6(E11.5,1x))
          endif
cccccccccccccccccccccccccccccccccccccccccc

          do 200 LM0=1,LM0N
          sum=0.d0
          do i=1,kkbeta
          sum=sum+dabs(qfuncLM0(LM0,i))
          enddo
          if(sum.lt.1.D-10) goto 200

          do i=1,kkbetaL
          vr_tmp(i)=qfuncLM0(LM0,i)
          enddo

cccccccccccccccccccccccccccccccccccccc 
          do 202 iq=1,mnq
          g1=(iq-1)*qmx_tmp/(mnq-1.d0)

          s=0.d0
          do i=2,kkbeta-1
          x=r(i)*g1
          if(x.gt.0.01d0) then
          aj0=dsin(x)/x
          else
          aj0=1.d0-x**2/6.d0+x**4/120.d0-x**6/5040.d0
          endif
          s=s+aj0*r(i)**2*vr_tmp(i)*(r(i+1)-r(i-1))/2
          enddo

          s=s*4*pi
          qi_tmp(iq)=g1

ccccccccccccc  different mask function
c          y=g1/qmx_tmp*pi/2
c          amask=1.d0-dsin(y)**20
ccccccccccccccccccccccccccc
          amask=1.d0
          if(g1.gt.qmx_tmp*0.8) amask=0.d0
cccccccccccccccccccccccc
ccccc those few masks (inluding amask=1) gives very similar results in practice. 
ccccc here we will use the 0.8 cutoff mask, since it has some possible theoretical meaning

          vq_tmp(iq)=s*amask

202       enddo
cccccccccccccccccccccccccccccccccccccc 
          if(LM0.eq.1.and.inode.eq.1) then
          open(10,file="qfuncLM0.q")
          do iq=1,mnq
          write(10,400) qi_tmp(iq),vq_tmp(iq),vq_tmp(iq)*qi_tmp(iq)**2
          enddo
          close(10)
          endif

          do 203 i=1,kkbetaL

          s=0.d0
          do iq=1,mnq

          x=qi_tmp(iq)*r(i)
          if(x.gt.0.01d0) then
          aj0=dsin(x)/x
          else
          aj0=1.d0-x**2/6.d0+x**4/120.d0-x**6/5040.d0
          endif
          s=s+aj0*qi_tmp(iq)**2*vq_tmp(iq)*
     &          (qi_tmp(iq+1)-qi_tmp(iq-1))/2

          enddo
          
          qfuncLM0(LM0,i)=s/(2*pi**2)

          y=r(i)/rcut_q2t*pi/2
          amask=1.d0-dsin(y)**10
          qfuncLM0(LM0,i)=qfuncLM0(LM0,i)*amask

203       continue

200       continue      ! diff LM0

ccccccccccccccccccccccccccccccccccccccccccccc
          if(inode.eq.1) then
          open(10, file="qfunLM0.new")
          rewind(10)
          do i=1,kkbetaL
          write(10,400) r(i),qfuncLM0(1,i),qfuncLM0(2,i),
     &  qfuncLM0(3,i),qfuncLM0(4,i),qfuncLM0(5,i)
          enddo
          close(10)
          endif

          return
          end subroutine qfuncLM0_maskr
cccccccccccccccccccccccccccccccccccccccccccccc

          end

          
