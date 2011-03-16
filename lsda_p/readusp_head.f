       subroutine readusp_head(filename,iiatom,nref)

c      subroutine rwps(iread,title,zps,zvps,exftps,nvalps,irel,
c     +  z,zv,exfact,nvales,mesh,etot,nnlzps,wwnlps,eeps,
c     +  nnlz,wwnl,ee,ncores,keyps,ifpcor,rpcor,rinner,rc,nbeta,
c     +  kkbeta,lll,eee,iptype,npf,ptryc,
c     +  beta,ddd0,ddd,qqq,qfunc,qfcoef,
c     +  rcloc,vloc0,vloc,snl,lloc,eloc,rsatom,rspsco,flname,
c     +  ifqopt,nqf,qtryc,nang,nbl0,nbl,r,rab,iver,idmy,
c     +  idim1,idim2,idim3,idim4,idim5,idim6,idim8)

c
c----------------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)

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
     +qfcoef(idim8,idim6,idim3,idim3),vloc(idim1),vloc0(idim1)
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
      character*20 filename
c
      integer stderr
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
        open( unit = iops , file = filename , status = 'old' ,
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
        iiatom=zps*1.001
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
        
        nref=0
        do 100 j=1,nbeta
          read (iops) lll(j),eee(j),(beta(i,j),i=1,kkbeta)
        nref=nref+lll(j)*2+1
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
        
        return
        end

