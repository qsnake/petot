      subroutine vionic(ispp,itype,icorr,ifcore,zsh,rsh,
     1 lmax,nr,a,b,r,rab,nameat,ncore,znuc,
     2 cdd,cdu,cdc,viod,viou)
c
c  Vionic sets up the ionic potential.
c  Note that viod/viou is the ionic potential times r.
c
c  njtj ***  major modifications  ***
c    If a potential does not exist, it is approximated
c    by an existing potential.
c    A nonspin or spin-polarized pseudo test, uses the 
c    down(nonspin generation), weighted average(spin-
c    polarized), or averaged(relativistic) potentials.
c    A relativistic pseudo test, must use relativistic
c    generated potentials.  The Schroedinger equation is
c    used to integrate a relativistic pseudo test,
c    not the Dirac equation. 
c  njtj  ***  major modifications  ***
c
c  njtj
c  ###  Cray conversions  
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
c  jlm  version 5.60
c
      implicit double precision (a-h,o-z)
c
      parameter (zero=0.D0)
Cray      parameter (zero=0.0)
c
      character*1 ispp
      character*2 icorr,icorrt,nameat,namet
      character*3 irel
      character*4 nicore
      character*10 iray(6),ititle(7)

      dimension r(nr),rab(nr),cdd(nr),cdu(nr),cdc(nr),
     1 viod(lmax,nr),viou(lmax,nr),npd(5),npu(5)
c
c  2*znuc part
c
      ifcore = 0
      if (itype .lt. 4) then
        do 10 i=1,lmax
          do 12 j=1,nr
            viod(i,j) = -2*znuc
            viou(i,j) = -2*znuc
 12       continue
 10     continue 
      else
c
c  read pseudopotentials from tape1
c
        rewind 1
        read(1) namet,icorrt,irel,nicore,(iray(i),i=1,6),
     1   (ititle(i),i=1,7),npotd,npotu,nrm,a,b,zion
        if(nicore.eq.'fcec'.or.nicore.eq.'pcec') ifcore = 1
        if(nicore.eq.'fche'.or.nicore.eq.'pche') ifcore = 2
        nr = nrm+1
        read(1) (r(i),i=2,nr)
        r(1) = zero             
c
c   down potentials (or average relativistic potentials)
c
c njtj  ***  major start  ***
c   if a potential does not exist, it is replaced by the
c   next existing lower angular momentum potential or 
c   the next existing higher if no lower exist.
c        
        do 15 i=1,lmax
          npd(i)=0
 15     continue
        do 20 i=1,npotd
          read(1) loi,(viod(loi+1,j),j=2,nr)
          viod(loi+1,1) = zero
          npd(loi+1)=1
 20     continue
        if (npd(1) .eq. 0) then
          do 25 i=2,lmax
            if (npd(i) .gt. 0) then
              do 24 j=1,nr
                viod(1,j)=viod(i,j)
 24           continue
              goto 30
            endif
 25       continue
        endif                   
 30     do 33 i=2,lmax
          if (npd(i) .eq. 0) then
            do 32 j=1,nr
              viod(i,j)=viod(i-1,j)
 32         continue
          endif
 33     continue    
c
c   up potentials (or spin orbit potentials)
c
        if (npotu .le. 0) goto 49
        do 35 i=1,lmax
          npu(i)=0
 35     continue
        do 37 i=1,npotu
          read(1) loi,(viou(loi+1,j),j=2,nr)
          viou(loi+1,1) = zero
          npu(loi+1)=1
 37     continue
        if (npu(1) .eq. 0) then
          do 38 i=2,lmax
            if (npu(i) .gt. 0) then
              do 39 j=1,nr
                viou(1,j)=viou(i,j)
 39           continue
              goto 40
            endif
 38       continue
        endif                   
 40     do 45 i=2,lmax
          if (npu(i) .eq. 0) then
            do 43 j=1,nr
              viou(i,j)=viou(i-1,j)
 43         continue
          endif
 45     continue    
c
c  njtj  ***  major end  ***
c
c
c  core and valence charges
c
 49     read(1) (cdc(i),i=2,nr)
        cdc(1) = zero
c
c  replace valence charge on tape(valence charge modify)
c    
        if (itype .eq. 6) then
          write(1) (cdd(i)+cdu(i),i=2,nr)
          return
        endif
        read(1) (cdd(i),i=2,nr)
        cdd(1) = zero
c
c  njtj  ***   major start  ***
c   distribute charge as up and down charge
c   generate radial intergration grid
c   set up potentials equal to down potentials for 
c   spin-polarized pseudo test of nonspin and relativistic
c   generated potentails.  Construct spin-orbit potentials
c   from relativistic sum and difference potentials and
c   change ispp='r' to ispp=' '.
c
        do 50 i=1,nr
          rab(i) = (r(i)+a)*b
          cdd(i) = cdd(i)/2
          cdu(i) = cdd(i)
 50     continue
        if (ispp .eq. 's' .and. irel .ne. 'isp') then
          do 51 i=1,lmax
            do 52 j=1,nr
              viou(i,j) = viod(i,j)
 52         continue
 51       continue
        endif
        if (ispp .eq. 'r') then
          ispp=' '
          if (irel .ne. 'rel') then
            write(6,130)irel
 130  format(//,'Pseudopotentail is not relativistic!!!!',/
     1 ' setting up potentials equal to down!!!',//) 
            do 53 i=1,lmax
              do 54 j=1,nr
                viou(i,j) = viod(i,j)
 54           continue
 53         continue
          else
            do 57 j=1,nr
              viou(1,j)=viod(1,j)
 57         continue
            do 58 i=2,lmax
              do 56 j=1,nr
                vsum=viod(i,j)
                vdiff=viou(i,j)
                viod(i,j)=vsum-i*vdiff/2
                viou(i,j)=vsum+(i-1)*vdiff/2
 56           continue
 58         continue
          endif
        endif
c
c   njtj  ***  major end   ***
c
c
c   printout
c
        write(6,60) namet,icorrt,irel,nicore,(iray(i),i=1,6),
     1   (ititle(i),i=1,7)
 60   format(//,1x,a2,2x,a2,2x,a3,2x,a4,
     1 '  pseudopotential read from tape',
     2 /,1x,2a10,5x,4a10,/,1x,7a10,//)
        if (nameat .ne. namet) write(6,70) nameat,namet 
 70   format(' input element ',a2,
     1 ' not equal to element on tape ',a2,//)
        if (icorr .ne. icorrt) write(6,80) icorr,icorrt
 80   format(' input correlation ',a2,
     1 ' not equal to correlation from tape ',a2,//)
        write(6,90) r(2),nr,r(nr)
cjlm
 90   format(' radial grid parameters',//,
     1 ' r(1) = .0 , r(2) =',e9.3,' , ... , r(',i4,') =',
     2 f8.3,//)
cjlm
      endif
c
c   add potential from shell charge
c
      if (abs(zsh) .gt. 0.e-5) then
        do 110 i=1,lmax
          do 120 j=1,nr
            if (r(j) .ge. rsh) then
              viod(i,j) = viod(i,j) - 2*zsh
              viou(i,j) = viou(i,j) - 2*zsh
            else
              viod(i,j) = viod(i,j) - 2*zsh*r(i)/rsh
              viou(i,j) = viou(i,j) - 2*zsh*r(i)/rsh 
            endif
 120      continue
 110    continue
       endif
       return
       end
