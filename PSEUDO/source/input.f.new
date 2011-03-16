      subroutine input(itype,ikerk,icorr,ispp,zsh,rsh,
     1 nr,a,b,r,rab,nameat,norb,ncore,no,lo,so,zo,
     2 znuc,zel,evi,nval_orig)
c
c  subroutine to read input parameters
c
c  njtj ***  modifications  ***
c    The input and output variables passed have been changed.
c    There are five new pseudopotential generation options
c    The input variables znuc,zsh,rsh,rmax,aa,bb are 
c    compared to a small positive value - eliminates
c    floating point comparisions errors(zero is 
c    not always zero).
c  njtj ***  modifications  ***
c
c  lcb
c    modified for GGA by Carlos Balbas,   January 97.
c  lcb
c
c  jlm   version 5.61
c
c  njtj
c  ###  Cray conversions  
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      parameter (one=1.D0,zero=0.D0,pfive=0.5D0)
Cray      parameter (one=1.0,zero=0.0,pfive=0.5)

      character*1 ispp
      character*2 type,icorr,nameat
      character*3 name,kerker
      character*10 iray(5)
      character*10 ititle
c
c  dimension of transfered data
c
      dimension r(nr),rab(nr),no(norb),lo(norb),so(norb),
     1 zo(norb),evi(norb)
c
c  dimensions of data used in routine
c
      dimension nc(15),lc(15),nomin(5)
c
c  data for orbitals, 1s,2s,2p,3s,3p,3d,4s,4p,4d,5s,5p,4f,5d,6s,6p
c
      data nc /1,2,2,3,3,3,4,4,4,5,5,4,5,6,6/
      data lc /0,0,1,0,1,2,0,1,2,0,1,3,2,0,1/ 
c
      do 5 i=1,5
        nomin(i)=10
 5    continue
      do 6 i=1,norb
        no(i)=0
        lo(i)=0
        so(i)=zero
        zo(i)=zero
        evi(i)=zero
 6    continue
c
c  read the type of calculation and title card
c   itype =
c   ae = 0 all electron calculation
c   pg = 1 pseudopotential generation w/o core correction
c   pe = 2 pseudopotential generation w/  core correction exchange
c   ph = 3 pseudopotential generation w/  core correction hartree/exc
c   pt = 4 pseudopotential test
c   pm = 5 pseudopotential test + valence charge modify
c

      read(5,*) type,ititle


 10   format(3x,a2,5a10)
c
c  if type = ' ' , no more data, program ends
c
      if (type .eq. 'ae') then
        itype=0
      elseif (type .eq. 'pg') then
        itype=1
      elseif (type .eq. 'pe') then
        itype=2
      elseif (type .eq. 'ph') then
        itype=3
      elseif (type .eq. 'pt') then
        itype=4
      elseif (type .eq. 'pm') then
        itype=5
      else
        itype=-1
        return
      endif
c
c  njtj  ***  major modification  start  ***
c  There are seven ways to generate the pseudopotential :
c    kerker = van Vanderbilt
c    kerker = tam Troullier and Martins
c    kerker = ker (yes) Kerker
c    kerker = hsc (no)  Hamann Schluter and Chiang
c    kerker = min (oth) datafile made for minimization
c    kerker = bhs Bachelet, Hamann and Schluter
c    kerker = tm2 Improved Troullier and Martins
c
      if (itype.gt.0) then
        read(5,*)kerker
   11 format(8x,a3)
        if(kerker .eq. 'tm2' .or. kerker .eq. 'TM2') then
          ikerk = 6
        elseif(kerker .eq. 'bhs' .or. kerker .eq. 'BHS') then
          ikerk = 5 
        elseif(kerker .eq. 'oth' .or. kerker .eq. 'OTH' .or.
     1   kerker .eq. 'min' .or. kerker .eq. 'MIN') then
          ikerk = 4
        elseif (kerker .eq. 'van' .or. kerker .eq.'VAN') then
          ikerk = 3
        elseif (kerker .eq. 'tbk' .or. kerker .eq. 'TBK'
     1   .or. kerker .eq. 'tam' .or. kerker .eq. 'TAM') then
          ikerk = 2
        elseif (kerker .eq. 'yes' .or. kerker .eq. 'YES' .or.
     1   kerker .eq. 'ker' .or. kerker .eq. 'KER') then
          ikerk = 1
        elseif (kerker .eq. 'no ' .or. kerker .eq. ' no' .or.
     1   kerker .eq. 'NO ' .or. kerker .eq. ' NO' .or. kerker
     2   .eq. 'hsc' .or. kerker .eq. 'HSC') then
          ikerk = 0
        else
          write(6,1000)kerker
          call ext(150)
        endif
      endif                              
 1000 format(//,'error in input - kerker =',a3,' unknown')
c  njtj  ***  major modification end  ***
c
c   read element name and correlation type
c   ispp = ' ' - nonspin calculation
c   ispp = s  - spin polarized calculation
c   ispp = r  - relativistic calculation
c
      read(5,*) nameat,icorr,ispp
 15   format(3x,a2,3x,a2,a1)

      if (ispp .ne. 's' .and. ispp .ne. 'r') ispp=' '
      if (ispp .eq. 's' .and. icorr .eq. 'xa') ispp=' '
      if (ispp .eq. 's' .and. icorr .eq. 'wi') ispp=' '
      if (ispp .eq. 's' .and. icorr .eq. 'hl') ispp=' '
c
c  njtj   ***  major modification start  ***
c   Floating point comparison error modification.
c   Read the atomic number (nuclear charge),
c   shell charge and radius (added to the nuclear potential),
c   and radial grid parameters.
c
      read(5,*) znuc,zsh,rsh,rmax,aa,bb
 20   format(6f10.3)

      if (abs(znuc) .le. 0.00001) znuc=charge(nameat)
      if (itype .lt. 4) then
c
c   set up grid
c
        if (abs(rmax) .lt. 0.00001) rmax=120.0
        if (abs(aa) .lt. 0.00001) aa = 6.0  
        if (abs(bb) .lt. 0.00001) bb = 80.0 
c       bb = 40.0 standard
c       bb = 80.0 2 * standard
c       bb = 120.0 3 * standard
        a = exp(-aa)/znuc
        b = 1/bb
        do 30 i=1,nr
          if (i .eq. nr) then
            write(6,50)
            call ext(100)
          endif
          r(i) = a*(exp(b*(i-1))-1)
          rab(i) = (r(i)+a)*b
          if (r(i) .gt. rmax) goto 60
 30     continue
 60     nr = i-1
      endif
 50   format(/,' error in input - arraylimits',
     1 ' for radial array exceeded',/)          
c  njtj  ***  major modification end  ***
c
c   read the number of core and valence orbitals
c

      read(5,*) ncore,nval
      nval_orig=nval
 70   format(2i5,4f10.3)
      if (ncore .gt. 15) then
        write(6,1010)
        call ext(101)
      endif          
 1010 format(//,'error in input - max number of core orbitals',
     1 'is 15')
c
c   compute occupation numbers and orbital energies for the core
c  
      zcore = zero
      if (ncore .eq. 0) goto 85
      sc = zero
      if (ispp .ne. ' ') sc=-pfive
      norb = 0
      do 80 i=1,ncore
        do 80 j=1,2
          if (ispp .eq. ' ' .and. j .eq. 2) goto 80
          norb = norb + 1
          no(norb) = nc(i)
          lo(norb) = lc(i)
          so(norb) = sc
          zo(norb) = 2*lo(norb)+1
          if (ispp .eq. ' ') zo(norb) = 2*zo(norb)
          if (ispp .eq. 'r') zo(norb) = 2*(lo(norb)+sc)+1
          zcore = zcore + zo(norb)
          if (abs(zo(norb)) .lt. 0.1) norb=norb-1
          if (ispp .ne. ' ') sc=-sc
 80   continue
      ncore = norb
c
c   for the valence orbitals
c
 85   if (itype .ge. 4) ncore =0
      norb = ncore
      zval = zero
      if (nval .eq. 0) goto 105
      do 90 i=1,nval
ccc        read(5,*) ni,li,zd,zu,evd
        read(5,*) ni,li,zd,zu
        evd=0.d0           ! changed by LW

        si = zero
        if (ispp .ne. ' ') si=pfive
        do 90 j=1,2
          if (ispp .eq. ' ' .and. j .eq. 2) goto 90
          norb = norb + 1
          if (ispp .ne. ' ') si=-si
          no(norb) = ni
          lo(norb) = li
          so(norb) = si
          zo(norb) = zd+zu 
          if (zo(norb) .eq. zero) evi(norb)=evd
          if (ispp .eq. 's') then
            if (si .lt. 0.1) then
              zo(norb) = zd
            else
              zo(norb) = zu
            endif
          elseif (ispp .eq. 'r') then
            zo(norb)=zo(norb)*(2*(li+si)+1)/(4*li+2)
          endif
          zval = zval + zo(norb)
          if (ispp .eq. 'r' .and. li+si .lt. zero) norb=norb-1
          if (norb .eq. 0) goto 90
          if (nomin(lo(norb)+1) .gt. no(norb))
     1     nomin(lo(norb)+1)=no(norb)
 90   continue
c
c   abort if two orbitals are equal
c
      nval = norb - ncore
      do 100 i=1,norb
        do 100 j=1,norb
          if (i .le. j) goto 100
          if (no(i) .ne. no(j)) goto 100
          if (lo(i) .ne. lo(j)) goto 100
          if (abs(so(i)-so(j)) .gt. 0.001) goto 100
          write(6,1020)i
          call ext(110+i)
 100  continue          
 1020 format(//,'error in input - orbital ',i2,
     1 'is already occupied')
c
c   reduce n quantum number if pseudoatom
c
      if (itype .ge. 4) then
        do 103 i=1,nval
          no(i) = no(i)-nomin(lo(i)+1)+lo(i)+1
 103    continue
      endif
 105  zion = znuc - zcore - zval
      zel = zval
      if (itype .lt. 4) then
        zel=zel+zcore
      else
        znuc=znuc-zcore
      endif
c
c   find jobname and date and printout, zedate is a machine dependent
c   routine
c
      if (icorr.eq.'pb') then
         iray(1)='atom-GGA96'
      elseif (icorr.eq.'pw') then
         iray(1)='atom-LDApw' 
      elseif (icorr.eq.'ca') then
         iray(1)='atom-LDAca'
      else
         stop 'unrecognized correlation choice'
      endif
      call zedate(iray(2))
c
c   printout
c
      write(6,110) iray(1),iray(2),ititle
 110  format(1x,a10,a10,5x,5a10,/,21('*'),/)
      if (itype .eq. 0) then
        write(6,120) nameat
      elseif (itype .lt. 4) then
        write(6,121) nameat
      elseif (itype .eq. 4) then
        write(6,124) nameat
      elseif (itype .eq. 5) then
        write(6,125) nameat
      endif
 120  format(1x,a2,' all electron calculation ',/,1x,27('-'),/)
 121  format(1x,a2,' pseudopotential generation',/,1x,29('-'),/)
 124  format(1x,a2,' pseudopotential test',/,1x,23('-'),/)
 125  format(1x,a2,' pseudo test + charge mod ',/,1x,27('-'),/)
      if (ispp .eq. 'r') then
        write(6,150)
 150  format(' r e l a t i v i s t i c ! !',/)
        name = '   ' 
      elseif (ispp .eq. ' ') then
        name = 'non'
      else
        name = '   '
      endif
      write(6,160) icorr,name
 160  format(' correlation = ',a2,3x,a3,'spin-polarized',/)
      write(6,170) znuc,ncore,nval,zel,zion
 170  format(' nuclear charge             =',f10.6,/,
     1       ' number of core orbitals    =',i3,/,
     2       ' number of valence orbitals =',i3,/,
     3       ' electronic charge          =',f10.6,/,
     4       ' ionic charge               =',f10.6,//)
      if (zsh .gt. 0.00001) write(6,175) zsh,rsh
 175  format(' shell charge =',f6.2,' at radius =',f6.2,//)
      write(6,180)
 180  format(' input data for orbitals',//,
     1 '  i    n    l    s     j     occ',/)
      xji = zero
      do 200 i=1,norb
        if (ispp .eq. 'r') xji = lo(i) + so(i)
        write(6,190) i,no(i),lo(i),so(i),xji,zo(i)
 190  format(1x,i2,2i5,2f6.1,f10.4)
 200  continue
      if (itype .lt. 4) write(6,210) r(2),nr,r(nr),aa,bb
 210  format(//,' radial grid parameters',//,
     1 ' r(1) = .0 , r(2) =',e9.3,' , ... , r(',i4,') =',f8.3,
     2 /,' a =',f7.3,'  b =',f8.3,/)
      return
      end
