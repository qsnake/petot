      subroutine warn(i,xarg,iarg)
c     
c     subroutine evaluates and prints warning messages
c     $Id: warn.x,v 1.3 89/04/24 15:16:32 sverre Exp $
c
c     $Log:	warn.x,v $
c     Revision 1.3  89/04/24  15:16:32  sverre
c     checked in with -k by jimb at 89.05.02.18.03.39.
c     
c     Revision 1.3  89/04/24  15:16:32  sverre
c     Added warning 411 for subroutine setst.
c     
c     Revision 1.2  89/04/14  20:58:30  sverre
c     Renamed and added warnings for intpnt.
c     
c     Revision 1.1  89/03/11  19:09:51  sverre
c     Initial revision
c     
      implicit DOUBLE PRECISION(a-h,o-z)
c     
c     the warning number i is of the form xxy,
c     xx identifies the calling subroutine
c     y is the sequence number within the routine.
c     The following numbers have been assigned.
c     
c     1    main
c     3    tpage
c     5    crstl
c     7    paramt
c     9    gspace
c     11   struct
c     13   pseudo
c     15   ewald
c     17   screen
c     19   intpnt
c     21   matrix
c     22   matrix.iter
c     23   flevel
c     25   mtlden
c     27   velect
c     29   etotal
c     31   mixer
c     32   saveta
c     33   forces
c     34   dipole
c     35   denplt
c     37   decomp
c     39   savein
c     41   setst
c     101  symchk
c     102  intsub
c     103  msolve
c     105  charge
c     107  excorr
c     109  dencub
c     111  dentet
c     113  warn
c     
c     force static storage
c     

      save nwarn                    ! added by LWW
      data nwarn /0/                ! added by LWW
      common /warn1/ nout,narg
c
c     rcs id string - allows use of ident to identify binaries
c
      character rcsid*50
      rcsid = '$RCSfile: warn.x,v $$Revision: 1.3 $'
c     
c     initialize or get total number of warnings
c     
      if (i .eq. 0) then
         if (iarg .le. 0) then
            nout = -iarg
            nwarn = 0
         else
            iarg = nwarn
         end if
         return
      end if
      nwarn = nwarn + 1
      write(nout,1) i
    1 format(/,1x,10('*'),2x,'warning',2x,10('*'),2x,'(',i4,' )')
c     
c     decode warning number
c     
      i10 = i/10
      if (i10 .ge. 100) goto 3
      irm = i - 10*i10
      if (irm .eq. 0) goto 1130
      if (i10 .le. 0) goto 1130
      if (i10 .le. 10)
     +     goto (10,1130,30,1130,50,1130,70,1130,90,1130), i10
      i10 = i10 - 10
      if (i10 .le. 10)
     +     goto (110,1130,130,1130,150,1130,170,1130,190,1130), i10
      i10 = i10 - 10
      if (i10 .le. 10)
     +     goto (210,220,230,1130,250,1130,270,1130,290,1130), i10
      i10 = i10 - 10
      if (i10 .le. 10)
     +     goto (310,320,330,340,350,1130,370,1130,390,1130), i10
      i10 = i10 - 10
      if (i10 .le. 2)
     +     goto (410,1130), i10
      goto 1130
    3 i10 = i/10 - 100
      irm = i - 10*(i10+100)
      if (irm .eq. 0) goto 1130
      if (i10 .le. 0) goto 1130
      if (i10 .le. 10)
     +     goto (1010,1020,1030,1130,1050,1130,1070,1130,1090,1130), i10
      i10 = i10 - 10
      if (i10 .le. 3)
     +     goto (1110,1130,1130), i10
      goto 1130
c     
c     main
c     
   10 if (irm .gt. 3) goto 1130
      if (irm .eq. 1) write(nout,11)
   11 format(' improper program state for this operation')
      if (irm .eq. 2) write(nout,12)
   12 format(' complex structure factor - use complex program')
      if (irm .eq. 3) write(nout,13) iarg
   13 format(' arraylimit for energylevels exceeded - k point no',i3)
      stop
   30 goto 1130
c     
c     crstl
c     
   50 if (irm .gt. 6) goto 1130
      if (irm .eq. 1) write(nout,51) iarg
   51 format(' too many types of atoms - max is',i3)
      if (irm .eq. 2) write(nout,52) iarg
   52 format(' too many atoms of same type - max is',i4)
      if (irm .eq. 3) write(nout,53) iarg
   53 format(' unknown basis for symmetry operations - ibas is',i4)
      if (irm .eq. 4) write(nout,54) iarg
   54 format(' symmetry operation',i3,' does not preserve length')
      if (irm .eq. 5) write(nout,55) iarg
   55 format(' check symmetry operation',i3,' against basis')
      if (irm .eq. 6) write(nout,56)
   56 format(' fatal error in symmetry operations')
      if (irm .ne. 4 .and. irm .ne. 5) stop
      return
c     
c     paramt
c     
   70 if (irm .gt. 2) goto 1130
      if (irm .eq. 1) write(nout,71)
   71 format(' emax is set equal to gmax**2')
      if (irm .eq. 2) write(nout,72)
   72 format(' illegal input for metal / insulator / band occupation')
      if (irm .eq. 2) stop
      return
c     
c     gspace
c     
   90 if (irm .gt. 5) goto 1130
      if (irm .eq. 1) write(nout,91) xarg
   91 format(' gmax too large for arrays kx ...',
     +     ' - gmax reduced to',f10.6)
      if (irm .eq. 2) write(nout,92) xarg
   92 format(' gmax too large for arrays over stars',
     +     ' - gmax reduced to',f10.6)
      if (irm .eq. 3) write(nout,93) xarg,iarg
   93 format(' no match for prototype g-vector of length',f10.4,
     +     ' - symmetry operation',i3)
      if (irm .eq. 4) write(nout,94) iarg,xarg
   94 format(' phasefactor for g-vector no.',i5,
     +     ' has length',f10.6)
      if (irm .eq. 5) write(nout,95) xarg
   95 format(' gmax too large for index array indv',
     +     ' - gmax reduced to',f10.6)
      if (irm .eq. 3 .or. irm .eq. 4) stop
      return
c     
c     struct
c     
  110 if (irm .gt. 1) goto 1130
      write(nout,111) xarg,iarg
  111 format(' structure and phase factors do not agree',
     +     ' - max err =',e9.2,' - no err =',i5)
      if (xarg .gt. 1.0D-4) stop
      return
c     
c     pseudo
c     
  130 if (irm .gt. 6) goto 1130
      if (irm .eq. 1) write(nout,131) iarg
  131 format(' local potential too large for internal array',
     +     ' in pseudo - max is',i5)
      if (irm .eq. 2) write(nout,132) iarg
  132 format(' nonlocal arrays are too small - nqmax =',i4)
      if (irm .eq. 3) write(nout,133)
  133 format(' element name on tape does not agree with input')
      if (irm .eq. 4) write(nout,134)
  134 format(' correlation on tape does not agree with input')
      if (irm .eq. 5) write(nout,135) xarg
  135 format(' max g on tape is less than gmax - max g =',f6.2,
     +     ' - local potential')
      if (irm .eq. 6) write(nout,136) xarg
  136 format(' max g on tape is less than gmax and sqrt(emax)',
     +     ' - max g =',f6.2,' - non-local potential')
      if (irm .le. 2) stop
      return
c     
c     ewald
c     
  150 if (irm .gt. 0) goto 1130
      return
c     
c     screen
c     
  170 if (irm .gt. 4) goto 1130
      if (irm .eq. 1) write(nout,171)
  171 format(' unknown screening - atomic screening used')
      if (irm .eq. 2) write(nout,172) iarg
  172 format(' size of mixing matrix reduced to number of g stars',
     +     ' - size was',i4)
      if (irm .eq. 3) write(nout,173) iarg
  173 format(' screening with data from different structure',
     +     ' - g-vectors do not match - ierr =',i4)
      if (irm .eq. 4) write(nout,174)
  174 format(' tape10 empty - atomic screening used')
      return
c     
c     intpnt
c     
  190 if (irm .gt. 6) goto 1130
      if (irm .eq. 1) write(nout,191) iarg
  191 format(' internal array bound in intpnt exceeded',
     +     '  -  nmax =',i5)
      if (irm .eq. 2) write(nout,192) iarg
  192 format(' cannot express k point parameter vectors as',
     +     ' fractions of basis vectors - check vector #',i2)
      if (irm .eq. 3) write(nout,193) xarg
  193 format(' parameter vectors for k points are linearly dependent',
     +     '  -  xvol =',e9.2)
      if (irm .eq. 4) write(nout,194) iarg,xarg
  194 format(' k point grid is not translationally invariant',
     +     '  -  check vector #',i2,' - ',f10.3)
      if (irm .eq. 5) write(nout,195) iarg,xarg
  195 format (' k point grid generators fail symmetry check',
     +     '  -  check vector #',i2,' - ',f10.3)
      if (irm .eq. 6) write(nout,196) xarg
  196 format(' k-point sum is not zero or one - wsum-1.0 =',e9.2)
      if (irm .ne. 6) stop
      return
c     
c     matrix
c     
  210 if (irm .gt. 6) goto 1130
      if (irm .eq. 1) write(nout,211) xarg
  211 format(' new emax2 =',f10.6,' set in matrix -',
     +     ' emax2 was larger than gmax**2')
      if (irm .eq. 2) write(nout,212) xarg,iarg
  212 format(' new emax2 =',f10.6,' set in matrix -',
     +     ' matrix size =',i5,' - too large for arrays')
      if (irm .eq. 3) write(nout,213) xarg
  213 format(' new emax1 =',f10.6,' set in matrix -',
     +     ' emax1 was larger than emax2')
      if (irm .eq. 4) write(nout,214)
  214 format(' degenerate eigenvalue missing')
      if (irm .eq. 5) write(nout,215) iarg
  215 format(' size of eigenvector array in matrix exceeded',
     +     ', maxmv must be at least',i6)
      if (irm .eq. 6) write(nout,216) iarg
  216 format(' eigenvector',i3,' did not converge (tinvit)')
      if (irm .eq. 5) stop
      return
c     
c     matrix (iter special)
c     
  220 if (irm .gt. 6) goto 1130
      if (irm .eq. 1) write(nout,221) iarg
  221 format(' matrix array a too small - rayleigh ritz - maxr =',i4)
      if (irm .eq. 2) write(nout,222) iarg
  222 format(' matrix array a too small - submatrix - max1 =',i4)
      if (irm .eq. 3) write(nout,223) iarg
  223 format(' matrix array a too small - matrix - max2 =',i4)
      if (irm .eq. 4) write(nout,224) iarg
  224 format(' matrix size =',i4,
     +     ' too small - use different routine')
      if (irm .eq. 5) write(nout,225) xarg,iarg
  225 format(' new emax1 =',f10.6,' set in matrix -',
     +     ' submatrix size =',i5,' - too large for arrays')
      if (irm .eq. 6) write(nout,226) iarg
  226 format(' error in matrix read/write - unit no',i3)
      if (irm .le. 3 .or. irm .ge. 6) stop
      return
c     
c     flevel
c     
  230 if (irm .gt. 2) goto 1130
      if (irm .eq. 1) write(nout,231) xarg
  231 format(' number of electrons exceed number of bands',
     +     ' - excess charge =',e9.2)
      if (irm .eq. 2) write(nout,232) xarg
  232 format(' unable to compute fermi level - delta z =',e9.2)
      return
c     
c     mtlden
c     note that formats 251 and 252 are used by forces
c     
  250 if (irm .gt. 2) goto 1130
      if (irm .eq. 1) write(nout,251) iarg
  251 format(' no eigenvectors at k-point number',i4,
     +     ' - possible file error')
      if (irm .eq. 2) write(nout,252) iarg
  252 format(' missing eigenvector(s) at k-point number',i4)
      return
c     
c     velect
c     
  270 if (irm .gt. 2) goto 1130
      if (irm .eq. 1) write(nout,271) iarg
  271 format(' n(fft) reduced in velect - ierr =',i4)
      if (irm .eq. 2) write(nout,272) iarg,xarg
  272 format(' complex charge density in real space - ierr =',i5,
     +     ' - maxerr =',e9.2)
      if (irm .eq. 2 .and. xarg .gt. 1.0D-4) stop
      return
c     
c     etotal
c     
  290 if (irm .gt. 1) goto 1130
      write(nout,291) xarg
  291 format(' this system is not charge neutral',
     +     ' - excess el density =',e9.2)
      return
c
c     mixer
c
  310 goto 1130
c
c     saveta
c
  320 if (irm .gt. 1) goto 1130
      write(nout,321)
  321 format(' unable to open save file - saveta')
      stop

c     
c     forces
c     
  330 if (irm .gt. 3) goto 1130
      if (irm .eq. 1) write(nout,251) iarg
      if (irm .eq. 2) write(nout,252) iarg
      if (irm .eq. 3) write(nout,333) iarg
  333 format(' symmetry error (in forces) - operation',i3)
      return
c     
c     dipole
c     
  340 if (irm .gt. 1) goto 1130
      write(nout,341) xarg
  341 format(' nonorthogonal eigenvectors (in dipole) - ',e9.2)
      return
c     
c     denplt
c     
  350 if (irm .gt. 1) goto 1130
      write(nout,351) iarg
  351 format(' too many eigenvectors for array in denplt',
     +     ' - kpoint no',i4)
      stop
c     
c     decomp
c     
  370 if (irm .gt. 1) goto 1130
      write(nout,371) iarg
  371 format(' too many eigenvectors for array in decomp',
     +     ' - kpoint no',i4)
      stop
c
c     savein
c
  390 if (irm .gt. 1) goto 1130
      write(nout,391)
  391 format(' unable to open save file - savein')
      stop
c
c     setst
c
  410 if (irm .gt. 1) goto 1130
      write(nout,411)
  411 format(' invalid state flag (ignored)')
      return
c     
c     symchk
c     
 1010 if (irm .gt. 2) goto 1130
      if (irm .eq. 1) write(nout,1011) iarg
 1011 format(' the symmetry operations as input',
     +     ' do not form a group, check operation no',i4)
      if (irm .eq. 2) write(nout,1012)
 1012 format(' more than 48 symmetry operations')
      return
c     
c     intsub
c     
 1020 if (irm .gt. 1) goto 1130
      if (irm .eq. 1) write(nout,1021) iarg
 1021 format(' number of k points is greater than array bound',
     +     '  -  nrkmax =',i3)
      stop
c     
c     msolve
c     
 1030 if (irm .gt. 9) goto 1130
      if (irm .eq. 1) write(nout,1031)
 1031 format(' array a in msolve too small - no sub-matrix vectors')
      if (irm .eq. 2) write(nout,1032) iarg
 1032 format(' array a in msolve too small - blocksize - maxb =',i5)
      if (irm .eq. 3) write(nout,1033) iarg
 1033 format(' rayleigh-ritz basis linearly dependent - ierr =',i2)
      if (irm .eq. 4) write(nout,1034) iarg,xarg
 1034 format(' eigenvector',i3,' did not converge - resid =',e9.2)
      if (irm .eq. 5) write(nout,1035) xarg,iarg
 1035 format(' small norm =',e9.2,'  for eigenvector',i3)
      if (irm .eq. 6) write(nout,1036) iarg
 1036 format(' size of eigenvector array in eismat exceeded',
     +     ' - need space for',i3,' degenerate vectors')
      if (irm .eq. 7) write(nout,1037) iarg
 1037 format(' rayleigh-ritz basis too small - only',i4,
     +     ' linearly independent vectors')
      if (irm .eq. 8) write(nout,1038) iarg,xarg
 1038 format(' eigenvalue error detected on block boundary - redoing',
     +     i3,' vectors - xerr =',e9.2)
      if (irm .eq. 9) write(nout,1039) iarg,xarg
 1039 format(' converged eigenvectors not orthogonal - nvec =',i4,
     +     ' max overlap =',e9.2)
      if (irm .le. 2 .or. irm .eq. 6 .or. irm .eq. 7) stop
      return
 1050 goto 1130
c     
c     excorr
c     
 1070 if (irm .gt. 1) goto 1130
      write(nout,1071) iarg
 1071 format(1x,a2,' correlation unknown')
      stop
 1090 goto 1130
c     
c     dentet
c     
 1110 if (irm .gt. 1) goto 1130
      ryd = 13.60580
      if (irm .eq. 1) write(nout,1111) ryd*xarg
 1111 format(' 4 energies equal in dentet - e =',f7.3)
      return
c     
c     warn
c     
 1130 write(nout,1131) iarg,xarg
 1131 format(' warning unknown - iarg =',i5,' - xarg =',e9.2)
      return
      end
