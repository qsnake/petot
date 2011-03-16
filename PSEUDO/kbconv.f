      PROGRAM KBCONV 
c
c ********************************************************************
c *                                                                  *
c *  Program converts a nonlocal or spin-orbit pseudopotential       *
c *  to the Kleinman and Bylander form and does the `Fourier'        *
c *  like Bessel transforms.  See Phys. Rev. Lett. v.48 no.20        *
c *  page 1425 as ref. The local potential indicator                 *
c *  are stored at the end of the                                    *
c *  pseudokb.dat data file.  Note that the projector is stored      *
c *  as Phi(r)deltaV(r), not as rPhi(r)deltaV(r).                    *
c *                                                                  *
c *  Written by Norm Troullier while at the U. of MN                 *
c *  Copyright Norm Troullier and Jose Luis Martins                  *
c *  Version 1.25 Dated Nov. 12, 1990                                *
c *  Minor modifications by jlm, 2/2/96, 13/6/97                     *
c *                                                                  *
c ********************************************************************
c
      IMPLICIT REAL*8 (A-H,O-Z)
c
      PARAMETER (NRMAX=2000,LMAX=5)
c
      PARAMETER (ZERO=0.0D0,SMALL=0.5D-5,ONE=1.0D0,TWO=2.0D0)
c
      DIMENSION R(NRMAX),RAB(NRMAX),NO(2*LMAX),LO(2*LMAX),SO(LMAX),
     1 VIOD(LMAX,NRMAX),VIOU(LMAX,NRMAX),VID(NRMAX),CDD(NRMAX),
     2 CDC(NRMAX),EV(LMAX),Y(NRMAX),YP(NRMAX),YPP(NRMAX),W(NRMAX*3),
     3 CDU(NRMAX),VIU(NRMAX),S1(NRMAX),S2(NRMAX),ETOT(10),EVL(2),
     4 AR(NRMAX,LMAX),ANORM(LMAX),INORM(LMAX),VQL(NRMAX),ARD(NRMAX),
     5 EVI(5),EVD(5)
c
c     CHARACTER*1 ISPP
      CHARACTER*2 ICORR,NAMEAT 
      CHARACTER*3 IREL
      CHARACTER*4 NICORE 
      CHARACTER*10 IRAY(6),ITITLE(7),DATED
c
      DATA SO/5*ZERO/
      DATA EVI/5*ZERO/
      DATA EVD/5*ZERO/
      DO 1 I=1,NRMAX
        CDU(I)=ZERO
 1    CONTINUE
c  
c  Open the `kb.dat' input file - unit=2, the `kb.out'
C  output file - unit=3, and the `kbplot.dat' plotting 
c  file - unit=4.
c  Read in the angular momentum quantum number 
c  for the local potential part, the number of q-space
c  points and the seperation of the q-space points(a.u.).
c
c      OPEN (UNIT=2,FILE='kb.dat',FORM='FORMATTED',
c     1 STATUS='OLD')
      OPEN (UNIT=6,FILE='kb.out',FORM='FORMATTED',
     1 STATUS='unknown')
      OPEN (UNIT=4,FILE='kbplot.dat',FORM='FORMATTED',
     1 STATUS='unknown')
      
      open(unit=2, file='atom.input',form='formatted',
     &  status='old')
      rewind(2)
      read(2,*) 
      read(2,*) 
      read(2,*) 
      read(2,*) 
      read(2,*) in,nv
      do i=1,nv
      read(2,*)
      enddo
      read(2,*)
      read(2,*)
      read(2,*) LLOCAL
      close(2)

c      READ(2,2000)LLOCAL
c      READ(2,2001)NQL,DELQL
      NQL=1600
      DELQL=0.015d0
 2000 FORMAT(5X,I2)
 2001 FORMAT(10X,I4,10X,F10.8) 
      DO 2 J=1,5
c        READ(2,4001)EVI(J)
      EVI(J)=0.0d0
 2    CONTINUE
 4001 FORMAT(5X,F10.8)
c      CLOSE(UNIT=2)
c
c  Print out heading to `kb.out' file.
c
      CALL ZEDATE(DATED)
      WRITE(6,2002)DATED
 2002 FORMAT(2X,'Kleinman and Bylander pseudopotential conversion',
     1 /,'  and Bessel/Fourier transform program. Version 1.25',//,
     2 '  Ran on ',A10)
      WRITE(6,2003)LLOCAL
 2003 FORMAT(//,'  The l=',I2,' pseudopotentail is treated as the',
     1 /,'  local potential')
      WRITE(6,2004)NQL,DELQL
 2004 FORMAT(/,'  NQL=',I4,' and DELQL=',F10.8)  
c
c  Open and read in data from file `pseudo.dat'.
c  Open and output some of the same data to file `pseudokb.dat'
c  for use in the atomkb program.  Open and output some of
c  the same data to `potfourkb.dat' for use in the planewave
c  program.
c
      OPEN (UNIT=7,FILE='pseudo.dat',FORM='UNFORMATTED')
      READ(7) NAMEAT,ICORR,IREL,NICORE,(IRAY(I),I=1,6),
     1 (ITITLE(I),I=1,7),NORB,NUMU,NR,A,B,ZNUC
      NR=NR+1
      IZNUC=INT(ZNUC+0.1)
      READ(7) (R(I),I=2,NR)
      DO 5 I=1,NORB
        READ(7) LO(I),(VIOD(LO(I)+1,J),J=2,NR)
 5    CONTINUE
      DO 10 I=1,NUMU
        READ(7) LO(I+NORB),(VIOU(LO(I+NORB)+1,J),J=2,NR)
 10   CONTINUE
      READ(7) (CDC(I),I=2,NR)
      READ(7) (CDD(I),I=2,NR)
      CLOSE (UNIT=7)
c
c  Set up RAB integration grid
c
      R(1)=ZERO
      DO 15 I=1,NR
        RAB(I) = (R(I)+A)*B
 15   CONTINUE
c
c   Find the total amount of charge in the valence charge.
c   Fit cdd to splines and integrate.
c
      IF (NICORE .EQ. 'nc  ') THEN
        ITYPE=0
      ELSEIF(NICORE .EQ. 'pcec' .OR. NICORE .EQ. 'fcec') THEN
        ITYPE=1
      ELSE
        ITYPE=2 
      ENDIF
      Y(1) = ZERO
      DO 20 I=2,NR
        Y(I) = CDD(I)
 20   CONTINUE
      ISX = 0
      A1 = ZERO
      AN = ZERO
      B1 = ZERO
      BN = ZERO
      NRM = NR
      CALL SPLIFT(R,Y,YP,YPP,NRM,W,IERR,ISX,A1,B1,AN,BN)
      IF ( IERR .NE. 1) THEN
        WRITE(6,2005)IERR
        STOP 'SPLIFT'
      ENDIF
 2005 FORMAT(1X,'****** Error in splift ierr =',I2)
      XLO = ZERO 
      NUP = 1
      CALL SPLIQ(R,Y,YP,YPP,NRM,XLO,R(NR),NUP,TOTVEL,IERR)
      IF ( IERR .NE. 1 ) THEN
        WRITE(6,2006)IERR
        STOP 'SPLIQ'
      ENDIF
 2006 FORMAT(1X,'****** Error in spliq ierr =',I2)
c
c  Write out `pseudo.dat' info to `kb.out'.
c
      WRITE(6,2007)NAMEAT,ICORR,IREL,NICORE,TOTVEL,
     1 (IRAY(J),J=1,6),
     1 (ITITLE(J),J=1,7),A,B
 2007 FORMAT(/,A2,2X,A2,2X,A3,2X,A4,' pseudopotential with ',
     1 F10.6,' valence electrons',//,1X,6A10,4X,//,1X,7A10,//,
     2 ' a = ',F9.7,' b = ',F9.7,//)
c
c  Find el-el potential
c
      CALL VELECT(0,0,ICORR,' ',ITYPE,NR,R,RAB,TOTVEL,CDD,
     1 CDU,CDC,VID,VIU,ETOT,Y,YP,YPP,S1,S2,W)
c
c  Set up potentials
c
      C2 = -ONE/B**2
      C1 = -2*C2 + ONE/4
      DO 25 I=1,NORB
        LP=LO(I)+1
        NO(I)=LP
        LLP=LP*(LP-1)
c
c  Set up hamiltonian matrix for kinetic energy,
c  only the diagonal depends on the potential.
c
        Y(1)  = C1 / (R(2)+A)**2
        YP(1)  = ZERO
        YPP(1) = ZERO
        DO 30 K=3,NR
          Y(K-1)  = C1 / (R(K)+A)**2
          YP(K-1)  = C2 / ((R(K)+A)*(R(K-1)+A))
          YPP(K-1) = YP(K-1)**2
 30     CONTINUE
        DO 35 K=2,NR
          Y(K-1)=Y(K-1)+(VIOD(LP,K)+LLP/R(K))/R(K)+VID(K)
 35     CONTINUE
c
c  Diagonalize and find wave function and store in AR().
c               
        EPS = -ONE
        CALL TRIDIB(NR-1,EPS,Y,YP,YPP,BL,BU,1,1,E,IND,IERR,
     1   S1,S2) 
        EV(I)=E
        DO 40 J=2,NR
          Y(J)=(VIOD(LP,J)+LLP/R(J))/R(J)+VID(J)
 40     CONTINUE
        IFLAG = 0
        CALL DIFNRL(1,I,Y,AR(1,LP),YP,LMAX,NR,A,B,R,RAB,NORB,NO,
     1   LO,SO,ZNUC,VIOD,VIOU,VID,VIU,EV,IFLAG,YPP,S1,S2,EVI)
 25   CONTINUE
c
c  Create deltaV(r) for nonlocal parts, and find local
c  part of potential.
c                    
      IF (LLOCAL .GE. 0 ) THEN
        DO 48 J=2,NR
          VQL(J) = VIOD(LLOCAL+1,J)
 48     CONTINUE
      ELSEIF(LLOCAL .EQ. -1) THEN
        K1=NORB
        K2=NORB
        K3=NORB
        K4=NORB
        K5=NORB
        K6=NORB
        K7=NORB
        DO J=2,NR
          VQL(J) = VIOD(NORB,J)
          KI = NORB
          DO I=1,NORB-1
            IF (VIOD(I,J) .LT. VQL(J)) THEN
              VQL(J)=VIOD(I,J)
              KI=I
            ENDIF
          ENDDO
          K7=K6
          K6=K5
          K5=K4
          K4=K3
          K3=K2
          K2=K1
          K1=KI
          IF (K7 .NE. K2 .AND. J .GT. 9) THEN
            DO I=1,1
              VQL(J-I) = (VQL(J-I-1)/R(J-I-1)+VQL(J-I+1)/R(J-I+1))
     1         *0.5*R(J-I)
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      DO 50 I=1,NORB
        LP = LO(I) + 1
        ANORM(LP) = ZERO 
        A2NORM = ZERO
        DO 54 J=2,30
          VIOD(LP,J)=VIOD(LP,J)/R(J)-VQL(J)/R(J)
 54     CONTINUE
        DO 55 J=31,NR
          VIOD(LP,J)=(VIOD(LP,J)-VQL(J))/R(J)
 55     CONTINUE                             
c
c  Calculate the normalizing coeffcient for nonlocal parts, 
c  uses Bode's rule for integration.
c
        AR(1,LP) = ZERO
        AR(NR+1,LP) = ZERO
        AR(NR+2,LP) = ZERO
        AR(NR+3,LP) = ZERO
        AR(NR+4,LP) = ZERO 
        W(J+1) = ZERO 
        W(J+2) = ZERO
        W(J+3) = ZERO
        W(J+4) = ZERO
        DO 56 J=1,NR
          W(J) = AR(J,LP)*AR(J,LP)*RAB(J)*VIOD(LP,J)
 56     CONTINUE
        DO 60 J=1,NR,4        
          ANORM(LP)=ANORM(LP)+7*(W(J)+W(J+4))+
     1     32*(W(J+1)+W(J+3))+12*W(J+2)
 60     CONTINUE
        ANORM(LP)=2*ANORM(LP)/45
C
CJLM
C
        IF (ABS(ANORM(LP)) .LT. 1.0E-12 .AND. 
     1                   LLOCAL .NE. LO(I)) THEN
          PRINT*,' WARNING ANORM IS SMALL',ANORM(LP)
C
CJLM
C
        ENDIF
        IF ( ANORM(LP) .LT. ZERO) THEN
          INORM(LP) = -1
          ANORM(LP) = SQRT(-ANORM(LP))
        ELSEIF(ANORM(LP) .GT. ZERO) THEN
          INORM(LP) = 1   
          ANORM(LP) = SQRT(ANORM(LP))
        ENDIF
        IF (LLOCAL .EQ. LO(I)) THEN
          INORM(LP) = 0
          ANORM(LP) = ONE
        ENDIF
        IF (INORM(LP) .NE. 0) THEN
          WMAX = ZERO
          DO 59 J=1,NR
            WMAX=MAX(WMAX,ABS(W(J)))
 59       CONTINUE
          DO 61 J=1,NR
            Y(J) = W(J)/WMAX
 61       CONTINUE
        ENDIF
        CALL PLOTKB(NR,Y,R,INORM(I),LO(I),'t')
        DO 57 J=1,NR
          W(J) = AR(J,LP)*AR(J,LP)*RAB(J)*ABS(VIOD(LP,J))
 57     CONTINUE
        DO 58 J=1,NR,4        
          A2NORM=A2NORM+7*(W(J)+W(J+4))+32*(W(J+1)+W(J+3))+12*W(J+2)
 58     CONTINUE
        A2NORM=2*A2NORM/45
        IF ( A2NORM .NE. ZERO) THEN
          RATIO = ANORM(LP)*ANORM(LP)/A2NORM
          WRITE(6,2021)LP-1,INORM(LP)*RATIO
        ENDIF
 2021 FORMAT(1X,'Ratio for l=',I1,' is ',F10.7)
c
c  Calculate projector, NOTE stored as Phi(r)deltaV(r).
c
        IF (ANORM(LP) .EQ. ZERO) THEN
          DO 64 J=2,NR
            VIOD(LP,J)=ZERO
 64       CONTINUE
        ELSE
          DO 65 J=2,NR
            VIOD(LP,J)=(AR(J,LP)/R(J))*VIOD(LP,J)/ANORM(LP)
 65       CONTINUE
        ENDIF
 50   CONTINUE
c
      DO 67 J=2,NR
        VQL(J) = VQL(J)/R(J)
 67   CONTINUE 
c
c   Find the 1st and 2nd eigenvalues for all angular 
c   momentum using just the local potential.  Print
c   out results, this is used to find if any ghost 
c   states exist for the potential.  See (preprint)paper
c   of Gonze, Kackell, and Scheffler, Phy. Rev. B. 
c
      WRITE(6,301)
 301  FORMAT(//,'  Ghost State Existence Data',/
     1 ' l     0,node-eigen   1,node-eigen    True-eigen   inorm',/)
      DO 300 I=1,NORB
c
c  Set up the potential.
c   
        IF (LLOCAL .EQ. LO(I)) GOTO 300
        EVT=EV(I)
        IF (EV(I) .GE. 0.0) EV(I)=-SMALL
        LP=LO(I)+1
        LLP = LP * (LP-1)
        STORE = VIOD(I,2)
        VIOD(I,2) = VQL(2)*R(2)
        DO 318 J=2,NR
          Y(J) = VQL(J) + LLP/R(J)**2 + VID(J)
 318    CONTINUE
c
c  Call the integration routine.
c
        NOT=NO(I)
        DO 345 JJ=0,1
          DO 317 J=1,NR
            ARD(J)=ZERO
 317      CONTINUE
          NO(I)=NOT+JJ
          IFLAG = 0
          CALL DIFNRL(1,I,Y,ARD,YP,LMAX,NR,A,B,R,RAB,NORB,NO,
     1     LO,SO,ZNUC,VIOD,VIOU,VID,VIU,EV,IFLAG,YPP,S1,S2,EVD)
           EVL(JJ+1)=EV(I)
 345    CONTINUE
c
        VIOD(I,2) = STORE
        EV(I) = EVT
        NO(I) = NOT
c
        WRITE(6,355)LO(I),EVL(1),EVL(2),EV(I),INORM(I)
 355  FORMAT(1X,I1,5X,3(F12.9,4X),I2)
        IF (INORM(I) .LT. 0 ) THEN
          IF (EV(I) .GT. EVL(1)) THEN
            WRITE(6,302)
            IF (EVI(I) .NE. ZERO) WRITE(6,903)
          ENDIF
        ELSEIF (INORM(I) .GT. 0 ) THEN
          IF (EV(I) .GT. EVL(2)) THEN
            WRITE(6,302)
            IF (EVI(I) .NE. ZERO) WRITE(6,903)
          ENDIF
        ENDIF
 302  FORMAT(1X,' WARNING GHOST STATE WILL BE PRESENT!!!!!!')
 903  FORMAT(1X,' WARNING: GHOST STATE MAY BE DUE TO NON-EIGENVALUE')
 300  CONTINUE
c
c  Open and write heading to `pseudokb.dat' file.
c  Store data back into file `pseudokb.dat'.
c
      OPEN (UNIT=8,FILE='pseudokb.dat',FORM='UNFORMATTED')
      WRITE(8) NAMEAT,ICORR,IREL,NICORE,(IRAY(I),I=1,6),
     1 (ITITLE(I),I=1,7),NORB,NUMU,NR-1,A,B,ZNUC
      WRITE(8) (R(I),I=2,NR)
      DO 70 I=1,NORB
        WRITE(8) LO(I),(VIOD(I,J),J=2,NR)
 70   CONTINUE
      DO 75 I=1,NUMU
        WRITE(8) LO(I+NORB),(VIOU(LO(I+NORB)+1,J),J=2,NR)
 75   CONTINUE
      WRITE(8) (CDC(I),I=2,NR)
      WRITE(8) (CDD(I),I=2,NR)
      WRITE(8) (VQL(J),J=2,NR)
      WRITE(8) NORB
      DO 80 I=1,NORB
        WRITE(8)INORM(LO(I)+1),ANORM(LO(I)+1)
 80   CONTINUE
      CLOSE(UNIT=8)
c
c   Send operator data to plotkb for writing of data to 
c   the `kbplot.dat' file.
c
      DO 85 I=1,NORB
        DO 90 J=2,NR
          Y(J) = VIOD(LO(I)+1,J)
 90     CONTINUE
       CALL PLOTKB(NR,Y,R,INORM(LO(I)+1),LO(I),'r')
 85   CONTINUE
c
c  Print out the local and K & B potentials
c
      WRITE(6,2008)NORB,NR
 2008 FORMAT(//,' Pseudo potentials in real space',//,
     1 ' number of potentials =',I2,/,
     2 ' number of radial grid points =',I4,//)
      IF (NORB .EQ. 1) THEN
        WRITE(6,2009)
      ELSEIF (NORB .EQ. 2) THEN
        WRITE(6,2010)
      ELSEIF (NORB .EQ. 3) THEN
        WRITE(6,2011)
      ELSEIF (NORB .EQ. 4) THEN
        WRITE(6,2012)          
      ELSE
        WRITE(6,3013)
      ENDIF
 2009 FORMAT('Index   r(I)    Local    S pot.   Core c  Val c',/)
 2010 FORMAT('Index   r(I)    Local    S pot.   P pot.   Core c',
     1 '  Val c',/)
 2011 FORMAT('Index   r(I)    Local    S pot.   P pot.   D pot.',
     1 '  Core c  Val c',/)
 2012 FORMAT('Index   r(I)    Local    S pot.   P pot.   D pot.',
     1 '  F pot.  Core c  Val c',/)
 3013 FORMAT('Index   r(I)    Local    S pot.   P pot.   D pot.',
     1 '   F pot.  G pot.  Core c  Val c',/)
      DO 95 J=2,NR
        WRITE(6,2013) J,R(J),VQL(J),(VIOD(I,J),I=1,NORB),CDC(J),CDD(J)
 95   CONTINUE
 2013 FORMAT(I4,1X,F9.5,6F9.4,2F6.4)
c
c  Fourier transform local potential, first for q=0 and
c  then for rest q>0.  Use Bode's rule to integrate.
c
      PI4 = 16*ATAN(ONE)
      VQL0 = ZERO
      W(1) =ZERO
      VT = 2*ZNUC 
      LOC = LLOCAL + 1
      DO 100 K=2,NR
        W(K) = RAB(K)*R(K)*(R(K)*VQL(K)+VT)
 100  CONTINUE
      DO 105 K=NR+1,NR+4
        W(K) = ZERO
 105  CONTINUE
      DO 110 K=1,NR,4
        VQL0 = VQL0 + 7*W(K) + 32*W(K+1) + 12*W(K+2) +
     1   32*W(K+3) + 7*W(K+4)
 110  CONTINUE
      VQL0 = 2 * VQL0 * PI4 / 45
c
      Y(1) = ZERO
      Y(NR+1) = ZERO
      Y(NR+2) = ZERO
      Y(NR+3) = ZERO
      Y(NR+4) = ZERO 
      DO 111 K=2,NR
        W(K) = W(K)/R(K)
 111  CONTINUE
      YP(1) = ZERO
      DO 112 J=2,NQL+1
        YP(J) = DELQL * (J-1)
 112  CONTINUE
      DO 115 J=1,NQL
        VQL(J) = ZERO
        DO 120 K=2,NR
          Y(K) = SIN(YP(J+1)*R(K))*W(K)
 120    CONTINUE
        DO 125 K=1,NR,4
          VQL(J) = VQL(J) + 7*Y(K) + 32*Y(K+1) + 12*Y(K+2) +
     1     32*Y(K+3) + 7*Y(K+4)
 125    CONTINUE
        VQL(J) = PI4 * (2 * VQL(J) / 45 - VT/YP(J+1))/YP(J+1)
 115  CONTINUE
c
c  Do loop over K&B projector potentials
c
      DO 130 I=1,NORB
        DO 131 K=NR,2,-1
          IF (VIOD(I,K) .EQ. ZERO) THEN
            NRM=K
          ELSE
            GOTO 134
          ENDIF
 131    CONTINUE
 134    DO 135 K=2,NRM
          W(K) = RAB(K)*R(K)*R(K)*VIOD(I,K)
 135    CONTINUE
        Y(1) = ZERO
        DO 136 K=NRM+1,NRM+7
          Y(K) = ZERO
 136    CONTINUE
        CRNORM = 7*SQRT((2*LO(I)+1)*PI4)/17280
        DO 145 J=1,NQL+1
          DO 150 K=2,NRM
            Y(K) = SBESSJ(I-1,YP(J)*R(K))*W(K)
 150      CONTINUE
c
c  Due to the high number of occilations in the intagrand,
c  an eight point Newton-Cotes intagration method is used.
c  See  Abramowitz and Stegun Eq. 25.4.17
c
          AR(J,I) = ZERO
          DO 155 K=1,NRM,7
            AR(J,I) = AR(J,I)+751*(Y(K)+Y(K+7))+3577*(Y(K+1)+Y(K+6))+
     1       1323*(Y(K+2)+Y(K+5))+2989*(Y(K+3)+Y(K+4))
 155      CONTINUE
          AR(J,I) = CRNORM * AR(J,I) 
 145    CONTINUE
 130  CONTINUE
c
c  Fourtier transform core charge density, ignore if
c  itype is equal to 0.
c      
      IF (ITYPE .EQ. 0) THEN
        DO 160 J=2,NQL+1
          CDC(J) = ZERO
 160    CONTINUE
      ELSE 
        DO 165 J=2,NR
          W(J) = RAB(J)*CDC(J)/R(J)
 165    CONTINUE
        DO 170 K=2,NQL+1
          CDC(K) = ZERO
          DO 175 J=1,NR
            Y(J) = W(J)*SIN(YP(K)*R(J))
 175      CONTINUE
          DO 180 J=1,NR,4
            CDC(K) = CDC(K) + 7*Y(J) + 32*Y(J+1) + 12*Y(J+2) +
     1       32*Y(J+3) + 7*Y(J+4)
 180      CONTINUE
          CDC(K) = 2*CDC(K)/45/YP(K)
 170    CONTINUE 
      ENDIF
c
c  Fourier transform the valence charge density.
c
      DO 185 J=2,NR
        W(J) = RAB(J)*CDD(J)/R(J)
 185  CONTINUE
      DO 190 K=2,NQL+1
        CDD(K) = ZERO
        DO 195 J=1,NR
          Y(J) = W(J)*SIN(YP(K)*R(J))
 195    CONTINUE
        DO 200 J=1,NR,4
          CDD(K) = CDD(K) + 7*Y(J) + 32*Y(J+1) + 12*Y(J+2) +
     1     32*Y(J+3) + 7*Y(J+4)
 200    CONTINUE
        CDD(K) = 2*CDD(K)/45/YP(K)
 190  CONTINUE
c
c  Send K&B potentials to plotkb
c
      DO 205 J=1,NORB
        IF (J .NE. LOC) THEN
          CALL PLOTKB(NQL,AR(1,J),YP,INORM(J),LO(J),'q')
        ENDIF
 205  CONTINUE 
      CLOSE(UNIT=4)
c
c  Open and write heading to `potfourkb.dat' file.
c  Write the transforms to file `potfourkb.dat'
c
      OPEN (UNIT=9,FILE='potfourkb.dat',FORM='UNFORMATTED',
     1 RECL=16000)
      WRITE(9) NAMEAT,ICORR,IREL,NICORE,(IRAY(I),I=1,6),
     1 (ITITLE(I),I=1,7)
      WRITE(9) IZNUC,NQL,DELQL,VQL0
      WRITE(9) NORB,(INORM(I),I=1,NORB)
      WRITE(9) (VQL(J),J=1,NQL)
      DO 210 I=1,NORB
        WRITE(9) (AR(J,I),J=1,NQL+1)
 210  CONTINUE
      WRITE(9) (CDC(J),J=2,NQL+1)
      WRITE(9) (CDD(J),J=2,NQL+1)
      CLOSE(UNIT=9)
c
c  Print out info to `kb.out' file.
c        
      WRITE(6,2014)NQL,VQL0
 2014 FORMAT(//,' Pseudo potentials in fourier space',//,
     1 ' Number of q-space points = ',I4,/,
     2 ' V(q=0) = ',f10.5,/)
      WRITE(6,2020)NORB,(INORM(I),I=1,NORB)
 2020 FORMAT(1X,'Number of potentials = ',I1,/
     1 ' Normalization indexs ',5(I2,5X))
      IF (NORB .EQ. 1) THEN
        WRITE(6,2015)
      ELSEIF (NORB .EQ. 2) THEN
        WRITE(6,2016)
      ELSEIF (NORB .EQ. 3) THEN
        WRITE(6,2017)
      ELSEIF (NORB .EQ. 4) THEN
        WRITE(6,2018)          
      ELSE
        WRITE(6,3018)
      ENDIF
 2015 FORMAT(/,' I   q(I)       Local      S pot.   Core c    Val c',/)                      
 2016 FORMAT(/,' I   q(I)       Local      S pot.   P pot.   Core c ',
     1 '   Val c',/)                      
 2017 FORMAT(/,' I   q(I)       Local      S pot.   P pot.   D pot. ',
     1 '   Core c   Val c',/)                                         
 2018 FORMAT(/,' I   q(I)       Local      S pot.   P pot.   D pot. ',
     1 '   F pot.   Core c   Val c',/)
 3018 FORMAT(/,' I   q(I)       Local      S pot.   P pot.   D pot. ',
     1 '    F pot.   G pot.   Core c   Val c',/)
      WRITE(6,2030)(AR(1,I),I=1,NORB)
 2030 FORMAT(2X,'0',3X,'0.000',14X,5F9.5)
      DO 215 J=1,NQL
        WRITE(6,2019)J,YP(J+1),VQL(J),(AR(J+1,I),I=1,NORB),
     1   CDC(J+1),CDD(J+1)
 215  CONTINUE
 2019 FORMAT(I4,F7.3,F14.5,7F9.5)
      CLOSE(UNIT=6)
c
      STOP
      END

      subroutine plotkb(nr,vd,r,inorm,lo,rq)
c
c ***********************************************************
c *                                                         *
c *  This is a plotting routine for kbconv.  Prints         *
c *  out the K & B operator or transform.                   *
c *                                                         *
c ***********************************************************
c  
c  njtj
c  ###  Cray conversions  
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter statement
c  ###    to single precision statement.
c  ###  Cray conversions
c  njtj
c  
      implicit double precision (a-h,o-z)
c
      parameter (zero=0.D0,pzf=0.05D0,small=0.0001D0)
Cray      parameter (zero=0.0,pzf=0.05)
c
      character*1 rq
      character*3 marker
c
      dimension vd(nr),r(nr)
c
c  Step size of 0.05 is adjustable as seen fit to give 
c  a reasonalble plot.
c
      if (inorm .eq. 0) then
        return
      elseif (rq .eq. 'r' .or. rq .eq. 't') then
        nrm = nr - 120
        nst = 2
        step=small
      elseif (rq .eq. 'q') then
        nrm = nr
        nst = 1
        step=zero
      endif
      do 150,j=nst,nrm
        if (r(j) .ge. step) then
          write(4,6000)r(j),vd(j)
          step=step+pzf
        endif
 150  continue             
      if (rq .eq. 'r') then
        if (lo .eq. 0) then
          marker='psr'
        elseif (lo .eq. 1) then
          marker='ppr'
        elseif (lo .eq. 2) then
          marker='pdr'
        elseif (lo .eq. 3) then
          marker='pfr'         
        elseif (lo .eq. 4) then 
          marker='pgr'
        endif
      elseif (rq .eq. 'q') then
        if (lo .eq. 0) then
          marker='psq'
        elseif (lo .eq. 1) then
          marker='ppq'
        elseif (lo .eq. 2) then
          marker='pdq'
        elseif (lo .eq. 3) then
          marker='pfq'         
        elseif (lo .eq. 4) then
          marker='pgq'
        endif
      else
        if (lo .eq. 0) then
          marker='pst'
        elseif (lo .eq. 1) then
          marker='ppt'
        elseif (lo .eq. 2) then
          marker='pdt'
        elseif (lo .eq. 3) then
          marker='pft'         
        elseif (lo .eq. 4) then
          marker='pgt'
        endif
      endif
      write(4,6001)marker
      return
c
c  Format statements
c
 6000 format(1x,f7.4,3x,f10.6)
 6001 format(1x,'marker ',a3)
      end
