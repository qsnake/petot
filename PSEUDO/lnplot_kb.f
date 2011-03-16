      PROGRAM LNPLOTKB 
c
c    This program used with lnpt, will construct a d(ln[R])/dr
c  vs Energy plot for a K & B operator pseudopotential.  
c  General warnings:  This program is not bullet proof and the
c  user should carefully check their results.  Some of the
c  more probable errors are:
c  1)  The pseudopotential file has fewer potentials then the 
c  datafile, ie. lmax of the pseudofile is less then lmax of the 
c  datafile.
c  2)  The configuration of the datafile is not the same as the
c  pseudofile.  In this case either create a datafile with the
c  same configuration or use the atom program to modify the
c  charge density of the pseudofile('pm' option).
c  3)  The datafile and pseudfile are not the same method,
c  ie. relativistic, spinpolarized, nonspinpolarized.  Since this
c  routine was designed for nonspin test, a relativistic spin are 
c  reduced to nonspin.  A spinpolarized should only be done with
c  a both up and down spin occupations equal(nonspin) for the 
c  datafile.  This does not mean a pseudofile with spinpolarized
c  cannot be used, but that a spinpolarized pseduofile shopuld be
c  modified with the atom program ('pm' option) into a 
c  nonspinpolarized configuration.
c  

cccccc
c mmga
cccccc


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

cccccc
c mmga
c
      PARAMETER(MAXPTS=2001,MTOT=10005)
c      PARAMETER(MAXPTS=800,MTOT=4000)
c
      DIMENSION R(2000),RAB(2000),NO(40),LO(40),SO(40),
     1 VIOD(5,2000),VIOU(5,2000),VID(2000),VIU(2000),EV(40),
     2 CDD(2000),CDC(2000),CDU(2000),DUM2(2000),ZO(40)
c
      DIMENSION DRORYT(MAXPTS,5),DRORYP(MAXPTS,5),DUMU(MAXPTS,5),
     1 DRORXT(MAXPTS),V(2000),DUM(2000),AG(2000),AGP(2000),
     2 EIGEN(5),EIGLN(5),DREG(40),ZTOT(5),WK1(2000),DUMD(MAXPTS,5),
     3 WK2(2000),WK3(2000),WK4(2000),WK5(2000),WKB(6000),IS(5),
     4 INORM(5),ANORM(5),VQL(2000),VS(2000)

c      DIMENSION R(1000),RAB(1000),NO(40),LO(40),SO(40),
c     1 VIOD(5,1000),VIOU(5,1000),VID(1000),VIU(1000),EV(40),
c     2 CDD(1000),CDC(1000),CDU(1000),DUM2(1000),ZO(40)
cc
c      DIMENSION DRORYT(MAXPTS,5),DRORYP(MAXPTS,5),DUMU(MAXPTS,5),
c     1 DRORXT(MAXPTS),V(1000),DUM(1000),AG(1000),AGP(1000),
c     2 EIGEN(5),EIGLN(5),DREG(40),ZTOT(5),WK1(1000),DUMD(MAXPTS,5),
c     3 WK2(1000),WK3(1000),WK4(1000),WK5(1000),WKB(6000),IS(5),
c     4 INORM(5),ANORM(5),VQL(1000),VS(1000)
c                                                      
      CHARACTER*1  ISPP,OPT,PAS
c      CHARACTER*1 MARK,ISPP
      CHARACTER*2 ICORR,NAMEAT
      CHARACTER*3 IREL
      CHARACTER*4 NICORE,C_VAL(5) 
c      CHARACTER*4 NICORE 
      CHARACTER*10 IRAY(6),ITITLE(7)
      CHARACTER*15 DATAFILE,PSEUDO
c      CHARACTER*15 DATAFILE,PSEUDO,HEADING
      CHARACTER*17 FNAME,GNAME
c
      DATA DUM2/2000*0.D0/
c      DATA DUM2/1000*0.D0/
      DATA ZTOT/5*0.D0/
      DATA EIGEN/5*0.0D0/
      DATA EIGLN/5*0.0D0/
      DATA CDU/2000*0.D0/
c      DATA CDU/1000*0.D0/
      DATA ((DRORYT(I,J),I=1,2001),J=1,5)/MTOT*0.D0/
c      DATA ((DRORYT(I,J),I=1,800),J=1,5)/MTOT*0.D0/
      DATA IS/5*0/
      DATA C_VAL /'log0','log1','log2','log3','log4'/
c
c mmga
cccccc
c mmga
c
      CALL SPACE1
      WRITE(*,*)' Enter pseudopotential file:'
      WRITE(*,*)'(default value = pseudokb.dat)'
      READ(*,6002) PSEUDO
      IF (PSEUDO.eq.' ') PSEUDO='pseudokb.dat'
      CALL SPACE1

      WRITE(*,*)' Enter data file:'
      WRITE(*,*)'(default value = datafile.dat)'
      READ(*,6002) DATAFILE
      IF (DATAFILE.eq.' ') DATAFILE='datafile.dat'
      CALL SPACE1

      WRITE(*,*)' Enter distance (in a.u.) at which the'
      WRITE(*,*)' logaritmic derivatives are calculated'
      READ(*,*)RPOINT
      CALL SPACE2

c      OPEN (UNIT=4,FILE='temp1',STATUS='OLD')
c        READ (4,6002)DATAFILE
c        READ (4,6002)PSEUDO
c      CLOSE (UNIT=4,STATUS='DELETE') 
c      OPEN (UNIT=15,FILE='temp2',STATUS='OLD')
c        READ (15,6002)HEADING
c      CLOSE (UNIT=15,STATUS='DELETE')
c      OPEN (UNIT=16,FILE='temp3',STATUS='OLD')
c        READ (16,6003)MARK
c      CLOSE (UNIT=16,STATUS='DELETE') 
c      OPEN (UNIT=17,FILE='temp4',STATUS='OLD')
c        READ (17,6004)RPOINT
c      CLOSE (UNIT=17,STATUS='DELETE')
 6001 FORMAT(A2)
 6002 FORMAT(A15)
 6003 FORMAT(A1)                     
 6004 FORMAT(F5.3)
c
c mmga
cccccc
c
c  Open and read in data from file DATAFILE
c
      OPEN (UNIT=7,FILE=DATAFILE,STATUS='OLD',
     1 FORM='UNFORMATTED',RECL=8000)
        READ(7)ITYPE,ICORR,ISPP,NR,A,B
        READ(7)(R(I),I=1,NR) 
        READ(7)(RAB(I),I=1,NR)
        READ(7)LMAX,NAMEAT,NORB,NCORE 
        READ(7)(NO(I),I=1,NORB)
        READ(7)(LO(I),I=1,NORB) 
        READ(7)(SO(I),I=1,NORB)
        READ(7)(ZO(I),I=1,NORB)     
        READ(7)ZNUC,DUM(1)
        READ(7)(DUM(I),I=1,NR)
        DO 1,J=1,LMAX
          READ(7)(VIOD(J,I),I=1,NR)
 1      CONTINUE
        DO 2,J=1,LMAX
          READ(7)(VIOU(J,I),I=1,NR)
 2      CONTINUE
        READ(7)(VID(I),I=1,NR)
        READ(7)(VIU(I),I=1,NR)
        READ(7)(EV(I),I=1,NORB)
      CLOSE (UNIT=7)               
c
c  Find eigen values and range for plotting. 
c          
      ZVAL=0.D0
      DO 10,I=NCORE+1,NORB
        ZVAL=ZVAL+ZO(I)
 10   CONTINUE         
cccccc
c mmga
c
      EVHI=10.0D0
      EVLOW=-10.0D0
c      EVHI=1.0D0
c      EVLOW=-2.0D0

      WRITE(*,*) '... CALCULATING LOGARITMIC DERIVATIVES ... '
      CALL SPACE1

c mmga
c mmga You can modify the range of energies changing
c mmga the values of EVLOW and EVHI in this file, but ...
c mmga !!! EVHI minus EVLOW musn't be greater than 20 Ry !!!
c mmga

      NUMEV=0
      DO 11 ENERGY=EVLOW,EVHI,0.01D0
c      DO 11 ENERGY=EVLOW,EVHI+0.01D0,0.02D0 
c
c mmga
cccccc
        NUMEV=NUMEV+1
        DRORXT(NUMEV)=ENERGY
 11   CONTINUE
      NUMTOTAL=NUMEV
c
c  Find the point closest to RPOINT
c  
      DO 12 I=1,NR
        IF (RPOINT .LT. R(I)) THEN
          NPOINT=I
          GOTO 14
        ENDIF
 12   CONTINUE
 14   CONTINUE
c
c  Find ln[] vs ev for the true potential.
c
      DO 15 I=NCORE+1,NORB
        LP=LO(I)+1
        IF (ISPP .EQ. 'r' ) THEN
          IF (SO(I) .LT. 0.1D0) THEN
            DO 35 J=2,NR
              V(J)=VIOD(LP,J)/R(J)+VID(J)
 35         CONTINUE
          ELSE
            DO 40 J=2,NR
              V(J)=VIOU(LP,J)/R(J)+VIU(J)
 40         CONTINUE
          ENDIF
        ELSEIF (ISPP .EQ. 's') THEN 
          IF (SO(I) .LT. 0.1D0 ) THEN
            DO 25 J=2,NR
              V(J)=VIOD(LP,J)/R(J)+VID(J)+LP*(LP-1)/R(J)**2
 25         CONTINUE
          ELSE
            DO 30 J=2,NR
              V(J)=VIOU(LP,J)/R(J)+VIU(J)+LP*(LP-1)/R(J)**2
 30         CONTINUE
          ENDIF
        ELSE
          DO 20 J=2,NR
            V(J)=VIOD(LP,J)/R(J)+VID(J)+LP*(LP-1)/R(J)**2
 20       CONTINUE
        ENDIF 
        NUMEV=0
c mmga
        DO 45 ENERGY=EVLOW,EVHI,0.01D0
c        DO 45 ENERGY=EVLOW,EVHI+0.01D0,0.02D0 
c mmga
          NUMEV=NUMEV+1
          IF (ISPP .EQ. 'r') THEN
            CALL DIFRLE(I,V,NR,R,RAB,NORB,NO,LO,SO,ZNUC,
     1       VIOD,VIOU,VID,VIU,NPOINT,ENERGY,
     2       DUMMY)
            IF (SO(I) .GT. 0.1) THEN
              DRORYT(NUMEV,LP)=DRORYT(NUMEV,LP)
     1         +DUMMY*LP/(2*LO(I)+1)
            ELSE
              DRORYT(NUMEV,LP)=DRORYT(NUMEV,LP)
     1         +DUMMY*LO(I)/(2*LO(I)+1)
            ENDIF
          ELSE 
            CALL DIFNRE(I,V,NR,A,B,R,RAB,NORB,LO,SO,
     1       ZNUC,VIOD,VIOU,VID,VIU,NPOINT,ENERGY,
     2       DUMMY)
            IF (ISPP .EQ. ' ') THEN
              DRORYT(NUMEV,LP)=DUMMY
            ELSE
              IF( SO(I) .GT. 0.1D0) THEN
                DUMU(NUMEV,LP)=DUMMY
              ELSE
                DUMD(NUMEV,LP)=DUMMY
              ENDIF
            ENDIF
          ENDIF
 45     CONTINUE
 15   CONTINUE
c
c  Store eigen values and find R'/R(eigen) 
c  for later transfer to plotit
c
      DO 70 I=NCORE+1,NORB
        LP=LO(I)+1
        IF (ISPP .EQ. 'r' ) THEN
          IF (SO(I) .LT. 0.1D0) THEN
            DO 85 J=2,NR
              V(J)=VIOD(LP,J)/R(J)+VID(J)
 85         CONTINUE
          ELSE
            DO 90 J=2,NR
              V(J)=VIOU(LP,J)/R(J)+VIU(J)
 90         CONTINUE
          ENDIF
        ELSEIF (ISPP .EQ. 's') THEN 
          IF (SO(I) .LT. 0.1D0 ) THEN
            DO 95 J=2,NR
              V(J)=VIOD(LP,J)/R(J)+VID(J)+LP*(LP-1)/R(J)**2
 95         CONTINUE
          ELSE
            DO 100 J=2,NR
              V(J)=VIOU(LP,J)/R(J)+VIU(J)+LP*(LP-1)/R(J)**2
 100        CONTINUE
          ENDIF
        ELSE
          DO 105 J=2,NR
            V(J)=VIOD(LP,J)/R(J)+VID(J)+LP*(LP-1)/R(J)**2
 105      CONTINUE
        ENDIF 
        IF (ISPP .EQ. 'r') THEN
          CALL DIFRLE(I,V,NR,R,RAB,NORB,NO,LO,SO,ZNUC,
     1     VIOD,VIOU,VID,VIU,NPOINT,EV(I),DREG(I))
        ELSE 
          CALL DIFNRE(I,V,NR,A,B,R,RAB,NORB,LO,SO,
     1     ZNUC,VIOD,VIOU,VID,VIU,NPOINT,EV(I),DREG(I))
        ENDIF
 70   CONTINUE
      LMAX2=0
      DO 75 I=NCORE+1,NORB
        LP=LO(I)+1
        IF (LP .GT. LMAX2) LMAX2=LP
        ZTOT(LP)=ZO(I)+ZTOT(LP)
        IF(ZTOT(LP) .NE. 0.D0 .AND. ISPP .NE. 's') THEN
          EIGEN(LP)=ZO(I)*EV(I)+EIGEN(LP)
          EIGLN(LP)=DREG(I)*ZO(I)+EIGLN(LP)
        ELSEIF(ZTOT(LP) .EQ. 0.D0 .AND. ISPP .NE. 's') THEN
          IF (ISPP .EQ. ' ') THEN
            EIGEN(LP)=EV(I)
            EIGLN(LP)=DREG(I)
          ELSEIF(ISPP .EQ. 'r') THEN
            IF (SO(I) .GT. 0.1) THEN
              EIGEN(LP)=EV(I)*(LP)/(2*LO(I)+1)+EIGEN(LP)
              EIGLN(LP)=DREG(I)*(LP)/(2*LO(I)+1)+EIGLN(LP) 
            ELSE
              EIGEN(LP)=EV(I)*LO(I)/(2*LO(I)+1)+EIGEN(LP)
              EIGLN(LP)=DREG(I)*LO(I)/(2*LO(I)+1)+EIGLN(LP)
            ENDIF
          ENDIF
       ELSEIF (ISPP .EQ. 's') THEN
         IF (IS(LP) .EQ. 0) THEN
           IF (ZO(I) .EQ. 0.D0) THEN
             EIGEN(LP)=EV(I)/2
             EIGLN(LP)=DREG(I)/2
             IS(LP)=2
             IF (SO(I) .GT. 0.1D0) THEN
               DO 71 JI=1,NUMTOTAL
                 DRORYT(JI,LP)=DUMU(JI,LP)/2
 71            CONTINUE
             ELSE
               DO 72 JI=1,NUMTOTAL
                 DRORYT(JI,LP)=DUMD(JI,LP)/2
 72            CONTINUE
             ENDIF
           ELSE
             EIGEN(LP)=EV(I)*ZO(I)
             EIGLN(LP)=DREG(I)*ZO(I)
             IS(LP)=1
             IF (SO(I) .GT. 0.1D0) THEN
               DO 73 JI=1,NUMTOTAL
                 DRORYT(JI,LP)=DUMU(JI,LP)*ZO(I)
 73            CONTINUE
             ELSE
               DO 74 JI=1,NUMTOTAL
                 DRORYT(JI,LP)=DUMD(JI,LP)*ZO(I)
 74            CONTINUE
             ENDIF
           ENDIF
         ELSEIF (IS(LP) .EQ. 1) THEN
           EIGEN(LP)=EIGEN(LP)+EV(I)*ZO(I)
           EIGLN(LP)=EIGLN(LP)+DREG(I)*ZO(I)
           IF (SO(I) .GT. 0.1D0) THEN
             DO 76 JI=1,NUMTOTAL
               DRORYT(JI,LP)=(DRORYT(JI,LP)+DUMU(JI,LP)*ZO(I))/ZTOT(LP)
 76          CONTINUE
           ELSE
             DO 77 JI=1,NUMTOTAL
               DRORYT(JI,LP)=(DRORYT(JI,LP)+DUMD(JI,LP)*ZO(I))/ZTOT(LP)
 77          CONTINUE
           ENDIF
         ELSEIF (IS(LP) .EQ. 2) THEN
           IF(ZO(I) .EQ. 0.D0) THEN
             EIGEN(LP)=EIGEN(LP)+EV(I)/2
             EIGLN(LP)=EIGLN(LP)+DREG(I)/2
             IF (SO(I) .GT. 0.1D0) THEN
               DO 78 JI=1,NUMTOTAL
                 DRORYT(JI,LP)=DRORYT(JI,LP)+DUMU(JI,LP)/2
 78            CONTINUE
             ELSE
               DO 79 JI=1,NUMTOTAL
                 DRORYT(JI,LP)=DRORYT(JI,LP)+DUMD(JI,LP)/2
 79            CONTINUE
             ENDIF
           ELSE
             EIGEN(LP)=EV(I)*ZO(I)
             EIGLN(LP)=DREG(I)*ZO(I)
             IF (SO(I) .GT. 0.1D0) THEN
               DO 81 JI=1,NUMTOTAL
                 DRORYT(JI,LP)=DUMU(JI,LP)
 81            CONTINUE
             ELSE
               DO 82 JI=1,NUMTOTAL
                 DRORYT(JI,LP)=DUMD(JI,LP)
 82            CONTINUE
             ENDIF
           ENDIF
         ENDIF
       ENDIF
 75   CONTINUE
      DO 80 I=1,LMAX2
        IF (ZTOT(I) .NE. 0.D0) THEN
          EIGEN(I)=EIGEN(I)/ZTOT(I)
          EIGLN(I)=EIGLN(I)/ZTOT(I)
        ENDIF
 80   CONTINUE
c
c  Open and read in data from file PSEUDO
c
      OPEN (UNIT=8,FILE=PSEUDO,FORM='UNFORMATTED',
     1 RECL=9000)
        READ(8) NAMEAT,ICORR,IREL,NICORE,(IRAY(I),I=1,6),
     1   (ITITLE(I),I=1,7),NORB,NUMU,NR,DUM(1),DUM(1),ZNUC
        NR=NR+1
        READ(8) (R(I),I=2,NR)
        DO 388,I=1,NORB
          READ(8) LO(I),(VIOD(I,J),J=2,NR)
          SO(I)=0.D0
 388    CONTINUE
        DO 390,I=1,NUMU
          READ(8) IDUM,(VIOU(I,J),J=2,NR)
 390    CONTINUE
        READ(8) (CDC(I),I=2,NR)
        READ(8) (CDD(I),I=2,NR)
        READ(8) (VQL(J),J=2,NR)
        READ(8) NORBS
        DO 392 I=1,NORBS
          READ(8) INORM(I),ANORM(I)
 392    CONTINUE
        NVAL=NORB
c
c  Find el-el potential
c
      CLOSE (UNIT=8)
      IF (NICORE .EQ. 'nc  ') THEN
        ITYPE=1
      ELSEIF(NICORE .EQ. 'pcec' .OR. NICORE .EQ. 'fcec') THEN
        ITYPE=2
      ELSE
        ITYPE=3 
      ENDIF
      CALL VELECT(0,1,ICORR,' ',ITYPE-1,NR,R,RAB,ZVAL,
     1 CDD,CDU,CDC,VID,VIU,DUM,WK1,WK2,WK3,WK4,WK5,WKB)
c
c  Set up potentials
c
      DO 50 I=1,NORB
        LP=LO(I)+1
        LLP=LP*(LP-1)
        AG(1) = 0.D0
        AG(NPOINT+1) = 0.D0
        AG(NPOINT+2) = 0.D0
        AG(NPOINT+3) = 0.D0
        AG(NPOINT+4) = 0.D0
        DO 101 J=2,NPOINT
          AG(J) = VIOD(LP,J)*RAB(J)*R(J)
 101    CONTINUE
        DO 55 J=2,NPOINT+4
          VS(J)=VQL(J)+VID(J)+LLP/R(J)**2 
          V(J)=VS(J)
 55     CONTINUE
        CALL DIFKBY(I,V,NR,A,B,R,RAB,NORB,LO,SO,
     1   ZNUC,VIOD,VIOU,VID,VIU,NPOINT,EVLOW,
     2   DUMEV,WK1)
        NUMEV=0
c mmga
        DO 60 ENERGY=EVLOW,EVHI,0.01D0
c        DO 60 ENERGY=EVLOW,EVHI+0.01D0,0.02D0 
c mmga
          NUMEV=NUMEV+1
c mmga
          N_ERROR=0
          ALPHA=10.
c mmga
 125      PROWAV =0.D0
          DO 110  J=1,NPOINT,4
            PROWAV = PROWAV + 7*(AG(J)*WK1(J)+AG(J+4)*WK1(J+4))+
     1       32*(AG(J+1)*WK1(J+1)+AG(J+3)*WK1(J+3))+
     2       12*AG(J+2)*WK1(J+2)
 110      CONTINUE
          PROWAV = 2*PROWAV/45
          DO 130 J=2,NPOINT+4
            AGP(J) = INORM(LP)*PROWAV*VIOD(LP,J)*R(J)/WK1(J)
 130      CONTINUE
          DO 135 J=2,NPOINT+4
            V(J) = VS(J) + AGP(J)
 135      CONTINUE
          CALL DIFKBY(I,V,NR,A,B,R,RAB,NORB,LO,SO,
     1     ZNUC,VIOD,VIOU,VID,VIU,NPOINT,ENERGY,
     2     DRORYP(NUMEV,I),WK2)
          ERR = 0.D0
          DO 140 J=2,NPOINT 
            ERR = MAX(ABS(WK2(J)-WK1(J)),ERR)
 140      CONTINUE 
          IF (ERR .GT. 1.D-4) THEN
            DO 145 J = 2,NPOINT+4
c mmga
              WK1(J) = ((20.-ALPHA)*WK2(J) + 
     .                   ALPHA*WK1(J))/20. + 1.0D-28
c              WK1(J) = (19*WK2(J) + WK1(J))/20 + 1.0D-28
 145        CONTINUE
            N_ERROR=N_ERROR+1
            IF (N_ERROR.EQ.1000) THEN 
             N_ERROR=0
             ALPHA=ALPHA-0.2
              IF (ALPHA.LE.0.) THEN 
               WRITE(*,*) 'CAN NOT REACH SELF-CONSISTENCY AT:',
     .                     ENERGY,' (Ry)'
               STOP
              END IF
            END IF
c mmga
c            print*,lp,energy,err
            GOTO 125
          ENDIF
 60     CONTINUE
 50   CONTINUE
cccccc
c mmga
c
      WRITE(*,*) '... WRITING TRUE DERIVATIVES IN log_true.gp ...'
      CALL SPACE1

      OPEN (UNIT=18,FILE='log_true.gp',STATUS='UNKNOWN',
     .      FORM='FORMATTED')
      rewind(18)
       WRITE(18,2003) (C_VAL(I),I=1,NVAL)
       DO 2000 I=1,NUMEV
        WRITE(18,2002) DRORXT(I),(DRORYT(I,k),k=1,NVAL)
 2000  CONTINUE
      CLOSE(18)

      WRITE(*,*) '... WRITING PSEUDO DERIVATIVES IN log_ppkb.gp ...'
      CALL SPACE1

      OPEN (UNIT=18,FILE='log_ppkb.gp',STATUS='UNKNOWN',
     .      FORM='FORMATTED')
       rewind(18)
       WRITE(18,2003) (C_VAL(I),I=1,NVAL)
       DO 2001 I=1,NUMEV
        WRITE(18,2002) DRORXT(I),(DRORYP(I,k),k=1,NVAL)
 2001  CONTINUE
      CLOSE(18)


      CALL SPACE2
      WRITE(*,*)'!!! QUIT OTHER GNUPLOTS      !!!'
      WRITE(*,*)'!!! THEN ANSWER THE QUESTION !!!'


      CALL SPACE2


      FNAME='command_logkb.gp'
      GNAME='command_logkb.ps'

      OPEN (UNIT=18,FILE=FNAME,STATUS='UNKNOWN',FORM='FORMATTED')

      X_MIN = EVLOW
      X_MAX = EVHI
      Y_MIN =-20.
      Y_MAX = 20.

      DO 2004 I=1,NVAL

 2014 WRITE(18,2006)

 2015 WRITE(18,2005) RPOINT

      WRITE(18,2007)Y_MIN,Y_MAX,X_MIN,X_MAX

      WRITE(18,2008)I+1,I-1
      WRITE(18,2009)I+1,I-1

      WRITE(18,2011)

2004  CONTINUE

      CLOSE(18)


c mmga
c mmga *** FORMATS ***
c mmga

 2002 FORMAT(1pe13.6,20(2x,e13.6))
 2003 FORMAT('#',4x,'E(Ry)',10x,5(a4,11x))

 2006 FORMAT('set term x11')

 2005 FORMAT(36Hset title 'LOGARITMIC DERIVATIVES AT,
     .f5.2,8H (a.u.)',/,
     .'set nolabel',/,
     .'set key',/,
     .'set notime',/,
     .'set noparametric',/,
     .'set size',/,
     .20Hset format xy '%.1f',/,
     .'set nogrid',/,
     .24Hset xlabel 'Energy (Ry)',/,
     .39Hset ylabel 'Log. Derivatives (1/au)' -2)

 2007 FORMAT('set yrange [',f7.1,':',f7.1,']',/,
     .'set xrange [',f7.1,':',f7.1,']')

 2008 FORMAT('plot ',13H'log_true.gp',' using 1:',i1,
     .        26H ti 'true functions for l=,i1,1H',' w points 1,',1H\)
 2009 FORMAT(13H'log_ppkb.gp',' using 1:',i1,
     .        34H ti 'K & B pseudo functions for l=,i1,1H',
     .          ' w points 2 ')

 2011 FORMAT('pause -1')

 2018 FORMAT('... SAVING IN FILE ',a17,' ...')

 2019 FORMAT('set term postscript color ',11H'Helvetica',/,
     .12Hset output ',a17,1H')

c
c mmga
cccccc

c
c  Inputting is finished, now set up plot.
c      open(unit=12,file='kbdata.dat',status='new',
c     1 form='formatted',RECL=10000)
c      mpoint=201
c      write(12,*) mpoint
c      do 1012 i=1,3
c        do 1014 j=1,mpoint
c          write(12,*)drorxt(j),droryt(j,i)
c 1014   continue
c 1012 continue
c      do 1013 i=1,3
c        do 1015 j=1,mpoint
c          write(12,*)drorxt(j),droryp(j,i)
c 1015   continue
c 1013 continue
c      close(unit=12)
c
c Specify the Apollo screen for output.  Hit return
c to clear the plot
c
c 300  IF (MARK .EQ. 'A' .OR. MARK .EQ. 'B') THEN
c        CALL APOLLO('drct')
c        CALL PLOTIT(HEADING,EVLOW,EVHI,R(NPOINT),NVAL,
c     1   NUMEV,DRORYT,DRORXT,DRORYP,EIGEN,EIGLN,1)
c      ENDIF
c
c Specify POSTSCRIPT output for the laser printer.
c Resulting output file std00001.dat is printed with
c prf -s //whereever -transparent std00001.dat
c
c CAUTION: This file collects all plots sent to it 
c since it was last removed.  Before you run this 
c program to get hard copy I suggest you remove 
c std00001.dat or you may get lots of old plots
c                     
c      IF (MARK .EQ. 'P' .OR. MARK .EQ. 'B') THEN
c        CALL POSTSCRIPT(7.99,10.78,0.0139)
c        CALL PLOTIT(HEADING,EVLOW,EVHI,R(NPOINT),NVAL,
c     1   NUMEV,DRORYT,DRORXT,DRORYP,EIGEN,EIGLN,0)
c      ENDIF
c
      STOP
      END



c      SUBROUTINE PLOTIT(HEADING,EVLOW,EVHI,RPOINT,NVAL,NUMEV,
c     1 DRORYT,DRORXT,DRORYP,EIGEN,EIGLN,IP)
cc
c      IMPLICIT REAL*8 (B-H,O-Z)
cc
c      PARAMETER(MAXPTS=800)
cc
c      DIMENSION DRORYT(MAXPTS,5),DRORXT(MAXPTS),DRORYP(MAXPTS,5),
c     1 ADUMMY(MAXPTS),ADUMMX(MAXPTS),EIGEN(5),EIGLN(5)
cc                        
c      CHARACTER*3 EIGPL(5)
c      CHARACTER*15 HEADING
c      CHARACTER*28 THEAD
cc                          
c      DATA EIGPL/'l=0','l=1','l=2','l=3','l=4'/
cc
c      WRITE(THEAD,100)HEADING,'Rc=',RPOINT,' (au)'
c 100  FORMAT(A15,A3,F5.3,A5)
cc
c      CALL RESET('all')
cc
cc  Plot 
c
c      CALL SETDEV(0,0)
c      CALL NOBRDR
c      CALL AREA2D(6.0,8.0)
c      CALL HEADIN('d(Ln(R))/dr vs. Energy - K & B',-30,2.,2)
c      CALL COMPLX
c      CALL HEADIN(THEAD,28,1.5,2)
c      CALL FRAME
c      CALL NOCHEK
c      CALL GRACE(0.0)
c      CALL YNAME('Logarithmic derivatives (1/au)',30)
c      CALL XNAME('Energy (Ry)',11)
c      CALL YAXANG(0.0)
c      ALOW=SNGL(EVLOW)
c      AHI=SNGL(EVHI)
c      CALL GRAF(ALOW,1.0,AHI,-4.0,1.0,4.0)
c      CALL SPLINE
c      DO 30 K=1,NVAL
c        ADUMMY(0)=SNGL(DRORYT(1,K))
c        IFLAG=0
c        J=1
c        DO 5,I=1,NUMEV
c          ADUMMX(J)=SNGL(DRORXT(I))
c          ADUMMY(J)=SNGL(DRORYT(I,K))
c          IF (IFLAG .EQ. -1) THEN
c            IF (ADUMMY(1) .GT. 1.0) THEN
c              IFLAG=0
c              J=2    
c            ENDIF
c          ELSE
c            IF (ADUMMY(J) .GT. 4.0 .AND. ADUMMY(J-1) .GT. 4.0) THEN
c              ADUMMY(J-1)=ADUMMY(J)
c              ADUMMX(J-1)=ADUMMX(J)
c            ELSE
c              IF( ADUMMY(J) .LT. -4.0 ) THEN
c                IF (J .GT. 5) CALL CURVE(ADUMMX,ADUMMY,J,0)
c                IFLAG=-1 
c                J=1 
c              ELSE
c                J=J+1
c              ENDIF
c            ENDIF
c          ENDIF
c 5      CONTINUE
c        IF (J .GT. 6) CALL CURVE(ADUMMX,ADUMMY,J-1,0)
c 30   CONTINUE
c      CALL DASH
c      DO 40 K=1,NVAL
c        ADUMMY(0)=SNGL(DRORYP(1,K))
c        IFLAG=0
c        J=1
c        DO 45,I=1,NUMEV
c          ADUMMX(J)=SNGL(DRORXT(I))
c          ADUMMY(J)=SNGL(DRORYP(I,K))
c          IF (IFLAG .EQ. -1) THEN
c            IF (ADUMMY(1) .GT. 1.0) THEN
c              IFLAG=0
c              J=2    
c            ENDIF
c          ELSE
c            IF (ADUMMY(J) .GT. 4.0 .AND. ADUMMY(J-1) .GT. 4.0) THEN
c              ADUMMY(J-1)=ADUMMY(J)
c              ADUMMX(J-1)=ADUMMX(J)
c            ELSE
c              IF( ADUMMY(J) .LT. -4.0 ) THEN
c                IF (J .GT. 5) CALL CURVE(ADUMMX,ADUMMY,J,0)
c                IFLAG=-1 
c                J=1 
c              ELSE
c                J=J+1
c              ENDIF
c            ENDIF
c          ENDIF
c 45     CONTINUE
c        IF (J .GT. 6) CALL CURVE(ADUMMX,ADUMMY,J-1,0)
c 40   CONTINUE
c      CALL BASALF('L/CSCRIPT')
c      CALL MSHIFT(0.05,0.02)
c      DO 50 I=1,NVAL
c        CALL RLMESS(EIGPL(I),3,SNGL(EIGEN(I)),SNGL(EIGLN(I)))
c 50   CONTINUE      
c      CALL HEIGHT(0.38)
c      CALL MSHIFT(-0.08,-0.45)
c      CALL BASALF('MATHEMMATIC')
c      DO 55 I=1,NVAL
c        CALL RLMESS('U',1,SNGL(EIGEN(I)),SNGL(EIGLN(I)))
c 55   CONTINUE
c      CALL ENDPL(0)
cc
c      CALL DONEPL
c      RETURN
c      END
c
      SUBROUTINE DIFNRE(IORB,V,NR,A,B,R,RAB,NORB,LO,SO,
     1 ZNUC,VIOD,VIOU,VID,VIU,NPOINT,ENV,DROR)
c
      IMPLICIT REAL*8 (A-H,O-Z)
c
c  DIFNRL INTEGRATES THE SCHROEDINGER EQUATION
c  IT FINDS THE [ R'/R ] VS ENERGY
c
cccccc
c mmga
c
      DIMENSION V(NR),R(NR),RAB(NR),LO(NORB),
     1 SO(NORB),VIOD(4,NR),VIOU(4,NR),VID(NR),VIU(NR),
     2 A1(7),RR(7),COE(7),AR(2000),BR(2000)

c      DIMENSION V(NR),R(NR),RAB(NR),LO(NORB),
c     1 SO(NORB),VIOD(4,NR),VIOU(4,NR),VID(NR),VIU(NR),
c     2 A1(7),RR(7),COE(7),AR(1000),BR(1000)
c
c mmga
cccccc
c
c  INTEGRATION COEFFICIENTS
c
      ABC1 = 1901.D0/720.D0
      ABC2 = -1387.D0/360.D0
      ABC3 = 109.D0/30.D0
      ABC4 = -637.D0/360.D0
      ABC5 = 251.D0/720.D0
      AMC0 = 251.D0/720.D0
      AMC1 = 323.D0/360.D0
      AMC2 = -11.D0/30.D0
      AMC3 = 53.D0/360.D0
      AMC4 = -19.D0/720.D0
c
c  DETERMINE EFFECTIVE CHARGE AND VZERO FOR STARTUP OF
c  OUTWARD INTEGRATION
c  AR = R**(L+1) * (1 + AA R + BB R**2 + ... )
c  AA = -ZNUC / LP     BB = (-2 ZNUC AA + V(0) - E)/(4 L + 6)
c
      LP = LO(IORB)+1
      ZEFF = 0.D0
      IF (SO(IORB) .LT. 0.1D0 .AND. VIOD(LP,2) .LT. -0.1D0) ZEFF=ZNUC
      IF (SO(IORB) .GT. 0.1D0 .AND. VIOU(LP,2) .LT. -0.1D0) ZEFF=ZNUC
      AA = -ZEFF/LP
      VZERO = -2*ZEFF*AA
      IF (SO(IORB) .LT. 0.1D0 .AND. ZEFF .EQ. 0.D0)
     1 VZERO=VIOD(LP,2)/R(2)
      IF (SO(IORB) .GT. 0.1D0 .AND. ZEFF .EQ. 0.D0)
     1 VZERO=VIOU(LP,2)/R(2)
      IF (SO(IORB) .LT. 0.1D0) VZERO=VZERO+VID(2)
      IF (SO(IORB) .GT. 0.1D0) VZERO=VZERO+VIU(2)
      VAR0 = 0.D0
      IF (LO(IORB) .EQ. 0) VAR0=-2*ZEFF
      IF (LO(IORB) .EQ. 1) VAR0=2.D0
c
c  OUTWARD INTEGRATION FROM 1 TO NPOINT+3
c
      BB = (VZERO-ENV)/(4*LP+2)
      AR(1) = 0.D0
      BR(1) = 0.D0
      IF (LO(IORB) .EQ. 0) BR(1) = B*A
      DO 35 J=2,5
        AR(J) = R(J)**LP * (1+(AA+BB*R(J))*R(J))
        BR(J) = RAB(J)*R(J)**LO(IORB)
     1   * (LP+(AA*(LP+1)+BB*(LP+2)*R(J))*R(J))
 35   CONTINUE
      FA5 = BR(1)
      FB5 = B*BR(1) + RAB(1)*RAB(1)*VAR0
      FA4 = BR(2)
      FB4 = B*BR(2) + RAB(2)*RAB(2)*(V(2)-ENV)*AR(2)
      FA3 = BR(3)
      FB3 = B*BR(3) + RAB(3)*RAB(3)*(V(3)-ENV)*AR(3)
      FA2 = BR(4)
      FB2 = B*BR(4) + RAB(4)*RAB(4)*(V(4)-ENV)*AR(4)
      FA1 = BR(5)
      FB1 = B*BR(5) + RAB(5)*RAB(5)*(V(5)-ENV)*AR(5)
c
c  INTERGRATION LOOP
c
      DO 40 J=6,NPOINT+3
c
c  PREDICTOR (ADAMS-BASHFORTH)
c
        ARP = AR(J-1) + ABC1*FA1+ABC2*FA2+ABC3*FA3+ABC4*FA4+ABC5*FA5
        BRP = BR(J-1) + ABC1*FB1+ABC2*FB2+ABC3*FB3+ABC4*FB4+ABC5*FB5
        FA0 = BRP
        FB0 = B*BRP + RAB(J)*RAB(J)*(V(J)-ENV)*ARP
c
c  CORRECTOR (ADAMS-MOULTON)
c
        ARC = AR(J-1) + AMC0*FA0+AMC1*FA1+AMC2*FA2+AMC3*FA3+AMC4*FA4
        BRC = BR(J-1) + AMC0*FB0+AMC1*FB1+AMC2*FB2+AMC3*FB3+AMC4*FB4
        FA5 = FA4
        FB5 = FB4
        FA4 = FA3
        FB4 = FB3
        FA3 = FA2
        FB3 = FB2
        FA2 = FA1
        FB2 = FB1
        FA1 = BRC
        FB1 = B*BRC + RAB(J)*RAB(J)*(V(J)-ENV)*ARC
        AR(J) = ARC + AMC0*(FA1-FA0)
        BR(J) = BRC + AMC0*(FB1-FB0)
        FA1 = BR(J)
        FB1 = B*BR(J) + RAB(J)*RAB(J)*(V(J)-ENV)*AR(J)
 40   CONTINUE
c
c  END OUTWARD INTEGRATION
c
      A1(1)=AR(NPOINT-3)/R(NPOINT-3)
      A1(2)=AR(NPOINT-2)/R(NPOINT-2)
      A1(3)=AR(NPOINT-1)/R(NPOINT-1)
      A1(4)=AR(NPOINT)/R(NPOINT)
      A1(5)=AR(NPOINT+1)/R(NPOINT+1)
      A1(6)=AR(NPOINT+2)/R(NPOINT+2)
      A1(7)=AR(NPOINT+3)/R(NPOINT+3)
      RR(1)=R(NPOINT-3)-R(NPOINT)
      RR(2)=R(NPOINT-2)-R(NPOINT)
      RR(3)=R(NPOINT-1)-R(NPOINT)
      RR(4)=0.D0
      RR(5)=R(NPOINT+1)-R(NPOINT)
      RR(6)=R(NPOINT+2)-R(NPOINT)
      RR(7)=R(NPOINT+3)-R(NPOINT)
      CALL POLCOE(RR,A1,7,COE)
      IF (ABS(A1(4)) .LT. 1.D-10) THEN
        DROR=1.D10
      ELSE
        DROR=COE(2)/A1(4)
      ENDIF
      RETURN
      END
       SUBROUTINE DIFRLE(IORB,V,NR,R,RAB,NORB,NO,LO,SO,
     1  ZNUC,VIOD,VIOU,VID,VIU,NPOINT,EVE,DROR)
c
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
c  DIFREL INTEGRATES THE RELATIVISTIC DIRAC EQUATION
c  IT FINDS THE [ R'/R ] vs ENERGY.
c
cccccc
c mmga
c
      DIMENSION V(NR),AR(2000),BR(2000),R(NR),RAB(NR),NO(NORB),
     1 LO(NORB),SO(NORB),VIOD(4,NR),VIOU(4,NR),VID(NR),VIU(NR),
     2 RR(7),AA(7),COE(7)

c      DIMENSION V(NR),AR(1000),BR(1000),R(NR),RAB(NR),NO(NORB),
c     1 LO(NORB),SO(NORB),VIOD(4,NR),VIOU(4,NR),VID(NR),VIU(NR),
c     2 RR(7),AA(7),COE(7)
c
c mmga
cccccc
c
       AI = 2*137.04D0
       AZ = ZNUC/(2*AI)
       KA = LO(IORB)+1
       IF (SO(IORB) .LT. 0.1D0 .AND. LO(IORB) .NE. 0) KA=-LO(IORB)
c
c  INTEGRATION COEFFICIENTS
c
       ABC1 = 1901.D0/720.D0
       ABC2 = -1387.D0/360.D0
       ABC3 = 109.D0/30.D0
       ABC4 = -637.D0/360.D0
       ABC5 = 251.D0/720.D0
       AMC0 = 251.D0/720.D0
       AMC1 = 323.D0/360.D0
       AMC2 = -11.D0/30.D0
       AMC3 = 53.D0/360.D0
       AMC4 = -19.D0/720.D0
c
c  DETERMINE EFFECTIVE CHARGE AND VZERO FOR STARTUP OF THE
c  OUTWARD INTEGRATION
c   AR = R**S * (1  + A1 R + A2 R**2 + ... )
c   BR = R**S * (B0 + B1 R + B2 R**2 + ... )
c   S = SQRT (KA**2 - AZ**2)    B0 = - AZ / (S + KA)
c   AN = (AZ (V0 - E) A(N-1) - (S + N + KA) (V0 - E - AI**2) B(N-1))
c    / (N AI (2 S + N))
c   BN = ((V0 - E) A(N-1) - 2 ZNUC AN ) / ( AI (S + N + KA))
c
       S = SQRT(KA*KA-AZ*AZ)
       IF (KA .GT. 0) B0 = -AZ/(S+KA)
       IF (KA .LE. 0) B0 = (S-KA)/AZ
       IF (SO(IORB) .LT. 0.1D0) VZERO=VID(2)
       IF (SO(IORB) .GT. 0.1D0) VZERO=VIU(2)
c
c  OUTWARD INTEGRATION FROM 1 TO NPOINT+3
c
       A1 = (AZ*(VZERO-EVE)-(S+1+KA)*(VZERO-EVE-AI**2)*B0)/ (AI*(2*S+1))
       B1 = ((VZERO-EVE)-2*ZNUC*A1) / (AI*(S+1+KA))
       A2 = (AZ*(VZERO-EVE)*A1-(S+2+KA)*(VZERO-EVE-AI**2)*B1)
     1    / (2*AI*(2*S+2))
       B2 = ((VZERO-EVE)*A1-2*ZNUC*A2) / (AI*(S+2+KA))
       AR(1) = 0.D0
       BR(1) = 0.D0
       DO 35 J=2,5
         AR(J) = R(J)**S * (1 +(A1+A2*R(J))*R(J))
         BR(J) = R(J)**S * (B0+(B1+B2*R(J))*R(J))
 35    CONTINUE
       FA5 = 0.D0
       FB5 = 0.D0
       FA4 = RAB(2)*(+KA*AR(2)/R(2)+(EVE-V(2)+AI*AI)*BR(2)/AI)
       FB4 = RAB(2)*(-KA*BR(2)/R(2)-(EVE-V(2))*AR(2)/AI)
       FA3 = RAB(3)*(+KA*AR(3)/R(3)+(EVE-V(3)+AI*AI)*BR(3)/AI)
       FB3 = RAB(3)*(-KA*BR(3)/R(3)-(EVE-V(3))*AR(3)/AI)
       FA2 = RAB(4)*(+KA*AR(4)/R(4)+(EVE-V(4)+AI*AI)*BR(4)/AI)
       FB2 = RAB(4)*(-KA*BR(4)/R(4)-(EVE-V(4))*AR(4)/AI)
       FA1 = RAB(5)*(+KA*AR(5)/R(5)+(EVE-V(5)+AI*AI)*BR(5)/AI)
       FB1 = RAB(5)*(-KA*BR(5)/R(5)-(EVE-V(5))*AR(5)/AI)
c
c  INTERGRATION LOOP
c
       DO 40 J=6,NPOINT+3
c
c  PREDICTOR (ADAMS-BASHFORTH)
c
         ARP = AR(J-1) + ABC1*FA1+ABC2*FA2+ABC3*FA3+ABC4*FA4+ABC5*FA5
         BRP = BR(J-1) + ABC1*FB1+ABC2*FB2+ABC3*FB3+ABC4*FB4+ABC5*FB5
         FA0 = RAB(J)*(+KA*ARP/R(J)+(EVE-V(J)+AI*AI)*BRP/AI)
         FB0 = RAB(J)*(-KA*BRP/R(J)-(EVE-V(J))*ARP/AI)
c
c  CORRECTOR (ADAMS-MOULTON)
c
         ARC = AR(J-1) + AMC0*FA0+AMC1*FA1+AMC2*FA2+AMC3*FA3+AMC4*FA4
         BRC = BR(J-1) + AMC0*FB0+AMC1*FB1+AMC2*FB2+AMC3*FB3+AMC4*FB4
         FA5 = FA4
         FB5 = FB4
         FA4 = FA3
         FB4 = FB3
         FA3 = FA2
         FB3 = FB2
         FA2 = FA1
         FB2 = FB1
         FA1 = RAB(J)*(+KA*ARC/R(J)+(EVE-V(J)+AI*AI)*BRC/AI)
         FB1 = RAB(J)*(-KA*BRC/R(J)-(EVE-V(J))*ARC/AI)
         AR(J) = ARC + AMC0*(FA1-FA0)
         BR(J) = BRC + AMC0*(FB1-FB0)
         FA1 = RAB(J)*(+KA*AR(J)/R(J)+(EVE-V(J)+AI*AI)*BR(J)/AI)
         FB1 = RAB(J)*(-KA*BR(J)/R(J)-(EVE-V(J))*AR(J)/AI)
 40    CONTINUE
c
c  END OUTWARD INTEGRATION
c
       AA(1)=SQRT(AR(NPOINT-3)**2+BR(NPOINT-3)**2)/R(NPOINT-3)
       AA(2)=SQRT(AR(NPOINT-2)**2+BR(NPOINT-2)**2)/R(NPOINT-2)
       AA(3)=SQRT(AR(NPOINT-1)**2+BR(NPOINT-1)**2)/R(NPOINT-1)
       AA(4)=SQRT(AR(NPOINT)**2+BR(NPOINT)**2)/R(NPOINT)
       AA(5)=SQRT(AR(NPOINT+1)**2+BR(NPOINT+1)**2)/R(NPOINT+1)
       AA(6)=SQRT(AR(NPOINT+2)**2+BR(NPOINT+2)**2)/R(NPOINT+2)
       AA(7)=SQRT(AR(NPOINT+3)**2+BR(NPOINT+3)**2)/R(NPOINT+3)
       RR(1)=R(NPOINT-3)-R(NPOINT)
       RR(2)=R(NPOINT-2)-R(NPOINT)
       RR(3)=R(NPOINT-1)-R(NPOINT)
       RR(4)=0.D0
       RR(5)=R(NPOINT+1)-R(NPOINT)
       RR(6)=R(NPOINT+2)-R(NPOINT)
       RR(7)=R(NPOINT+3)-R(NPOINT)
       CALL POLCOE(RR,AA,7,COE)
       IF (AA(4) .LT. 1.D-10) THEN
         DROR=1.D10
       ELSE
         DROR=COE(2)/AA(4)
       ENDIF
       RETURN
       END
c
c
c
      SUBROUTINE DIFKBY(IORB,V,NR,A,B,R,RAB,NORB,LO,SO,
     1 ZNUC,VIOD,VIOU,VID,VIU,NPOINT,ENV,DROR,AR)
c
      IMPLICIT REAL*8 (A-H,O-Z)
c
c  DIFKBY INTEGRATES THE SCHROEDINGER EQUATION
c  IT FINDS THE [ R'/R ] VS ENERGY AND RETURNS
c  THE WAVEFUNCTION NORMALIZED TO ONE OUT TO RC(NPOINTS +4)
c
cccccc
c mmga
c
      DIMENSION V(NR),R(NR),RAB(NR),LO(NORB),
     1 SO(NORB),VIOD(4,NR),VIOU(4,NR),VID(NR),VIU(NR),
     2 A1(7),RR(7),COE(7),AR(2000),BR(2000),TEMP(2000)

c      DIMENSION V(NR),R(NR),RAB(NR),LO(NORB),
c     1 SO(NORB),VIOD(4,NR),VIOU(4,NR),VID(NR),VIU(NR),
c     2 A1(7),RR(7),COE(7),AR(1000),BR(1000),TEMP(1000)
c
c mmga
cccccc
c
c  INTEGRATION COEFFICIENTS
c
      ABC1 = 1901.D0/720.D0
      ABC2 = -1387.D0/360.D0
      ABC3 = 109.D0/30.D0
      ABC4 = -637.D0/360.D0
      ABC5 = 251.D0/720.D0
      AMC0 = 251.D0/720.D0
      AMC1 = 323.D0/360.D0
      AMC2 = -11.D0/30.D0
      AMC3 = 53.D0/360.D0
      AMC4 = -19.D0/720.D0
c
c  DETERMINE EFFECTIVE CHARGE AND VZERO FOR STARTUP OF
c  OUTWARD INTEGRATION
c  AR = R**(L+1) * (1 + AA R + BB R**2 + ... )
c  AA = -ZNUC / LP     BB = (-2 ZNUC AA + V(0) - E)/(4 L + 6)
c
      LP = LO(IORB)+1
      VZERO=V(2)-LP*(LP-1)/R(2)**2
      VAR0 = 0.D0
      IF (LO(IORB) .EQ. 1) VAR0=2.D0
c
c  OUTWARD INTEGRATION FROM 1 TO NPOINT+4
c
      BB = (VZERO-ENV)/(4*LP+2)
      AR(1) = 0.D0
      BR(1) = 0.D0
      IF (LO(IORB) .EQ. 0) BR(1) = B*A
      DO 35 J=2,5
        AR(J) = R(J)**LP * (1+(BB*R(J))*R(J))
        BR(J) = RAB(J)*R(J)**LO(IORB)
     1   * (LP+(BB*(LP+2)*R(J))*R(J))
 35   CONTINUE
      FA5 = BR(1)
      FB5 = B*BR(1) + RAB(1)*RAB(1)*VAR0
      FA4 = BR(2)
      FB4 = B*BR(2) + RAB(2)*RAB(2)*(V(2)-ENV)*AR(2)
      FA3 = BR(3)
      FB3 = B*BR(3) + RAB(3)*RAB(3)*(V(3)-ENV)*AR(3)
      FA2 = BR(4)
      FB2 = B*BR(4) + RAB(4)*RAB(4)*(V(4)-ENV)*AR(4)
      FA1 = BR(5)
      FB1 = B*BR(5) + RAB(5)*RAB(5)*(V(5)-ENV)*AR(5)
c
c  INTERGRATION LOOP
c
      DO 40 J=6,NPOINT+4
c
c  PREDICTOR (ADAMS-BASHFORTH)
c
        ARP = AR(J-1) + ABC1*FA1+ABC2*FA2+ABC3*FA3+ABC4*FA4+ABC5*FA5
        BRP = BR(J-1) + ABC1*FB1+ABC2*FB2+ABC3*FB3+ABC4*FB4+ABC5*FB5
        FA0 = BRP
        FB0 = B*BRP + RAB(J)*RAB(J)*(V(J)-ENV)*ARP
c
c  CORRECTOR (ADAMS-MOULTON)
c
        ARC = AR(J-1) + AMC0*FA0+AMC1*FA1+AMC2*FA2+AMC3*FA3+AMC4*FA4
        BRC = BR(J-1) + AMC0*FB0+AMC1*FB1+AMC2*FB2+AMC3*FB3+AMC4*FB4
        FA5 = FA4
        FB5 = FB4
        FA4 = FA3
        FB4 = FB3
        FA3 = FA2
        FB3 = FB2
        FA2 = FA1
        FB2 = FB1
        FA1 = BRC
        FB1 = B*BRC + RAB(J)*RAB(J)*(V(J)-ENV)*ARC
        AR(J) = ARC + AMC0*(FA1-FA0)
        BR(J) = BRC + AMC0*(FB1-FB0)
        FA1 = BR(J)
        FB1 = B*BR(J) + RAB(J)*RAB(J)*(V(J)-ENV)*AR(J)
 40   CONTINUE
c
c  END OUTWARD INTEGRATION
c
      AMAXW = 0.D0
      DO 61 J=1,NPOINT+4
        AMAXW = MAX(ABS(AR(J)),AMAXW)
 61   CONTINUE
      DO 62 J=1,NPOINT+4 
        AR(J) = AR(J)/AMAXW
 62   CONTINUE            
c
c   Normalize AR to one. R(0) to Rc(Npoints+3)
c
      DO 60 J=1,NPOINT+4
        TEMP(J) = RAB(J)*AR(J)*AR(J)
 60   CONTINUE
      ARNORM = 0.D0
      DO 65 J=1,NPOINT,4
        ARNORM = ARNORM + 7*(TEMP(J)+TEMP(J+4)) +
     1   32*(TEMP(J+1)+TEMP(J+3)) + 12*TEMP(J+2)
 65   CONTINUE
      ARNORM = 2 * ARNORM / 45
      ARNORM = SQRT(ARNORM)
      DO 70 J=1,NPOINT+4
        AR(J) = AR(J)/ARNORM
 70   CONTINUE
c
c  Find Phi'/Phi
c
      A1(1)=AR(NPOINT-3)/R(NPOINT-3)
      A1(2)=AR(NPOINT-2)/R(NPOINT-2)
      A1(3)=AR(NPOINT-1)/R(NPOINT-1)
      A1(4)=AR(NPOINT)/R(NPOINT)
      A1(5)=AR(NPOINT+1)/R(NPOINT+1)
      A1(6)=AR(NPOINT+2)/R(NPOINT+2)
      A1(7)=AR(NPOINT+3)/R(NPOINT+3)
      RR(1)=R(NPOINT-3)-R(NPOINT)
      RR(2)=R(NPOINT-2)-R(NPOINT)
      RR(3)=R(NPOINT-1)-R(NPOINT)
      RR(4)=0.D0
      RR(5)=R(NPOINT+1)-R(NPOINT)
      RR(6)=R(NPOINT+2)-R(NPOINT)
      RR(7)=R(NPOINT+3)-R(NPOINT)
      CALL POLCOE(RR,A1,7,COE)
      IF (ABS(A1(4)) .LT. 1.D-10) THEN
        DROR=1.D10
      ELSE
        DROR=COE(2)/A1(4)
      ENDIF
      RETURN
      END
c
c      subroutine postscript(val1, val2, val3) 
cc
cc This subroutine must be included to have postscript written
cc to a file.  (See /usr/local/disspla/dis10p0/doc/pscript)
cc
c      dimension ibuf(16)
c      ibuf(1)=16
c      call iomgr(ibuf,-101)
c      ibuf(1)=5
c      call iomgr(ibuf,-102)
c      ibuf(1)=0
c      call iomgr(ibuf,-104)
c      ibuf(1)=0
c      call iomgr(ibuf,-105)
c      call pscrpt(val1, val2, val3)
c      return
c      end

cccccc
c mmga
c
c mmga
c
        SUBROUTINE SPACE1

        WRITE(*,*)

        RETURN
        END
c
c mmga
c
        SUBROUTINE SPACE2

        WRITE(*,*)
        WRITE(*,*)

        RETURN
        END
c
c mmga
cccccc
