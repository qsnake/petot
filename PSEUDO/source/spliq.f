       SUBROUTINE SPLIQ(X,Y,YP,YPP,N,XLO,XUP,NUP,ANS,IERR)
C 
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  
C  NJTJ
C  ###  CRAY CONVERSIONS  
C  ###    1)Comment out implicit double precision.
C  ###  CRAY CONVERSIONS  
C  NJTJ
C
       DIMENSION X(N),Y(N),YP(N),YPP(N),XUP(NUP),ANS(NUP)
C
C     SANDIA MATHEMATICAL PROGRAM LIBRARY
C     APPLIED MATHEMATICS DIVISION 2613
C     SANDIA LABORATORIES
C     ALBUQUERQUE, NEW MEXICO  87185
C     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C                    ISSUED BY SANDIA LABORATORIES
C  *                   A PRIME CONTRACTOR TO THE
C  *                UNITED STATES DEPARTMENT OF ENERGY
C  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * *
C  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE
C  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE
C  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,
C  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES
C  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL
C  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR
C  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS
C  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
C  * OWNED RIGHTS.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS
C  * PART IS SAND77-1441.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     THIS ROUTINE WAS WRITTEN BY M. K. GORDON
C
C     ABSTRACT
C
C     SUBROUTINE SPLIQ INTEGRATES A CUBIC SPLINE (GENERATED BY
C     SPLIFT, SMOO, ETC.) ON THE INTERVALS (XLO,XUP(I)), WHERE XUP
C     IS A SEQUENCE OF UPPER LIMITS ON THE INTERVALS OF INTEGRATION.
C     THE ONLY RESTRICTIONS ON XLO AND XUP(*) ARE
C                XLO .LT. XUP(1),
C                XUP(I) .LE. XUP(I+1)   FOR EACH I .
C     ENDPOINTS BEYOND THE SPAN OF ABSCISSAS ARE ALLOWED.
C     THE SPLINE OVER THE INTERVAL (X(I),X(I+1)) IS REGARDED
C     AS A CUBIC POLYNOMIAL EXPANDED ABOUT X(I) AND IS INTEGRATED
C     ANALYTICALLY.
C
C     DESCRIPTION OF ARGUMENTS
C         THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST
C         E.G.  X(N), Y(N), YP(N), YPP(N), XUP(NUP), ANS(NUP)
C
C      --INPUT--
C
C        X    - ARRAY OF ABSCISSAS (IN INCREASING ORDER) THAT DEFINE TH
C               SPLINE.  USUALLY X IS THE SAME AS X IN SPLIFT OR SMOO.
C        Y    - ARRAY OF ORDINATES THAT DEFINE THE SPLINE.  USUALLY Y I
C               THE SAME AS Y IN SPLIFT OR AS R IN SMOO.
C        YP   - ARRAY OF FIRST DERIVATIVES OF THE SPLINE AT ABSCISSAS.
C               USUALLY YP IS THE SAME AS YP IN SPLIFT OR R1 IN SMOO.
C        YPP  - ARRAY OF SECOND DERIVATIVES THAT DEFINE THE SPLINE.
C               USUALLY YPP IS THE SAME AS YPP IN SPLIFT OR R2 IN SMOO.
C        N    - THE NUMBER OF DATA POINTS THAT DEFINE THE SPLINE.
C        XLO  - LEFT ENDPOINT OF INTEGRATION INTERVALS.
C        XUP  - RIGHT ENDPOINT OR ARRAY OF RIGHT ENDPOINTS OF
C               INTEGRATION INTERVALS IN ASCENDING ORDER.
C        NUP  - THE NUMBER OF RIGHT ENDPOINTS.  IF NUP IS GREATER THAN
C               1, THEN XUP AND ANS MUST BE DIMENSIONED AT LEAST NUP.
C
C      --OUTPUT--
C
C        ANS -- ARRAY OF INTEGRAL VALUES, THAT IS,
C               ANS(I) = INTEGRAL FROM XLO TO XUP(I)
C        IERR -- ERROR STATUS
C                = 1 INTEGRATION SUCCESSFUL
C                = 2 IMPROPER INPUT - N.LT.4 OR NUP.LT.1
C                = 3 IMPROPER INPUT - ABSCISSAS NOT IN
C                        STRICTLY ASCENDING ORDER
C                = 4 IMPROPER INPUT - RIGHT ENDPOINTS XUP NOT
C                        IN ASCENDING ORDER
C                = 5 IMPROPER INPUT - XLO.GT.XUP(1)
C                = 6 INTEGRATION SUCCESSFUL BUT AT LEAST ONE ENDPOINT
C                        NOT WITHIN SPAN OF ABSCISSAS
C              ** NOTE.  ERRCHK PROCESSES DIAGNOSTICS FOR CODES 2,3,4,5
C
C   CHECK FOR IMPROPER INPUT
C
       IERR = 2
       IF(N .LT. 4  .OR.  NUP .LT. 1) THEN 
         RETURN
       ENDIF
       NM1 = N-1
       NM2 = N-2
       IERR = 3
       DO 2 I = 1,NM1
         IF(X(I) .GE. X(I+1)) THEN
           RETURN
         ENDIF
 2     CONTINUE
       IF(NUP .NE. 1) THEN
         IERR = 4
         DO 3 I = 2,NUP
           IF(XUP(I-1) .GT. XUP(I)) THEN
             RETURN
           ENDIF
 3       CONTINUE
       ENDIF
       IERR = 5
       IF(XLO .GT. XUP(1)) THEN
         RETURN
       ENDIF
       IERR = 1
       IF(XLO .LT. X(1)  .OR.  XUP(NUP) .GT. X(N)) IERR = 6
C
C   LOCATE XLO IN INTERVAL (X(I),X(I+1))
C
       DO 10 I = 1,NM2
         IF(XLO .LT. X(I+1)) GO TO 20
 10      CONTINUE
       I = NM1
 20    HLO = XLO-X(I)
       HLO2 = HLO*HLO
       HI = X(I+1)-X(I)
       HI2 = HI*HI
       DO 30 J = 1,NUP
         IF(XUP(J) .GT. X(I+1)  .AND.  XLO .LT. X(NM1)) GO TO 40
C
C   COMPUTE SPECIAL CASES OF XUP IN INTERVAL WITH XLO
C
         HUP = XUP(J)-X(I)
         HSUM = HUP+HLO
         HDIFF = HUP-HLO
         HUP2 = HUP*HUP
         SUM = (YPP(I+1)-YPP(I))*HSUM*HDIFF*(HUP2+HLO2)/(24*HI)
         SUM = SUM + YPP(I)*HDIFF*(HUP2+HLO*HUP+HLO2)/6
         SUM = SUM + YP(I)*HDIFF*HSUM/2
         SUM = SUM + Y(I)*HDIFF
 30    ANS(J) = SUM
       RETURN
C
C   COMPUTE INTEGRAL BETWEEN XLO AND X(I+1) AS FOUR TERMS IN TAYLOR
C   POLYNOMIAL AND ADVANCE I TO I+1
C
 40    HDIFF = HI-HLO
       HSUM = HI+HLO
       SUM0 = Y(I)*HDIFF
       SUM1 = YP(I)*HDIFF*HSUM
       SUM2 = YPP(I)*HDIFF*(HI2+HI*HLO+HLO2)
       SUM3 = (YPP(I+1)-YPP(I))*HDIFF*HSUM*(HI2+HLO2)/HI
       I = I+1
C
C   LOCATE EACH XUP(M) IN INTERVAL (X(I),X(I+1))
C
       DO 80 M = J,NUP
 50      IF(XUP(M) .LT. X(I+1)  .OR.  I .EQ. NM1) GO TO 60
C
C   AUGMENT INTEGRAL BETWEEN ABSCISSAS TO INCLUDE INTERVAL
C   (X(I),X(I+1)) AND ADVANCE I TO I+1
C
         HI = X(I+1)-X(I)
         HI2 = HI*HI
         HI3 = HI2*HI
         SUM0 = SUM0 + Y(I)*HI
         SUM1 = SUM1 + YP(I)*HI2
         SUM2 = SUM2 + YPP(I)*HI3
         SUM3 = SUM3 + (YPP(I+1)-YPP(I))*HI3
         I = I+1
         GO TO 50
C
C   INTEGRAL BETWEEN X(I) AND XUP(M) IS ZERO
C
 60      IF(XUP(M) .NE. X(I)) THEN
C
C   COMPUTE INTEGRAL BETWEEN X(I) AND XUP(M) AND EVALUATE
C   TAYLOR POLYNOMIAL IN REVERSE ORDER
C
           HUP = XUP(M)-X(I)
           HUP2 = HUP*HUP
           HUP3 = HUP2*HUP
           HUP4 = HUP3*HUP
           HI = X(I+1)-X(I)
           PSUM0 = Y(I)*HUP
           PSUM1 = YP(I)*HUP2
           PSUM2 = YPP(I)*HUP3
           PSUM3 = (YPP(I+1)-YPP(I))*HUP4/HI
           SUM = (SUM3+PSUM3)/24 + (SUM2+PSUM2)/6
           SUM = SUM + (SUM1+PSUM1)/2
           SUM = SUM + (SUM0+PSUM0)
         ELSE
           SUM = ((SUM3/24 + SUM2/6) + SUM1/2) + SUM0
         ENDIF
 80    ANS(M) = SUM
       RETURN
       END
