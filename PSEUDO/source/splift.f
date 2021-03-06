       SUBROUTINE SPLIFT (X,Y,YP,YPP,N,W,IERR,ISX,A1,B1,AN,BN)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER (FOUR=4.D0)
CRAY      PARAMETER (FOUR=4.0)
C  
C  NJTJ
C  ###  CRAY CONVERSIONS  
C  ###    1)Comment out the implicit double precision.
C  ###    2)Switch double precision parameter 
C  ###      to single precision parameter
C  ###  CRAY CONVERSIONS  
C  NJTJ
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
C     WRITTEN BY RONDALL E. JONES
C
C     ABSTRACT
C         SPLIFT FITS AN INTERPOLATING CUBIC SPLINE TO THE N DATA POINT
C         GIVEN IN X AND Y AND RETURNS THE FIRST AND SECOND DERIVATIVES
C         IN YP AND YPP.  THE RESULTING SPLINE (DEFINED BY X, Y, AND
C         YPP) AND ITS FIRST AND SECOND DERIVATIVES MAY THEN BE
C         EVALUATED USING SPLINT.  THE SPLINE MAY BE INTEGRATED USING
C         SPLIQ.  FOR A SMOOTHING SPLINE FIT SEE SUBROUTINE SMOO.
C
C     DESCRIPTION OF ARGUMENTS
C         THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST
C         E.G.   X(N), Y(N), YP(N), YPP(N), W(3N)
C
C       --INPUT--
C
C         X    - ARRAY OF ABSCISSAS OF DATA (IN INCREASING ORDER)
C         Y    - ARRAY OF ORDINATES OF DATA
C         N    - THE NUMBER OF DATA POINTS.  THE ARRAYS X, Y, YP, AND
C                YPP MUST BE DIMENSIONED AT LEAST N.  (N .GE. 4)
C         ISX  - MUST BE ZERO ON THE INITIAL CALL TO SPLIFT.
C                IF A SPLINE IS TO BE FITTED TO A SECOND SET OF DATA
C                THAT HAS THE SAME SET OF ABSCISSAS AS A PREVIOUS SET,
C                AND IF THE CONTENTS OF W HAVE NOT BEEN CHANGED SINCE
C                THAT PREVIOUS FIT WAS COMPUTED, THEN ISX MAY BE
C                SET TO ONE FOR FASTER EXECUTION.
C         A1,B1,AN,BN - SPECIFY THE END CONDITIONS FOR THE SPLINE WHICH
C                ARE EXPRESSED AS CONSTRAINTS ON THE SECOND DERIVATIVE
C                OF THE SPLINE AT THE END POINTS (SEE YPP).
C                THE END CONDITION CONSTRAINTS ARE
C                        YPP(1) = A1*YPP(2) + B1
C                AND
C                        YPP(N) = AN*YPP(N-1) + BN
C                WHERE
C                        ABS(A1).LT. 1.0  AND  ABS(AN).LT. 1.0.
C
C                THE SMOOTHEST SPLINE (I.E., LEAST INTEGRAL OF SQUARE
C                OF SECOND DERIVATIVE) IS OBTAINED BY A1=B1=AN=BN=0.
C                IN THIS CASE THERE IS AN INFLECTION AT X(1) AND X(N).
C                IF THE DATA IS TO BE EXTRAPOLATED (SAY, BY USING SPLIN
C                TO EVALUATE THE SPLINE OUTSIDE THE RANGE X(1) TO X(N))
C                THEN TAKING A1=AN=0.5 AND B1=BN=0 MAY YIELD BETTER
C                RESULTS.  IN THIS CASE THERE IS AN INFLECTION
C                AT X(1) - (X(2)-X(1)) AND AT X(N) + (X(N)-X(N-1)).
C                IN THE MORE GENERAL CASE OF A1=AN=A  AND B1=BN=0,
C                THERE IS AN INFLECTION AT X(1) - (X(2)-X(1))*A/(1.0-A)
C                AND AT X(N) + (X(N)-X(N-1))*A/(1.0-A).
C
C                A SPLINE THAT HAS A GIVEN FIRST DERIVATIVE YP1 AT X(1)
C                AND YPN AT Y(N) MAY BE DEFINED BY USING THE
C                FOLLOWING CONDITIONS.
C
C                A1=-0.5
C
C                B1= 3.0*((Y(2)-Y(1))/(X(2)-X(1))-YP1)/(X(2)-X(1))
C
C                AN=-0.5
C
C                BN=-3.0*((Y(N)-Y(N-1))/(X(N)-X(N-1))-YPN)/(X(N)-X(N-1)
C
C       --OUTPUT--
C
C         YP   - ARRAY OF FIRST DERIVATIVES OF SPLINE (AT THE X(I))
C         YPP  - ARRAY OF SECOND DERIVATIVES OF SPLINE (AT THE X(I))
C         IERR - A STATUS CODE
C              --NORMAL CODE
C                 1 MEANS THAT THE REQUESTED SPLINE WAS COMPUTED.
C              --ABNORMAL CODES
C                 2 MEANS THAT N, THE NUMBER OF POINTS, WAS .LT. 4.
C                 3 MEANS THE ABSCISSAS WERE NOT STRICTLY INCREASING.
C
C       --WORK--
C
C         W    - ARRAY OF WORKING STORAGE DIMENSIONED AT LEAST 3N.
       DIMENSION X(N),Y(N),YP(N),YPP(N),W(N,3)
C
       IF (N.LT.4) THEN
         IERR = 2
         RETURN
       ENDIF
       NM1  = N-1
       NM2  = N-2
       IF (ISX.GT.0) GO TO 40
       DO 5 I=2,N
         IF (X(I)-X(I-1) .LE. 0) THEN
           IERR = 3
           RETURN
         ENDIF
 5     CONTINUE
C
C     DEFINE THE TRIDIAGONAL MATRIX
C
       W(1,3) = X(2)-X(1)
       DO 10 I=2,NM1
         W(I,2) = W(I-1,3)
         W(I,3) = X(I+1)-X(I)
 10      W(I,1) = 2*(W(I,2)+W(I,3))
       W(1,1) = FOUR
       W(1,3) =-4*A1
       W(N,1) = FOUR
       W(N,2) =-4*AN
C
C     L U DECOMPOSITION
C
       DO 30 I=2,N
         W(I-1,3) = W(I-1,3)/W(I-1,1)
 30    W(I,1) = W(I,1) - W(I,2)*W(I-1,3)
C
C     DEFINE *CONSTANT* VECTOR
C
 40   YPP(1) = 4*B1
      DOLD = (Y(2)-Y(1))/W(2,2)
      DO 50 I=2,NM2
        DNEW   = (Y(I+1) - Y(I))/W(I+1,2)
        YPP(I) = 6*(DNEW - DOLD)
        YP(I)  = DOLD
 50   DOLD = DNEW
      DNEW = (Y(N)-Y(N-1))/(X(N)-X(N-1))
      YPP(NM1) = 6*(DNEW - DOLD)
      YPP(N) = 4*BN
      YP(NM1)= DOLD
      YP(N) = DNEW
C
C     FORWARD SUBSTITUTION
C
      YPP(1) = YPP(1)/W(1,1)
      DO 60 I=2,N
 60   YPP(I) = (YPP(I) - W(I,2)*YPP(I-1))/W(I,1)
C
C     BACKWARD SUBSTITUTION
C
       DO 70 J=1,NM1
         I = N-J
   70 YPP(I) = YPP(I) - W(I,3)*YPP(I+1)
C
C     COMPUTE FIRST DERIVATIVES
C
      YP(1) = (Y(2)-Y(1))/(X(2)-X(1)) - (X(2)-X(1))*(2*YPP(1)
     1  + YPP(2))/6
      DO 80 I=2,NM1
 80   YP(I) = YP(I) + W(I,2)*(YPP(I-1) + 2*YPP(I))/6
      YP(N) = YP(N) + (X(N)-X(NM1))*(YPP(NM1) + 2*YPP(N))/6
C
      IERR = 1
      RETURN
      END
