      SUBROUTINE TINVIT(NM,N,D,E,E2,M,W,IND,Z,
     X                  IERR,RV1,RV2,RV3,RV4,RV6)
C
c  njtj
c  ###  Cray conversions  
c  ###    1)Switch double precision to real.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
C
      INTEGER I,J,M,N,P,Q,R,S,II,IP,JJ,NM,ITS,TAG,IERR,GROUP
      INTEGER IND(M)
      DOUBLE PRECISION D(N),E(N),E2(N),W(M),Z(NM,M),
     X       RV1(N),RV2(N),RV3(N),RV4(N),RV6(N)
Cray      REAL D(N),E(N),E2(N),W(M),Z(NM,M),
Cray     X       RV1(N),RV2(N),RV3(N),RV4(N),RV6(N)
      DOUBLE PRECISION U,V,UK,XU,X0,X1,EPS2,EPS3,EPS4,
     X       NORM,ORDER,MACHEP
Cray      REAL U,V,UK,XU,X0,X1,EPS2,EPS3,EPS4,NORM,ORDER,MACHEP
C
      DOUBLE PRECISION ZERO,ONE,ONEM3,TWO
      PARAMETER(ZERO=0.D0,ONE=1.D0,ONEM3=1.D-3,TWO=2.D0)
Cray      PARAMETER(ZERO=0.0,ONE=1.0,ONEM3=1.E-3,TWO=2.0)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE INVERSE ITERATION TECH-
C     NIQUE IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVECTORS OF A TRIDIAGONAL
C     SYMMETRIC MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES,
C     USING INVERSE ITERATION.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E,
C          WITH ZEROS CORRESPONDING TO NEGLIGIBLE ELEMENTS OF E.
C          E(I) IS CONSIDERED NEGLIGIBLE IF IT IS NOT LARGER THAN
C          THE PRODUCT OF THE RELATIVE MACHINE PRECISION AND THE SUM
C          OF THE MAGNITUDES OF D(I) AND D(I-1).  E2(1) MUST CONTAIN
C          0.0 IF THE EIGENVALUES ARE IN ASCENDING ORDER, OR 2.0
C          IF THE EIGENVALUES ARE IN DESCENDING ORDER.  IF  BISECT,
C          TRIDIB, OR  IMTQLV  HAS BEEN USED TO FIND THE EIGENVALUES,
C          THEIR OUTPUT E2 ARRAY IS EXACTLY WHAT IS EXPECTED HERE,
C
C        M IS THE NUMBER OF SPECIFIED EIGENVALUES,
C
C        W CONTAINS THE M EIGENVALUES IN ASCENDING OR DESCENDING ORDER,
C
C        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
C          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
C          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
C          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.
C
C     ON OUTPUT-
C
C        ALL INPUT ARRAYS ARE UNALTERED,
C
C        Z CONTAINS THE ASSOCIATED SET OF ORTHONORMAL EIGENVECTORS.
C          ANY VECTOR WHICH FAILS TO CONVERGE IS SET TO ZERO,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          -R         IF THE EIGENVECTOR CORRESPONDING TO THE R-TH
C                     EIGENVALUE FAILS TO CONVERGE IN 5 ITERATIONS,
C
C        RV1, RV2, RV3, RV4, AND RV6 ARE TEMPORARY STORAGE ARRAYS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C
C                **********
      MACHEP = TWO**(-40)
C
      IERR = 0
      IF (M .EQ. 0) GO TO 1001
      TAG = 0
      ORDER = ONE - E2(1)
      Q = 0
C     ********** ESTABLISH AND PROCESS NEXT SUBMATRIX **********
  100 P = Q + 1
C
      DO 120 Q = P, N
         IF (Q .EQ. N) GO TO 140
         IF (E2(Q+1) .EQ. ZERO) GO TO 140
  120 CONTINUE
C     ********** FIND VECTORS BY INVERSE ITERATION **********
  140 TAG = TAG + 1
      S = 0
C
      DO 920 R = 1, M
         IF (IND(R) .NE. TAG) GO TO 920
         ITS = 1
         X1 = W(R)
         IF (S .NE. 0) GO TO 510
C     ********** CHECK FOR ISOLATED ROOT **********
         XU = ONE
         IF (P .NE. Q) GO TO 490
         RV6(P) = ONE
         GO TO 870
  490    NORM = ABS(D(P))
         IP = P + 1
C
         DO 500 I = IP, Q
  500    NORM = NORM + ABS(D(I)) + ABS(E(I))
C     ********** EPS2 IS THE CRITERION FOR GROUPING,
C                EPS3 REPLACES ZERO PIVOTS AND EQUAL
C                ROOTS ARE MODIFIED BY EPS3,
C                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW **********
         EPS2 = ONEM3 * NORM
         EPS3 = MACHEP * NORM
         UK = REAL(Q-P+1)
         EPS4 = UK * EPS3
         UK = EPS4 / SQRT(UK)
         S = P
  505    GROUP = 0
         GO TO 520
C     ********** LOOK FOR CLOSE OR COINCIDENT ROOTS **********
  510    IF (ABS(X1-X0) .GE. EPS2) GO TO 505
         GROUP = GROUP + 1
         IF (ORDER * (X1 - X0) .LE. ZERO) X1 = X0 + ORDER * EPS3
C     ********** ELIMINATION WITH INTERCHANGES AND
C                INITIALIZATION OF VECTOR **********
  520    V = ZERO
C
         DO 580 I = P, Q
            RV6(I) = UK
            IF (I .EQ. P) GO TO 560
            IF (ABS(E(I)) .LT. ABS(U)) GO TO 540
C     ********** WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF
C                E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY **********
            XU = U / E(I)
            RV4(I) = XU
            RV1(I-1) = E(I)
            RV2(I-1) = D(I) - X1
            RV3(I-1) = ZERO
            IF (I .NE. Q) RV3(I-1) = E(I+1)
            U = V - XU * RV2(I-1)
            V = -XU * RV3(I-1)
            GO TO 580
  540       XU = E(I) / U
            RV4(I) = XU
            RV1(I-1) = U
            RV2(I-1) = V
            RV3(I-1) = ZERO
  560       U = D(I) - X1 - XU * V
            IF (I .NE. Q) V = E(I+1)
  580    CONTINUE
C
         IF (U .EQ. ZERO) U = EPS3
         RV1(Q) = U
         RV2(Q) = ZERO
         RV3(Q) = ZERO
C     ********** BACK SUBSTITUTION
C                FOR I=Q STEP -1 UNTIL P DO -- **********
  600    DO 620 II = P, Q
            I = P + Q - II
            RV6(I) = (RV6(I) - U * RV2(I) - V * RV3(I)) / RV1(I)
            V = U
            U = RV6(I)
  620    CONTINUE
C     ********** ORTHOGONALIZE WITH RESPECT TO PREVIOUS
C                MEMBERS OF GROUP **********
         IF (GROUP .EQ. 0) GO TO 700
         J = R
C
         DO 680 JJ = 1, GROUP
  630       J = J - 1
            IF (IND(J) .NE. TAG) GO TO 630
            XU = ZERO
C
            DO 640 I = P, Q
  640       XU = XU + RV6(I) * Z(I,J)
C
            DO 660 I = P, Q
  660       RV6(I) = RV6(I) - XU * Z(I,J)
C
  680    CONTINUE
C
  700    NORM = ZERO
C
         DO 720 I = P, Q
  720    NORM = NORM + ABS(RV6(I))
C
         IF (NORM .GE. ONE) GO TO 840
C     ********** FORWARD SUBSTITUTION **********
         IF (ITS .EQ. 5) GO TO 830
         IF (NORM .NE. ZERO) GO TO 740
         RV6(S) = EPS4
         S = S + 1
         IF (S .GT. Q) S = P
         GO TO 780
  740    XU = EPS4 / NORM
C
         DO 760 I = P, Q
  760    RV6(I) = RV6(I) * XU
C     ********** ELIMINATION OPERATIONS ON NEXT VECTOR
C                ITERATE **********
  780    DO 820 I = IP, Q
            U = RV6(I)
C     ********** IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE
C                WAS PERFORMED EARLIER IN THE
C                TRIANGULARIZATION PROCESS **********
            IF (RV1(I-1) .NE. E(I)) GO TO 800
            U = RV6(I-1)
            RV6(I-1) = RV6(I)
  800       RV6(I) = U - RV4(I) * RV6(I-1)
  820    CONTINUE
C
         ITS = ITS + 1
         GO TO 600
C     ********** SET ERROR -- NON-CONVERGED EIGENVECTOR **********
  830    IERR = -R
         XU = ZERO
         GO TO 870
C     ********** NORMALIZE SO THAT SUM OF SQUARES IS
C                1 AND EXPAND TO FULL ORDER **********
  840    U = ZERO
C
         DO 860 I = P, Q
  860    U = U + RV6(I)**2
C
         XU = ONE / SQRT(U)
C
  870    DO 880 I = 1, N
  880    Z(I,R) = ZERO
C
         DO 900 I = P, Q
  900    Z(I,R) = RV6(I) * XU
C
         X0 = X1
  920 CONTINUE
C
      IF (Q .LT. N) GO TO 100
 1001 RETURN
C     ********** LAST CARD OF TINVIT **********
      END
