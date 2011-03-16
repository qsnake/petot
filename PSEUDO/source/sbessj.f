       REAL*8 FUNCTION SBESSJ(N,X)
       IMPLICIT REAL*8 (A-H,O-Z)
C
C  
C
       PARAMETER(ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,ZERO=0.0D0)
       PARAMETER(FIVE = 5.0D0,TEN = 10.0D0,FOURTN = 14.0D0)
C
C      SPHERICAL BESSEL FUNCTION OF THE FIRST KIND
C
       IF(ABS(X) .GT. 0.001) THEN
         SB0 = SIN(X)/X
       ELSE
         X2 = X*X/TWO
         SB0 = ONE - (X2/THREE)*(ONE - X2/TEN)
       ENDIF
       IF(N .EQ. 0) THEN
         SBESSJ = SB0
       ELSE
         IF(ABS(X) .GT. 0.001) THEN
           SB1 = (SIN(X)/X - COS(X)) / X
         ELSE
           X2 = X*X/TWO
           SB1 = (X/THREE)*(ONE - (X2/FIVE)*(ONE - X2/FOURTN))
         ENDIF
         IF(N .EQ. 1) THEN
           SBESSJ = SB1
         ELSEIF(X .EQ. ZERO) THEN
           SBESSJ = ZERO
         ELSE
           BY = SB1
           BYM = SB0
           UX = ONE / X
           DO 10 J=1,N-1
             BYP = REAL(2*J+1)*UX*BY - BYM
             BYM = BY
             BY = BYP
 10        CONTINUE
           SBESSJ = BY
         ENDIF
       ENDIF
       RETURN
       END
