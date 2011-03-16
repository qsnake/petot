       SUBROUTINE ZEDATE(BDATE)
C   
C   Gets the data (DAY-MONTH-YEAR)
C   DEC Alpha version
C
       CHARACTER*10 BDATE
       CHARACTER*9 CDATE
       CHARACTER*1 NULL
       NULL = ' '
       CALL DATE(CDATE)
       BDATE = CDATE // NULL
       RETURN
       END
C
       SUBROUTINE ZESEC(T)
C
C      Gets the computing time in seconds
C      DEC Alpha version
C
       REAL*8 T
       REAL*4 T1
       T1 = SECNDS(0.0)
       T = DBLE(T1)
       RETURN
       END

