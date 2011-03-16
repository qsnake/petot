      SUBROUTINE ZEDATE(BDATE)
C   
C   Gets the data (DAY-MONTH-YEAR)
C   Silicon Graphics version
C
       INTEGER IMON,IDAY,IYEAR
       CHARACTER*10 BDATE
       CHARACTER*3 MONTH(12)
       DATA MONTH /'JAN','FEB','MAR','APR','MAY','JUN',
     1             'JUL','AUG','SEP','OCT','NOV','DEC'/
       CALL IDATE(IMON,IDAY,IYEAR)
       WRITE(BDATE,100) IDAY,MONTH(IMON),IYEAR
 100   FORMAT(I2,'-',A3,'-',I2,' ')
       RETURN
       END
C
      SUBROUTINE ZESEC(T)
C
C     SILICON GRAPHICS VERSION
C     COMPUTING TIME IN SECONDS
C
      INTEGER IT0,IT1
      INTEGER TIME
      REAL*8 T
      SAVE IT0
      DATA IFIRST /0/
      EXTERNAL TIME
C
      IF(IFIRST .EQ. 0) THEN
        IT0 = TIME()
        IFIRST = 1
        T = 0.0
      ELSE      
        IT1 = TIME()
	T = REAL(IT1-IT0) 
      ENDIF
      RETURN
      END
