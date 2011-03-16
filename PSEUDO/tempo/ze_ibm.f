           SUBROUTINE ZEDATE(BDATE)
C
C   GETS THE DATE (DAY-MONTH-YEAR)
C   IBM version
C   Written by Pedro Gomes da Costa
C
           CHARACTER*10 BDATE
           CHARACTER*9 ADATE,UNXCOM*26
           WRITE(UNXCOM,'(a)') 'date +"%d-%h-%y" > zzf.tmp'
           call system(UNXCOM)
           open(31,file='zzf.tmp',status='old')
           read(31,'(a)')ADATE
           close(31,status='delete')
           WRITE(BDATE,100) ADATE
100        FORMAT(A9,' ')
           RETURN
           END
C
          SUBROUTINE ZESEC(TIM)
C
C   GETS CPU TIME IN SECONDS
C   IBM version
C   Written by Pedro Gomes da Costa
C
           CHARACTER*1 s1,s2
           CHARACTER*20 UNXCOM
           DOUBLE PRECISION TIM
           WRITE(UNXCOM,'(a)') 'date +"%T" > zzf.tmp'
           call system(UNXCOM)
           open(31,file='zzf.tmp',status='old')
           read(31,200)i,s1,j,s2,k
200        format(i2,a1,i2,a1,i2)
           close(31,status='delete')
           TIM=float(i*3600+j*60+k)
           RETURN
           END
