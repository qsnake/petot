      SUBROUTINE DJACOBI(A,N,NP,D,V,NROT)
      implicit double precision(a-h,o-z)
      PARAMETER (NMAX=510)
      DIMENSION A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)
      DO 12 IP=1,N
        DO 11 IQ=1,N
          V(IP,IQ)=0.d0
11      CONTINUE
        V(IP,IP)=1.d0
12    CONTINUE
      DO 13 IP=1,N
        B(IP)=A(IP,IP)
        D(IP)=B(IP)
        Z(IP)=0.d0
13    CONTINUE
      NROT=0
      DO 24 I=1,50
        SM=0.d0
        DO 15 IP=1,N-1
          DO 14 IQ=IP+1,N
            SM=SM+dABS(A(IP,IQ))
14        CONTINUE
15      CONTINUE
        IF(SM.EQ.0.d0)RETURN
        IF(I.LT.4)THEN
          TRESH=0.2d0*SM/N**2
        ELSE
          TRESH=0.d0
        ENDIF
        DO 22 IP=1,N-1
          DO 21 IQ=IP+1,N
            G=100.d0*dABS(A(IP,IQ))
            IF((I.GT.4).AND.(dABS(D(IP))+G.EQ.dABS(D(IP)))
     *         .AND.(dABS(D(IQ))+G.EQ.dABS(D(IQ))))THEN
              A(IP,IQ)=0.d0
            ELSE IF(dABS(A(IP,IQ)).GT.TRESH)THEN
              H=D(IQ)-D(IP)
              IF(dABS(H)+G.EQ.dABS(H))THEN
                T=A(IP,IQ)/H
              ELSE
                THETA=0.5d0*H/A(IP,IQ)
                T=1.d0/(dABS(THETA)+dSQRT(1.d0+THETA**2))
                IF(THETA.LT.0.)T=-T
              ENDIF
              C=1.d0/dSQRT(1+T**2)
              S=T*C
              TAU=S/(1.d0+C)
              H=T*A(IP,IQ)
              Z(IP)=Z(IP)-H
              Z(IQ)=Z(IQ)+H
              D(IP)=D(IP)-H
              D(IQ)=D(IQ)+H
              A(IP,IQ)=0.d0
              DO 16 J=1,IP-1
                G=A(J,IP)
                H=A(J,IQ)
                A(J,IP)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
16            CONTINUE
              DO 17 J=IP+1,IQ-1
                G=A(IP,J)
                H=A(J,IQ)
                A(IP,J)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
17            CONTINUE
              DO 18 J=IQ+1,N
                G=A(IP,J)
                H=A(IQ,J)
                A(IP,J)=G-S*(H+G*TAU)
                A(IQ,J)=H+S*(G-H*TAU)
18            CONTINUE
              DO 19 J=1,N
                G=V(J,IP)
                H=V(J,IQ)
                V(J,IP)=G-S*(H+G*TAU)
                V(J,IQ)=H+S*(G-H*TAU)
19            CONTINUE
              NROT=NROT+1
            ENDIF
21        CONTINUE
22      CONTINUE
        DO 23 IP=1,N
          B(IP)=B(IP)+Z(IP)
          D(IP)=B(IP)
          Z(IP)=0.d0
23      CONTINUE
24    CONTINUE
      PAUSE '50 iterations should never happen'
      RETURN
      END


