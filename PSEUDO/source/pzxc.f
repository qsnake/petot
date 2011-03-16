      SUBROUTINE PZXC( IREL, NSP, DS, EX, EC, VX, VC )

C *****************************************************************
C  Perdew-Zunger parameterization of Ceperley-Alder exchange and 
C  correlation. Ref: Perdew & Zunger, Phys. Rev. B 23 5075 (1981).
C  Adapted by J.M.Soler from routine velect of Froyen's 
C    pseudopotential generation program. Madrid, Jan'97. Version 0.5.
C **** Input *****************************************************
C INTEGER IREL    : relativistic-exchange switch (0=no, 1=yes)
C INTEGER NSP     : spin-polarizations (1=>unpolarized, 2=>polarized)
C REAL*8  DS(NSP) : total (nsp=1) or spin (nsp=2) electron density
C **** Output *****************************************************
C REAL*8  EX      : exchange energy density
C REAL*8  EC      : correlation energy density
C REAL*8  VX(NSP) : (spin-dependent) exchange potential
C REAL*8  VC(NSP) : (spin-dependent) correlation potential
C **** Units *******************************************************
C Densities in electrons/Bohr**3
C Energies in Hartrees
C *****************************************************************

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION DS(NSP), VX(NSP), VC(NSP)

       PARAMETER (ZERO=0.D0,ONE=1.D0,PFIVE=.5D0,OPF=1.5D0,PNN=.99D0)
       PARAMETER (PTHREE=0.3D0,PSEVF=0.75D0,C0504=0.0504D0) 
       PARAMETER (C0254=0.0254D0,C014=0.014D0,C0406=0.0406D0)
       PARAMETER (C15P9=15.9D0,C0666=0.0666D0,C11P4=11.4D0)
       PARAMETER (C045=0.045D0,C7P8=7.8D0,C88=0.88D0,C20P59=20.592D0)
       PARAMETER (C3P52=3.52D0,C0311=0.0311D0,C0014=0.0014D0)
       PARAMETER (C0538=0.0538D0,C0096=0.0096D0,C096=0.096D0)
       PARAMETER (C0622=0.0622D0,C004=0.004D0,C0232=0.0232D0)
       PARAMETER (C1686=0.1686D0,C1P398=1.3981D0,C2611=0.2611D0)
       PARAMETER (C2846=0.2846D0,C1P053=1.0529D0,C3334=0.3334D0)
Cray       PARAMETER (ZERO=0.0,ONE=1.0,PFIVE=0.5,OPF=1.5,PNN=0.99)
Cray       PARAMETER (PTHREE=0.3,PSEVF=0.75,C0504=0.0504) 
Cray       PARAMETER (C0254=0.0254,C014=0.014,C0406=0.0406)
Cray       PARAMETER (C15P9=15.9,C0666=0.0666,C11P4=11.4)
Cray       PARAMETER (C045=0.045,C7P8=7.8,C88=0.88,C20P59=20.592)
Cray       PARAMETER (C3P52=3.52,C0311=0.0311,C0014=0.0014)
Cray       PARAMETER (C0538=0.0538,C0096=0.0096,C096=0.096)
Cray       PARAMETER (C0622=0.0622,C004=0.004,C0232=0.0232)
Cray       PARAMETER (C1686=0.1686,C1P398=1.3981,C2611=0.2611)
Cray       PARAMETER (C2846=0.2846,C1P053=1.0529,C3334=0.3334)

C    Ceperly-Alder 'ca' constants. Internal energies in Rydbergs.
       PARAMETER (CON1=1.D0/6, CON2=0.008D0/3, CON3=0.3502D0/3) 
       PARAMETER (CON4=0.0504D0/3, CON5=0.0028D0/3, CON6=0.1925D0/3)
       PARAMETER (CON7=0.0206D0/3, CON8=9.7867D0/6, CON9=1.0444D0/3)
       PARAMETER (CON10=7.3703D0/6, CON11=1.3336D0/3)
Cray       PARAMETER (CON1=1.0/6, CON2=0.008/3, CON3=0.3502/3) 
Cray       PARAMETER (CON4=0.0504/3, CON5=0.0028/3, CON6=0.1925/3)
Cray       PARAMETER (CON7=0.0206/3, CON8=9.7867/6, CON9=1.0444/3)
Cray       PARAMETER (CON10=7.3703/6, CON11=1.3336/3) 

C      X-alpha parameter:
       PARAMETER ( ALP = 2.D0 / 3.D0 )

C      Other variables converted into parameters by J.M.Soler
       PARAMETER ( PI   = 3.14159265358979312D0 )
       PARAMETER ( HALF = 0.5D0 ) 
       PARAMETER ( TRD  = 1.D0 / 3.D0 ) 
       PARAMETER ( FTRD = 4.D0 / 3.D0 )
       PARAMETER ( TFTM = 0.51984209978974638D0 )
       PARAMETER ( A0   = 0.52106176119784808D0 )
       PARAMETER ( CRS  = 0.620350490899400087D0 )
       PARAMETER ( CXP  = - 3.D0 * ALP / (PI*A0) )
       PARAMETER ( CXF  = 1.25992104989487319D0 )

C      Find density and polarization
       IF (NSP .EQ. 2) THEN
         D = DS(1) + DS(2)
         IF (D .LE. ZERO) THEN
           EX = ZERO
           EC = ZERO
           VX(1) = ZERO
           VX(2) = ZERO
           VC(1) = ZERO
           VC(2) = ZERO
           RETURN
         ENDIF
         Z = (DS(1) - DS(2)) / D
         FZ = ((1+Z)**FTRD+(1-Z)**FTRD-2)/TFTM
         FZP = FTRD*((1+Z)**TRD-(1-Z)**TRD)/TFTM 
       ELSE
         D = DS(1)
         IF (D .LE. ZERO) THEN
           EX = ZERO
           EC = ZERO
           VX(1) = ZERO
           VC(1) = ZERO
           RETURN
         ENDIF
         Z = ZERO
         FZ = ZERO
         FZP = ZERO
       ENDIF
       RS = CRS / D**TRD

C      Exchange
       VXP = CXP / RS
       EXP = 0.75D0 * VXP
       IF (IREL .EQ. 1) THEN
         BETA = C014/RS
         SB = SQRT(1+BETA*BETA)
         ALB = LOG(BETA+SB)
         VXP = VXP * (-PFIVE + OPF * ALB / (BETA*SB))
         EXP = EXP *(ONE-OPF*((BETA*SB-ALB)/BETA**2)**2) 
       ENDIF
       VXF = CXF * VXP
       EXF = CXF * EXP

C      Correlation 
       IF (RS .GT. ONE) THEN  
         SQRS=SQRT(RS)
         TE = ONE+CON10*SQRS+CON11*RS
         BE = ONE+C1P053*SQRS+C3334*RS
         ECP = -C2846/BE
         VCP = ECP*TE/BE
         TE = ONE+CON8*SQRS+CON9*RS
         BE = ONE+C1P398*SQRS+C2611*RS
         ECF = -C1686/BE
         VCF = ECF*TE/BE
       ELSE
         RSLOG=LOG(RS)
         ECP=(C0622+C004*RS)*RSLOG-C096-C0232*RS
         VCP=(C0622+CON2*RS)*RSLOG-CON3-CON4*RS
         ECF=(C0311+C0014*RS)*RSLOG-C0538-C0096*RS
         VCF=(C0311+CON5*RS)*RSLOG-CON6-CON7*RS
       ENDIF

C      Find up and down potentials
       IF (NSP .EQ. 2) THEN
         EX    = EXP + FZ*(EXF-EXP)
         EC    = ECP + FZ*(ECF-ECP)
         VX(1) = VXP + FZ*(VXF-VXP) + (1-Z)*FZP*(EXF-EXP)
         VX(2) = VXP + FZ*(VXF-VXP) - (1+Z)*FZP*(EXF-EXP)
         VC(1) = VCP + FZ*(VCF-VCP) + (1-Z)*FZP*(ECF-ECP)
         VC(2) = VCP + FZ*(VCF-VCP) - (1+Z)*FZP*(ECF-ECP)
       ELSE
         EX    = EXP
         EC    = ECP
         VX(1) = VXP
         VC(1) = VCP
       ENDIF

C      Change from Rydbergs to Hartrees
       EX = HALF * EX
       EC = HALF * EC
       DO 10 ISP = 1,NSP
         VX(ISP) = HALF * VX(ISP)
         VC(ISP) = HALF * VC(ISP)
   10  CONTINUE
      END
