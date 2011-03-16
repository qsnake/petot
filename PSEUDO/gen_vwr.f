      PROGRAM GEN_VWR
c
c ********************************************************************
c *                                                                  *
c *  Program reads out the pseudopotentials from pseudo.dat01
c *  and regenerate the wavefunctions using DIFNRL, then
c *  writes out them in vwr.atom file for the use of PEtot.
c *  
c *  Written by Lin-Wang Wang at NERSC, 02/27/2001
c *  
c *                                                                  *
c ********************************************************************
c
      IMPLICIT REAL*8 (A-H,O-Z)
c
      PARAMETER (NRMAX=2000,LMAX=5)
c
      PARAMETER (ZERO=0.0D0,SMALL=0.5D-5,ONE=1.0D0,TWO=2.0D0)
c
      DIMENSION R(NRMAX),RAB(NRMAX),NO(2*LMAX),LO(2*LMAX),SO(LMAX),
     1 VIOD(LMAX,NRMAX),VIOU(LMAX,NRMAX),VID(NRMAX),CDD(NRMAX),
     2 CDC(NRMAX),EV(LMAX),Y(NRMAX),YP(NRMAX),YPP(NRMAX),W(NRMAX*3),
     3 CDU(NRMAX),VIU(NRMAX),S1(NRMAX),S2(NRMAX),ETOT(10),EVL(2),
     4 AR(NRMAX,LMAX),ANORM(LMAX),INORM(LMAX),VQL(NRMAX),ARD(NRMAX),
     5 EVI(5),EVD(5)

       integer ispd(5)

c
c     CHARACTER*1 ISPP
      CHARACTER*2 ICORR,NAMEAT 
      CHARACTER*3 IREL
      CHARACTER*4 NICORE 
      CHARACTER*10 IRAY(6),ITITLE(7),DATED

      character*2 type_tmp,icorr_tmp,nameat_tmp
      character*3 name_tmp,kerker_tmp
      character*1 ispp_tmp
      real*8 zd_tmp(5),zu_tmp(5)
      integer n_tmp(5),l_tmp(5)
c
      DATA SO/5*ZERO/
      DATA EVI/5*ZERO/
      DATA EVD/5*ZERO/
      DO 1 I=1,NRMAX
        CDU(I)=ZERO
 1    CONTINUE
c  
      
      open(unit=2, file='atom.input',form='formatted',
     &  status='old')
      rewind(2)
      read(2,*) type_tmp 
      read(2,*) kerker_tmp
      read(2,*) nameat_tmp,icorr_tmp,ispp_tmp
      read(2,*) 
      read(2,*) ncore_tmp,nval_tmp
      do i=1,nval_tmp
      read(2,*) n_tmp(i),l_tmp(i),zd_tmp(i),zu_tmp(i)
      enddo
      read(2,*) rcut_s_tmp,rcut_p_tmp,rcut_d_tmp
      read(2,*)
      read(2,*) LLOCAL, s_occ,p_occ,d_occ
      close(2)

c
c
c  Open and read in data from file `pseudo.dat01'.
c
      OPEN (UNIT=7,FILE='pseudo.dat01',FORM='UNFORMATTED')
      READ(7) NAMEAT,ICORR,IREL,NICORE,(IRAY(I),I=1,6),
     1 (ITITLE(I),I=1,7),NORB,NUMU,NR,A,B,ZNUC
      write(6,*) "norb,numu,nr,A,B,Znuc",norb,numu,nr,A,B,Znuc
      WRITE(6,*) "nameat,icorr,irel,nicore ",
     & NAMEAT,ICORR,IREL,NICORE

      NR=NR+1
      IZNUC=INT(ZNUC+0.1)
      READ(7) (R(I),I=2,NR)
      DO 5 I=1,NORB
        READ(7) LO(I),(VIOD(LO(I)+1,J),J=2,NR)
 5    CONTINUE
      DO 10 I=1,NUMU
        READ(7) LO(I+NORB),(VIOU(LO(I+NORB)+1,J),J=2,NR)
 10   CONTINUE
      READ(7) (CDC(I),I=2,NR)
      READ(7) (CDD(I),I=2,NR)
      CLOSE (UNIT=7)

c
c  Set up RAB integration grid
c
      R(1)=ZERO
      DO 15 I=1,NR
        RAB(I) = (R(I)+A)*B
 15   CONTINUE
c
c   Find the total amount of charge in the valence charge.
c   Fit cdd to splines and integrate.
c
      IF (NICORE .EQ. 'nc  ') THEN
        ITYPE=0
      ELSEIF(NICORE .EQ. 'pcec' .OR. NICORE .EQ. 'fcec') THEN
        ITYPE=1
      ELSE
        ITYPE=2 
      ENDIF
      Y(1) = ZERO
      DO 20 I=2,NR
        Y(I) = CDD(I)
 20   CONTINUE
      ISX = 0
      A1 = ZERO
      AN = ZERO
      B1 = ZERO
      BN = ZERO
      NRM = NR
      CALL SPLIFT(R,Y,YP,YPP,NRM,W,IERR,ISX,A1,B1,AN,BN)
      IF ( IERR .NE. 1) THEN
        WRITE(6,2005)IERR
        STOP 'SPLIFT'
      ENDIF
 2005 FORMAT(1X,'****** Error in splift ierr =',I2)
      XLO = ZERO 
      NUP = 1
      CALL SPLIQ(R,Y,YP,YPP,NRM,XLO,R(NR),NUP,TOTVEL,IERR)
      IF ( IERR .NE. 1 ) THEN
        WRITE(6,2006)IERR
        STOP 'SPLIQ'
      ENDIF
 2006 FORMAT(1X,'****** Error in spliq ierr =',I2)
c
c  Write out `pseudo.dat' info to `kb.out'.
c
      WRITE(6,2007)NAMEAT,ICORR,IREL,NICORE,TOTVEL,
     1 (IRAY(J),J=1,6),
     1 (ITITLE(J),J=1,7),A,B
 2007 FORMAT(/,A2,2X,A2,2X,A3,2X,A4,' pseudopotential with ',
     1 F10.6,' valence electrons',//,1X,6A10,4X,//,1X,7A10,//,
     2 ' a = ',F9.7,' b = ',F9.7,//)
c
c  Find el-el potential
c
      CALL VELECT(0,0,ICORR,' ',ITYPE,NR,R,RAB,TOTVEL,CDD,
     1 CDU,CDC,VID,VIU,ETOT,Y,YP,YPP,S1,S2,W)
c
c  Set up potentials
c
      C2 = -ONE/B**2
      C1 = -2*C2 + ONE/4
      DO 25 I=1,NORB
        LP=LO(I)+1
        NO(I)=LP
        LLP=LP*(LP-1)
c
c  Set up hamiltonian matrix for kinetic energy,
c  only the diagonal depends on the potential.
c
        Y(1)  = C1 / (R(2)+A)**2
        YP(1)  = ZERO
        YPP(1) = ZERO
        DO 30 K=3,NR
          Y(K-1)  = C1 / (R(K)+A)**2
          YP(K-1)  = C2 / ((R(K)+A)*(R(K-1)+A))
          YPP(K-1) = YP(K-1)**2
 30     CONTINUE
        DO 35 K=2,NR
          Y(K-1)=Y(K-1)+(VIOD(LP,K)+LLP/R(K))/R(K)+VID(K)
 35     CONTINUE
c
c  Diagonalize and find wave function and store in AR().
c               
        EPS = -ONE
        CALL TRIDIB(NR-1,EPS,Y,YP,YPP,BL,BU,1,1,E,IND,IERR,
     1   S1,S2) 
        EV(I)=E
        DO 40 J=2,NR
          Y(J)=(VIOD(LP,J)+LLP/R(J))/R(J)+VID(J)
 40     CONTINUE
        IFLAG = 0
        CALL DIFNRL(1,I,Y,AR(1,LP),YP,LMAX,NR,A,B,R,RAB,NORB,NO,
     1   LO,SO,ZNUC,VIOD,VIOU,VID,VIU,EV,IFLAG,YPP,S1,S2,EVI)
 25   CONTINUE

*******************************************************************
      iref_s=1      !flag for evaluting KB projector for s in PEtot
      iref_p=1
      iref_d=1
      iTB_s=1       !flag for TB initialization using s, not used in PEtot 
      iTB_p=1
      iTB_d=0

      if(llocal.eq.0) iref_s=0
      if(llocal.eq.1) iref_p=0
      if(llocal.eq.2) iref_d=0
   

       ispd(1)=0
       ispd(2)=0
       ispd(3)=0
       ispd(4)=0
       ispd(5)=0
       do i=1,norb
       ispd(LO(i)+1)=1
       enddo

       if(ispd(4).eq.1.or.ispd(5).eq.1) then
       write(6,*) "PETOT does not support pseudopot. for F and G orbits"
       write(6,*) "stop"
       stop
       endif

       
       if(ispd(llocal+1).eq.0) then
       write(6,*) "THE LOCAL L NOT CALCULATED AS ONE VALENCE ORBITAL"
       write(6,*) "RESET atom.input, STOP"
       stop
       endif
       
       do jj=1,3 
       if(ispd(jj).eq.0) then
         if(jj.eq.1) iref_s=0
         if(jj.eq.2) iref_p=0
         if(jj.eq.3) iref_d=0
       do i=2,nr
       viod(jj,i)=viod(llocal+1,i)
       viou(jj,i)=0.d0
       AR(i,jj)=1.d0*R(i)
       enddo
       endif
       enddo

       

cccccc When lsda is used, an averaged v_up and v_down is used to 
cccccc generate the wavefunction and used for the potential. 
cccccc There is no different references for up and down in the PEtot
cccccc calculation and in the final pseudopotential output here. 
       numcore=0
       numrel=0

       if(NICORE.eq."nc  ") numcore=0      ! no core charge
       if(NICORE.eq."pcec".or.NICORE.eq."fcec") numcore=1      ! with core charge, for exchange
       if(NICORE.eq."pche".or.NICORE.eq."fche") then           ! with core charge, for Hartree+XC
       numcore=1
       write(6,*) 
     &  "PEtot doesn't support core charge for Hatr.+exch+corr"
       write(6,*) "it supports core charge for exch+corr only"
       write(6,*) "change **ph** to **pe** in atom.input, stop"
       stop
       endif

      
       if(IREL.eq."nrl") numrel=0      ! non-relativisitc calculation
       if(IREL.eq."isp") numrel=0      ! spin, lsda, but no spin-orbit coupling
       if(IREL.eq."rel") numrel=1      ! relativisitic calc. include spin-orbit coupling



      iatom=charge(NAMEAT)+0.1

      if(NR.gt.800.and.R(800).gt.40.d0) NR=800

      if(nval_tmp.gt.3) then
      write(6,*) "WARNING: number of valence orbit.gt.3 !!" 
      nval_tmp=3
      endif

      pi=4*datan(1.d0)

      do i=2,NR
      if(dabs(AR(i,1)).lt.1.D-80) AR(i,1)=0.d0
      if(dabs(AR(i,2)).lt.1.D-80) AR(i,2)=0.d0
      if(dabs(AR(i,3)).lt.1.D-80) AR(i,3)=0.d0
      if(dabs(CDC(i)).lt.1.D-80) CDC(i)=0.d0
      CDC(i)=CDC(i)/(4*pi*R(i)**2)
      enddo



       open(10,file="vwr."//NAMEAT)
       rewind(10)

       write(10,900) NR-1,numcore,iatom,znuc,llocal+1,
     & s_occ,p_occ,d_occ,numrel

c       write(10,901) iref_s,iref_p,iref_d,iTB_s,iTB_p,iTB_d,
c     & type_tmp,kerker_tmp,icorr_tmp,ispp_tmp,
c     & rcut_s_tmp,rcut_p_tmp,rcut_d_tmp,
c     & ncore_tmp,nval_tmp,
c     & ((n_tmp(i),l_tmp(i),zd_tmp(i),zu_tmp(i)),i=1,nval_tmp)

      if(nval_tmp.eq.1) then
       write(10,901) iref_s,iref_p,iref_d,iTB_s,iTB_p,iTB_d,
     & type_tmp,kerker_tmp,icorr_tmp,ispp_tmp,
     & rcut_s_tmp,rcut_p_tmp,rcut_d_tmp,
     & ncore_tmp,nval_tmp,
     & n_tmp(1),l_tmp(1),zd_tmp(1),zu_tmp(1)
        endif
 
       if(nval_tmp.eq.2) then
       write(10,901) iref_s,iref_p,iref_d,iTB_s,iTB_p,iTB_d,
     & type_tmp,kerker_tmp,icorr_tmp,ispp_tmp,
     & rcut_s_tmp,rcut_p_tmp,rcut_d_tmp,
     & ncore_tmp,nval_tmp,
     & n_tmp(1),l_tmp(1),zd_tmp(1),zu_tmp(1),
     & n_tmp(2),l_tmp(2),zd_tmp(2),zu_tmp(2)
        endif

        if(nval_tmp.eq.3) then
       write(10,901) iref_s,iref_p,iref_d,iTB_s,iTB_p,iTB_d,
     & type_tmp,kerker_tmp,icorr_tmp,ispp_tmp,
     & rcut_s_tmp,rcut_p_tmp,rcut_d_tmp,
     & ncore_tmp,nval_tmp,
     & n_tmp(1),l_tmp(1),zd_tmp(1),zu_tmp(1),
     & n_tmp(2),l_tmp(2),zd_tmp(2),zu_tmp(2),
     & n_tmp(3),l_tmp(3),zd_tmp(3),zu_tmp(3)
        endif


900   format(i4,", ",i1,",",i3,", ",f4.1,", ",i1,", ",
     & f5.2,", ",f5.2,", ",f5.2,", ",i1,
     & "  |nrr,icor,iatom,z,spd_loc, occ_s,occ_p,occ_d,iso")    

901   format(3(i1,1x),2x,3(i1,1x),
     &  "  |<-iref_s,p,d,iTB_s,p,d; inform->|  ",
     &  a2,1x,a3,1x,a2,1x,a1,";rc=",3(f3.1,1x),
     &  ";nc,nv=",i1,1x,i1,";n,l,od,ou=",
     &   3(i1,1x,i1,1x,f4.1,1x,f4.1,":"))
     

       if(numcore.eq.0.and.numrel.eq.0) then
       do i=2,NR
       write(10,1000) R(i),VIOD(1,i)/R(i)/2,VIOD(2,i)/R(i)/2,
     & VIOD(3,i)/R(i)/2,AR(i,1)/R(i),AR(i,2)/R(i),AR(i,3)/R(i)
       enddo
       write(10,*) "#r, vs(Hart.),vp,vd,phi_s/sqr(4pi),phi_p/sqr4pi,"//
     &   "phi_d/sqr4pi"
       endif

       if(numcore.eq.1.and.numrel.eq.0) then
       do i=2,NR
       write(10,1001) R(i),VIOD(1,i)/R(i)/2,VIOD(2,i)/R(i)/2,
     & VIOD(3,i)/R(i)/2,AR(i,1)/R(i),AR(i,2)/R(i),AR(i,3)/R(i),
     & CDC(i)
       enddo
       write(10,*) "#r, vs(Hart.),vp,vd,phi_s/sqr(4pi),phi_p/sqr4pi,"//
     &   "phi_d/sqr4pi,rho_core"

       endif

       if(numcore.eq.0.and.numrel.eq.1) then
       do i=2,NR
       write(10,1002) R(i),VIOD(1,i)/R(i)/2,VIOD(2,i)/R(i)/2,
     & VIOD(3,i)/R(i)/2,AR(i,1)/R(i),AR(i,2)/R(i),AR(i,3)/R(i),
     & -VIOU(2,i)/R(i)/2,-VIOU(3,i)/R(i)/2
       enddo
       write(10,*) "#r, vs(Hart.),vp,vd,phi_s/sqr(4pi),phi_p/sqr4pi,"//
     &   "phi_d/sqr4pi,v_(p+1/2)-v_(p-1/2),v_(d+1/2)-v_(d-1/2)"

       endif

       if(numcore.eq.1.and.numrel.eq.1) then
       do i=2,NR
       write(10,1003) R(i),VIOD(1,i)/R(i)/2,VIOD(2,i)/R(i)/2,
     & VIOD(3,i)/R(i)/2,AR(i,1)/R(i),AR(i,2)/R(i),AR(i,3)/R(i),
     & CDC(i),-VIOU(2,i)/R(i)/2,-VIOU(3,i)/R(i)/2
       enddo
       write(10,*) "#r, vs(Hart.),vp,vd,phi_s/sqr(4pi),"//
     & "phi_p/sqr4pi,phi_d/sqr4pi,rho_core,v_(p+1/2)-v_(p-1/2),"//
     &    "v_(d+1/2)-v_(d-1/2)"
       endif

       close(10)


1000   format(7(E14.8,1x))
1001   format(8(E14.8,1x))
1002   format(9(E14.8,1x))
1003   format(10(E14.8,1x))

      STOP
      END

