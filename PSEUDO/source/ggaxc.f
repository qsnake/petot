      SUBROUTINE GGAXC( AUTHOR, IREL, NSPIN, D, GD,
     .                  EPSX, EPSC, DEXDD, DECDD, DEXDGD, DECDGD )

C Finds the exchange and correlation energies at a point, and their
C derivatives with respect to density and density gradient, in the
C Generalized Gradient Correction approximation.
C Lengths in Bohr, energies in Hartrees
C Written by L.C.Balbas and J.M.Soler, Dec'96. Version 0.5.

      IMPLICIT          NONE
      CHARACTER*(*)     AUTHOR
      INTEGER           IREL, NSPIN
      DOUBLE PRECISION  D(NSPIN), DECDD, DECDGD, DEXDD, DEXDGD,
     .                  EPSC, EPSX, GD(3,NSPIN)

      IF (AUTHOR.EQ.'PBE' .OR. AUTHOR.EQ.'pbe') THEN
        CALL PBEXC( IREL, NSPIN, D, GD,
     .              EPSX, EPSC, DEXDD, DECDD, DEXDGD, DECDGD )
      ELSE
        WRITE(6,*) 'GGAXC: Unknown author ', AUTHOR
        STOP
      ENDIF
      END
