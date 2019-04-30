C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   LIMITS FOR COS(THETA(K,P)) FOR KPS-KP PARAMETRIZATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSLZTS(ZMIN,ZMAX)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSPSPC/ IPHSPC
C
      IPHSPC=0
      CALP = (PS-PEL*COSE)/PQ
      IF (DABS(CALP).GT.1D0) THEN
        IPHSPC=1
        RETURN
      ENDIF
      ATAU = DELTA/OMEGA*( PQH*EH + PH*EQH)
      ALAM = PQH*EEL + PH*EQ
      AMU = PQH*PS*DSQRT((1D0-CALP)*(1D0+CALP))
      AXI = -PH*PQ+PQH*(PQ-PS*CALP)
C-----
      IF (AMU.GT.0D0) THEN
C
      F2U = (ATAU + AXI - ALAM)/AMU
      F2O = (ATAU - AXI - ALAM)/AMU
C
      IF ( (F2U.LT.0D0).AND.(F2O.LT.0D0) ) THEN
         ZMIN = -1D0
         ZMAX =  1D0
         RETURN
      ENDIF
C
      EZ = F2U*F2O
      IF ( EZ.LE.0D0 ) THEN
         AZ = AXI*AXI + AMU*AMU
         BZ = 2D0*AXI*(ALAM - ATAU)
         CZ = (ALAM-ATAU)*(ALAM-ATAU) - AMU*AMU
         DZ = BZ*BZ - 4D0*AZ*CZ
         IF(DZ.LT.0D0)DZ=0D0
         IF ((AXI/AMU).LT.0D0) THEN
            ZMIN = -1D0
            IF (BZ.GT.0D0) THEN
               ZMAX = 2D0*CZ/(-BZ - DSQRT(DZ))
            ELSE
               ZMAX = (- BZ + DSQRT(DZ))/2D0/AZ
            ENDIF
            RETURN
         ELSE
            ZMAX = 1D0
            IF (BZ.GT.0D0) THEN
               ZMIN = (- BZ - DSQRT(DZ))/2D0/AZ
            ELSE
               ZMIN = 2D0*CZ/(-BZ + DSQRT(DZ))
            ENDIF
            RETURN
         ENDIF
      ENDIF
C
      IF ( (F2U.GT.0D0).AND.(F2O.GT.0D0) ) THEN
         AZ = AXI*AXI + AMU*AMU
         BZ = 2D0*AXI*(ALAM - ATAU)
         CZ = (ALAM-ATAU)*(ALAM-ATAU) - AMU*AMU
         DZ = BZ*BZ - 4D0*AZ*CZ
         IF (DZ.LT.0D0) THEN
            ZMIN = 0D0
            ZMAX = 0D0
            RETURN
         ELSE
            IF (BZ.GT.0D0) THEN
               ZMIN = (- BZ - DSQRT(DZ))/2D0/AZ
               ZMAX = 2D0*CZ/(-BZ - DSQRT(DZ))
               RETURN
            ELSE
               ZMAX = (- BZ + DSQRT(DZ))/2D0/AZ
               ZMIN = 2D0*CZ/(-BZ + DSQRT(DZ))
               RETURN
            ENDIF
         ENDIF
      ENDIF
      ENDIF
C
C-----
      IF (AMU.EQ.0D0) THEN
        IF(AXI.GT.0D0) THEN
           ZMAX = 1D0
           ZMIN = DMAX1(-1D0,(ATAU-ALAM)/AXI)
           IF (ZMIN.GT.1D0) ZMIN = 1D0
           RETURN
         ELSEIF(AXI.LT.0D0) THEN
           ZMIN =-1D0
           ZMAX = DMIN1(1D0,(ATAU-ALAM)/AXI)
           IF (ZMAX.LT.-1D0) ZMAX =-1D0
           RETURN
         ELSE
           IF (ALAM.GE.ATAU) THEN
              ZMIN = -1D0
              ZMAX = 1D0
              RETURN
            ELSE
              IPHSPC=1
            ENDIF
         ENDIF
      ENDIF
      IF (AMU.LT.0D0)THEN
          IPHSPC=1
      ENDIF
      END
