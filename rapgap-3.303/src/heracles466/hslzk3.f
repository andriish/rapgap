C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   LIMITS FOR COS(THETA(K,Q)) FOR KQ-PARAMETRIZATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSLZK3(ZMIN,ZMAX)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSCMS1/ COSQ,SINQ
      COMMON /HSPSPC/ IPHSPC
C
      IPHSPC=0
      DEPS = DELTA/OMEGA*( PQH*EH + PH*EQH)
      ALAM = PQH*EEL + PH*EQ
      TAU  = PQH*( PS*COSQ - PQ ) + PH*PQ
      F2U = -TAU - ALAM + DEPS
      F2O =  TAU - ALAM + DEPS
C
      IF ( (F2U.LT.0D0).AND.(F2O.LT.0D0) ) THEN
         ZMIN = -1D0
         ZMAX =  1D0
         RETURN
      ENDIF
C
      EZ = F2U*F2O
      IF ( EZ.LE.0D0 ) THEN
         PPP4 = PQH*PQH*PS*PS*SINQ*SINQ
         DZ = 4D0*PPP4*( PPP4 - EZ )
         AZ = TAU*TAU + PPP4
         BZ = 2D0*TAU*(DEPS - ALAM)
         CZ = (DEPS-ALAM)*(DEPS-ALAM) - PPP4
         IF (TAU.GT.0D0) THEN
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
         PPP4 = PQH*PQH*PS*PS*SINQ*SINQ
         DZ = 4D0*PPP4*( PPP4 - EZ )
         IF (DZ.LT.0D0) THEN
            IPHSPC=1
            RETURN
         ELSE
         AZ = TAU*TAU + PPP4
         BZ =  2D0*TAU*(DEPS - ALAM)
         CZ = (DEPS - ALAM)*(DEPS - ALAM) - PPP4
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
      END
