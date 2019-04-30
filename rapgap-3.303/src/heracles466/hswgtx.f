C
C***********************************************************************
C
C***********************************************************************
C
C   GIVE EXTERNAL WEIGHT TO EVENT
C
      SUBROUTINE HSWGTX(X,Y,IACPT)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSWGTC/ IWEIGS
      LOGICAL LFIRST
      DATA LFIRST/.TRUE./
 
      IF (LFIRST) THEN
        Q20=100D0
        LFIRST=.FALSE.
      ENDIF
 
      IACPT=1
      WI=1D0
      IF (IWEIGS.EQ.1) THEN
        WI=X
      ELSEIF (IWEIGS.EQ.2) THEN
        Q2=X*Y*SP
        IF (Q2.LT.Q20) THEN
          WI=Q2/Q20
          ELSE
          WI=1D0
        ENDIF
      ELSE
        RETURN
      ENDIF
      IF (HSRNDM(-1).LE.WI) RETURN
      IACPT=0
      RETURN
      END
